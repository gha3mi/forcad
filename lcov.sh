#!/usr/bin/env bash
#
# Build and run the GFortran debug tests with coverage instrumentation, then
# generate line, function, and branch reports for production sources only.
#
# Output:
#   coverage/lcov/index.html
#
# Requirements:
#   fpm, gfortran, lcov, genhtml
#
# Usage:
#   ./lcov.sh

set -euo pipefail

readonly ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly OUTDIR="$ROOT/coverage/lcov"
readonly INFO="$OUTDIR/lcov.info"

for tool in fpm gfortran lcov genhtml; do
    command -v "$tool" >/dev/null 2>&1 || {
        printf 'Required tool not found: %s\n' "$tool" >&2
        exit 127
    }
done

cd "$ROOT"

# Keep downloaded dependencies so coverage remains reproducible offline.
fpm clean --skip --verbose
fpm test \
    --compiler gfortran \
    --profile debug \
    --flag "-DFOR_DEBUG --coverage -O0 -g" \
    --link-flag "--coverage" \
    --verbose

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

# GFortran-generated type copy/final symbols are not ForCAD source procedures.
lcov \
    --capture \
    --base-directory "$ROOT" \
    --directory "$ROOT/build" \
    --include "$ROOT/src/*" \
    --exclude '*/dependencies/*' \
    --exclude '*/test/*' \
    --erase-functions '.*_MOD___(copy|final).*' \
    --rc external=0 \
    --rc function_coverage=1 \
    --rc branch_coverage=1 \
    --rc geninfo_unexecuted_blocks=1 \
    --rc geninfo_adjust_testname=1 \
    --no-checksum \
    --output-file "$INFO" \
    --quiet

lcov --summary "$INFO" --rc branch_coverage=1

# LCOV 2.0-1 can misformat `lcov --list`; print the tracefile totals directly.
awk '
    function rate(hit, found) { return found > 0 ? 100.0*hit/found : 0.0 }
    BEGIN {
        printf "%-34s %14s %14s %14s\n", "Source", "Lines", "Functions", "Branches"
    }
    /^SF:/  { file=$0; sub(/^SF:.*\//, "", file) }
    /^LF:/  { lf=substr($0,4)+0 }
    /^LH:/  { lh=substr($0,4)+0 }
    /^FNF:/ { fnf=substr($0,5)+0 }
    /^FNH:/ { fnh=substr($0,5)+0 }
    /^BRF:/ { brf=substr($0,5)+0 }
    /^BRH:/ { brh=substr($0,5)+0 }
    /^end_of_record/ {
        printf "%-34s %6.1f%% %5d %6.1f%% %5d %6.1f%% %5d\n", \
            file, rate(lh,lf), lf, rate(fnh,fnf), fnf, rate(brh,brf), brf
    }
' "$INFO"

genhtml \
    "$INFO" \
    --title "ForCAD coverage report" \
    --header-title "ForCAD coverage" \
    --rc genhtml_dark_mode=0 \
    --rc legend=1 \
    --rc genhtml_hierarchical=1 \
    --rc genhtml_sort=1 \
    --rc function_coverage=1 \
    --rc branch_coverage=1 \
    --rc genhtml_highlight=1 \
    --frames \
    --output-directory "$OUTDIR" \
    --quiet

printf 'HTML report: %s/index.html\n' "$OUTDIR"
