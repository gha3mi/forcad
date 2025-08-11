#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# release.sh — Automate changelog generation, version bumping, tagging, and release.
#
# Part of the setup-fortran-conda project: https://github.com/fortran-lang/setup-fortran-conda
#
# This script streamlines the release process by:
#   - Analyzing Git commit messages and merged pull requests (PRs) to determine the next semantic version bump (major, minor, or patch).
#   - Generating a well-structured CHANGELOG.md section using PR titles, commit messages, and contributor names.
#   - Updating the VERSION file with the new semantic version.
#   - Committing the above changes (CHANGELOG.md and VERSION) with meaningful commit messages.
#   - Optionally pushing changes to the repository and tagging the release (`vX.Y.Z`).
#   - Creating a GitHub Release with the generated changelog as release notes.
#   - Tagging the latest release with both the new version tag and a floating `latest` tag.
#
# It supports two optional modes:
#   --dry-run  : Preview all steps without writing files or pushing anything.
#   --local    : Commit changes locally, but skip pushing, tagging, and publishing a GitHub release.
#
# Requires:
#   - GitHub CLI (`gh`) — must be installed and authenticated.
#   - `jq` — used for parsing GitHub API responses.
#
# Author: Seyed Ali Ghasemi
# License: MIT
# ------------------------------------------------------------------------------
set -euo pipefail

print_help() {
    cat <<EOF
Usage: ./release.sh [OPTIONS]

Automatically generate a semantic version bump, update CHANGELOG.md and VERSION,
create a Git tag, and publish a GitHub Release based on recent commits and merged PRs.

Options:
  --dry-run     Preview the version bump and generated changelog.
                No files will be written or pushed.

  --local       Perform all steps locally: update CHANGELOG.md and VERSION,
                commit changes, but do NOT push, tag, or create a GitHub Release.

  --help        Show this help message and exit.

Requirements:
  - GitHub CLI ('gh') must be authenticated and installed
  - 'jq' must be installed

Example:
  ./release.sh --dry-run
EOF
}

DRY_RUN=false
LOCAL_ONLY=false
INCLUDE_PR_COMMITS=false

while [[ $# -gt 0 ]]; do
    case "$1" in
    --dry-run) DRY_RUN=true ;;
    --local) LOCAL_ONLY=true ;;
    --help)
        print_help
        exit 0
        ;;
    *)
        echo "❌ Unknown option: $1"
        print_help
        exit 1
        ;;
    esac
    shift
done

RELEASE_IGNORE_PATTERNS=(
    "^docs: update CHANGELOG"
    "^chore: update VERSION"
    "^Merge pull request"
)

should_ignore() {
    local msg="$1"
    for pattern in "${RELEASE_IGNORE_PATTERNS[@]}"; do
        if [[ "$msg" =~ $pattern ]]; then
            return 0
        fi
    done
    return 1
}

echo "🧪 Dry run mode: $DRY_RUN"
echo "🔧 Local-only mode: $LOCAL_ONLY"

cd "$(git rev-parse --show-toplevel)"

if ! command -v gh >/dev/null 2>&1; then
    echo "❌ GitHub CLI (gh) not found. Please install it: https://cli.github.com/"
    exit 1
fi

if ! gh auth status &>/dev/null; then
    echo "❌ GitHub CLI is not authenticated. Run: gh auth login"
    exit 1
fi

repo=$(gh repo view --json nameWithOwner -q .nameWithOwner)
previous_tag=$(git tag --sort=-v:refname | grep '^v' | head -n1)
previous_tag=${previous_tag:-v0.0.0}
echo "📌 Previous tag: $previous_tag"

commits=$(git log "$previous_tag"..HEAD --pretty=format:"%s")
if [[ -z "$commits" ]]; then
    echo "✅ No commits since last release."
    exit 0
fi

bump="patch"
if echo "$commits" | grep -qE "^feat(\(.+\))?!?: "; then bump="minor"; fi
if echo "$commits" | grep -qE "BREAKING CHANGE|!: "; then bump="major"; fi

IFS=. read -r major minor patch <<<"${previous_tag#v}"
case "$bump" in
major)
    major=$((major + 1))
    minor=0
    patch=0
    ;;
minor)
    minor=$((minor + 1))
    patch=0
    ;;
patch) patch=$((patch + 1)) ;;
esac
next_tag="v$major.$minor.$patch"
echo "🆕 Next version: $next_tag"

echo "🔍 Fetching commits between $previous_tag and HEAD..."
compare=$(gh api "/repos/$repo/compare/$previous_tag...HEAD")
commit_shas=$(echo "$compare" | jq -r '.commits[].sha')
all_commit_shas=()
while IFS= read -r sha; do all_commit_shas+=("$sha"); done <<<"$commit_shas"

echo "🔍 Fetching recent closed PRs from GitHub API..."
if ! raw_prs=$(gh api "/repos/$repo/pulls?state=closed&per_page=100" 2>&1); then
    echo "❌ Failed to fetch PRs from GitHub:"
    echo "$raw_prs"
    exit 1
fi

if ! echo "$raw_prs" | jq empty >/dev/null 2>&1; then
    echo "❌ GitHub API returned invalid JSON:"
    echo "$raw_prs"
    exit 1
fi
mapfile -t all_prs < <(echo "$raw_prs" | jq -c '.[]')
echo "ℹ️ Retrieved ${#all_prs[@]} PRs"

features=()
fixes=()
breaking=()
others=()
contributors=()

classify() {
    local msg="$1" line="$2"
    if [[ "$msg" =~ (!:|BREAKING CHANGE) ]]; then
        breaking+=("$line")
    elif [[ "$msg" =~ ^feat(\(.+\))?!?:\  || "$msg" =~ ^feat:\  ]]; then
        features+=("$line")
    elif [[ "$msg" =~ ^fix(\(.+\))?!?:\  || "$msg" =~ ^fix:\  ]]; then
        fixes+=("$line")
    else
        others+=("$line")
    fi
}

matched_prs=()
pr_commit_shas=()
for pr_json in "${all_prs[@]}"; do
    merged_at=$(echo "$pr_json" | jq -r '.merged_at')
    [[ "$merged_at" == "null" ]] && continue

    pr_number=$(echo "$pr_json" | jq -r '.number')
    mapfile -t pr_commits < <(gh api "/repos/$repo/pulls/$pr_number/commits" | jq -r '.[].sha')

    pr_matches_range=false
    for sha in "${pr_commits[@]}"; do
        if [[ " ${all_commit_shas[*]} " =~ " $sha " ]]; then
            pr_matches_range=true
            break
        fi
    done
    [[ "$pr_matches_range" == false ]] && continue

    title=$(echo "$pr_json" | jq -r '.title')
    number=$(echo "$pr_json" | jq -r '.number')
    author=$(echo "$pr_json" | jq -r '.user.login')

    if should_ignore "$title"; then
        echo "⏭️  Ignored PR: $title"
        continue
    fi

    matched_prs+=("$pr_json")
    pr_commit_shas+=("${pr_commits[@]}")

    line="* $title ([#${number}](https://github.com/$repo/pull/${number})) by [@$author](https://github.com/$author)"
    classify "$title" "$line"
    contributors+=("$author")

    if [[ "${INCLUDE_PR_COMMITS:-true}" == true ]]; then
        for sha in "${pr_commits[@]}"; do
            if [[ " ${all_commit_shas[*]} " =~ " $sha " ]]; then
                msg=$(git show -s --format="%s" "$sha")
                if should_ignore "$msg"; then continue; fi
                commit_author=$(gh api "/repos/$repo/commits/$sha" --jq '.author.login' 2>/dev/null || echo "")
                author_name=${commit_author:-$(git show -s --format="%an" "$sha")}
                short_sha=$(git rev-parse --short "$sha")
                subline="* [#${number}] $msg ([${short_sha}](https://github.com/$repo/commit/$sha)) by [@$author_name](https://github.com/$author_name)"
                classify "$msg" "$subline"
                contributors+=("$author_name")
            fi
        done
    fi
done

for sha in "${all_commit_shas[@]}"; do
    [[ " ${pr_commit_shas[*]} " =~ " $sha " ]] && continue
    msg=$(git show -s --format="%s" "$sha")
    if should_ignore "$msg"; then
        echo "⏭️  Ignored commit: $msg"
        continue
    fi

    commit_author=$(gh api "/repos/$repo/commits/$sha" --jq '.author.login' 2>/dev/null || echo "")
    author=${commit_author:-$(git show -s --format="%an" "$sha")}
    short_sha=$(git rev-parse --short "$sha")
    line="* $msg ([${short_sha}](https://github.com/$repo/commit/$sha)) by [@$author](https://github.com/$author)"

    classify "$msg" "$line"
    contributors+=("$author")
done

readarray -t unique_contributors < <(printf "%s\n" "${contributors[@]}" | grep -v '^null$' | sort -u)

today=$(date +%F)
url="https://github.com/$repo/compare/$previous_tag...$next_tag"

changelog_header="## [$next_tag]($url) - $today\n\n"
changelog_body=""
changelog_contrib=""

section() {
    local title="$1"
    shift
    local items=("$@")
    if [[ ${#items[@]} -gt 0 ]]; then
        changelog_body+=$'\n'"### $title"$'\n\n'
        for item in "${items[@]}"; do
            changelog_body+="$item"$'\n'
        done
    fi
}

section "Features" "${features[@]}"
section "Fixes" "${fixes[@]}"
section "Breaking Changes" "${breaking[@]}"
section "Others" "${others[@]}"

if [[ ${#unique_contributors[@]} -gt 0 ]]; then
    changelog_contrib+="\n\n### Contributors\n"
    for c in "${unique_contributors[@]}"; do
        changelog_contrib+="- [@$c](https://github.com/$c)\n"
    done
    changelog_contrib+="\n"
fi

changelog_tail="\n\nFull Changelog: [$previous_tag...$next_tag]($url)\n"
changelog_md="${changelog_header}${changelog_body}${changelog_contrib}${changelog_tail}"
changelog_release="${changelog_header}${changelog_body}${changelog_tail}"

echo -e "\n📝 Generated CHANGELOG:\n"
echo -e "$changelog_md"

# DRY RUN: Skip all changes
if [[ "$DRY_RUN" == true ]]; then
    echo "🧪 Dry run complete — no files written, no commits made, no tags created."
    exit 0
fi

# ✅ STEP 1: Confirm writing CHANGELOG.md
if [[ -f CHANGELOG.md && $(grep -c "## \[$next_tag\]" CHANGELOG.md) -gt 0 ]]; then
    echo "⚠️ Changelog for $next_tag already exists — skipping update."
else
    echo
    read -rp "📄 Write and commit CHANGELOG.md? [y/N] " confirm
    confirm="${confirm,,}"
    if [[ "$confirm" == "y" || "$confirm" == "yes" ]]; then
        echo -e "$changelog_md\n$(cat CHANGELOG.md 2>/dev/null)" >.CHANGELOG.md && mv .CHANGELOG.md CHANGELOG.md
        git add CHANGELOG.md
        git commit -m "docs: update CHANGELOG for $next_tag"
    else
        echo "⏭️ Skipping CHANGELOG.md update."
    fi
fi

# ✅ STEP 2: Confirm writing VERSION
echo
read -rp "📝 Write and commit VERSION file? [y/N] " confirm
confirm="${confirm,,}"
if [[ "$confirm" == "y" || "$confirm" == "yes" ]]; then
    echo "${next_tag#v}" >VERSION
    git add VERSION
    git commit -m "chore: update VERSION to ${next_tag}"
else
    echo "⏭️ Skipping VERSION file update."
fi

# ✅ STEP 3: Confirm pushing commits
if [[ "$LOCAL_ONLY" == false ]]; then
    echo
    read -rp "📤 Push commits to GitHub? [y/N] " confirm
    confirm="${confirm,,}"
    if [[ "$confirm" == "y" || "$confirm" == "yes" ]]; then
        git push
    else
        echo "⏭️ Skipping git push."
    fi
fi

# ✅ STEP 4: Confirm tagging release
if [[ "$LOCAL_ONLY" == false ]]; then
    echo
    read -rp "🏷️  Create and push Git tag $next_tag? [y/N] " confirm
    confirm="${confirm,,}"
    if [[ "$confirm" == "y" || "$confirm" == "yes" ]]; then
        if ! git rev-parse "$next_tag" >/dev/null 2>&1; then
            git tag "$next_tag"
            git push origin "$next_tag"
        else
            echo "⚠️ Tag $next_tag already exists — skipping."
        fi
        git tag -f latest "$next_tag"
        git push origin -f latest
    else
        echo "⏭️ Skipping Git tag and latest update."
    fi
fi

# ✅ STEP 5: Confirm creating GitHub Release
if [[ "$LOCAL_ONLY" == false ]]; then
    echo
    read -rp "🚀 Create GitHub Release $next_tag? [y/N] " confirm
    confirm="${confirm,,}"
    if [[ "$confirm" == "y" || "$confirm" == "yes" ]]; then
        gh release create "$next_tag" --title "$next_tag" --notes-file <(echo -e "$changelog_release")
        echo "✅ GitHub Release $next_tag created."
    else
        echo "⏭️ Skipping GitHub Release."
    fi
fi

echo "🎉 Done."