title: Build environments

Pixi selects the compiler and fpm environment. The core uses fpm 0.14 and the
compiler flags in the checkout's `fpm.toml`. Studio uses the GNU environment
and qt-fortran, with its own manifest at `fpm-studio.toml`.

| Pixi environment | Compiler | Available on |
| --- | --- | --- |
| `gnu` | GNU Fortran | Linux, macOS (Apple Silicon and Intel), Windows |
| `intel` | Intel ifx | Linux, Windows |
| `flang` | LLVM Flang | Linux, Windows |
| `nvhpc` | nvfortran from an installed NVIDIA HPC SDK | Linux |
| `studio` | GNU Fortran with Qt and a C++ toolchain | Linux, macOS (Apple Silicon and Intel), Windows |

Run `pixi run -e gnu test --profile debug` to test with GNU Fortran. Substitute `intel` or `flang` on a supported platform. The shorter `pixi run build` and `pixi run test` commands select GNU automatically. To build the Studio Debug executable in `build/studio-install/bin`, run `pixi run build-studio debug`. See [Launch ForCAD Studio](../how-to/install.html) for platform prerequisites.

[TOC]

## Dependencies

Pixi's [environment manifest](https://github.com/gha3mi/forcad/blob/dev/pixi.toml) specifies compiler and tool versions. The [library manifest](https://github.com/gha3mi/forcad/blob/dev/fpm.toml) lists the Fortran dependencies, including those required only by particular examples. The [Studio manifest](https://github.com/gha3mi/forcad/blob/dev/fpm-studio.toml) specifies the Qt bindings.

## Working precision

ForCAD requests 15 decimal digits by default. Use the library's `rk` kind in
your code. Set at most one macro for a build.

| Macro | Kind selection |
| --- | --- |
| `REAL32` | `selected_real_kind(6)` |
| `REAL64` (default) | `selected_real_kind(15)` |
| `REALXDP` | `selected_real_kind(18)` |
| `REAL128` | `selected_real_kind(33)` |

The compiler must support the requested precision.

## Generate documentation locally

Run `pixi run docs` from the checkout root. Open `docs/ford/index.html` to read the handbook and public API reference. The [guide and API index](../../search.html) works directly from disk. FORD accepts LaTeX in the guide text and Fortran comments. The browser loads MathJax for equations and Bootstrap for API navigation from a CDN.

To use full-text search locally, serve the built site from the checkout root. Open <http://127.0.0.1:8000/> in your browser. The server uses only the local machine.

```sh
pixi run -e docs python -m http.server 8000 --bind 127.0.0.1 --directory docs/ford
```

The search index is generated locally by [Pagefind](https://pagefind.app/). Its [license](pagefind-license.txt) is included with the handbook.
