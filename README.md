[![GitHub](https://img.shields.io/badge/GitHub-ForCAD-blue.svg?style=social&logo=github)](https://github.com/gha3mi/forcad)
[![Version](https://img.shields.io/github/v/tag/gha3mi/forcad?label=version&sort=semver)](https://github.com/gha3mi/forcad/releases)
[![Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://gha3mi.github.io/forcad/)
[![ForCAD CI](https://github.com/gha3mi/forcad/actions/workflows/CI-CD.yml/badge.svg?branch=main)](https://github.com/gha3mi/forcad/actions/workflows/CI-CD.yml)
[![codecov](https://codecov.io/gh/gha3mi/forcad/branch/main/graph/badge.svg?token=N69NG6C86I)](https://codecov.io/gh/gha3mi/forcad)
[![License](https://img.shields.io/github/license/gha3mi/forcad?color=green)](https://github.com/gha3mi/forcad/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/778032800.svg)](https://zenodo.org/doi/10.5281/zenodo.10904447)

![](logo/logo.png)

**ForCAD**: A parallel Fortran library for geometric modeling using NURBS (Non-Uniform Rational B-Splines).

ForCAD supports **B-Spline**, **NURBS**, **Bezier** and **Rational Bezier** curves, surfaces and volumes.

## Table of Contents
- [Table of Contents](#table-of-contents)
- [Main Features](#main-features)
- [Examples](#examples)
- [Installation](#installation)
  - [Using Pixi](#using-pixi)
  - [Compiler environments](#compiler-environments)
  - [Using fpm](#using-fpm)
    - [Running Examples with fpm](#running-examples-with-fpm)
    - [Using ForCAD as a fpm Dependency](#using-forcad-as-a-fpm-dependency)
- [Configuration](#configuration)
  - [Do Concurrent Support](#do-concurrent-support)
  - [Precision Configuration](#precision-configuration)
- [Dependencies](#dependencies)
- [CI Status](#ci-status)
- [API documentation](#api-documentation)
- [Contributing](#contributing)
- [Citation](#citation)
- [References](#references)


## Main Features

- Parallelized using `do concurrent`.
- Create NURBS objects by specifying control points, weights and knots.
- Refine NURBS objects by inserting or removing knots and elevating degree.
- Compute analytical basis functions and their first and second derivatives for NURBS and B-Spline objects.
- Generation of IGA-compatible element connectivity and shape functions.
- Obtain visualized elements connectivity and coordinates for geometry and control geometry.
- Mesh insertion into a NURBS object.
- Export NURBS objects to VTK files for visualization.
- Export of NURBS curves and surfaces to IGES format (volumes currently not supported).
- Includes predefined NURBS shapes: Circle, Half Circle, Tetragon, Hexahedron, 2D Ring, Half 2D Ring, 3D Ring, Half 3D Ring, C-shapes.
- Rotate and translate NURBS objects.
- Least squares fitting for NURBS curves, surfaces and volumes.
- Numerical integration of: NURBS curve length, NURBS surface area and NURBS volume.

## Examples

Below are some sample outputs from the `examples` directory:

<img src="https://github.com/gha3mi/forcad/raw/main/vtk/1.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/2.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/3.png" width="30%">

<img src="https://github.com/gha3mi/forcad/raw/main/vtk/4.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/5.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/6.png" width="30%">

<img src="https://github.com/gha3mi/forcad/raw/main/vtk/example_bend_pipe90.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/example_bend_pipe270.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/example_bend_pipe360.png" width="30%">

<img src="https://github.com/gha3mi/forcad/raw/main/vtk/example_helix_pipe.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/example_twist_taper.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/example_toroidal_pipe.png" width="30%">

<img src="https://github.com/gha3mi/forcad/raw/main/vtk/example_plate_hole_4_3d.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/example_plate_hole_1_3d.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/example_plate_hole_2_2d.png" width="30%">

<img src="https://github.com/gha3mi/forcad/raw/main/ppm/example_ppm1.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/ppm/example_ppm2.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/ppm/example_ppm3.png" width="30%">

<img src="https://github.com/gha3mi/forcad/raw/main/vtk/lsq_fit_bspline_3d.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/poisson_iga_solver_2d.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/poisson_iga_solver_3d.png" width="30%">

## Installation

### Using Pixi

Install [Pixi](https://pixi.sh/) and run these commands from the checkout.

```sh
pixi run build --profile release
pixi run test --profile debug
pixi run start
```

Pixi installs fpm 0.14 and GNU Fortran by default. The core builds without Qt.
ForCAD Studio uses GNU Fortran and qt-fortran. Pixi copies
`fpm-studio.toml` into `build/studio/` before building. Run
`pixi run build-studio debug` to build a Debug executable in
`build/studio-install/bin`. The `start` task builds and launches Release.

macOS builds need the Xcode command-line tools. On Windows, ifx, Flang and
ForCAD Studio need Visual Studio Build Tools with the C++ workload.

### Compiler environments

Choose the compiler with `-e` and the build profile with `--profile`.

| Environment | Compiler | Platform |
| --- | --- | --- |
| `default` | GNU Fortran | Linux, macOS (Apple Silicon and Intel), Windows |
| `intel` | Intel ifx | Linux, Windows |
| `flang` | LLVM Flang | Linux, Windows |
| `nvhpc` | NVIDIA nvfortran | Linux, with NVIDIA HPC SDK already on the path |

```sh
pixi run -e intel build --profile release
pixi run -e flang test --profile debug
pixi run lint
```

Pixi provisions the compiler and tools for each environment. The `nvhpc`
environment uses the compiler from the NVIDIA HPC SDK and installs fpm.
Compiler flags for Debug and Release live in `fpm.toml`.

### Using fpm

#### Running Examples with fpm

```shell
fpm run --example <file name excluding the .f90 extension> --compiler gfortran --profile release --flag "-ftree-parallelize-loops=8 -march=native"
```
Examples that export VTK write their files to the `vtk` directory. Open those files in ParaView or another VTK viewer.


#### Using ForCAD as a fpm Dependency

If you want to use ForCAD as a dependency in your own fpm project,
you can easily include it by adding the following line to your `fpm.toml` file:

```toml
[dependencies]
forcad = {git="https://github.com/gha3mi/forcad.git"}
```

## Configuration

### Do Concurrent Support

Compiler flags for enabling `do concurrent` parallelism:

| Compiler      | Flag(s)                                             |
| ------------- | --------------------------------------------------- |
| `gfortran`    | `-fopenmp -ftree-parallelize-loops=n`               |
| `ifx`         | `-qopenmp -fopenmp-target-do-concurrent`            |
| `nvfortran`   | `-stdpar=multicore,gpu -Minfo=stdpar,accel`         |
| `flang(-new)` | `-fopenmp -fdo-concurrent-to-openmp=[host\|device]` |
| `lfortran`    | ?                                                   |

*flang-new 21.1.0: warning: Mapping `do concurrent` to OpenMP is still experimental.*

Compiler flags can be passed to fpm using the `--flag` option, for example:

```shell
fpm build --flag "-stdpar=multicore,gpu -Minfo=stdpar,accel"
```

Store recurring compiler flags in the profiles in `fpm.toml`.

### Precision Configuration

The library uses **double precision** (`real64`) by default for all real-valued computations. To change the precision, you can define one of the following preprocessor flags during compilation:

| Preprocessor Flag    | Fortran Kind             | Description               |
| -------------------- | ------------------------ | ------------------------- |
| `REAL32`             | `selected_real_kind(6)`  | Single precision          |
| `REAL64` *(default)* | `selected_real_kind(15)` | Double precision          |
| `REALXDP`            | `selected_real_kind(18)` | Extended double precision |
| `REAL128`            | `selected_real_kind(33)` | Quadruple precision       |

**Note**: The examples `example_ppm1.f90`, `example_ppm2.f90` and `example_ppm3.f90` use the `ForColormap` library, which only supports `REAL64` precision.

Example: Building with double precision

```bash
fpm build --profile release --flag "-DREAL64"
```

## Dependencies

```mermaid
flowchart LR
  N1[forcad]
  N2[forIGES]
  N3[fordebug]
  N4[forunittest]
  N5[forimage]
  N6[forcolormap]
  N7[fortime]
  N8[FACE]

  click N1 href "https://github.com/gha3mi/forcad" "A parallel Fortran library for geometric modeling using NURBS"
  click N2 href "https://github.com/rweed/forIGES" "Modern Fortran Library for Reading and Writing IGES CAD Files"
  click N3 href "https://github.com/gha3mi/fordebug" "A Fortran library for handling errors, warnings and info messages with debugging support in pure procedures"
  click N4 href "https://github.com/gha3mi/forunittest" "A Fortran library for unit testing"
  click N5 href "https://github.com/gha3mi/forimage" "A Fortran library for processing and editing images and managing colors"
  click N6 href "https://github.com/vmagnin/forcolormap" "A Fortran colormap library"
  click N7 href "https://github.com/gha3mi/fortime" "A Fortran library for measuring elapsed time, DATE_AND_TIME time, CPU time, OMP time and MPI time"
  click N8 href "https://github.com/szaghi/FACE" "Fortran Ansi Colors (and Styles) Environment"

  %% core dependencies
  N1 --> N2
  N1 --> N3
  N3 --> N7
  N7 --> N8

  %% example dependency (dashed line)
  N1 -.->|example| N6
  N6 --> N5

  %% test dependency (dotted line)
  N1 ...->|test| N4
  N4 --> N8
```
Graph generated with [fpm-deps](https://github.com/ivan-pi/fpm-deps), modified to include example and test dependencies.

## CI Status

CI runs the core tests with fpm in Debug and Release and checks that Studio
renders Home and the bundled project views.

| Platform | Core compilers | Studio smoke check |
| --- | --- | --- |
| Linux | GNU Fortran, Intel ifx, LLVM Flang | GNU Fortran |
| macOS (Apple Silicon and Intel) | GNU Fortran | GNU Fortran |
| Windows | GNU Fortran, Intel ifx, LLVM Flang | GNU Fortran |

## API documentation

The most up-to-date API documentation for the master branch is available
[here](https://gha3mi.github.io/forcad/).
To generate the API documentation for ForCAD using
[ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following
command:

```shell
pixi run docs
```

## Contributing

To contribute to ForCAD, please review the [CONTRIBUTING.md](https://github.com/gha3mi/forcad/blob/main/CONTRIBUTING.md).

## Citation

If you use ForCAD in your research, please cite it as follows:


```bibtex
@software{seyed_ali_ghasemi_2025_10904447,
  author       = {Ghasemi, S. A.},
  title        = {gha3mi/ForCAD},
  year         = {2025},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.10904447},
  url          = {https://doi.org/10.5281/zenodo.10904447}
}
```

## References

- Piegl, L., & Tiller, W. (1995). The NURBS Book. In Monographs in Visual Communications. Springer Berlin Heidelberg. [https://doi.org/10.1007/978-3-642-97385-7](https://doi.org/10.1007/978-3-642-97385-7)

- An Introduction to NURBS. (2001). Elsevier. [https://doi.org/10.1016/b978-1-55860-669-2.x5000-3](https://doi.org/10.1016/b978-1-55860-669-2.x5000-3)

- Ahrens, James, Geveci, Berk, Law, Charles, ParaView: An End-User Tool for Large Data Visualization, Visualization Handbook, Elsevier, 2005, ISBN-13: 9780123875822
