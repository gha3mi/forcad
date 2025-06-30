[![GitHub](https://img.shields.io/badge/GitHub-ForCAD-blue.svg?style=social&logo=github)](https://github.com/gha3mi/forcad)
[![Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://gha3mi.github.io/forcad/)
[![fpm](https://github.com/gha3mi/forcad/actions/workflows/fpm.yml/badge.svg)](https://github.com/gha3mi/forcad/actions/workflows/fpm.yml)
[![doc](https://github.com/gha3mi/forcad/actions/workflows/doc.yml/badge.svg)](https://github.com/gha3mi/forcad/actions/workflows/doc.yml) 
[![codecov](https://codecov.io/gh/gha3mi/forcad/branch/dev/graph/badge.svg?token=N69NG6C86I)](https://codecov.io/gh/gha3mi/forcad)
[![License](https://img.shields.io/github/license/gha3mi/forcad?color=green)](https://github.com/gha3mi/forcad/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/778032800.svg)](https://zenodo.org/doi/10.5281/zenodo.10904447)

![](logo/logo.png)

**ForCAD**: A parallel Fortran library for geometric modeling using NURBS (Non-Uniform Rational B-Splines).

ForCAD supports **B-Spline**, **NURBS**, **Bezier**, and **Rational Bezier** curves, surfaces, and volumes.

## Main Features

- Parallelized using `OpenMP` and `do concurrent`.
- Create NURBS objects by specifying control points, weights and knots.
- Refine NURBS objects by inserting or removing knots and elevating degree.
- Compute basis functions and derivatives of NURBS objects.
- Obtain IGA elements connectivity.
- Obtain visualized elements connectivity and coordinates for geometry and control geometry.
- Mesh insertion into a NURBS object.
- Export NURBS objects to VTK files for visualization.
- Export of NURBS curves and surfaces to IGES format (volumes currently not supported).
- Includes predefined NURBS shapes: Circle, Half Circle, Tetragon, Hexahedron, 2D Ring, Half 2D Ring, 3D Ring, Half 3D Ring, C-shapes.
- Rotate and translate NURBS objects.
- Visualization using provided python PyVista scripts.

## Examples

<img src="https://github.com/gha3mi/forcad/raw/main/vtk/1.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/2.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/3.png" width="30%">

<img src="https://github.com/gha3mi/forcad/raw/main/vtk/4.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/5.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/vtk/6.png" width="30%">

<img src="https://github.com/gha3mi/forcad/raw/main/ppm/example_ppm1.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/ppm/example_ppm2.png" width="30%"> <img src="https://github.com/gha3mi/forcad/raw/main/ppm/example_ppm3.png" width="30%">

## Installation

### Requirements

- A Fortran compiler, such as [GNU Fortran](https://gcc.gnu.org/fortran/) (`gfortran`), [Intel Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html) (`ifx`) or [NVIDIA HPC SDK Fortran compiler](https://developer.nvidia.com/hpc-sdk) (`nvfortran`).
- The Fortran Package Manager [fpm](https://fpm.fortran-lang.org/).
- Optional: [PyVista](https://pyvista.org/) (Recommended) or [ParaView](https://www.paraview.org/) for visualization.

### Clone the repository

Clone the ForCAD repository from GitHub:

```shell
git clone https://github.com/gha3mi/forcad.git
cd forcad
```

### Install PyVista (Optional)

To install PyVista, run the following command:

```shell
pip install pyvista
```

### Running Examples with fpm

```shell
fpm run --example <file name excluding the .f90 extension>
```
After executing the examples, `.vtk` files will be generated in the `vtk` directory. To visualize these files, a `show()` method is provided which utilizes PyVista. Alternatively, other visualization tools like ParaView can also be used.

### Using ForCAD as a fpm Dependency

If you want to use ForCAD as a dependency in your own fpm project,
you can easily include it by adding the following line to your `fpm.toml` file:

```toml
[dependencies]
forcad = {git="https://github.com/gha3mi/forcad.git"}
```

### Precision Configuration

The library uses **double precision** (`real64`) by default for all real-valued computations. To change the precision, you can define one of the following preprocessor flags during compilation:

| Preprocessor Flag    | Fortran Kind             | Description               |
|----------------------|--------------------------|---------------------------|
| `REAL32`             | `selected_real_kind(6)`  | Single precision          |   
| `REAL64` *(default)* | `selected_real_kind(15)` | Double precision          |
| `REALXDP`            | `selected_real_kind(18)` | Extended double precision |
| `REAL128`            | `selected_real_kind(33)` | Quadruple precision       |

**Note**: The examples `example_ppm1.f90`, `example_ppm2.f90`, and `example_ppm3.f90` use the `ForColormap` library, which only supports `REAL64` precision.

#### Example: Building with double precision

```bash
fpm build --profile release --flag "-DREAL64"
```

## Status

<!-- STATUS:setup-fortran-conda:START -->
<!-- STATUS:setup-fortran-conda:END -->

## API documentation

The most up-to-date API documentation for the master branch is available
[here](https://gha3mi.github.io/forcad/).
To generate the API documentation for ForCAD using
[ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following
command:

```shell
ford README.md
```

## Roadmap

For a detailed roadmap outlining upcoming features and enhancements, please refer to [ROADMAP.md](https://github.com/gha3mi/forcad/blob/main/ROADMAP.md).

## Contributing

To contribute to ForCAD, please review the [CONTRIBUTING.md](https://github.com/gha3mi/forcad/blob/main/CONTRIBUTING.md).

## Citation

If you use ForCAD in your research, please cite it as follows:


```bibtex
@software{seyed_ali_ghasemi_2024_10904447,
  author       = {Ghasemi, S. A.},
  title        = {gha3mi/ForCAD},
  year         = 2024,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.10904447},
  url          = {https://doi.org/10.5281/zenodo.10904447}
}
```

## References

- Piegl, L., & Tiller, W. (1995). The NURBS Book. In Monographs in Visual Communications. Springer Berlin Heidelberg. [https://doi.org/10.1007/978-3-642-97385-7](https://doi.org/10.1007/978-3-642-97385-7)

- An Introduction to NURBS. (2001). Elsevier. [https://doi.org/10.1016/b978-1-55860-669-2.x5000-3](https://doi.org/10.1016/b978-1-55860-669-2.x5000-3)

- Sullivan et al., (2019). PyVista: 3D plotting and mesh analysis through a streamlined interface for the Visualization Toolkit (VTK). Journal of Open Source Software, 4(37), 1450, https://doi.org/10.21105/joss.01450

- Ahrens, James, Geveci, Berk, Law, Charles, ParaView: An End-User Tool for Large Data Visualization, Visualization Handbook, Elsevier, 2005, ISBN-13: 9780123875822