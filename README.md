[![GitHub](https://img.shields.io/badge/GitHub-ForCAD-blue.svg?style=social&logo=github)](https://github.com/gha3mi/forcad)
[![Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://gha3mi.github.io/forcad/)
[![fpm](https://github.com/gha3mi/forcad/actions/workflows/fpm.yml/badge.svg)](https://github.com/gha3mi/forcad/actions/workflows/fpm.yml)
[![doc](https://github.com/gha3mi/forcad/actions/workflows/doc.yml/badge.svg)](https://github.com/gha3mi/forcad/actions/workflows/doc.yml) 
[![License](https://img.shields.io/github/license/gha3mi/forcad?color=green)](https://github.com/gha3mi/forcad/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/778032800.svg)](https://zenodo.org/doi/10.5281/zenodo.10904447)

**ForCAD**: A Fortran library for Geometric Modeling using NURBS (Non-Uniform Rational B-Splines).

ForCAD supports **B-Spline**, **NURBS**, **Bezier**, and **Rational Bezier** curves, surfaces, and volumes.

<img alt="example_bezier" src="https://github.com/gha3mi/forcad/raw/main/vtk/example_bezier.png" width="750">

## Key Features

- Create NURBS objects by specifying control points, weights and knots.
- Refine NURBS objects by inserting or removing knots and elevating degree.
- Compute basis functions and derivatives of NURBS objects.
- Obtain IGA elements connectivity.
- Obtain visualized elements connectivity and coordinates for geometry and control geometry.
- Mesh insertion into a NURBS object.
- Export NURBS objects to VTK files for visualization.
- Includes predefined NURBS shapes: Circle, Tetragon, Hexahedron.
- Rotate and translate NURBS objects.

## Installation

### Reuirements

- A Fortran compiler, such as [GNU Fortran](https://gcc.gnu.org/fortran/) (`gfortran`), [Intel Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html) (`ifx/ifort`) or [NVIDIA HPC SDK Fortran compiler](https://developer.nvidia.com/hpc-sdk) (`nvfortran`).
- The Fortran Package Manager [fpm](https://fpm.fortran-lang.org/).

### Clone the repository

Clone the ForCAD repository from GitHub:

```shell
git clone https://github.com/gha3mi/forcad.git
cd forcad
```

### Running Examples with fpm

```shell
fpm run --example <file name excluding the .f90 extension>
```
Once the examples have been executed, `.vtk` files will be generated within the `vtk` directory. These files can then be visualized using tools such as [ParaView](https://www.paraview.org/).

### Using ForCAD as a fpm Dependency

If you want to use ForCAD as a dependency in your own fpm project,
you can easily include it by adding the following line to your `fpm.toml` file:

```toml
[dependencies]
forcad = {git="https://github.com/gha3mi/forcad.git"}
```

## API documentation

The most up-to-date API documentation for the master branch is available
[here](https://gha3mi.github.io/forcad/).
To generate the API documentation for ForCAD using
[ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following
command:

```shell
ford ford.yml
```

## Roadmap

For a detailed roadmap outlining upcoming features and enhancements, please refer to [ROADMAP.md](https://github.com/gha3mi/forcad/blob/main/ROADMAP.md).

## Contributing

Contributions to ForCAD are welcome!

- If you find any issues or would like to suggest improvements, please open an issue.
- If you've implemented new features, fixed bugs, or enhanced existing functionality, please consider submitting a pull request (PR).
- Please share your examples by submitting a pull request (PR).

## Citation

If you use ForCAD in your research, please cite it as follows:


```bibtex
@software{seyed_ali_ghasemi_2024_10938054,
  author       = {Seyed Ali Ghasemi},
  title        = {gha3mi/ForCAD},
  year         = 2024,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.10938054},
  url          = {https://doi.org/10.5281/zenodo.10938054}
}
```

## References

- Piegl, L., & Tiller, W. (1995). The NURBS Book. In Monographs in Visual Communications. Springer Berlin Heidelberg. [https://doi.org/10.1007/978-3-642-97385-7](https://doi.org/10.1007/978-3-642-97385-7)

- An Introduction to NURBS. (2001). Elsevier. [https://doi.org/10.1016/b978-1-55860-669-2.x5000-3](https://doi.org/10.1016/b978-1-55860-669-2.x5000-3)