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
- Includes predefined NURBS shapes.
- Rotate and translate NURBS objects.

## Predefined NURBS Shapes

- Curves: Circle
- Surfaces: Tetragon
- Volumes: Hexahedron

## fpm dependency

If you want to use ForCAD as a dependency in your own fpm project,
you can easily include it by adding the following line to your `fpm.toml` file:

```toml
[dependencies]
forcad = {git="https://github.com/gha3mi/forcad.git"}
```

## How to run examples

To get started, follow these steps:

**Clone the repository:**

Clone the ForCAD repository from GitHub:

```shell
git clone https://github.com/gha3mi/forcad.git
cd forcad
```

### Using fpm


```shell
fpm run --example <file name excluding the .f90 extension>
```
Once the examples have been executed, `.vtk` files will be generated within the `vtk` directory. These files can then be visualized using tools such as [ParaView](https://www.paraview.org/).

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

This roadmap outlines upcoming features and enhancements for ForCAD. Contributions and suggestions are welcome!

- v0.2.0:
    - [x] Add `insert_knots()` method for curves, surfaces and volumes.
    - [x] Add `elevate_degree()` method for curves, surfaces and volumes.
    - [x] Add `derivative()` method for curves, surfaces and volumes.

- v0.3.0:
    - [x] Add `remove_knots()` method for curves, surfaces and volumes.
    - [ ] Add `reduce_degree()` method for curves, surfaces and volumes.
    - [x] Add `put_to_nurbs()` method for volumes.
    - [x] Get IGA elements connectivity.
    - [x] Add predefined shapes: `Circle`, `Tetragon`, `Hexahedron`.
    - [x] Add `rotate_Xc` and `rotate_Xg` methods for curves, surfaces and volumes.
    - [x] Add `translate_Xc` and `translate_Xg` methods for curves, surfaces and volumes.
    - [x] Add basic unit tests.
    - [x] Add simple examples.

- Future Tasks:
    - [ ] Add `morph()` method for morphing box.
    - [ ] Design a logo.
    - [ ] Add support binary `vtk` files.
    - [ ] Export to `IGES` format.
    - [ ] Add support for multiple patches.
    - [ ] Add extraction of piecewise Bezier objects from NURBS.
    - [ ] Add more unit tests.
    - [ ] Add more examples.
    - [ ] Add more predefined shapes.

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