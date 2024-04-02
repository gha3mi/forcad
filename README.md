[![GitHub](https://img.shields.io/badge/GitHub-ForCAD-blue.svg?style=social&logo=github)](https://github.com/gha3mi/forcad)
[![Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://gha3mi.github.io/forcad/)
[![fpm](https://github.com/gha3mi/forcad/actions/workflows/fpm.yml/badge.svg)](https://github.com/gha3mi/forcad/actions/workflows/fpm.yml)
[![doc](https://github.com/gha3mi/forcad/actions/workflows/doc.yml/badge.svg)](https://github.com/gha3mi/forcad/actions/workflows/doc.yml) 
[![License](https://img.shields.io/github/license/gha3mi/forcad?color=green)](https://github.com/gha3mi/forcad/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/778032800.svg)](https://zenodo.org/doi/10.5281/zenodo.10904447)

**ForCAD**: A Fortran library for Geometric Modeling.

ForCAD supports **B-Spline**, **NURBS**, **Bezier**, and **Rational Bezier** curves, surfaces, and volumes.

<img alt="example_bezier" src="https://github.com/gha3mi/forcad/raw/main/vtk/example_bezier.png" width="750">

## fpm dependency

If you want to use `ForCAD` as a dependency in your own fpm project,
you can easily include it by adding the following line to your `fpm.toml` file:

```toml
[dependencies]
forcad = {git="https://github.com/gha3mi/forcad.git"}
```

## How to run examples

To get started, follow these steps:

**Clone the repository:**

Clone the `ForCAD` repository from GitHub:

```shell
git clone https://github.com/gha3mi/forcad.git
```

```shell
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
To generate the API documentation for `ForCAD` using
[ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following
command:

```shell
ford ford.yml
```

## ToDo

- [ ] Add method `insert_knot()` for B-splines and NURBS.
- [ ] Add method `remove_knot()` for B-splines and NURBS.
- [ ] Add method `elevate_degree()` for Bezier, Rational Bezier, B-splines, and NURBS.
- [x] Add method `elevate_degree()` for Bezier and Rational Bezier curves.
- [ ] Add method `reduce_degree()` for Bezier, Rational Bezier, B-splines, and NURBS.
- [ ] Add method `derivative()`.
- [ ] Add support for multiple patches.
- [ ] Implement extraction of piecewise Bezier objects from NURBS.

## Contributing

Contributions to `ForCAD` are welcome!
If you find any issues or would like to suggest improvements, please open an issue.

## References

- [The NURBS Book](https://doi.org/10.1007/978-3-642-97385-7) by Les Piegl, Wayne Tiller
- [An Introduction to NURBS](https://doi.org/10.1016/B978-1-55860-669-2.X5000-3) by David F. Rogers