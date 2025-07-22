## [v0.10.0](https://github.com/gha3mi/forcad/compare/v0.9.0...v0.10.0) - 2025-07-22


### Features

* feat: add CMake ([e38acd0b7](https://github.com/gha3mi/forcad/commit/e38acd0b7371897a2fdbd1e96d5702e8c9b25f4f)) by [@gha3mi](https://github.com/gha3mi)

### Fixes

* fix: use VERSION file and add description ([8bc700155](https://github.com/gha3mi/forcad/commit/8bc700155af83d8b8fde8657b6a9f64d15f36678)) by [@gha3mi](https://github.com/gha3mi)
* fix: update comments ([357d6b14f](https://github.com/gha3mi/forcad/commit/357d6b14fd351af2f91ae70bd594c2648382b6dc)) by [@gha3mi](https://github.com/gha3mi)

### Others

* Update README.md status table [ci skip] (#40) ([19349cccb](https://github.com/gha3mi/forcad/commit/19349cccb4002fc7e455c6a33d32be38d8a0b78c)) by [@gha3mi](https://github.com/gha3mi)


### Contributors
- [@gha3mi](https://github.com/gha3mi)



Full Changelog: [v0.9.0...v0.10.0](https://github.com/gha3mi/forcad/compare/v0.9.0...v0.10.0)

## [v0.9.0](https://github.com/gha3mi/forcad/compare/v0.8.0...v0.9.0) - 2025-07-21


### Features

* feat: implement new CI/CD workflow using setup-fortran-conda ([c97459817](https://github.com/gha3mi/forcad/commit/c97459817a647d5294cbb759de46d56362718ac8)) by [@gha3mi](https://github.com/gha3mi)
* feat: update CI/CD workflow and add fpm.rsp configuration for multiple platforms ([03c80752f](https://github.com/gha3mi/forcad/commit/03c80752f3603abb9e832201a6d8be57259c6fef)) by [@gha3mi](https://github.com/gha3mi)
* feat: add release.sh automation script ([9158e6f9b](https://github.com/gha3mi/forcad/commit/9158e6f9bc80157a7d93323cc880b5332f059fa5)) by [@gha3mi](https://github.com/gha3mi)
* feat: optimize B-spline basis and derivatives ([7c749cbd2](https://github.com/gha3mi/forcad/commit/7c749cbd29833ca2606ed7ffa9d09657ef12d091)) by [@gha3mi](https://github.com/gha3mi)

### Fixes

* fix: update condition for README.md status table update job ([ebd39f89c](https://github.com/gha3mi/forcad/commit/ebd39f89cec6df1d0dab781557614d0fa8791b1b)) by [@gha3mi](https://github.com/gha3mi)
* fix: remove extra-packages from CI/CD workflow for all platforms ([c153fe1f2](https://github.com/gha3mi/forcad/commit/c153fe1f236cb0026bdc88d03808d5d99d18c626)) by [@gha3mi](https://github.com/gha3mi)
* fix(CI): install llvm-openmp for flang compiler ([82db44e36](https://github.com/gha3mi/forcad/commit/82db44e36fd89e2ec7f066c422bd9df58bafcc36)) by [@gha3mi](https://github.com/gha3mi)
* fix: add optimization flags for gfortran ([14a6c7cbd](https://github.com/gha3mi/forcad/commit/14a6c7cbd7c00135d48e17fb554f406f0ae041ad)) by [@gha3mi](https://github.com/gha3mi)
* fix: use setup-fortran-conda for codecov workflow ([d12299fbb](https://github.com/gha3mi/forcad/commit/d12299fbbf5d68c01fe392fee386452aaeac9d46)) by [@gha3mi](https://github.com/gha3mi)
* fix: update citation year in README.md to 2025 ([8392b81f2](https://github.com/gha3mi/forcad/commit/8392b81f2c6a5e4bd57cb19510903ae27af13928)) by [@gha3mi](https://github.com/gha3mi)

### Others

* update README.md status table ([#36](https://github.com/gha3mi/forcad/pull/36)) by [@gha3mi](https://github.com/gha3mi)
* refactor: use setup-fortran-conda ([#34](https://github.com/gha3mi/forcad/pull/34)) by [@gha3mi](https://github.com/gha3mi)
* refactor: move dependencies to example sections in fpm.toml ([#33](https://github.com/gha3mi/forcad/pull/33)) by [@gha3mi](https://github.com/gha3mi)
* Update README.md status table [ci skip] (#37) ([60eef8c54](https://github.com/gha3mi/forcad/commit/60eef8c544838b9d97b8c133fd6267fe92549917)) by [@gha3mi](https://github.com/gha3mi)
* Update README.md status table [ci skip] (#38) ([f76e6d2dd](https://github.com/gha3mi/forcad/commit/f76e6d2ddefeaaf83846ba0e5e33ed48522f4a4b)) by [@gha3mi](https://github.com/gha3mi)


### Contributors
- [@gha3mi](https://github.com/gha3mi)



Full Changelog: [v0.8.0...v0.9.0](https://github.com/gha3mi/forcad/compare/v0.8.0...v0.9.0)

# Changelog

All notable changes to [ForCAD](https://github.com/gha3mi/forcad) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),

## [0.8.0]

### Added

* New module procedure interfaces.

### Changed

* Replaced OpenMP loops with `do concurrent`.
* Moved external procedures into their corresponding modules.
* Replaced several allocatable arrays with fixed-size arrays.
* Updated the Fortitude ignore list.
* Renamed variable `dim` to `d`.

### Fixed

* Corrected variable usage in NURBS derivative calls.
* Removed unused procedures from module imports.
* Updated derivative calls in NURBS tests.
* Removed unused variable.
* Fixed line truncation issues.

## [0.7.0]

### Features

* Add initial environment configuration in `environment.yml`
* Add pre-commit configuration for Fortitude hooks
* Add `extra.fortitude.check` section to `fpm.toml` for improved static checks
* Add IGES export functionality for NURBS curves and surfaces
* Add VSCode configuration for running example programs
* Add type-bound procedure `export_Xth()`

### Bug Fixes

* Update ignore list in `extra.fortitude.check`
* Explicitly add `implicit none` to prevent implicit typing
* Correct filename of Dependabot configuration (`dependabot.yml`)
* Remove redundant `isinf/isnan` calls from `forcad_utils`
* Fix IGES file extensions
* Match evaluation points (`Xt`) to knot vector domain
* Resolve import conflict by importing only `timer` from `fortime`

### Documentation

* Update roadmap with GUI implementation suggestion using OpenGL
* Refactor and update Ford documentation
* Mark binary VTK support as completed in `ROADMAP.md`
* Replace outdated references (`ifx/ifort`) with `ifx` in `README.md`

### Chores (CI/CD & Dependency Improvements)

* Add Dependabot configuration for automatic GitHub Actions updates
* Update Codecov action to v5
* Update example PPM files to consistently use `wp` kind real numbers
* Improve the `export_vtk_legacy` function
* Exclude Intel-classic compilers from CI tests on Windows/macOS
* Include Intel compilers in GitHub Actions CI (`fpm.yml`)
* Explicitly update fpm dependencies
* Remove trailing whitespace from source files
* Include export functionality in unit tests

## [0.6.1]

### Changed
- Improved `nearest_point2()` procedures.
- Updated and improved unit tests.
- Fixed bugs in `set4` and `compute_dTgc_bspline_2d_vector` in `nurbs_surface`.

## [0.6.0]

### Added

- Added generic `get_nc()` method to `nurbs_surface` and `nurbs_volume` derived types.
- Added `det`, `inv`, `dyad` and `gauss_leg` procedures in the `forcad_utils` module.
- Added procedure `set1a` to the `nurbs_curve` derived type.
- Added procedure `set4` to `nurbs_curve`,  `nurbs_area` and `nurbs_volume` derived types.
- Added optional input variable `elem` to `derivative_scalar` procedures.
- Added `ansatz` procedures to compute shape functions, derivatives of shape functions and (dV, dA, dL).
- Added `cmp_length()` to compute the length of a NURBS curve.
- Added `cmp_area()` to compute the area of a NURBS surface.
- Added `cmp_volume()` to compute the volume of a NURBS volume.
- Added examples for `cmp_length()`, `cmp_area()`, and `cmp_volume()`.

## [0.5.1]

### Changed

- Fixed initial guess in `nearest_point2()` procedures

## [0.5.0]

### Added

- Added `CHANGELOG.md` file.
- Added `CONTRIBUTING.md` file.
- Add NURBS surface to PPM conversion examples.
- Utilized ForUnitTest for testing.
- Added Support for `OpenMP` and `do concurrent`.
- Implemented memory cleanup in the examples and tests.
- Added cleanup for colormap type in the ppm examples.
- Added screenshots to the README file.
- Included nvidia compiler in the CI.
- Added missing allocation for `Tgc` in the `compute_Xg_nurbs_1d` subroutine.
- Added `nearest_point()` procedures to the `nurbs_curve`, `nurbs_surface` and `nurbs_volume`.
- Added examples for the `nearest_point()` method.
- Added `Xt` to the `nurbs_surface` and `nurbs_volume` derived types.
- Added codecoverage to the CI.
- Added allocate(Tgci) and allocate(dTgci) to `compute_Tgc_nurbs_*d()` and `compute_dTgc_nurbs_*d()` subroutines.
- Added `cmp_elemFace_Xc_vis()` to extract element connectivity of the faces of control geometry.
- Added `cmp_elemFace_Xg_vis()` to extract element connectivity of the faces of visualization geometry points.
- Added `cmp_elemFace()` to extract element connectivity of faces.
- Added `cmp_degreeFace()` to extract the degrees of faces.
- Updated `example_volume_1` to use `cmp_elemFace()` and `cmp_degreeFace()`.
- Added `cmp_Xg()` to evaluate the geometry points.
- Added generic method `derivative2()` to compute the second derivative of a NURBS objects.
- Added `nearest_point2()` to compute the nearest point on a NURBS object using optimization.
- Added Interfaces.
- Added new tests: fdm_curve.f90, fdm_surface.f90, fdm_volume.f90.
- Updated nearest_point_* examples to use the new `nearest_point2()` method.

### Changed

- Updated `README.md` file.
- Updated tests to use ForUnitTest.
- Added ParaView to the list of References in the README file.
- Used `matmul` instead of `dot_product` in the `put_to_nurbs` subroutine.
- Removed unused variables in the `put_to_nurbs` subroutine.
- Removed unused variables in the `example_ppm3.f90` file.
- Renamed program name in the example `put_to_nurbs.f90` file.
- Parallelized calculation of distance in the `nearest_point()` method.
- Updated `ROADMAP.md` file.
- Updated unit tests for `nurbs_curve`.
- Improved `elemConn_C0` to support higher order elements and check for the number of elements.
- Renamed test files to match the module name.
- Used `do concurrent` within `kron` and `basis_bspline` functions.
- Converted the `basis_bspline_der` function to a subroutine.
- Fixed NURBS derivative calculations.
- Made `basis` and `derivative` generic methods.
- Fixed initial guess in `nearest_point2()` procedures

### Removed

- Excluded macOS from CI due to a problem with fpm.
- Removed `intel-classic` from CI (`ifort` is deprecated).
