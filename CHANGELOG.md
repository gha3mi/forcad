# Changelog

All notable changes to [ForCAD](https://github.com/gha3mi/forcad) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),

## [Unreleased]

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

### Removed

- Excluded macOS from CI due to a problem with fpm.