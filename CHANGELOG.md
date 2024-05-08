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

### Changed

- Updated `README.md` file.
- Updated tests to use ForUnitTest.
- Added ParaView to the list of References in the README file.
- Used `matmul` instead of `dot_product` in the `put_to_nurbs` subroutine.

### Removed

- Excluded macOS from CI due to a problem with fpm.
- Excluded NVIDIA from CI due to a problem with ForColormap or a bug with the NVIDIA compiler.
