## [v0.14.0](https://github.com/gha3mi/forcad/compare/v0.13.0...v0.14.0) - 2025-08-11


### Features

* [#67] feat: add kron, linspace and eye ([b91d98082](https://github.com/gha3mi/forcad/commit/b91d980828dc2cbba0cb1e92d68f85dc1e0b9db1)) by [@gha3mi](https://github.com/gha3mi)
* [#67] feat: add default case for matrix inversion ([5bad62d05](https://github.com/gha3mi/forcad/commit/5bad62d05174cd354ac0108f880e494484b38944)) by [@gha3mi](https://github.com/gha3mi)
* [#67] feat: add option to disable PyVista visualization via preprocessor flag ([c381bace3](https://github.com/gha3mi/forcad/commit/c381bace392b32ef7c0008ea1443f7b543b5aa9a)) by [@gha3mi](https://github.com/gha3mi)
* [#67] feat: add unit tests for forcad_utils module ([e78f50091](https://github.com/gha3mi/forcad/commit/e78f500912239ace0d25d7c84e636384e7fbd9b8)) by [@gha3mi](https://github.com/gha3mi)
* [#67] feat: add optional number of Gauss points ([3538af3ad](https://github.com/gha3mi/forcad/commit/3538af3ad6c720aeb2f5d70a9912119ee8e316e1)) by [@gha3mi](https://github.com/gha3mi)
* [#67] feat: update compute_dTgc calls to handle optional elem parameter ([ee9126887](https://github.com/gha3mi/forcad/commit/ee91268870a84274a94988fc5d047bf3d6bffa5d)) by [@gha3mi](https://github.com/gha3mi)
* [#67] feat: make compute_Tgc, compute_dTgc public ([a899fcb34](https://github.com/gha3mi/forcad/commit/a899fcb34b47ecbbfeee33ca0785be7b65f3336f)) by [@gha3mi](https://github.com/gha3mi)
* [#67] feat: add cmp_elem_Xth to NURBS objects ([8e84ac22e](https://github.com/gha3mi/forcad/commit/8e84ac22e6520c6d664928f07ac347ac3f0872cf)) by [@gha3mi](https://github.com/gha3mi)
* [#67] feat: replace local element connection handling with cmp_elem_Xth ([d2e2147fe](https://github.com/gha3mi/forcad/commit/d2e2147fec7ab26806753d7553ab8880d2d024f3)) by [@gha3mi](https://github.com/gha3mi)
* [#67] feat: add export_Xth_in_Xg for surface and volume ([859791c42](https://github.com/gha3mi/forcad/commit/859791c4276594f6b374f840f840cc788efc5dc8)) by [@gha3mi](https://github.com/gha3mi)

### Fixes

* [#67] fix: condition in knot multiplicity calculation ([6bf3cb8eb](https://github.com/gha3mi/forcad/commit/6bf3cb8ebd613424f43ef0bdbc626a469bf6c42a)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: test30, 40, 42 ([0de619e68](https://github.com/gha3mi/forcad/commit/0de619e68e3ad1a4efd7b7781ff67eb70618791d)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: use do loops for nvcompiler ([0bd396e9e](https://github.com/gha3mi/forcad/commit/0bd396e9e143355379403c1d9001fd83fb28e076)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: NOSHOW_PYVISTA preprocessor ([6ba5e6770](https://github.com/gha3mi/forcad/commit/6ba5e67708277bd84ecb316adbfcefca58124a24)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: warnings for temporary arrays ([9e18ca99d](https://github.com/gha3mi/forcad/commit/9e18ca99ddd1b4a83d53685fb658df43490b0ceb)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: improve lsq_bspline_* ([3c69d15d1](https://github.com/gha3mi/forcad/commit/3c69d15d12ec79902f562c619671c9afd926f1b4)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: remove unused variable ([963f6a1d6](https://github.com/gha3mi/forcad/commit/963f6a1d63d31f23d8b9889c75344902020b958e)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: avoid temporary array ([baf3644cf](https://github.com/gha3mi/forcad/commit/baf3644cf2a2f2fa615878b5c2307277608aeb87)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: cmp_volume add missing ngauss_ ([261ee14db](https://github.com/gha3mi/forcad/commit/261ee14db6a063e141109eba3c4f463796fb632b)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: remove unused variables from test program ([e66b98579](https://github.com/gha3mi/forcad/commit/e66b98579e883d1b449ae03da65b658cb84730b8)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: avoid temporary array warnings ([f82ebdb85](https://github.com/gha3mi/forcad/commit/f82ebdb8500cb07208da0e544687f27f0dbd1cfd)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: print format for L2 error norm ([4b2012f70](https://github.com/gha3mi/forcad/commit/4b2012f70a3a31209f4bcc79006583c5aaa72d7d)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: update preprocessor condition ([a2f8b7e20](https://github.com/gha3mi/forcad/commit/a2f8b7e20e12c339197d08ee786452dfe1ecbcc4)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: change output array declarations to use explicit size ([660445ad0](https://github.com/gha3mi/forcad/commit/660445ad0aebaaf08290442d8542c2cfe1a12f73)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: update preprocessor condition for gfortran ([23686d9ec](https://github.com/gha3mi/forcad/commit/23686d9ec1174b8b7c8864d401f2ba8d77943833)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: correct error message ([05ea02219](https://github.com/gha3mi/forcad/commit/05ea02219f331332aad0637373eeb250aa033704)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: add local variables to do concurrent loops ([19e4d95e1](https://github.com/gha3mi/forcad/commit/19e4d95e19447f5eeaa550250711c78f6f73c4f7)) by [@gha3mi](https://github.com/gha3mi)
* [#67] fix: improve nearest_point* methods ([2223670ca](https://github.com/gha3mi/forcad/commit/2223670cae35ba30e03985529a70808238f952bd)) by [@gha3mi](https://github.com/gha3mi)
* fix: remove unused variables ([b612525d3](https://github.com/gha3mi/forcad/commit/b612525d32f48bc00b61b9204964fbfe65dbf5a1)) by [@gha3mi](https://github.com/gha3mi)
* fix: replace minloc with manual implementation due to nvfoortran bug ([393623c88](https://github.com/gha3mi/forcad/commit/393623c88f101a634a4c3b7a1d07ee3052d8f0c5)) by [@gha3mi](https://github.com/gha3mi)

### Others

* Feature Additions, Bug Fixes, and Compiler Compatibility Improvements ([#67](https://github.com/gha3mi/forcad/pull/67)) by [@gha3mi](https://github.com/gha3mi)
* [#67] refactor: do concurrent for connectivity/coordinates ([e6dddd778](https://github.com/gha3mi/forcad/commit/e6dddd7780a58be027f917fbf990ffd91d39a260)) by [@gha3mi](https://github.com/gha3mi)
* [#67] refactor: remove test for remove_knots ([2e84096c4](https://github.com/gha3mi/forcad/commit/2e84096c45777f94874153c1d8daf71231e465cc)) by [@gha3mi](https://github.com/gha3mi)
* chore(vscode): exclude vtk, iges, ppm from fortls ([1f510d18f](https://github.com/gha3mi/forcad/commit/1f510d18fafbccfefee277361143abe815138dd0)) by [@gha3mi](https://github.com/gha3mi)
* refactor: update module usage in lsqe_* and poisson_* examples ([8dc8d51ba](https://github.com/gha3mi/forcad/commit/8dc8d51bac2d539602d887f8cc047e2af159f46c)) by [@gha3mi](https://github.com/gha3mi)
* chore: migrate ROADMAP.md to issues [skip ci] ([ca550e2c7](https://github.com/gha3mi/forcad/commit/ca550e2c7a9c2e5dfa93add3d601e1a18deab8b8)) by [@gha3mi](https://github.com/gha3mi)
* chore: clean up conda environment file [skip ci] ([18486abad](https://github.com/gha3mi/forcad/commit/18486abad0a49fae9851455f9c8831ef85a0bb54)) by [@gha3mi](https://github.com/gha3mi)
* chore: update release.sh [skip ci] ([020962ee3](https://github.com/gha3mi/forcad/commit/020962ee34709802241d6ed9018d931d355fec6c)) by [@gha3mi](https://github.com/gha3mi)


### Contributors
- [@gha3mi](https://github.com/gha3mi)



Full Changelog: [v0.13.0...v0.14.0](https://github.com/gha3mi/forcad/compare/v0.13.0...v0.14.0)

## [v0.13.0](https://github.com/gha3mi/forcad/compare/v0.12.0...v0.13.0) - 2025-08-01


### Features

* feat: optimize memory management and enhance VTK export (#48) ([440022ead](https://github.com/gha3mi/forcad/commit/440022ead99b47e0510a4bc27a0ca4708c7fb3a7)) by [@gha3mi](https://github.com/gha3mi)
* feat: add 2D and 3D Poisson IGA solver examples and update README (#49) ([55652d145](https://github.com/gha3mi/forcad/commit/55652d1454d7e541724d720f0b6bcae40ff016dd)) by [@gha3mi](https://github.com/gha3mi)

### Fixes

* fix(ci): use separate build directories for test, static and shared builds ([b6e5e4d53](https://github.com/gha3mi/forcad/commit/b6e5e4d53e0a455917d07d2ce54e74ea348c3c6f)) by [@gha3mi](https://github.com/gha3mi)
* fix: update cmake compiler flags ([e571c0af9](https://github.com/gha3mi/forcad/commit/e571c0af9a11aec556ab3dc887da5f8aa9b4cf05)) by [@gha3mi](https://github.com/gha3mi)
* fix: ensure all fpm and cmake tests run regardless of previous failures ([2a6187e29](https://github.com/gha3mi/forcad/commit/2a6187e2925d17d542c5fe34bda9a846ce8b5f55)) by [@gha3mi](https://github.com/gha3mi)

### Others

* chore: add workflow_dispatch trigger to CI/CD configuration ([50e7b81b4](https://github.com/gha3mi/forcad/commit/50e7b81b4cbfdfef20f42bf99c726c101b42ae5f)) by [@gha3mi](https://github.com/gha3mi)
* chore: fix codecov workflow to trigger on pushes to main branch ([6354d885c](https://github.com/gha3mi/forcad/commit/6354d885c556834c216a67ac5a4b80691eae159d)) by [@gha3mi](https://github.com/gha3mi)
* Update README.md status table [ci skip] (#47) ([f219070fd](https://github.com/gha3mi/forcad/commit/f219070fd036a170ea8851bb8657ce2b47203863)) by [@gha3mi](https://github.com/gha3mi)


### Contributors
- [@gha3mi](https://github.com/gha3mi)



Full Changelog: [v0.12.0...v0.13.0](https://github.com/gha3mi/forcad/compare/v0.12.0...v0.13.0)

## [v0.12.0](https://github.com/gha3mi/forcad/compare/v0.11.0...v0.12.0) - 2025-07-29


### Features

* feat: enhance CMake configuration with install and uninstall targets ([391e2d03a](https://github.com/gha3mi/forcad/commit/391e2d03ab53743ae54f9009f00135b182828a88)) by [@gha3mi](https://github.com/gha3mi)
* feat: change functions to pure elemental in forcad_utils module ([cbd80377a](https://github.com/gha3mi/forcad/commit/cbd80377a43d6f3b07d1ef73ed3ab17a498f263f)) by [@gha3mi](https://github.com/gha3mi)

### Others

* docs: update README [skip ci] ([6c94b0b3c](https://github.com/gha3mi/forcad/commit/6c94b0b3c89d5dbacec6d0eec4f2c94682ad672f)) by [@gha3mi](https://github.com/gha3mi)


### Contributors
- [@gha3mi](https://github.com/gha3mi)



Full Changelog: [v0.11.0...v0.12.0](https://github.com/gha3mi/forcad/compare/v0.11.0...v0.12.0)

## [v0.11.0](https://github.com/gha3mi/forcad/compare/v0.10.1...v0.11.0) - 2025-07-28


### Features

* feat: B-spline fitting; ensure nvfortran multicore/gpu compat (#42) ([336b1d484](https://github.com/gha3mi/forcad/commit/336b1d4843d69b21f06d84af2afcdff010422340)) by [@gha3mi](https://github.com/gha3mi)

### Others

* chore: exclude examples and logo dirs, update README badges (#41) ([3e9bd91a3](https://github.com/gha3mi/forcad/commit/3e9bd91a3bd407e6855657e8310ead2c07184deb)) by [@gha3mi](https://github.com/gha3mi)
* docs: update README to include LLVM Flang as a supported Fortran compiler ([2f0a536c6](https://github.com/gha3mi/forcad/commit/2f0a536c6d34d668323e6fcc09a132a8eec96a7b)) by [@gha3mi](https://github.com/gha3mi)
* docs: correct typos and enhance clarity in CONTRIBUTING.md ([5ffb8b4c2](https://github.com/gha3mi/forcad/commit/5ffb8b4c22139110a4cea60b5bb77f71d460bd28)) by [@gha3mi](https://github.com/gha3mi)
* Update README.md status table [ci skip] (#43) ([289a72798](https://github.com/gha3mi/forcad/commit/289a727989489e7158655043dd6d62a261378a27)) by [@gha3mi](https://github.com/gha3mi)
* chore: add homepage field to fpm.toml [skip ci] ([2f8aec3e8](https://github.com/gha3mi/forcad/commit/2f8aec3e88c757a181d34722bf6888ad2db9a983)) by [@gha3mi](https://github.com/gha3mi)


### Contributors
- [@gha3mi](https://github.com/gha3mi)



Full Changelog: [v0.10.1...v0.11.0](https://github.com/gha3mi/forcad/compare/v0.10.1...v0.11.0)

## [v0.10.1](https://github.com/gha3mi/forcad/compare/v0.10.0...v0.10.1) - 2025-07-22


### Fixes

* fix: avoid real equality/inequality comparison (gfortran warning) ([60348eb43](https://github.com/gha3mi/forcad/commit/60348eb43a8c8faed58c8025f85b1d86de90e65a)) by [@gha3mi](https://github.com/gha3mi)

### Others

* refactor: restrict module usage to only necessary components ([33967129b](https://github.com/gha3mi/forcad/commit/33967129b50aa200f6e329185f461952321c96ae)) by [@gha3mi](https://github.com/gha3mi)
* refactor: add gfortran flag -Wno-maybe-uninitialized to fpm.rsp file ([4437f005f](https://github.com/gha3mi/forcad/commit/4437f005faa6e40c3d6f27e47124e1f3f9956f94)) by [@gha3mi](https://github.com/gha3mi)


### Contributors
- [@gha3mi](https://github.com/gha3mi)



Full Changelog: [v0.10.0...v0.10.1](https://github.com/gha3mi/forcad/compare/v0.10.0...v0.10.1)

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
