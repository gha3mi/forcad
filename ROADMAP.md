This roadmap outlines upcoming features and enhancements for ForCAD. Contributions and suggestions are welcome!

- v0.2.0:
    - [x] Add `insert_knots()` method for curves, surfaces and volumes.
    - [x] Add `elevate_degree()` method for curves, surfaces and volumes.
    - [x] Add `derivative()` method for curves, surfaces and volumes.

- v0.3.0:
    - [x] Add `remove_knots()` method for curves, surfaces and volumes.
    - [x] Add `put_to_nurbs()` method for volumes.
    - [x] Get IGA elements connectivity.
    - [x] Add predefined shapes: `Circle`, `Tetragon`, `Hexahedron`.
    - [x] Add `rotate_Xc` and `rotate_Xg` methods for curves, surfaces and volumes.
    - [x] Add `translate_Xc` and `translate_Xg` methods for curves, surfaces and volumes.
    - [x] Add basic unit tests.
    - [x] Add simple examples.

- v0.4.0:
    - [x] Visualization using PyVista.

- v0.4.1:
    - [x] Add new predefined shapes: `half circle`, `2D ring`, `half 2D ring`, `3D ring`, `half 3D ring`, `c-shapes`.
    - [x] Add new examples.
    - [x] Add generic procedures `get_Xc()`, `get_Xg()` and `get_Wc()`
    - [x] Add `cmp_degree()`.
    - [x] Design a logo with ForCAD.
    - [x] Improvements and bug fixes.

- Additional Ideas:
    - [ ] Add `reduce_degree()` method for curves, surfaces and volumes.
    - [ ] Add `morph()` method for morphing box.
    - [ ] Add support binary `vtk` files.
    - [ ] Export to `IGES` format.
    - [ ] Add support for multiple patches.
    - [ ] Add extraction of piecewise Bezier objects from NURBS.
    - [ ] Add more unit tests.
    - [ ] Add more examples.
    - [ ] Add more predefined shapes.
    - [ ] Improve PyVista scripts.
    - [ ] Use ForCAD for IGA simulations.
    - [ ] Use ForCAD for FEM simulations.
    - [ ] Get surfaces and edges from volumes.
    - [ ] Get edges from surfaces.
    - [ ] Add extrude method for surfaces.
    - [ ] Improve documentation.
    - [ ] Add offset method.