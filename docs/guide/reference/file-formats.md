title: File formats

| Format | Public Fortran API | ForCAD Studio |
| --- | --- | --- |
| FCAD (`.fcad`) · editable project | — | **Project → Open** opens another tab. **Import FCAD** adds to the active project. **Save** writes geometry, assemblies, generated meshes, and project settings. |
| IGES (`.igs`, `.iges`) · NURBS geometry | [`export_iges` on a curve or surface](../fortran/export-curve-iges.html) writes one patch. The writer uses an inch unit header. | **Project → Import IGES** reads supported B-spline curves (entity 126) and surfaces (entity 128). **Export IGES** writes the selected patch and uses the project length unit when one is set. |
| Legacy VTK (`.vtk`) · sampled geometry or mesh | [`create` and `export_Xg`](../fortran/export-volume-vtk.html) write a sampled patch in physical coordinates. | **Project → Export VTK** samples geometry. **Mesh → Export mesh (VTK)** writes the nodes and cells of a generated mesh. |
| Fortran (`.f90`) · reconstruction program | The exported program calls the public API. | **Project → Export Fortran** writes a program to rebuild supported geometry, multipatches, and meshes. |

[Choose an export in Studio](../how-to/exchange-files.html#choose-an-export) for instructions on saving each representation.
