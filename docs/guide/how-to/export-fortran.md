title: Export a project as Fortran
category: Projects and exchange

Export the active project as a standalone `.f90` program that reconstructs its supported patches, multipatch connections, and generated meshes.

1. Activate the project tab to export. If it contains a multipatch with an invalid interface, resolve the issue before exporting.
2. Choose **Project → Export Fortran** and save the file.
3. Open the source file and locate `program rebuild_forcad_project`. The program builds the patches first, then their assemblies and meshes. It uses the ForCAD library and can run independently of Studio.

To compile the generated program, use it as `app/main.f90` in an fpm project that depends on ForCAD. The [first Fortran tutorial](../tutorials/first-curve.html#add-the-dependency) shows the manifest.
