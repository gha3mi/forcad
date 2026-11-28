title: Export sampled geometry as VTK
category: Projects and exchange

Choose this export to view a NURBS patch or assembly in a VTK viewer. ForCAD Studio samples its geometry at export time.

1. Select a patch under **Scene → Geometry**, or an assembly under **Scene → Multipatches**.
2. Choose **Project → Export VTK**. Enter a filename ending in `.vtk` and save.
3. Open the file in a VTK viewer. The sampling resolution is independent of any mesh generated for that object.

<figure class="doc-figure">
  <a href="../../media/vtk-geometry-selection.webp"><img src="../../media/vtk-geometry-selection.webp" width="264" height="194" alt="In the Scene tree, Block is selected under Volumes. Mesh 1 · Block is listed separately under Meshes."></a>
  <figcaption>Select Block to export the sampled volume. To keep the 24 cells of its generated mesh, select Mesh 1 · Block and use the mesh export.</figcaption>
</figure>

[Export a generated mesh](export-mesh-vtk.html) to retain its nodes and cells. **Fortran equivalent.** [Sample a volume and write it to VTK](../fortran/export-volume-vtk.html) with explicit resolution in each direction.
