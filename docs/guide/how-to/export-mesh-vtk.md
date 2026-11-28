title: Export a generated mesh as VTK
category: Projects and exchange

Use this export to preserve the nodes and cell connectivity of a mesh you have already generated.

1. Select the mesh under **Scene → Meshes**. If its row shows **Needs update**, remesh its source first.
2. Choose **Mesh → Export mesh (VTK)**, or right-click the mesh row and choose **Export mesh (VTK)**. Enter a `.vtk` filename and save.
3. Open the file in a VTK viewer. Compare its point and cell counts with the counts in the Studio status bar for the selected mesh.

<figure class="doc-figure">
  <a href="../../media/vtk-mesh-command.webp"><img src="../../media/vtk-mesh-command.webp" width="264" height="301" alt="The context menu for Mesh 1 · Block has Export mesh (VTK) highlighted."></a>
  <figcaption>Each of the example block’s 24 hexahedra refers to eight of its 60 exported points.</figcaption>
</figure>

[Generate the block mesh](mesh-a-volume.html) or [export the sampled shape of its source geometry](export-geometry-vtk.html).
