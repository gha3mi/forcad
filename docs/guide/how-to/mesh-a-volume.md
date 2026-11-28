title: Mesh a single volume
category: Meshing

Generate a structured grid of 24 hexahedra on one volume. To reproduce the example, create a unit cube with **Create → Box**.

1. Select the volume under **Scene → Geometry → Volumes** and choose **Mesh → Generate mesh…**.
2. Choose **Cell counts**. Set **U cells** to **4**, **V cells** to **3** and **W cells** to **2**, then choose **Generate**.
3. Select the new row under **Scene → Meshes**. The status bar reports **24 cells** and **60 nodes**. To change the counts later, choose **Mesh → Remesh…**. Remeshing replaces that mesh.

<div class="doc-gallery">
  <figure class="doc-figure">
    <a href="../../media/volume-mesh-counts.webp">
      <picture>
        <source media="(max-width: 32rem)" srcset="../../media/volume-mesh-counts-detail.webp" width="342" height="65">
        <img src="../../media/volume-mesh-counts.webp" width="448" height="244" alt="The U, V and W cell counts are set to 4, 3 and 2.">
      </picture>
    </a>
    <figcaption>Enter counts along the patch’s parametric directions.</figcaption>
  </figure>
  <figure class="doc-figure">
    <a href="../../media/volume-mesh-grid.webp"><img src="../../media/volume-mesh-grid.webp" width="383" height="429" alt="Exterior of the generated block mesh. Its top face shows a four-by-three grid, and the side faces show two layers."></a>
    <figcaption>The top surface has 12 quadrilateral faces. Two layers fill the volume.</figcaption>
  </figure>
</div>

[Export these nodes and cells as VTK](export-mesh-vtk.html). To export a fresh sampling of the geometry instead, use [Export sampled geometry as VTK](export-geometry-vtk.html).

**Fortran equivalent.** [Sample a volume into hexahedra](../fortran/sample-volume-mesh.html) uses \(5\times4\times3\) sample points for the same cell counts.
