title: Inspect a layer of volume mesh cells
category: Meshing

1. Select a volume mesh under **Scene → Meshes**. If its eye is crossed out, click the eye to show it.
2. On the **Mesh** tab, enable **Meshes only** and choose **Cell layer**. Select **U**, **V**, or **W** and step through the layer numbers.
3. Enable **Pick mesh cells** and click a cell. The status bar identifies its cell, patch, and face.

<figure class="doc-figure" markdown="1">

[![A cell on the exposed middle ring of the White Rook hexahedral mesh is highlighted in magenta](../../media/mesh-cell-layer.webp)](../../media/mesh-cell-layer.webp)

<figcaption>White Rook at W layer 3 of 6. Each patch contributes a layer along its own W direction. The magenta cell is on the middle ring.</figcaption>
</figure>

Choose **Exterior** to restore the outside faces or **Full lattice** to see interior edges. Turn off **Meshes only** to show the geometry again.
