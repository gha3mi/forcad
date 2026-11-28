title: Mesh a multipatch volume
category: Meshing

Start with a volume assembly that has [joined faces](join-volume-patches.html). A conforming hex mesh gives hexahedra on either side of a joined face common node IDs.

## Set the sizing

1. Select the assembly under **Scene → Multipatches**. Choose **Multipatch → Check interfaces** and read the result in the status bar.
2. Choose **Mesh → Generate mesh…**, or **Remesh…** if the source already has a mesh.
3. Choose the sizing method below. For manual counts, both mapped directions on each joined face must match. A swap in the face mapping can make U on one member correspond to V on another.
4. Enable **Conforming hex mesh**. Choose **Generate** or **Remesh**, then select the mesh under **Scene → Meshes** to check its cell and node counts. Remeshing replaces the existing mesh for this assembly.

| Sizing method | Input |
| --- | --- |
| **Shared counts** | Enter one U, V, and W count for all members. |
| **Per patch** | Enter U, V, and W cell counts for each member independently. |
| **Target size** | Enter a positive **Target cell size** in the displayed length unit. Enable **Synchronize joined directions** and choose **Preview counts** to inspect the estimates for each member. |

<figure class="doc-figure" markdown="1">

[![Remesh dialog for White Rook with target size 1.5, synchronized joined directions, a per-patch count table, and Conforming hex mesh selected](../../media/mesh-counts.webp)](../../media/mesh-counts.webp)

<figcaption>The preview gives the rook’s first core patch 8 × 8 × 1 cells and its first sector 7 × 8 × 1. Counts can differ between members while their joined face directions still match.</figcaption>
</figure>

If generation fails, check the reported interface and the counts along its two mapped directions. [How multipatch meshes share nodes](../explanation/multipatch.html)

## Inspect the interior

[Show a layer of volume mesh cells](inspect-mesh-interior.html) to inspect the finished mesh.
