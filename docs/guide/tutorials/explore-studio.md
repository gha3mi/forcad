title: Explore a project in Studio

Open the bundled Chess project and follow one point from a volume's parameter
cube to its position on a chess piece. This exercise leaves the model intact.

## Open Chess

1. [Install and launch Studio](../how-to/install.html#launch-studio).
2. On Home, select **Chess** under **Samples**. The board and pieces appear in
   a project tab. The **Scene** tree lists the geometry.
3. In the viewport, drag with the right mouse button to orbit, drag with the
   middle mouse button to pan, and scroll to zoom.

<figure class="doc-figure" markdown="1">

![ForCAD Studio Home screen with New project, Open project, Chess, and Head](../../media/studio-home.webp)

<figcaption>Home offers a blank project, a file picker, and the Chess and Head samples.</figcaption>
</figure>

## Inspect a volume patch

1. In **Scene**, find **Geometry → Volumes**. If the branch is collapsed,
   click its arrow. Select one named patch inside it and press **Shift+F**
   to frame that patch.
2. Select **Patch → Parameter panels**. You should now see **Parameter Space
   (U–V–W)** and **Basis Functions**.
3. In **Parameter Space**, select **Probe** and drag inside the cube. Open
   **Inspector** to see the normalized U, V, W inputs and the evaluated
   position in physical space. To place the probe exactly halfway along each
   direction, set all three inputs to `0.5`.
4. Select **W** at the top of **Parameter Space**. The plot title changes to
   **Basis Functions (W)**. Switch the plot to **First derivative**, then back
   to **Basis functions**, to compare the two at the same probe position.
5. Press **F** in the 3D viewport to frame the whole scene again.

<figure class="doc-figure" markdown="1">

![Selected Chess volume patch with its parameter cube, basis plot, and Inspector](../../media/studio-parameters.webp)

<figcaption>The parameter cube shows U, V, and W. The 3D viewport shows where the selected patch sits in the model.</figcaption>
</figure>

Next, [edit a patch](../how-to/edit-a-patch.html) or use the
[Studio controls reference](../reference/studio.html).
