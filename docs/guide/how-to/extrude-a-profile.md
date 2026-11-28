title: Extrude a profile
category: Construction

Extrude a curve to make a surface, or extrude a surface to make a volume. The displacement gives both the direction and distance from the original profile.

1. Select a curve or surface under **Scene → Geometry**. To try the example shown here, create a **Rectangle** on the **Create** tab, name it **Base** and keep the default dimensions.
2. Choose **Create → Extrude**. Check **Input type** and **Profile**, then enter a name for the result.
3. Set the X, Y and Z components of **Displacement**. For the rectangle in the XY plane, enter **0**, **0**, **0.5** to make a block half a unit high.
4. Choose **Extrude**. The block appears under **Scene → Geometry → Volumes**, and the original rectangle remains under **Surfaces**.

<figure class="doc-figure" markdown="1">

[![Extrude dialog with the Base surface selected, the result named Block, and displacement X 0, Y 0, Z 0.5](../../media/extrude-profile.webp)](../../media/extrude-profile.webp)

<figcaption>Check the result type at the upper right before creating the patch.</figcaption>
</figure>

**Fortran equivalent.** [Extrude the same surface with the public API](../fortran/extrude-a-surface.html) and check the resulting volume.
