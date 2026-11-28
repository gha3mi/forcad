title: Split a patch
category: Patch editing

Divide a curve, surface, or volume into two independently editable patches at one parameter value.

1. Select one patch under **Scene → Geometry**. If it belongs to a multipatch, Studio requires you to detach it from the assembly before splitting it.
2. Choose **Patch → Split patch**. Select a parameter direction and enter an interior value in **Split at**. For example, for U in `[0, 1]`, enter `0.5` to divide the interval in half.
3. Choose **Split**. In **Scene**, the lower piece keeps the original name and the upper piece is a separate object. **Undo** restores the original patch.

<figure class="doc-figure" markdown="1">

[![Split patch dialog for Block. U is selected, its current range is 0 to 1, and Split at is 0.5.](../../media/split-volume.webp)](../../media/split-volume.webp)

<figcaption>For this unit cube, the two volumes retain U ranges [0, 0.5] and [0.5, 1].</figcaption>
</figure>

**Fortran equivalent.** [Split the unit cube](../fortran/split-a-volume.html) and check both volumes at their shared face. To add an interior span boundary while keeping one patch, [insert a knot](insert-knot.html).
