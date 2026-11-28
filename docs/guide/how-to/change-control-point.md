title: Move a control point
category: Patch editing

Use the coordinate table when you know the exact position a control point should occupy.

1. Select one patch under **Scene → Geometry** and open **Patch → Points + weights** (**Ctrl+E**).
2. Locate the point by its one-based index. A curve shows U, a surface shows U and V, and a volume shows U, V, and W.
3. Double-click the **X**, **Y**, or **Z** cell you want to change. Enter the new coordinate in the unit shown in the column heading.
4. Choose **Apply** and inspect the changed shape. **Undo** restores the previous coordinates.

<figure class="doc-figure" markdown="1">

[![Eight controls of a unit block. The Z coordinate of control 2, 2, 2 is selected and reads 1.25.](../../media/control-point-editor.webp)](../../media/control-point-editor.webp)

<figcaption>The last row shows the corner after its Z coordinate changed from 1 to 1.25. Its weight is still 1.</figcaption>
</figure>

**Fortran equivalent.** [Move the same corner by its control index](../fortran/move-volume-control.html). For a point shared across joined patches, use the [multipatch control editor](edit-shared-control.html).
