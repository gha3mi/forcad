title: Remove a knot
category: Patch editing

Select a single patch that has an interior knot in the direction you want to simplify.

1. Choose **Patch → Remove knot** and select the direction.
2. In **Knot value**, choose an existing interior value. Set **Copies** to the number of occurrences to remove. The dialog shows the current multiplicity.
3. Choose **Remove**. Check the status bar for the accepted count. ForCAD tests whether the geometry can be preserved and may accept fewer copies than requested.

<figure class="doc-figure">
  <a href="../../media/remove-knot.webp"><img src="../../media/remove-knot.webp" width="472" height="274" alt="Remove knot on U. Knot value is 0.25, Copies is 2, and Current multiplicity is 3."></a>
  <figcaption>Here the knot occurs three times. Asking for two copies to be removed leaves at least one occurrence.</figcaption>
</figure>

If the direction has no interior knots, the list shows **No interior knots**. To inspect a different direction, change the direction in the dialog.

**Fortran equivalent.** [Remove a knot from a curve and check how many copies were removed](../fortran/remove-a-knot.html).
