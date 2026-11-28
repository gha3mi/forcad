title: Subdivide the knot spans of a patch
category: Patch editing

Use **Refine patch** to divide every span along a direction by the same amount.

1. Select one patch and open **Patch → Refine patch**.
2. Select a direction, or **All** to apply the same settings along all available directions.
3. Enter **Divisions per span**. For example, enter **2** to split each span in two. Leave **Degree increase** at **0** to keep the current degree.
4. Choose **Refine patch**. Open **Patch → Parameter panels** and inspect **Knots** to see the new divisions. **Undo** restores the previous representation.

<figure class="doc-figure">
  <a href="../../media/refine-patch.webp">
    <picture>
      <source media="(max-width: 32rem)" srcset="../../media/refine-patch-detail.webp" width="273" height="75">
      <img src="../../media/refine-patch.webp" width="549" height="177" alt="Refine patch with U selected, two divisions per span and a degree increase of one.">
    </picture>
  </a>
  <figcaption>These settings raise the U degree by one, then halve each U span. V keeps its original spans.</figcaption>
</figure>

For a single new knot at a specific parameter, see [Insert a knot](insert-knot.html).

**Fortran equivalent.** [Subdivide the knot spans of a patch](../fortran/refine-a-patch.html).
