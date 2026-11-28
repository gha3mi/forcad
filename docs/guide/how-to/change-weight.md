title: Change a control weight
category: Patch editing

A positive weight changes a control point’s influence on the patch while keeping its coordinates fixed.

1. Select a patch under **Scene → Geometry**. Choose **Patch → Points + weights** (**Ctrl+E**).
2. Find the control’s index, double-click its **Weight** cell and enter a positive value. On the [three-control quadratic curve](../tutorials/first-curve.html), set the weight at U index **3** to **2**.
3. Choose **Apply**. Reopen the table to check the stored value. To check the shape, open **Inspector → Selection and probe**. At \(u=0.5\), the example curve now passes through \((1.2,0.8)\) in the XY plane.

<figure class="doc-figure">
  <a href="../../media/curve-control-weights.webp">
    <picture>
      <source media="(max-width: 32rem)" srcset="../../media/curve-control-weight-detail.webp" width="172" height="128">
      <img src="../../media/curve-control-weights.webp" width="604" height="128" alt="The three weights are 1, 1 and 2. The third control’s weight cell is selected.">
    </picture>
  </a>
  <figcaption>The first two weights stay at 1. Increasing the third pulls the curve toward its endpoint.</figcaption>
</figure>

**Fortran equivalent.** [Change the same curve weight](../fortran/change-curve-weight.html). For a \(G^1\)-connected assembly, [edit the shared weight](edit-shared-control.html) to preview the affected controls.
