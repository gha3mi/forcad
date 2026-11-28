title: Patch editing commands
category: Patch editing

Select a curve, surface, or volume in **Scene**, then open **Patch**.

<figure class="doc-figure" markdown="1">

<a href="../../media/patch-ribbon.webp">
  <picture>
    <source media="(max-width: 32rem)" srcset="../../media/patch-ribbon-detail.webp" width="363" height="108">
    <img src="../../media/patch-ribbon.webp" width="643" height="108" alt="Control and Parameters groups on the Patch ribbon, with Points + weights, Set domain, and Parameter panels">
  </picture>
</a>

<figcaption>Choose <strong>Points + weights</strong> to edit the control lattice or <strong>Parameter panels</strong> to inspect the knots and basis functions.</figcaption>
</figure>

## Set a control point or weight

In **Points + weights**, [move a control point](change-control-point.html) or [change its weight](change-weight.html).

## Edit a shared control

Select the assembly under **Scene → Multipatches**, then [edit its shared control](edit-shared-control.html) to preview dependent changes at \(G^1\) joins.

## Refine or reparameterize

- **Knots** — [Insert a knot](insert-knot.html) or [remove a knot](remove-knot.html).
- **Degree & spans** — [Elevate the degree](elevate-degree.html), [reduce the degree](reduce-degree.html), or [subdivide the spans](refine-patch.html).
- **Parameters** — [Set the domain](set-domain.html) to change the bounds of one parameter direction.
- **Topology** — [Split the patch](split-a-patch.html) at an interior parameter value.
