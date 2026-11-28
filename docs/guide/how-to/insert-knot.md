title: Insert a knot
category: Patch editing

Insert a knot to refine the representation of a single patch while keeping its shape.

1. Select a patch and choose **Patch → Insert knot**.
2. Choose **U**, **V**, or **W** from the directions available for the patch. Enter a **Knot value** strictly inside the interval displayed for that direction. Use the patch's parameter units.
3. Set **Copies** to the number of occurrences to add. For a new knot, start with one. The dialog shows the current multiplicity and the allowed limit.
4. Choose **Insert**. The main status bar reports how many copies were inserted and along which direction. To see the updated vector, open **Patch → Parameter panels**, then choose **Knots** in **Parameter Space**.

<figure class="doc-figure" markdown="1">

[![Insert knot dialog with the parameter direction, Knot value, and Copies fields](../../media/insert-knot.webp)](../../media/insert-knot.webp)

<figcaption>At 0.5 on the U axis, the dialog shows zero existing copies and a limit of three. Adding two copies gives this knot multiplicity two.</figcaption>
</figure>

To propagate refinement across joined patches, select their assembly and use **Multipatch → Insert knot**.

**Fortran equivalent.** [Insert a knot in a curve and check its midpoint](../fortran/insert-a-knot.html).
