title: Loft between two sections
category: Construction

Connect two circles of different radii to make a tapered surface. The order of the circles determines which end is the start of the loft.

1. Use **Create → Circle** twice. Name the first circle **Base**, set **Center** to `(0, 0, 0)` and **Radius** to `2`. Name the second **Rim**, set **Center** to `(0, 0, 3)` and **Radius** to `1`.
2. Choose **Create → From profiles → Loft** and select **Curve** as the **Input type**. Choose a circle in **Section** and click **Add**. Repeat for any circle missing from **Order**. Studio already lists the curve selected when you opened the dialog.
3. Use **Move up** or **Move down** to put **Base** first and **Rim** second. Leave **Loft degree** at `0` (Automatic) and choose **Loft**.

The new surface meets both circles. At half its height, the radius is 1.5 units. For other profiles, use sections with matching degrees, knot vectors, control counts and coordinate dimensions.

<figure class="doc-figure" markdown="1">

[![Loft dialog showing Base as section 1, Rim as section 2, and Loft degree 0 for automatic](../../media/loft-section-order.webp)](../../media/loft-section-order.webp)

<figcaption>With Base first and Rim second, the surface starts on the larger circle and ends on the smaller one.</figcaption>
</figure>

<span id="in-fortran"></span>
**Fortran equivalent.** [Loft a surface between two curves](../fortran/loft-sections.html).
