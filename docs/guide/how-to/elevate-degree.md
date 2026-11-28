title: Elevate the degree of a patch
category: Patch editing

Raise the polynomial degree in one direction while preserving the patch geometry.

1. Select a single patch and choose **Patch → Elevate degree**.
2. Choose the direction. Enter **Degrees to add** and compare **Current degree** with **Requested degree**. To change a quadratic direction to cubic, add one.
3. Choose **Elevate**. Read the result in the main status bar. Use **Undo** to restore the previous representation.

<figure class="doc-figure">
  <a href="../../media/elevate-degree.webp">
    <picture>
      <source media="(max-width: 32rem)" srcset="../../media/elevate-degree-detail.webp" width="352" height="58">
      <img src="../../media/elevate-degree.webp" width="472" height="290" alt="Elevate degree on U, with Degrees to add set to 1 and a preview from degree 2 to 3.">
    </picture>
  </a>
  <figcaption>The preview shows the degree along U. V keeps its current degree.</figcaption>
</figure>

To elevate across connected directions, select the assembly and use **Multipatch → Elevate degree**.

**Fortran equivalent.** [Raise the U degree of a rational surface](../fortran/elevate-degree.html).
