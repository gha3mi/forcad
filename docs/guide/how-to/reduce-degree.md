title: Reduce the degree of a patch
category: Patch editing

Use degree reduction to simplify a patch within a permitted physical error.

1. Select a single patch and choose **Patch → Reduce degree**.
2. Select the direction and enter **Degrees to remove**. Check **Requested degree** against **Current degree**.
3. Set **Error tolerance** in model units. This bounds the total permitted shape change.
4. Choose **Reduce**. Check the accepted reduction in the main status bar. If the tolerance admits fewer steps, review the result before attempting a further reduction.

<figure class="doc-figure">
  <a href="../../media/reduce-degree.webp">
    <picture>
      <source media="(max-width: 32rem)" srcset="../../media/reduce-degree-detail.webp" width="378" height="72">
      <img src="../../media/reduce-degree.webp" width="472" height="325" alt="Reduce degree on U. The preview shows 3 to 2 and the error tolerance is 1e-10 model units.">
    </picture>
  </a>
  <figcaption>The preview shows the requested degree. Check the status bar after applying to see how many steps the tolerance allowed.</figcaption>
</figure>

**Fortran equivalent.** [Reduce a curve degree and check the accepted steps and error bound](../fortran/reduce-degree.html).
