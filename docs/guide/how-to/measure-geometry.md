title: Measure selected geometry
category: Inspect and measure

1. Select a curve, surface or volume under **Scene → Geometry**. To measure just an edge or face, select that feature in the viewport.
2. Choose **Model → Measure**. The result is length for a curve or edge, area for a surface or face, and volume for a solid patch.
3. To read the result in a particular unit, close the result, choose **Project → Units**, set **Display unit** and measure again. Choose an explicit display unit to see its symbol beside the value.

<figure class="doc-figure">
  <a href="../../media/measure-plate.webp"><img src="../../media/measure-plate.webp" width="336" height="144" alt="Measurement dialog for Plate showing an area of 6 mm²."></a>
  <figcaption>This rectangular plate measures 2 mm by 3 mm. Switching the display unit to centimeters gives an area of 0.06 cm².</figcaption>
</figure>

**Fortran equivalent.** [Measure a patch](../fortran/measure-a-patch.html) gives the length, area and volume calls with runnable examples.
