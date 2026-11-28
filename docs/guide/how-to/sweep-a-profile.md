title: Sweep a profile along a curve
category: Construction

Sweep a tilted line around a circular path to make the side of an oblique cylinder. The line keeps its orientation as it travels.

1. Choose **Create → Circle**. Name it **Rail**, set **Center** to `(0, 0, 0)` and **Radius** to `2`, then choose **Create**.
2. Choose **Create → Line**. Name it **Edge**, set **Start** to `(2, 0, 0)` and the three end coordinates to `(3, 0, 1)`, then choose **Create**.
3. Choose **Create → From profiles → Sweep**. Select **Curve** as the **Input type**, **Edge** as the **Profile**, and **Rail** as the **Path**. Choose **Sweep**.

The upper end stays one unit farther along X and one unit higher along Z than the lower end throughout the turn. A surface profile swept along a curve produces a volume.

<figure class="doc-figure" markdown="1">

[![Oblique cylindrical surface with its lower circular rail and tilted line profile highlighted](../../media/sweep-circular-rail.webp)](../../media/sweep-circular-rail.webp)

<figcaption>The highlighted circle is the rail. The short highlighted line shows the profile joining the two offset circular edges.</figcaption>
</figure>

<span id="in-fortran"></span>
**Fortran equivalent.** [Sweep a curve along a path](../fortran/sweep-a-profile.html).
