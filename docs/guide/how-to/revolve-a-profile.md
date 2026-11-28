title: Revolve a profile
category: Construction

Make a quarter-cylinder wall by revolving a vertical line about the Z axis. The line is 2 units from the axis and 3 units long.

1. Choose **Create → Line**. Name it **Upright**. Set **Start** to `(2, 0, 0)` and **End X**, **End Y**, **End Z** to `2`, `0`, `3`. Choose **Create**.
2. Choose **Create → From profiles → Revolve**. Set **Input type** to **Curve** and choose **Upright** under **Profile**.
3. Set **Axis point** to `(0, 0, 0)`, **Axis direction** to `(0, 0, 1)`, and **Angle (°)** to `90`. Choose **Revolve**.

The new surface spans a quarter turn from +X toward +Y. Its upper edge ends at `(0, 2, 3)`. The original line remains in **Scene**.

<figure class="doc-figure" markdown="1">

[![The selected quarter-cylinder wall with two straight vertical sides and curved upper and lower edges](../../media/revolve-quarter-wall.webp)](../../media/revolve-quarter-wall.webp)

<figcaption>The straight sides are the initial and final positions of the line. The curved edges trace the paths of its endpoints.</figcaption>
</figure>

<span id="in-fortran"></span>
**Fortran equivalent.** [Revolve a curve into a surface](../fortran/revolve-a-profile.html).
