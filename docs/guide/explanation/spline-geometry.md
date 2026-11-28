title: NURBS geometry and parameter space

## Two coordinate systems

A patch maps parameters to a physical point. The **Parameter Space** panel shows one axis for a curve, a rectangle for a surface, and a box for a volume. Studio labels the axes **U**, **V**, and **W**. In formulas the corresponding coordinates are \(u\), \(v\), and \(w\). A curve in 3D still needs only \(u\) to locate a point on it.

The Inspector uses a normalized probe value \(\hat u\) between zero and one. If the active knot interval is \([a,b]\), the parameter passed to the patch is

\[
u=a+(b-a)\hat u.
\]

For example, a probe value of 0.5 on the interval from 2 to 6 evaluates the patch at 4. This is the midpoint of the parameter interval. The physical distance traveled along a curve also depends on its shape and parameterization.

## What a control point changes

A non-uniform rational B-spline (NURBS) curve blends its control points \(\mathbf{P}_i\) using the degree-\(p\) B-spline basis \(N_{i,p}\) and positive weights \(w_i\). With \(n_c\) control points, indexed from one as in Fortran,

\[
\mathbf{C}(u)=
\frac{\sum_{i=1}^{n_c}N_{i,p}(u)\,w_i\,\mathbf{P}_i}
     {\sum_{i=1}^{n_c}N_{i,p}(u)\,w_i}.
\]

The coefficient of each point in this sum is its rational basis value, \(R_i(u)=N_{i,p}(u)w_i / \sum_j N_{j,p}(u)w_j\). The Fortran binding `curve%basis` returns these coefficients in control-point order. [Evaluate them on a quarter circle](../fortran/evaluate-rational-basis.html).

Increasing one weight draws the curve toward its control point wherever that basis function contributes. Multiplying all weights by the same positive number leaves the curve unchanged. [Change a weight in Fortran](../fortran/change-curve-weight.html) shows the effect at the midpoint of a quadratic curve.

## Knots and smoothness

Distinct knots divide a parameter direction into spans. At an interior knot repeated \(s\) times, a degree-\(p\) spline has at least \(C^{p-s}\) continuity for \(1\le s\le p\). For a cubic curve, a simple knot guarantees continuity of the first two derivatives. A knot of multiplicity two guarantees continuity of the first derivative. Inserting a knot preserves the current geometry. The extra control points give later edits more freedom.

In Studio, **Basis Functions** follows the direction selected in **Parameter Space**. For a surface or volume, the other parameters stay at the probe position as the plot samples that direction. [Inspect the parameter cube and basis plot](../tutorials/explore-studio.html#inspect-a-volume-patch).

For the mathematical definitions, see [*The NURBS Book*](https://doi.org/10.1007/978-3-642-97385-7).
