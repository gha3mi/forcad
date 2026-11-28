title: Evaluate the rational basis of a curve
category: Evaluation

Find the contribution of each control point at the midpoint of a unit quarter circle. The three controls are \((1,0)\), \((1,1)\), and \((0,1)\). Giving the middle control a weight of \(1/\sqrt{2}\) makes this quadratic Bézier curve an exact circular arc.

```fortran
program quarter_circle_basis
    use forcad, only: rk, nurbs_curve
    implicit none
    type(nurbs_curve) :: arc
    real(rk) :: control(3, 2)
    real(rk), allocatable :: basis(:), point(:)

    control(1, :) = [1.0_rk, 0.0_rk]
    control(2, :) = [1.0_rk, 1.0_rk]
    control(3, :) = [0.0_rk, 1.0_rk]
    call arc%set(Xc=control, Wc=[1.0_rk, sqrt(0.5_rk), 1.0_rk])
    if (.not. arc%err%ok) error stop "Could not construct the arc"

    call arc%basis(Xt=0.5_rk, Tgc=basis)
    if (.not. arc%err%ok) error stop "Could not evaluate the basis"
    point = arc%cmp_Xg(0.5_rk)
    print '(3f9.6)', basis
    print '(2f9.6)', point
end program quarter_circle_basis
```

The first row prints `0.292893 0.414214 0.292893`. These are the rational basis values \(R_i(u)\), in control-point order. They sum to one. The second row prints the point `0.707107 0.707107`, at 45 degrees on the circle.

<figure class="doc-figure" markdown="1">

<a href="../../media/quarter-circle-basis.webp">
  <picture>
    <source media="(max-width: 40rem)" srcset="../../media/quarter-circle-basis-compact.webp" width="600" height="1120">
    <img src="../../media/quarter-circle-basis.webp" width="1200" height="560" alt="A unit quarter circle and its control polygon, with the three rational basis functions. At u = 0.5 the two outside functions agree and the middle function is largest.">
  </picture>
</a>

<figcaption>The dots mark the same parameter on the arc and on all three basis functions. The dashed segments join the control points.</figcaption>
</figure>

The [curve basis reference](../../type/nurbs_curve.html) gives the scalar and vector forms. In Studio, [open the parameter and basis panels](../tutorials/explore-studio.html#inspect-a-volume-patch) to follow a probe on a selected patch.
