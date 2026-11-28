title: Loft a surface between two curves
category: Construction

Construct a tapered surface between two circles. The radius changes from 2 to 1 over a height of 3.

Pass the sections to [`loft`](../../interface/loft.html) in their geometric order. With two sections the default degree along the loft is one.

```fortran
program tapered_wall
    use forcad, only: rk, nurbs_curve, nurbs_surface, loft
    implicit none
    type(nurbs_curve) :: sections(2)
    type(nurbs_surface) :: wall

    call sections(1)%set_circle(center=[0.0_rk, 0.0_rk, 0.0_rk], radius=2.0_rk)
    call sections(2)%set_circle(center=[0.0_rk, 0.0_rk, 3.0_rk], radius=1.0_rk)
    if (.not. all(sections%err%ok)) error stop "Could not create the sections"
    wall = loft(sections)
    if (.not. wall%err%ok) error stop 'Loft failed'
    print '(3f8.3)', wall%cmp_Xg(Xt=[0.0_rk, 0.5_rk])
end program
```

The program prints `(1.500, 0.000, 1.500)`, a point on the middle circle.

**Studio equivalent.** [Loft between two sections](../how-to/loft-sections.html).
