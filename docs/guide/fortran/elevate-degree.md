title: Elevate a surface degree
category: Patch editing

Raise the degree around a rational annulus from quadratic to cubic. U runs around the ring and V runs from its inner circle to its outer circle. In `elevate_degree(dir=1, t=1)`, `dir=1` selects U and `t=1` adds one degree.

```fortran
program elevate_surface_degree
    use forcad, only: nurbs_surface, rk
    implicit none
    type(nurbs_surface) :: ring
    real(rk), parameter :: probe(2) = [1._rk / 6._rk, 0.5_rk]
    real(rk), allocatable :: before(:), after(:)

    call ring%set_ring(center=[0._rk, 0._rk, 0._rk], radius1=1._rk, radius2=2._rk)
    if (.not. ring%err%ok) error stop "Could not create the annulus"
    before = ring%cmp_Xg(probe)

    call ring%elevate_degree(dir=1, t=1)
    if (.not. ring%err%ok) error stop "Degree elevation failed"
    after = ring%cmp_Xg(probe)

    print '(a,2(1x,i0))', "Degrees U V:", ring%get_degree()
    print '(a,2(1x,i0))', "Controls U V:", ring%get_nc()
    print '(a,3f8.3)', "Point:", after
    if (maxval(abs(after - before)) > 100 * epsilon(1._rk)) error stop "Probe moved"
end program
```

The degrees change from `2 1` to `3 1` and the control net grows from `7 × 2` to `10 × 2`. The point halfway across the ring at \(u=1/6\) still evaluates to `(0.750, 1.299, 0.000)` in model coordinates.

For a curve, use `curve%elevate_degree(t=1)`. A volume takes the same arguments as the surface, with `dir=3` available for W.

**Studio equivalent.** [Elevate the degree of a patch](../how-to/elevate-degree.html) with **U** selected and **Degrees to add** set to 1.
