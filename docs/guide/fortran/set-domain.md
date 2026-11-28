title: Change a parameter interval
category: Patch editing

Change the U interval of a unit cube from \([0,1]\) to \([2,6]\). The midpoint of that interval moves from \(u=0.5\) to \(u=4\). Evaluate the original and returned volumes at corresponding parameters to check that they give the same physical point.

```fortran
program change_parameter_interval
    use forcad, only: rk, nurbs_volume
    implicit none
    type(nurbs_volume) :: cube, rescaled
    real(rk), allocatable :: before(:), after(:)

    call cube%set_hexahedron(L=[1._rk, 1._rk, 1._rk], nc=[2, 2, 2])
    if (.not. cube%err%ok) error stop "Could not create the cube"
    rescaled = cube%reparameterize(dir=1, a=2._rk, b=6._rk)
    if (.not. rescaled%err%ok) error stop "Could not set the interval"

    before = cube%cmp_Xg([0.5_rk, 0.25_rk, 0.75_rk])
    after = rescaled%cmp_Xg([4._rk, 0.25_rk, 0.75_rk])
    print '(a,3f8.3)', "Original: ", before
    print '(a,3f8.3)', "Rescaled: ", after
    if (maxval(abs(after - before)) > 100 * epsilon(1._rk)) error stop "Probe moved"
end program
```

Both rows print `0.500 0.250 0.750`. The function returns a new patch. Pass finite bounds with `a < b`. A curve takes just the bounds, as in `curve%reparameterize(a=2._rk, b=6._rk)`. On a surface or volume, `dir` selects the parameter axis.

**Studio equivalent.** [Set a parameter domain](../how-to/set-domain.html) shows where to enter the bounds for one direction.
