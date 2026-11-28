title: Remove a knot from a curve
category: Patch editing

Use `remove_knots` to lower the multiplicity of an interior knot. `r` gives the number of copies to try. The result in `removed` tells you how many ForCAD accepted.

This quadratic curve already needs one knot at \(u=0.5\). Insert an extra copy there, then request removal of both copies. Evaluate the curve on each side of the knot and at the knot to check the edit.

```fortran
program remove_a_knot
    use forcad, only: rk, nurbs_curve
    implicit none
    type(nurbs_curve) :: curve
    real(rk), parameter :: knots(7) = [0._rk, 0._rk, 0._rk, 0.5_rk, 1._rk, 1._rk, 1._rk]
    real(rk), parameter :: probes(3) = [0.25_rk, 0.5_rk, 0.75_rk]
    real(rk) :: controls(4, 2)
    real(rk), allocatable :: before(:,:), after(:,:)
    integer, allocatable :: multiplicity(:), removed(:)

    controls(:,1) = [0._rk, 1._rk, 2._rk, 3._rk]
    controls(:,2) = [0._rk, 0._rk, 1._rk, 0._rk]
    call curve%set(knot=knots, Xc=controls, degree=2)
    call curve%insert_knots(Xth=[0.5_rk], r=[1])
    if (.not. curve%err%ok) error stop "Could not prepare the curve"

    before = curve%cmp_Xg(probes)
    multiplicity = curve%get_multiplicity()
    print '(a,i0)', "Multiplicity before: ", multiplicity(2)
    call curve%remove_knots(Xth=[0.5_rk], r=[2], removed=removed)
    if (.not. curve%err%ok) error stop "Invalid knot removal request"

    after = curve%cmp_Xg(probes)
    multiplicity = curve%get_multiplicity()
    print '(a,i0)', "Copies removed: ", removed(1)
    print '(a,i0)', "Multiplicity after: ", multiplicity(2)
    print '(a,es10.2)', "Largest sampled coordinate change: ", maxval(abs(after - before))
end program
```

The multiplicity falls from **2 to 1**. The largest coordinate change at the three probes is zero. One copy was removable with the default roundoff test. The remaining knot is needed to represent this curve. A valid request can return fewer copies than requested with `curve%err%ok` still true. Requesting more copies than the stored multiplicity is an invalid request.

Use an exact interior value from `get_knot()` when working with imported knots. For approximate removal, the optional `tol` sets a physical error budget in model units and `error_bound` reports a bound on the accepted change.

**Studio equivalent.** [Remove a knot](../how-to/remove-knot.html) shows the **Copies** field and where to read the accepted count.
