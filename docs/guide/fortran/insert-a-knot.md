title: Insert a knot in a curve
category: Patch editing

Insert one knot halfway along a quadratic arch. Check the control count and compare the physical midpoint before and after the edit.

```fortran
program insert_a_knot
    use forcad, only: rk, nurbs_curve
    implicit none
    type(nurbs_curve) :: curve
    real(rk) :: controls(3, 2)
    real(rk), allocatable :: before(:), after(:)

    controls(:, 1) = [0._rk, 1._rk, 2._rk]
    controls(:, 2) = [0._rk, 2._rk, 0._rk]
    call curve%set(Xc=controls)
    if (.not. curve%err%ok) error stop "Could not construct the curve"
    before = curve%cmp_Xg(0.5_rk)

    call curve%insert_knots(Xth=[0.5_rk], r=[1])
    if (.not. curve%err%ok) error stop "Could not insert the knot"
    after = curve%cmp_Xg(0.5_rk)
    print '(a,i0)', "Control points: ", curve%get_nc()
    print '(a,2f8.3)', "Midpoint: ", after
    if (maxval(abs(after - before)) > 100 * epsilon(1._rk)) error stop "Midpoint moved"
end program
```

The control count increases from three to four. The midpoint still prints as `1.000 1.000`. The array `Xth` supplies the knot values, and `r` specifies the number of copies to insert at each value. To undo an insertion programmatically, [request knot removal and check the accepted count](remove-a-knot.html).

For a surface or volume, the binding also takes `dir` to select the parameter direction. See the [surface](../../type/nurbs_surface.html) and [volume](../../type/nurbs_volume.html) declarations.

**Studio equivalent.** [Insert a knot](../how-to/insert-knot.html) in **Patch → Insert knot**.
