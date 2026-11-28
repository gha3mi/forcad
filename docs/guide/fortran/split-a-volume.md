title: Split a volume at a parameter
category: Patch editing

Split a unit cube halfway along U. The two output patches retain their original parameter intervals. The upper patch starts at \(u=0.5\), so use that value to evaluate either side of the cut.

```fortran
program split_a_volume
    use forcad, only: rk, nurbs_volume
    implicit none
    type(nurbs_volume) :: block, lower, upper
    real(rk), parameter :: seam(3) = [0.5_rk, 0.5_rk, 0.5_rk]
    real(rk) :: lower_volume, upper_volume

    call block%set_hexahedron(L=[1.0_rk, 1.0_rk, 1.0_rk], nc=[2, 2, 2])
    call block%split(dir=1, Xt=seam(1), left=lower, right=upper)
    if (.not. lower%err%ok .or. .not. upper%err%ok) error stop "Split failed"

    print '(a,3f8.3)', "Lower face ", lower%cmp_Xg(seam)
    print '(a,3f8.3)', "Upper face ", upper%cmp_Xg(seam)
    call lower%cmp_volume(lower_volume)
    call upper%cmp_volume(upper_volume)
    if (.not. lower%err%ok .or. .not. upper%err%ok) error stop "Measurement failed"
    print '(a,2f8.3)', "Volumes    ", lower_volume, upper_volume
end program
```

Both faces evaluate to `(0.500, 0.500, 0.500)` at the middle of the cut. The pieces each occupy `0.500` cubic model units. For a different cut, pass a strictly interior parameter on the chosen axis (`dir=1` for U, `2` for V, `3` for W).

**Studio equivalent.** [Split a patch](../how-to/split-a-patch.html) shows how to choose the same direction and parameter in the dialog.
