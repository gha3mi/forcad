title: Find the active control points on a knot span
category: Assembly and analysis

For isogeometric analysis, `cmp_elem()` returns one row per positive-length knot span. The entries are one-based indices into the patch's control points. A quadratic curve has three local basis functions per span.

This straight curve has four controls and two spans. The program pairs each span's control indices with their basis values at its midpoint. A one-point Gauss rule places the evaluation at the middle of each span.

```fortran
program active_controls
    use forcad, only: rk, nurbs_curve
    implicit none
    type(nurbs_curve) :: curve
    real(rk), parameter :: knots(7) = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
    real(rk) :: control(4, 1), weighted_length
    real(rk), allocatable :: basis(:), gradient(:,:)
    integer, allocatable :: spans(:,:), segments(:,:)
    integer :: i

    control(:,1) = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk]
    call curve%set(knot=knots, Xc=control, degree=2)
    if (.not. curve%err%ok) error stop "Could not create the curve"
    spans = curve%cmp_elem()
    do i = 1, size(spans, 1)
        call curve%ansatz(ie=i, ig=1, ngauss=1, Tgc=basis, dTgc_dXg=gradient, dL=weighted_length)
        if (.not. curve%err%ok) error stop "Could not evaluate the basis"
        print '(3i3,3f8.3)', spans(i,:), basis
    end do

    call curve%create(res=5)
    if (.not. curve%err%ok) error stop "Could not sample the curve"
    segments = curve%get_elem_Xg_vis()
    print '(2i3)', size(segments, 1), size(segments, 2)
end program
```

| Span | Midpoint | Control IDs | Basis values in the same order |
| --- | ---: | --- | --- |
| \([0,0.5]\) | 0.25 | 1, 2, 3 | 0.250, 0.625, 0.125 |
| \([0.5,1]\) | 0.75 | 2, 3, 4 | 0.125, 0.625, 0.250 |

Controls 2 and 3 contribute on both spans. The last line prints `4 2` because five sampled positions form four line segments with two endpoints each. Those endpoint indices refer to rows of `get_Xg()`. Changing `res` changes the sampled segments and leaves the two analysis elements and their active control IDs intact.
