title: Assemble the length sensitivity of a shared control point
category: Assembly and analysis

Join two straight curve patches at \((1,1)\), with outer endpoints at \((0,0)\) and \((2,0)\). Moving their common control point changes both segment lengths. Use `assemble_sensitivities` to add the two patch contributions to the gradient of the total length.

```fortran
program shared_control_sensitivity
    use forcad, only: rk, nurbs_curve, nurbs_multipatch_curve
    implicit none
    type(nurbs_curve) :: segment(2)
    type(nurbs_multipatch_curve) :: chain
    real(rk) :: xy(3, 2), controls(2, 2), local(4, 2)
    real(rk), allocatable :: basis(:,:), gradient(:,:,:), dL(:), dpoint(:,:,:), assembled(:,:)
    integer, allocatable :: map(:)
    integer :: i, id(2), joint

    xy(:,1) = [0.0_rk, 1.0_rk, 2.0_rk]
    xy(:,2) = [0.0_rk, 1.0_rk, 0.0_rk]
    do i = 1, 2
        controls = xy(i:i+1,:)
        call segment(i)%set(Xc=controls)
        call chain%add_patch(segment(i), id=id(i))
        if (.not. chain%err%ok) error stop "Could not register the segment"
        call segment(i)%ansatz_element(ie=1, Tgc=basis, dTgc_dXg=gradient, &
                                       dL=dL, dmeasure_dXc=dpoint)
        if (.not. segment(i)%err%ok) error stop "Could not evaluate the segment"
        local(2*i-1:2*i,:) = sum(dpoint, dim=3)
    end do

    call chain%connect(id(1), 'right', id(2), 'left', continuity=0)
    if (.not. chain%err%ok) error stop "Could not connect the segments"
    if (.not. chain%is_valid()) error stop "Invalid connection"
    call chain%assemble_sensitivities(local_dXc=local, dXc=assembled, map=map)
    if (size(map) /= 4) error stop "Could not assemble the gradient"
    joint = map(2)
    if (joint < 1 .or. joint /= map(3)) error stop "Expected a shared control"
    print '(a,2f10.6)', "Joint gradient (x, y): ", assembled(joint,:)
end program
```

Rows 2 and 3 of `local` refer to the same physical endpoint, stored once on each patch. Their map entries identify one independent control. The program prints `0.000000  1.414214`. The X derivatives cancel and the Y derivatives add. A geometric check is the total length when the joint is at height \(y\): \(L(y)=2\sqrt{1+y^2}\), so \(L'(1)=\sqrt{2}\).

For the quadrature terms on one rational patch, see [Differentiate arc length](curve-length-sensitivity.html). For a shared volume face, [differentiate both patch volumes](shared-face-volume-sensitivity.html). The [multipatch curve reference](../../type/nurbs_multipatch_curve.html) documents the gradient and map arguments.
