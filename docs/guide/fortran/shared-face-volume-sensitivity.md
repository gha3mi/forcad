title: Differentiate the volumes on either side of a shared face
category: Assembly and analysis

Moving an interior face can enlarge one volume and shrink its neighbor while the exterior stays fixed. This example computes how **each volume** responds when you change the physical \(x\) coordinate or the weight of the same NURBS control. It verifies each derivative by rebuilding and measuring both patches.

## Choose the shared control

Start with two boxes that fill \([0,3]\times[0,1]\times[0,1]\). The left box ends at \(x=1\), where its **U max** face meets the right box’s **U min** face. Represent each box as a triquadratic NURBS volume with \(3\times3\times3\) control points. The two face parameterizations have the same V and W directions.

Move the middle face control to \((1.5,0.5,0.5)\) and set its weight to \(1.25\) on both patches. Keep every other control fixed with unit weight. The face now bows into the right volume. With all nine face controls in the original plane, changing just the middle weight would leave the face planar and give a zero volume response.

<figure class="doc-figure" markdown="1">

[![Cross-section of the common face at z = 0.5. Both edits move its middle toward positive x while the endpoints stay at x = 1.](../../media/shared-face-volume-sensitivity.webp)](../../media/shared-face-volume-sensitivity.webp)

<figcaption>The face at z = 0.5. Starting from the baseline, move the control to x = 1.75 or increase its weight to 2. The curves show physical positions evaluated on the face. The control itself lies farther to the right.</figcaption>
</figure>

## Read the sensitivities

The baseline volumes are \(V_\mathrm{L}=1.06428469\) and \(V_\mathrm{R}=1.93571531\). The table compares ForCAD’s integrated derivatives with centered finite differences (FD). The coordinate column uses a displacement along positive \(x\), and the weight column uses a change in the ordinary, dimensionless weight \(w_c\).

| Volume | \(\partial V/\partial x_c\) | FD in \(x_c\) | \(\partial V/\partial w_c\) | FD in \(w_c\) |
| --- | ---: | ---: | ---: | ---: |
| Left | +0.12856938 | +0.12856938 | +0.03226242 | +0.03226242 |
| Right | −0.12856938 | −0.12856938 | −0.03226242 | −0.03226242 |
| Total | 0.00000000 | 0.00000000 | 0.00000000 | 0.00000000 |

For a small increase of 0.01 in the shared weight, the first-order prediction is a gain of about 0.000323 in the left volume and an equal loss in the right. The total volume stays at 3 because both patches share the deformed face and the outer boundary stays fixed.

## Run the check

1. **Accumulate each patch separately.** A quadratic patch in this example has one knot span in each direction. Sum the outputs of `ansatz_element` over its quadrature points. `dmeasure_dXc` gives coordinate derivatives. `dmeasure_dlogWc` gives derivatives with respect to \(\eta_c=\log w_c\), so divide by the baseline weight to obtain \(\partial V/\partial w_c=(\partial V/\partial\eta_c)/w_c\).
2. **Identify the shared variable.** With U varying fastest, the face control has local index 15 on the left and 13 on the right. Add the patch offsets and check that `assemble_sensitivities` maps both rows to one design variable. The two local entries give the separate volume responses. The assembled entry gives the response of their sum.
3. **Perturb both copies together.** For each variable \(s\in\{x_c,w_c\}\), measure \([V_p(s+h)-V_p(s-h)]/(2h)\) for each patch \(p\) with `cmp_volume`. Each trial rebuilds the connected geometry and changes only \(s\). The example uses the same eight-point Gauss rule per direction for the analytical and finite-difference calculations.

Use the [fpm project setup](../tutorials/first-curve.html#add-the-dependency), save the program as `app/main.f90`, and run `fpm run --profile release`.


```fortran
program shared_face_sensitivity
    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume
    implicit none
    integer, parameter :: face_control(2) = [15,13], rule(3) = [8,8,8]
    character(len=5), parameter :: label(3) = [character(len=5) :: 'Left', 'Right', 'Total']
    real(rk), parameter :: x0 = 1.5_rk, w0 = 1.25_rk, h = epsilon(1._rk)**0.25_rk
    type(nurbs_multipatch_volume) :: model
    type(nurbs_volume) :: patch
    integer, allocatable :: offset(:), map(:)
    integer :: p, first, last, row(2), joint
    real(rk), allocatable :: basis(:,:), grad(:,:,:), dV(:), dpoint(:,:,:), dlogw(:,:)
    real(rk), allocatable :: local_x(:,:), local_eta(:), global_x(:,:), global_eta(:)
    real(rk) :: volume(2), analytic(3,2), difference(3,2)

    call build_pair(x0, w0, model)
    offset = model%cmp_dof_offsets()
    allocate(local_x(offset(3),3), local_eta(offset(3)))
    do p = 1, 2
        patch = model%get_patch(p)
        call patch%ansatz_element(ie=1, Tgc=basis, dTgc_dXg=grad, dV=dV, ngauss=rule, &
                                  dmeasure_dXc=dpoint, dmeasure_dlogWc=dlogw, strict=.true.)
        if (.not. patch%err%ok) error stop "Volume quadrature failed"
        first = offset(p) + 1
        last = offset(p+1)
        volume(p) = sum(dV)
        local_x(first:last,:) = sum(dpoint, dim=3)
        local_eta(first:last) = sum(dlogw, dim=2)
    end do

    call model%assemble_sensitivities(local_dXc=local_x, local_dlogWc=local_eta, &
                                      dXc=global_x, dlogWc=global_eta, map=map)
    row = offset(1:2) + face_control
    if (.not. allocated(map)) error stop "Missing design variable map"
    if (size(map) /= offset(3)) error stop "Incomplete design variable map"
    joint = map(row(1))
    if (joint < 1 .or. joint /= map(row(2))) error stop "Expected one shared control"

    analytic(1:2,1) = local_x(row,1)
    analytic(1:2,2) = local_eta(row)/w0
    analytic(3,:) = [global_x(joint,1), global_eta(joint)/w0]

    difference(1:2,1) = (volumes_at(x0+h,w0) - volumes_at(x0-h,w0))/(2*h)
    difference(1:2,2) = (volumes_at(x0,w0+h) - volumes_at(x0,w0-h))/(2*h)
    difference(3,:) = sum(difference(1:2,:), dim=1)

    print '(a,2f12.8)', "Volumes (left, right): ", volume
    print '(a)', "                 dV/dx       FD(x)       dV/dw       FD(w)"
    do p = 1, 3
        print '(a5,8x,4f12.8)', label(p), analytic(p,1), difference(p,1), analytic(p,2), difference(p,2)
    end do
    if (.not. all(abs(analytic-difference) < 10*h**2)) error stop "Derivative check failed"
    if (abs(sum(volume)-3._rk) > 10*h**2) error stop "Incorrect total volume"
    if (any(abs(analytic(3,:)) > 10*h**2)) error stop "Exterior volume changed"

contains

    subroutine build_pair(x, weight, pair)
        real(rk), intent(in) :: x, weight
        type(nurbs_multipatch_volume), intent(out) :: pair
        type(nurbs_volume) :: member(2)
        integer :: i, id(2)
        do i = 1, 2
            call member(i)%set_hexahedron(L=[real(i,rk),1._rk,1._rk], nc=[3,3,3])
            if (i == 2) call member(i)%translate_Xc([1._rk,0._rk,0._rk])
            call member(i)%modify_Xc(X=x, num=face_control(i), dir=1)
            call member(i)%modify_Wc(W=weight, num=face_control(i))
            if (.not. member(i)%err%ok) error stop "Could not set the face control"
            call pair%add_patch(member(i), id=id(i))
            if (.not. pair%err%ok) error stop "Could not add the volume"
        end do
        call pair%connect(id(1), 'u_max', id(2), 'u_min', continuity=0)
        if (.not. pair%err%ok .or. .not. pair%is_valid()) error stop "Invalid face connection"
    end subroutine

    function volumes_at(x, weight) result(measured)
        real(rk), intent(in) :: x, weight
        real(rk) :: measured(2)
        type(nurbs_multipatch_volume) :: trial
        type(nurbs_volume) :: member
        integer :: i
        call build_pair(x, weight, trial)
        do i = 1, 2
            member = trial%get_patch(i)
            call member%cmp_volume(measured(i), ngauss=rule, strict=.true.)
            if (.not. member%err%ok) error stop "Trial measurement failed"
        end do
    end function
end program shared_face_sensitivity
```

For patches with several knot spans, accumulate each element’s derivatives into its patch control rows using the [element control IDs](control-ids-on-a-knot-span.html) before assembling across the interface. For the connection itself, see [Join two volume patches](join-volume-patches.html).
