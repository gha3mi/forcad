!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates a single volume with different nonuniform knots in each direction.
program volume_knot_repeated

    use forcad, only: rk, nurbs_volume

    implicit none

    type(nurbs_volume) :: volume
    real(rk), allocatable :: Xc(:,:), Wc(:), Tgc(:)
    real(rk) :: knot1(10), knot2(12), knot3(10), pt(3)
    real(rk) :: u, v, w, twist, scale, ca, sa, x0, y0, z0
    integer :: i, j, k, idx

    knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.16_rk, 0.36_rk, 0.36_rk, 0.70_rk, &
        1.0_rk, 1.0_rk, 1.0_rk]
    knot2 = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.18_rk, 0.42_rk, 0.42_rk, 0.73_rk, &
        1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
    knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.10_rk, 0.28_rk, 0.58_rk, 0.58_rk, &
        1.0_rk, 1.0_rk, 1.0_rk]

    allocate(Xc(392, 3), Wc(392))
    do concurrent (k = 1:7, j = 1:8, i = 1:7) local(idx, u, v, w, twist, scale, ca, sa, x0, y0, z0)
        idx = ((k - 1)*8 + j - 1)*7 + i
        u = real(i - 1, rk)/6.0_rk
        v = real(j - 1, rk)/7.0_rk
        w = real(k - 1, rk)/6.0_rk
        x0 = 3.0_rk*(u - 0.5_rk)
        y0 = 2.4_rk*(v - 0.5_rk)
        z0 = 2.8_rk*w
        twist = 0.55_rk*acos(-1.0_rk)*w + 0.12_rk*sin(2.0_rk*acos(-1.0_rk)*v)
        scale = 1.0_rk - 0.22_rk*w + 0.08_rk*sin(2.0_rk*acos(-1.0_rk)*u)*sin(acos(-1.0_rk)*v)
        ca = cos(twist)
        sa = sin(twist)
        Xc(idx,1) = scale*(ca*x0 - sa*y0)
        Xc(idx,2) = scale*(sa*x0 + ca*y0)
        Xc(idx,3) = z0 + 0.32_rk*sin(acos(-1.0_rk)*u)*sin(acos(-1.0_rk)*v) + &
            0.18_rk*cos(2.0_rk*acos(-1.0_rk)*w)
        Wc(idx) = 1.0_rk + 0.18_rk*exp(-4.0_rk*((u - 0.48_rk)**2 + (v - 0.52_rk)**2 + &
            (w - 0.50_rk)**2))
    end do

    pt = [0.44_rk, 0.57_rk, 0.49_rk]
    call volume%set(knot1, knot2, knot3, Xc, Wc=Wc, degree=[2, 3, 2])
    call volume%create(22, 24, 20)
    call volume%basis(pt, Tgc)

    write(*,"(a,*(1x,g0))") "complex volume degree:", volume%get_degree()
    write(*,"(a,*(1x,g0))") "complex volume nc    :", volume%get_nc()
    write(*,"(a,*(1x,g0))") "basis sum at point   :", sum(Tgc)
    write(*,"(a,*(1x,g0))") "knot1 multiplicity   :", volume%get_multiplicity(1)
    write(*,"(a,*(1x,g0))") "knot2 multiplicity   :", volume%get_multiplicity(2)
    write(*,"(a,*(1x,g0))") "knot3 multiplicity   :", volume%get_multiplicity(3)

    call volume%export_Xc("vtk/volume_knot_repeated_Xc.vtk")
    call volume%export_Xg("vtk/volume_knot_repeated_Xg.vtk")
    call volume%export_Xth_in_Xg("vtk/volume_knot_repeated_Xthg.vtk")
    call volume%show("vtk/volume_knot_repeated_Xc.vtk", "vtk/volume_knot_repeated_Xg.vtk", "vtk/volume_knot_repeated_Xthg.vtk")

    call volume%finalize()
    deallocate(Xc, Wc, Tgc)

end program volume_knot_repeated
