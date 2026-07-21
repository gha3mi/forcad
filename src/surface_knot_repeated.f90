!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates a single rational surface with anisotropic repeated knots.
program surface_knot_repeated

    use forcad, only: rk, nurbs_surface

    implicit none

    type(nurbs_surface) :: surface
    real(rk), allocatable :: Xc(:,:), Wc(:), Tgc(:)
    real(rk) :: knot1(13), knot2(10), pt(2)
    real(rk) :: u, v, x, y, r2
    integer :: i, j, idx

    knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
        0.12_rk, 0.30_rk, 0.30_rk, 0.55_rk, 0.82_rk, &
        1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
    knot2 = [0.0_rk, 0.0_rk, 0.0_rk, &
        0.20_rk, 0.42_rk, 0.42_rk, 0.75_rk, &
        1.0_rk, 1.0_rk, 1.0_rk]

    allocate(Xc(63, 3), Wc(63))
    do concurrent (j = 1:7, i = 1:9) local(idx, u, v, x, y, r2)
        idx = (j - 1)*9 + i
        u = real(i - 1, rk)/8.0_rk
        v = real(j - 1, rk)/6.0_rk
        x = 5.2_rk*(u - 0.5_rk) + 0.22_rk*sin(3.0_rk*acos(-1.0_rk)*v)
        y = 3.6_rk*(v - 0.5_rk) + 0.30_rk*sin(2.0_rk*acos(-1.0_rk)*u)
        r2 = x*x + 0.75_rk*y*y
        Xc(idx,1) = x*(1.0_rk + 0.06_rk*cos(2.0_rk*acos(-1.0_rk)*v))
        Xc(idx,2) = y*(1.0_rk + 0.08_rk*sin(2.0_rk*acos(-1.0_rk)*u))
        Xc(idx,3) = 0.85_rk*exp(-0.16_rk*r2) + 0.34_rk*sin(2.5_rk*acos(-1.0_rk)*u)* &
            cos(1.5_rk*acos(-1.0_rk)*v)
        Wc(idx) = 1.0_rk + 0.35_rk*exp(-2.0_rk*((u - 0.54_rk)**2 + (v - 0.48_rk)**2))
    end do

    pt = [0.47_rk, 0.58_rk]
    call surface%set(knot1, knot2, Xc, Wc=Wc, degree=[3, 2])
    call surface%create(70, 56)
    call surface%basis(pt, Tgc)

    write(*,"(a,*(1x,g0))") "complex surface degree:", surface%get_degree()
    write(*,"(a,*(1x,g0))") "complex surface nc    :", surface%get_nc()
    write(*,"(a,*(1x,g0))") "basis sum at point    :", sum(Tgc)
    write(*,"(a,*(1x,g0))") "knot1 multiplicity    :", surface%get_multiplicity(1)
    write(*,"(a,*(1x,g0))") "knot2 multiplicity    :", surface%get_multiplicity(2)

    call surface%export_Xc("vtk/surface_knot_repeated_Xc.vtk")
    call surface%export_Xg("vtk/surface_knot_repeated_Xg.vtk")
    call surface%export_iges("iges/surface_knot_repeated.iges")
    call surface%export_Xth_in_Xg("vtk/surface_knot_repeated_Xthg.vtk")
    call surface%show("vtk/surface_knot_repeated_Xc.vtk", "vtk/surface_knot_repeated_Xg.vtk", "vtk/surface_knot_repeated_Xthg.vtk")

    call surface%finalize()
    deallocate(Xc, Wc, Tgc)

end program surface_knot_repeated
