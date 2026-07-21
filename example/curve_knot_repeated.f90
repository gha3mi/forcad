!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates a single rational curve with nonuniform repeated knots.
program curve_knot_repeated

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: curve
    real(rk), allocatable :: Tgc(:)
    real(rk) :: Xc(11,3), Wc(11), knot(15)
    integer :: i

    knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
        0.08_rk, 0.18_rk, 0.18_rk, 0.35_rk, 0.62_rk, 0.62_rk, 0.82_rk, &
        1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    do i = 1, size(Xc, 1)
        Xc(i,1) = 0.55_rk*real(i - 1, rk)
        Xc(i,2) = 0.95_rk*sin(0.62_rk*real(i - 1, rk)) + &
            0.25_rk*cos(1.35_rk*real(i - 1, rk))
        Xc(i,3) = 0.22_rk*real(i - 1, rk) + 0.48_rk*cos(0.45_rk*real(i - 1, rk))
        Wc(i) = 1.0_rk + 0.38_rk*exp(-0.18_rk*(real(i - 6, rk))**2)
    end do

    Xc(5,2) = Xc(5,2) + 1.45_rk
    Xc(8,2) = Xc(8,2) - 1.10_rk
    Wc(5) = 1.85_rk
    Wc(8) = 0.62_rk

    call curve%set(knot, Xc, Wc)
    call curve%create(res=180)
    call curve%basis(0.47_rk, Tgc)

    write(*,"(a,*(1x,g0))") "complex curve degree:", curve%get_degree()
    write(*,"(a,*(1x,g0))") "complex curve nc    :", curve%get_nc()
    write(*,"(a,*(1x,g0))") "basis sum at 0.47  :", sum(Tgc)
    write(*,"(a,*(1x,g0))") "knot multiplicity  :", curve%get_multiplicity()
    write(*,"(a,*(1x,g0))") "knot continuity    :", curve%get_continuity()

    call curve%export_Xc("vtk/curve_knot_repeated_Xc.vtk")
    call curve%export_Xg("vtk/curve_knot_repeated_Xg.vtk")
    call curve%export_iges("iges/curve_knot_repeated.iges")
    call curve%export_Xth_in_Xg("vtk/curve_knot_repeated_Xthg.vtk")
    call curve%show("vtk/curve_knot_repeated_Xc.vtk", "vtk/curve_knot_repeated_Xg.vtk", "vtk/curve_knot_repeated_Xthg.vtk")

    call curve%finalize()
    deallocate(Tgc)

end program curve_knot_repeated
