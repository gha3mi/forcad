!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program curve_knot_clamped

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: curve
    real(rk) :: knot(9)
    real(rk), allocatable :: Xc(:,:), Wc(:)

    knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.25_rk, 0.55_rk, 0.80_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    allocate(Xc(6, 3), Wc(6))
    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [0.8_rk, 0.4_rk, 0.0_rk]
    Xc(3,:) = [1.7_rk, 1.2_rk, 0.2_rk]
    Xc(4,:) = [2.8_rk, 0.8_rk, 0.0_rk]
    Xc(5,:) = [3.6_rk, 0.2_rk, 0.3_rk]
    Xc(6,:) = [4.2_rk, 0.0_rk, 0.0_rk]
    Wc = [1.0_rk, 0.8_rk, 1.2_rk, 1.0_rk, 0.9_rk, 1.0_rk]

    call curve%set(knot, Xc, Wc, degree=2)
    call curve%create(120)
    call curve%export_Xc("vtk/curve_knot_clamped_Xc.vtk")
    call curve%export_Xg("vtk/curve_knot_clamped_Xg.vtk")
    call curve%export_Xth_in_Xg("vtk/curve_knot_clamped_Xthg.vtk")
    call curve%show("vtk/curve_knot_clamped_Xc.vtk", "vtk/curve_knot_clamped_Xg.vtk", "vtk/curve_knot_clamped_Xthg.vtk")

    write(*,"(a,i0)") "degree:", curve%get_degree()
    write(*,"(a,*(f6.2,1x))") "knot:", curve%get_knot()

    call curve%finalize()

end program curve_knot_clamped
