!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program curve_knot_unclamped

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: curve
    real(rk) :: knot(9)
    real(rk), allocatable :: Xc(:,:), Wc(:)

    knot = [-0.3_rk, 0.1_rk, 0.4_rk, 0.9_rk, 1.6_rk, 2.2_rk, 2.9_rk, 3.4_rk, 3.9_rk]

    allocate(Xc(6, 3), Wc(6))
    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [0.6_rk, 0.5_rk, 0.1_rk]
    Xc(3,:) = [1.4_rk, 0.9_rk, 0.3_rk]
    Xc(4,:) = [2.4_rk, 0.5_rk, 0.2_rk]
    Xc(5,:) = [3.1_rk, 0.1_rk, 0.4_rk]
    Xc(6,:) = [3.7_rk, 0.0_rk, 0.0_rk]
    Wc = [1.0_rk, 0.95_rk, 1.10_rk, 0.90_rk, 1.05_rk, 1.0_rk]

    call curve%set(knot, Xc, Wc, degree=2)
    call curve%create(120)
    call curve%export_Xc("vtk/curve_knot_unclamped_Xc.vtk")
    call curve%export_Xg("vtk/curve_knot_unclamped_Xg.vtk")
    call curve%export_Xth_in_Xg("vtk/curve_knot_unclamped_Xthg.vtk")
    call curve%show(&
        "vtk/curve_knot_unclamped_Xc.vtk", &
        "vtk/curve_knot_unclamped_Xg.vtk", &
        "vtk/curve_knot_unclamped_Xthg.vtk")

    write(*,"(a,2(f6.2,1x))") "active interval:", curve%get_knot(curve%get_degree() + 1), curve%get_knot(curve%get_nc() + 1)

    call curve%finalize()

end program curve_knot_unclamped
