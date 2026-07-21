!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Construct and visualize a rational Bezier curve.
program curve_basic

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: curve
    real(rk) :: Xc(4,3)
    real(rk), parameter :: Wc(4) = [1.0_rk, 0.75_rk, 1.25_rk, 1.0_rk]

    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [1.0_rk, 1.6_rk, 0.4_rk]
    Xc(3,:) = [2.2_rk,-0.8_rk, 1.0_rk]
    Xc(4,:) = [3.2_rk, 0.2_rk, 0.3_rk]

    call curve%set(&
        Xc = Xc,&
        Wc = Wc)
    call curve%create(res = 121)

    write(*,"(a,i0)") "Degree        : ", curve%get_degree()
    write(*,"(a,i0)") "Control points: ", curve%get_nc()
    write(*,"(a,l1)") "Rational      : ", curve%is_rational()

    call curve%export_Xc("vtk/curve_basic_Xc.vtk")
    call curve%export_Xg("vtk/curve_basic_Xg.vtk")
    call curve%export_Xth_in_Xg("vtk/curve_basic_Xth.vtk")
    call curve%show(&
        vtkfile_Xc        = "vtk/curve_basic_Xc.vtk",&
        vtkfile_Xg        = "vtk/curve_basic_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/curve_basic_Xth.vtk")

end program curve_basic
