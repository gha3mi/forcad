!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Construct a rational curve with a nonuniform periodic knot vector.
program curve_knot_periodic

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: curve
    real(rk) :: Xc(6,3)
    real(rk), parameter :: knot(9) = [&
        -2.2_rk, -1.6_rk, 0.0_rk, 0.7_rk, 1.8_rk, 2.4_rk, 4.0_rk, 4.7_rk, 5.8_rk]
    real(rk), parameter :: Wc(6) = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]

    Xc(1,:) = [ 1.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [ 0.0_rk, 1.0_rk, 0.2_rk]
    Xc(3,:) = [-1.0_rk, 0.0_rk, 0.0_rk]
    Xc(4,:) = [ 0.0_rk,-1.0_rk,-0.2_rk]
    Xc(5,:) = Xc(1,:)
    Xc(6,:) = Xc(2,:)

    call curve%set(&
        knot            = knot,&
        Xc              = Xc,&
        Wc              = Wc,&
        degree          = 2,&
        wrap_parameters = .true.)
    call curve%create(res = 161)

    write(*,"(a,i0)") "Degree            : ", curve%get_degree()
    write(*,"(a,i0)") "Control points     : ", curve%get_nc()
    write(*,"(a,l1)") "Parameter wrapping : ", curve%get_parameter_wrapping()

    call curve%export_Xc("vtk/curve_knot_periodic_Xc.vtk")
    call curve%export_Xg("vtk/curve_knot_periodic_Xg.vtk")
    call curve%export_Xth_in_Xg("vtk/curve_knot_periodic_Xth.vtk")
    call curve%show(&
        vtkfile_Xc        = "vtk/curve_knot_periodic_Xc.vtk",&
        vtkfile_Xg        = "vtk/curve_knot_periodic_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/curve_knot_periodic_Xth.vtk")

end program curve_knot_periodic
