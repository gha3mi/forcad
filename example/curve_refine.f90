!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Refine a rational curve while preserving its geometry.
program curve_refine

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: curve
    real(rk), parameter :: Xt(7) = [&
        0.0_rk, 0.15_rk, 0.32_rk, 0.50_rk, 0.71_rk, 0.88_rk, 1.0_rk]
    real(rk) :: Xc(5,3), reference(3,7), point(3), geometry_error
    integer :: i

    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [0.8_rk, 1.2_rk, 0.2_rk]
    Xc(3,:) = [1.8_rk, 0.4_rk, 0.8_rk]
    Xc(4,:) = [2.7_rk,-0.7_rk, 0.4_rk]
    Xc(5,:) = [3.5_rk, 0.2_rk, 0.0_rk]

    call curve%set(&
        degree = 3,&
        nc     = 5,&
        Xc     = Xc)
    do i = 1, size(Xt)
        reference(:,i) = curve%cmp_Xg(Xt(i))
    end do

    call curve%insert_knots(&
        Xth = [0.25_rk, 0.50_rk, 0.75_rk],&
        r   = [1, 1, 1])
    call curve%remove_knots(&
        Xth = [0.50_rk],&
        r   = [1])
    call curve%elevate_degree(t = 1)

    geometry_error = 0.0_rk
    do i = 1, size(Xt)
        point = curve%cmp_Xg(Xt(i))
        geometry_error = max(&
            geometry_error,&
            abs(point(1) - reference(1,i)),&
            abs(point(2) - reference(2,i)),&
            abs(point(3) - reference(3,i)))
    end do
    write(*,"(a,i0)") "Degree after refinement: ", curve%get_degree()
    write(*,"(a,i0)") "Control points         : ", curve%get_nc()
    write(*,"(a,es12.4)") "Geometry error         : ", geometry_error

    call curve%create(res = 141)
    call curve%export_Xc("vtk/curve_refine_Xc.vtk")
    call curve%export_Xg("vtk/curve_refine_Xg.vtk")
    call curve%export_Xth_in_Xg("vtk/curve_refine_Xth.vtk")
    call curve%show(&
        vtkfile_Xc        = "vtk/curve_refine_Xc.vtk",&
        vtkfile_Xg        = "vtk/curve_refine_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/curve_refine_Xth.vtk")

end program curve_refine
