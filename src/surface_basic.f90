!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Construct and visualize a rational Bezier surface.
program surface_basic

    use forcad, only: rk, nurbs_surface

    implicit none

    type(nurbs_surface) :: surface
    real(rk) :: Xc(9,3)
    real(rk), parameter :: Wc(9) = [&
        1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 0.70_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [1.0_rk, 0.0_rk, 0.3_rk]
    Xc(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
    Xc(4,:) = [0.0_rk, 1.0_rk, 0.2_rk]
    Xc(5,:) = [1.0_rk, 1.0_rk, 1.2_rk]
    Xc(6,:) = [2.0_rk, 1.0_rk, 0.3_rk]
    Xc(7,:) = [0.0_rk, 2.0_rk, 0.0_rk]
    Xc(8,:) = [1.0_rk, 2.0_rk, 0.4_rk]
    Xc(9,:) = [2.0_rk, 2.0_rk, 0.0_rk]

    call surface%set(&
        nc = [3, 3],&
        Xc = Xc,&
        Wc = Wc)
    call surface%create(&
        res1 = 51,&
        res2 = 51)

    write(*,"(a,2(i0,1x))") "Degrees       : ", surface%get_degree()
    write(*,"(a,2(i0,1x))") "Control points: ", surface%get_nc()
    write(*,"(a,l1)") "Rational      : ", surface%is_rational()

    call surface%export_Xc("vtk/surface_basic_Xc.vtk")
    call surface%export_Xg("vtk/surface_basic_Xg.vtk")
    call surface%export_Xth_in_Xg("vtk/surface_basic_Xth.vtk")
    call surface%show(&
        vtkfile_Xc        = "vtk/surface_basic_Xc.vtk",&
        vtkfile_Xg        = "vtk/surface_basic_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_basic_Xth.vtk")

end program surface_basic
