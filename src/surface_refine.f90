!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Refine a tensor-product surface directionally while preserving its geometry.
program surface_refine

    use forcad, only: rk, nurbs_surface

    implicit none

    type(nurbs_surface) :: surface
    real(rk), parameter :: Xt(2,6) = reshape([&
        0.0_rk, 0.0_rk, 0.18_rk, 0.76_rk, 0.35_rk, 0.24_rk,&
        0.62_rk, 0.58_rk, 0.83_rk, 0.91_rk, 1.0_rk, 1.0_rk], [2,6])
    real(rk) :: reference(3,6), point(3), geometry_error
    integer :: i

    call surface%set_tetragon(&
        L  = [3.0_rk, 2.0_rk],&
        nc = [4, 4])
    do i = 1, size(Xt, 2)
        reference(:,i) = surface%cmp_Xg(Xt(:,i))
    end do

    call surface%insert_knots(&
        dir = 1,&
        Xth = [0.25_rk, 0.50_rk, 0.75_rk],&
        r   = [1, 1, 1])
    call surface%insert_knots(&
        dir = 2,&
        Xth = [0.40_rk, 0.70_rk],&
        r   = [1, 1])
    call surface%remove_knots(&
        dir = 1,&
        Xth = [0.50_rk],&
        r   = [1])
    call surface%elevate_degree(&
        dir = 2,&
        t   = 1)

    geometry_error = 0.0_rk
    do i = 1, size(Xt, 2)
        point = surface%cmp_Xg(Xt(:,i))
        geometry_error = max(&
            geometry_error,&
            abs(point(1) - reference(1,i)),&
            abs(point(2) - reference(2,i)),&
            abs(point(3) - reference(3,i)))
    end do
    write(*,"(a,2(i0,1x))") "Degrees after refinement: ", surface%get_degree()
    write(*,"(a,2(i0,1x))") "Control points         : ", surface%get_nc()
    write(*,"(a,es12.4)") "Geometry error         : ", geometry_error

    call surface%create(&
        res1 = 61,&
        res2 = 49)
    call surface%export_Xc("vtk/surface_refine_Xc.vtk")
    call surface%export_Xg("vtk/surface_refine_Xg.vtk")
    call surface%export_Xth_in_Xg("vtk/surface_refine_Xth.vtk")
    call surface%show(&
        vtkfile_Xc        = "vtk/surface_refine_Xc.vtk",&
        vtkfile_Xg        = "vtk/surface_refine_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_refine_Xth.vtk")

end program surface_refine
