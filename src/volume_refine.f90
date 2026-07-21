!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Refine a tensor-product volume directionally while preserving its geometry.
program volume_refine

    use forcad, only: rk, nurbs_volume

    implicit none

    type(nurbs_volume) :: volume
    real(rk), parameter :: Xt(3,5) = reshape([&
        0.0_rk, 0.0_rk, 0.0_rk, 0.20_rk, 0.68_rk, 0.37_rk,&
        0.43_rk, 0.31_rk, 0.72_rk, 0.77_rk, 0.84_rk, 0.16_rk,&
        1.0_rk, 1.0_rk, 1.0_rk], [3,5])
    real(rk) :: reference(3,5), point(3), geometry_error
    integer :: i

    call volume%set_hexahedron(&
        L  = [3.0_rk, 2.0_rk, 1.5_rk],&
        nc = [3, 3, 3])
    do i = 1, size(Xt, 2)
        reference(:,i) = volume%cmp_Xg(Xt(:,i))
    end do

    call volume%insert_knots(&
        dir = 1,&
        Xth = [0.30_rk, 0.65_rk],&
        r   = [1, 1])
    call volume%insert_knots(&
        dir = 2,&
        Xth = [0.50_rk],&
        r   = [1])
    call volume%insert_knots(&
        dir = 3,&
        Xth = [0.25_rk, 0.75_rk],&
        r   = [1, 1])
    call volume%remove_knots(&
        dir = 3,&
        Xth = [0.25_rk],&
        r   = [1])
    call volume%elevate_degree(&
        dir = 1,&
        t   = 1)

    geometry_error = 0.0_rk
    do i = 1, size(Xt, 2)
        point = volume%cmp_Xg(Xt(:,i))
        geometry_error = max(&
            geometry_error,&
            abs(point(1) - reference(1,i)),&
            abs(point(2) - reference(2,i)),&
            abs(point(3) - reference(3,i)))
    end do
    write(*,"(a,3(i0,1x))") "Degrees after refinement: ", volume%get_degree()
    write(*,"(a,3(i0,1x))") "Control points         : ", volume%get_nc()
    write(*,"(a,es12.4)") "Geometry error         : ", geometry_error

    call volume%create(&
        res1 = 23,&
        res2 = 19,&
        res3 = 17)
    call volume%export_Xc("vtk/volume_refine_Xc.vtk")
    call volume%export_Xg("vtk/volume_refine_Xg.vtk")
    call volume%export_Xth_in_Xg("vtk/volume_refine_Xth.vtk")
    call volume%show(&
        vtkfile_Xc        = "vtk/volume_refine_Xc.vtk",&
        vtkfile_Xg        = "vtk/volume_refine_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_refine_Xth.vtk")

end program volume_refine
