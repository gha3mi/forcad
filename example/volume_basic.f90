!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Construct and visualize a trilinear B-spline volume.
program volume_basic

    use forcad, only: rk, nurbs_volume

    implicit none

    type(nurbs_volume) :: volume
    real(rk) :: Xc(8,3)

    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [2.0_rk, 0.0_rk, 0.0_rk]
    Xc(3,:) = [0.0_rk, 1.5_rk, 0.0_rk]
    Xc(4,:) = [2.0_rk, 1.5_rk, 0.2_rk]
    Xc(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
    Xc(6,:) = [2.0_rk, 0.0_rk, 1.2_rk]
    Xc(7,:) = [0.0_rk, 1.5_rk, 1.1_rk]
    Xc(8,:) = [2.0_rk, 1.5_rk, 1.5_rk]

    call volume%set(&
        nc = [2, 2, 2],&
        Xc = Xc)
    call volume%create(&
        res1 = 17,&
        res2 = 15,&
        res3 = 13)

    write(*,"(a,3(i0,1x))") "Degrees       : ", volume%get_degree()
    write(*,"(a,3(i0,1x))") "Control points: ", volume%get_nc()
    write(*,"(a,l1)") "Rational      : ", volume%is_rational()

    call volume%export_Xc("vtk/volume_basic_Xc.vtk")
    call volume%export_Xg("vtk/volume_basic_Xg.vtk")
    call volume%export_Xth_in_Xg("vtk/volume_basic_Xth.vtk")
    call volume%show(&
        vtkfile_Xc        = "vtk/volume_basic_Xc.vtk",&
        vtkfile_Xg        = "vtk/volume_basic_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_basic_Xth.vtk")

end program volume_basic
