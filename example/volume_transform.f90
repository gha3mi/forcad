!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Transform volume control points and regenerate the geometry.
program volume_transform

    use forcad, only: rk, nurbs_volume

    implicit none

    type(nurbs_volume) :: volume

    call volume%set_ring(&
        center  = [0.0_rk, 0.0_rk, 0.0_rk],&
        radius1 = 0.8_rk,&
        radius2 = 1.8_rk,&
        length  = 1.2_rk)
    call volume%rotate_Xc(&
        alpha = 18.0_rk,&
        beta  = 12.0_rk,&
        theta = 28.0_rk)
    call volume%translate_Xc([0.5_rk,-0.8_rk, 1.0_rk])
    call volume%create(&
        res1 = 49,&
        res2 = 15,&
        res3 = 13)

    call volume%export_Xc("vtk/volume_transform_Xc.vtk")
    call volume%export_Xg("vtk/volume_transform_Xg.vtk")
    call volume%export_Xth_in_Xg("vtk/volume_transform_Xth.vtk")
    call volume%show(&
        vtkfile_Xc        = "vtk/volume_transform_Xc.vtk",&
        vtkfile_Xg        = "vtk/volume_transform_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_transform_Xth.vtk")

end program volume_transform
