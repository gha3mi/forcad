!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Transform surface control points and regenerate the geometry.
program surface_transform

    use forcad, only: rk, nurbs_surface

    implicit none

    type(nurbs_surface) :: surface

    call surface%set_ring(&
        center  = [0.0_rk, 0.0_rk],&
        radius1 = 0.8_rk,&
        radius2 = 2.0_rk)
    call surface%rotate_Xc(&
        alpha = 25.0_rk,&
        beta  = 0.0_rk,&
        theta = 0.0_rk)
    call surface%translate_Xc([0.7_rk, 1.1_rk])
    call surface%create(&
        res1 = 81,&
        res2 = 25)

    call surface%export_Xc("vtk/surface_transform_Xc.vtk")
    call surface%export_Xg("vtk/surface_transform_Xg.vtk")
    call surface%export_Xth_in_Xg("vtk/surface_transform_Xth.vtk")
    call surface%show(&
        vtkfile_Xc        = "vtk/surface_transform_Xc.vtk",&
        vtkfile_Xg        = "vtk/surface_transform_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_transform_Xth.vtk")

end program surface_transform
