!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Transform curve control points and regenerate the geometry.
program curve_transform

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: curve

    call curve%set_circle(&
        center = [0.0_rk, 0.0_rk],&
        radius = 1.0_rk)
    call curve%rotate_Xc(&
        alpha = 35.0_rk,&
        beta  = 0.0_rk,&
        theta = 0.0_rk)
    call curve%translate_Xc([1.5_rk,-0.4_rk])
    call curve%create(res = 161)

    call curve%export_Xc("vtk/curve_transform_Xc.vtk")
    call curve%export_Xg("vtk/curve_transform_Xg.vtk")
    call curve%export_Xth_in_Xg("vtk/curve_transform_Xth.vtk")
    call curve%show(&
        vtkfile_Xc        = "vtk/curve_transform_Xc.vtk",&
        vtkfile_Xg        = "vtk/curve_transform_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/curve_transform_Xth.vtk")

end program curve_transform
