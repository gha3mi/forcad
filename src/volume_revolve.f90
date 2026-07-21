!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Revolve an exact annular surface into a 270-degree toroidal NURBS volume.
program volume_revolve

    use forcad, only: rk, nurbs_surface, nurbs_volume, revolve

    implicit none

    real(rk), parameter :: pi = acos(-1.0_rk)

    type(nurbs_surface) :: cross_section
    type(nurbs_volume) :: torus_sector

    call cross_section%set_ring(&
        center  = [0.0_rk,0.0_rk,0.0_rk],&
        radius1 = 0.22_rk,&
        radius2 = 0.48_rk)
    call cross_section%rotate_Xc(&
        alpha = 90.0_rk,&
        beta  = 0.0_rk,&
        theta = 0.0_rk)
    call cross_section%translate_Xc([2.2_rk,0.0_rk,0.0_rk])
    torus_sector = revolve(&
        surface        = cross_section,&
        axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
        axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
        angle          = 1.5_rk*pi)

    write(*,"(a,l1)") "Valid revolution       : ", torus_sector%err%ok
    write(*,"(a,l1)") "Rational geometry      : ", torus_sector%is_rational()
    write(*,"(a,3(i0,1x))") "Degrees                 : ", torus_sector%get_degree()
    write(*,"(a,3(i0,1x))") "Control points          : ", torus_sector%get_nc()

    call torus_sector%create(61, 13, 81)
    call torus_sector%export_Xc("vtk/volume_revolve_Xc.vtk")
    call torus_sector%export_Xg("vtk/volume_revolve_Xg.vtk")
    call torus_sector%export_Xth_in_Xg("vtk/volume_revolve_Xth.vtk", 13)
    call torus_sector%show(&
        vtkfile_Xc        = "vtk/volume_revolve_Xc.vtk",&
        vtkfile_Xg        = "vtk/volume_revolve_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_revolve_Xth.vtk")

    call cross_section%finalize()
    call torus_sector%finalize()

end program volume_revolve
