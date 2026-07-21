!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Revolve a cubic vessel profile into an exact rational NURBS surface.
program surface_revolve

    use forcad, only: rk, nurbs_curve, nurbs_surface, revolve

    implicit none

    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: knot(11) = [0.0_rk,0.0_rk,0.0_rk,0.0_rk, &
        0.25_rk,0.50_rk,0.75_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk]
    real(rk), parameter :: Xc(7,3) = reshape([0.55_rk,0.78_rk,1.18_rk,0.92_rk,1.32_rk,1.02_rk,0.68_rk, &
        0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk, &
        -2.0_rk,-1.55_rk,-0.90_rk,-0.10_rk,0.75_rk,1.50_rk,2.10_rk], [7,3])

    type(nurbs_curve) :: profile
    type(nurbs_surface) :: vessel

    call profile%set(&
        knot   = knot,&
        Xc     = Xc,&
        degree = 3)
    vessel = revolve(&
        curve          = profile,&
        axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
        axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
        angle          = 2.0_rk*pi)

    write(*,"(a,l1)") "Valid revolution       : ", vessel%err%ok
    write(*,"(a,l1)") "Rational geometry      : ", vessel%is_rational()
    write(*,"(a,2(i0,1x))") "Degrees                 : ", vessel%get_degree()
    write(*,"(a,2(i0,1x))") "Control points          : ", vessel%get_nc()

    call vessel%create(61, 81)
    call vessel%export_Xc("vtk/surface_revolve_Xc.vtk")
    call vessel%export_Xg("vtk/surface_revolve_Xg.vtk")
    call vessel%export_Xth_in_Xg("vtk/surface_revolve_Xth.vtk", 17)
    call vessel%show(&
        vtkfile_Xc        = "vtk/surface_revolve_Xc.vtk",&
        vtkfile_Xg        = "vtk/surface_revolve_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_revolve_Xth.vtk")

    call profile%finalize()
    call vessel%finalize()

end program surface_revolve
