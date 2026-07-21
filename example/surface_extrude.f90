!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Extrude an exact rational circle into a cylindrical NURBS surface.
program surface_extrude

    use forcad, only: rk, nurbs_curve, nurbs_surface, extrude

    implicit none

    type(nurbs_curve) :: profile
    type(nurbs_surface) :: cylinder

    call profile%set_circle(&
        center = [0.0_rk,0.0_rk,0.0_rk],&
        radius = 1.25_rk)
    cylinder = extrude(&
        curve  = profile,&
        vector = [0.0_rk,0.0_rk,3.5_rk])

    write(*,"(a,l1)") "Valid extrusion        : ", cylinder%err%ok
    write(*,"(a,l1)") "Rational geometry      : ", cylinder%is_rational()
    write(*,"(a,2(i0,1x))") "Degrees                 : ", cylinder%get_degree()
    write(*,"(a,2(i0,1x))") "Control points          : ", cylinder%get_nc()

    call cylinder%create(81, 31)
    call cylinder%export_Xc("vtk/surface_extrude_Xc.vtk")
    call cylinder%export_Xg("vtk/surface_extrude_Xg.vtk")
    call cylinder%export_Xth_in_Xg("vtk/surface_extrude_Xth.vtk", 17)
    call cylinder%show(&
        vtkfile_Xc        = "vtk/surface_extrude_Xc.vtk",&
        vtkfile_Xg        = "vtk/surface_extrude_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_extrude_Xth.vtk")

    call profile%finalize()
    call cylinder%finalize()

end program surface_extrude
