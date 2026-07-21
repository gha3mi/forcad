!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Extrude an exact rational ring surface into an oblique annular NURBS volume.
program volume_extrude

    use forcad, only: rk, nurbs_surface, nurbs_volume, extrude

    implicit none

    type(nurbs_surface) :: profile
    type(nurbs_volume) :: sleeve

    call profile%set_ring(&
        center  = [0.0_rk,0.0_rk,0.0_rk],&
        radius1 = 0.55_rk,&
        radius2 = 1.25_rk)
    sleeve = extrude(&
        surface = profile,&
        vector  = [0.45_rk,-0.25_rk,3.5_rk])

    write(*,"(a,l1)") "Valid extrusion        : ", sleeve%err%ok
    write(*,"(a,l1)") "Rational geometry      : ", sleeve%is_rational()
    write(*,"(a,3(i0,1x))") "Degrees                 : ", sleeve%get_degree()
    write(*,"(a,3(i0,1x))") "Control points          : ", sleeve%get_nc()

    call sleeve%create(81, 17, 41)
    call sleeve%export_Xc("vtk/volume_extrude_Xc.vtk")
    call sleeve%export_Xg("vtk/volume_extrude_Xg.vtk")
    call sleeve%export_Xth_in_Xg("vtk/volume_extrude_Xth.vtk", 13)
    call sleeve%show(&
        vtkfile_Xc        = "vtk/volume_extrude_Xc.vtk",&
        vtkfile_Xg        = "vtk/volume_extrude_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_extrude_Xth.vtk")

    call profile%finalize()
    call sleeve%finalize()

end program volume_extrude
