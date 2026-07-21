!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Sweep an exact rational circle by translation along a cubic NURBS spine.
program surface_sweep

    use forcad, only: rk, nurbs_curve, nurbs_surface, sweep

    implicit none

    real(rk), parameter :: knot(11) = [0.0_rk,0.0_rk,0.0_rk,0.0_rk, &
        0.25_rk,0.50_rk,0.75_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk]
    real(rk), parameter :: Xc(7,3) = reshape([0.0_rk,0.7_rk,1.6_rk,2.5_rk,3.4_rk,4.2_rk,5.0_rk, &
        0.0_rk,0.6_rk,-0.4_rk,0.8_rk,-0.7_rk,0.5_rk,0.0_rk, &
        0.0_rk,0.5_rk,1.2_rk,1.8_rk,2.7_rk,3.3_rk,4.0_rk], [7,3])

    type(nurbs_curve) :: profile, spine
    type(nurbs_surface) :: tube

    call profile%set_circle(&
        center = [0.0_rk,0.0_rk,0.0_rk],&
        radius = 0.32_rk)
    call spine%set(&
        knot   = knot,&
        Xc     = Xc,&
        degree = 3)
    tube = sweep(&
        profile = profile,&
        spine   = spine,&
        origin  = [0.0_rk,0.0_rk,0.0_rk])

    write(*,"(a,l1)") "Valid sweep            : ", tube%err%ok
    write(*,"(a,l1)") "Rational geometry      : ", tube%is_rational()
    write(*,"(a,2(i0,1x))") "Degrees                 : ", tube%get_degree()
    write(*,"(a,2(i0,1x))") "Control points          : ", tube%get_nc()

    call tube%create(81, 61)
    call tube%export_Xc("vtk/surface_sweep_Xc.vtk")
    call tube%export_Xg("vtk/surface_sweep_Xg.vtk")
    call tube%export_Xth_in_Xg("vtk/surface_sweep_Xth.vtk", 17)
    call tube%show(&
        vtkfile_Xc        = "vtk/surface_sweep_Xc.vtk",&
        vtkfile_Xg        = "vtk/surface_sweep_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_sweep_Xth.vtk")

    call profile%finalize()
    call spine%finalize()
    call tube%finalize()

end program surface_sweep
