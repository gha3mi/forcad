!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Sweep an exact ring surface by translation along a cubic NURBS spine.
program volume_sweep

    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume, sweep

    implicit none

    real(rk), parameter :: knot(11) = [0.0_rk,0.0_rk,0.0_rk,0.0_rk, &
        0.18_rk,0.46_rk,0.76_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk]
    real(rk), parameter :: Xc(7,3) = reshape([0.0_rk,0.5_rk,1.3_rk,2.2_rk,3.1_rk,3.9_rk,4.6_rk, &
        0.0_rk,0.4_rk,-0.3_rk,0.7_rk,-0.5_rk,0.3_rk,0.0_rk, &
        0.0_rk,0.8_rk,1.5_rk,2.4_rk,3.1_rk,4.0_rk,4.8_rk], [7,3])

    type(nurbs_curve) :: spine
    type(nurbs_surface) :: profile
    type(nurbs_volume) :: duct

    call profile%set_ring(&
        center  = [0.0_rk,0.0_rk,0.0_rk],&
        radius1 = 0.24_rk,&
        radius2 = 0.52_rk)
    call spine%set(&
        knot   = knot,&
        Xc     = Xc,&
        degree = 3)
    duct = sweep(&
        profile = profile,&
        spine   = spine,&
        origin  = [0.0_rk,0.0_rk,0.0_rk])

    write(*,"(a,l1)") "Valid sweep            : ", duct%err%ok
    write(*,"(a,l1)") "Rational geometry      : ", duct%is_rational()
    write(*,"(a,3(i0,1x))") "Degrees                 : ", duct%get_degree()
    write(*,"(a,3(i0,1x))") "Control points          : ", duct%get_nc()

    call duct%create(61, 13, 71)
    call duct%export_Xc("vtk/volume_sweep_Xc.vtk")
    call duct%export_Xg("vtk/volume_sweep_Xg.vtk")
    call duct%export_Xth_in_Xg("vtk/volume_sweep_Xth.vtk", 13)
    call duct%show(&
        vtkfile_Xc        = "vtk/volume_sweep_Xc.vtk",&
        vtkfile_Xg        = "vtk/volume_sweep_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_sweep_Xth.vtk")

    call profile%finalize()
    call spine%finalize()
    call duct%finalize()

end program volume_sweep
