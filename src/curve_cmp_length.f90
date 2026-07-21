!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program curve_cmp_length

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: shape
    real(rk) :: length
    real(rk) :: Xc(2,3)

    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [2.0_rk, 0.0_rk, 0.0_rk]

    call shape%set(&
        knot=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
        Xc=Xc)

    call shape%cmp_length(length)
    write(*,"(a,1x,es16.8)") "length:", length

    call shape%create(21)
    call shape%export_Xc("vtk/curve_cmp_length_Xc.vtk")
    call shape%export_Xg("vtk/curve_cmp_length_Xg.vtk")
    call shape%export_Xth_in_Xg("vtk/curve_cmp_length_Xth.vtk")
    call shape%show("vtk/curve_cmp_length_Xc.vtk", "vtk/curve_cmp_length_Xg.vtk", "vtk/curve_cmp_length_Xth.vtk")
end program curve_cmp_length
