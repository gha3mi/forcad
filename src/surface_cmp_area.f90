!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program surface_cmp_area

    use forcad, only: rk, nurbs_surface

    implicit none

    type(nurbs_surface) :: shape
    real(rk) :: area
    real(rk) :: Xc(4,3)

    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [2.0_rk, 0.0_rk, 0.0_rk]
    Xc(3,:) = [0.0_rk, 2.0_rk, 0.0_rk]
    Xc(4,:) = [2.0_rk, 2.0_rk, 0.0_rk]

    call shape%set(&
        knot1=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
        knot2=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
        Xc=Xc)

    call shape%cmp_area(area)
    write(*,"(a,1x,es16.8)") "area:", area

    call shape%create(21, 21)
    call shape%export_Xc("vtk/surface_cmp_area_Xc.vtk")
    call shape%export_Xg("vtk/surface_cmp_area_Xg.vtk")
    call shape%export_Xth_in_Xg("vtk/surface_cmp_area_Xth.vtk")
    call shape%show("vtk/surface_cmp_area_Xc.vtk", "vtk/surface_cmp_area_Xg.vtk", "vtk/surface_cmp_area_Xth.vtk")
end program surface_cmp_area
