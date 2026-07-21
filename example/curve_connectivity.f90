!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Compare curve IGA, control, geometry, and knot-span connectivities.
program curve_connectivity

    use forcad, only: rk, nurbs_curve

    implicit none

    real(rk), parameter :: center(3) = [0.0_rk, 0.0_rk, 0.0_rk]

    type(nurbs_curve) :: curve
    integer, allocatable :: elem(:,:), elem_Xc(:,:), elem_Xg(:,:), elem_Xth(:,:)

    call curve%set_circle(&
        center = center,&
        radius = 2.0_rk)
    call curve%create(res = 33)

    elem = curve%cmp_elem()
    elem_Xc = curve%cmp_elem_Xc_vis()
    elem_Xg = curve%cmp_elem_Xg_vis()
    elem_Xth = curve%cmp_elem_Xth()

    write(*,"(a,l1)") "Valid curve       : ", curve%err%ok
    write(*,"(a,2(i0,1x))") "IGA connectivity   : ", size(elem,1), size(elem,2)
    write(*,"(a,2(i0,1x))") "Control mesh       : ", size(elem_Xc,1), size(elem_Xc,2)
    write(*,"(a,2(i0,1x))") "Geometry mesh      : ", size(elem_Xg,1), size(elem_Xg,2)
    write(*,"(a,2(i0,1x))") "Knot-span mesh     : ", size(elem_Xth,1), size(elem_Xth,2)
    write(*,"(a,*(i0,1x))") "First IGA element  : ", elem(1,:)

    call curve%export_Xc("vtk/curve_connectivity_Xc.vtk")
    call curve%export_Xg("vtk/curve_connectivity_Xg.vtk")
    call curve%export_Xth_in_Xg("vtk/curve_connectivity_Xth.vtk")
    call curve%show(&
        vtkfile_Xc        = "vtk/curve_connectivity_Xc.vtk",&
        vtkfile_Xg        = "vtk/curve_connectivity_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/curve_connectivity_Xth.vtk")

    call curve%finalize()

end program curve_connectivity
