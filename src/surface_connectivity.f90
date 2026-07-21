!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Compare surface IGA, control, geometry, and knot-span connectivities.
program surface_connectivity

    use forcad, only: rk, nurbs_surface

    implicit none

    real(rk), parameter :: center(3) = [0.0_rk, 0.0_rk, 0.0_rk]

    type(nurbs_surface) :: surface
    integer, allocatable :: elem(:,:), elem_Xc(:,:), elem_Xg(:,:), elem_Xth(:,:)

    call surface%set_ring(&
        center  = center,&
        radius1 = 1.0_rk,&
        radius2 = 2.0_rk)
    call surface%create(&
        res1 = 33,&
        res2 = 9)

    elem = surface%cmp_elem()
    elem_Xc = surface%cmp_elem_Xc_vis()
    elem_Xg = surface%cmp_elem_Xg_vis()
    elem_Xth = surface%cmp_elem_Xth()

    write(*,"(a,l1)") "Valid surface     : ", surface%err%ok
    write(*,"(a,2(i0,1x))") "IGA connectivity   : ", size(elem,1), size(elem,2)
    write(*,"(a,2(i0,1x))") "Control mesh       : ", size(elem_Xc,1), size(elem_Xc,2)
    write(*,"(a,2(i0,1x))") "Geometry mesh      : ", size(elem_Xg,1), size(elem_Xg,2)
    write(*,"(a,2(i0,1x))") "Knot-span mesh     : ", size(elem_Xth,1), size(elem_Xth,2)
    write(*,"(a,*(i0,1x))") "First IGA element  : ", elem(1,:)

    call surface%export_Xc("vtk/surface_connectivity_Xc.vtk")
    call surface%export_Xg("vtk/surface_connectivity_Xg.vtk")
    call surface%export_Xth_in_Xg("vtk/surface_connectivity_Xth.vtk")
    call surface%show(&
        vtkfile_Xc        = "vtk/surface_connectivity_Xc.vtk",&
        vtkfile_Xg        = "vtk/surface_connectivity_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_connectivity_Xth.vtk")

    call surface%finalize()

end program surface_connectivity
