!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Compare volume IGA, control, geometry, and knot-span connectivities.
program volume_connectivity

    use forcad, only: rk, nurbs_volume

    implicit none

    real(rk), parameter :: center(3) = [0.0_rk, 0.0_rk, 0.0_rk]

    type(nurbs_volume) :: volume
    integer, allocatable :: elem(:,:), elem_Xc(:,:), elem_Xg(:,:), elem_Xth(:,:)

    call volume%set_ring(&
        center  = center,&
        radius1 = 1.0_rk,&
        radius2 = 2.0_rk,&
        length  = 1.0_rk)
    call volume%create(&
        res1 = 25,&
        res2 = 9,&
        res3 = 5)

    elem = volume%cmp_elem()
    elem_Xc = volume%cmp_elem_Xc_vis()
    elem_Xg = volume%cmp_elem_Xg_vis()
    elem_Xth = volume%cmp_elem_Xth()

    write(*,"(a,l1)") "Valid volume      : ", volume%err%ok
    write(*,"(a,2(i0,1x))") "IGA connectivity   : ", size(elem,1), size(elem,2)
    write(*,"(a,2(i0,1x))") "Control mesh       : ", size(elem_Xc,1), size(elem_Xc,2)
    write(*,"(a,2(i0,1x))") "Geometry mesh      : ", size(elem_Xg,1), size(elem_Xg,2)
    write(*,"(a,2(i0,1x))") "Knot-span mesh     : ", size(elem_Xth,1), size(elem_Xth,2)
    write(*,"(a,*(i0,1x))") "First IGA element  : ", elem(1,:)

    call volume%export_Xc("vtk/volume_connectivity_Xc.vtk")
    call volume%export_Xg("vtk/volume_connectivity_Xg.vtk")
    call volume%export_Xth_in_Xg("vtk/volume_connectivity_Xth.vtk")
    call volume%show(&
        vtkfile_Xc        = "vtk/volume_connectivity_Xc.vtk",&
        vtkfile_Xg        = "vtk/volume_connectivity_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_connectivity_Xth.vtk")

    call volume%finalize()

end program volume_connectivity
