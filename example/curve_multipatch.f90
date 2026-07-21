!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Construct two connected curve patches and inspect their shared topology.
program curve_multipatch

    use forcad, only: rk, nurbs_curve, nurbs_multipatch_curve

    implicit none

    type(nurbs_curve) :: left, right
    type(nurbs_multipatch_curve) :: patches
    real(rk) :: Xc(3,3)
    integer, allocatable :: dof_map(:), elem(:,:)
    integer :: left_id, right_id

    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [0.5_rk, 0.8_rk, 0.2_rk]
    Xc(3,:) = [1.0_rk, 0.0_rk, 0.4_rk]
    call left%set(Xc = Xc)

    Xc(1,:) = [1.0_rk, 0.0_rk, 0.4_rk]
    Xc(2,:) = [1.5_rk,-0.8_rk, 0.2_rk]
    Xc(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
    call right%set(Xc = Xc)

    call patches%add_patch(left, left_id)
    call patches%add_patch(right, right_id)
    call patches%connect(&
        patch_a    = left_id,&
        side_a     = "right",&
        patch_b    = right_id,&
        side_b     = "left",&
        continuity = 0)

    dof_map = patches%cmp_dof_map()
    elem = patches%cmp_elem(shared = .true.)

    write(*,"(a,i0)") "Patches           : ", patches%get_npatch()
    write(*,"(a,i0)") "Connections       : ", patches%get_nconnection()
    write(*,"(a,l1)") "Topology valid    : ", patches%is_valid()
    write(*,"(a,i0)") "Shared control DOF: ", maxval(dof_map)
    write(*,"(a,i0)") "Elements          : ", size(elem, 1)

    call patches%create(res = 81)
    call patches%export_Xc("vtk/curve_multipatch")
    call patches%export_Xg("vtk/curve_multipatch")
    call patches%export_Xth_in_Xg("vtk/curve_multipatch")
    call patches%show(&
        vtkfile_Xc        = "vtk/curve_multipatch_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/curve_multipatch_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/curve_multipatch_Xth_*.vtk")

end program curve_multipatch
