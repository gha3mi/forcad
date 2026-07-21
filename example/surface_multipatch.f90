!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates the public multi-patch API on two connected NURBS surfaces.
program surface_multipatch

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface

    implicit none

    type(nurbs_surface) :: left, right
    type(nurbs_multipatch_surface) :: patches
    integer, allocatable :: dof_map(:), elem(:,:), elem_patch(:)
    integer :: left_id, right_id

    call left%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
    call right%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
    call right%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

    call patches%add_patch(left, left_id)
    call patches%add_patch(right, right_id)
    call patches%connect(left_id, "u_max", right_id, "u_min", 0)

    dof_map = patches%cmp_dof_map()
    elem = patches%cmp_elem(shared=.true.)
    elem_patch = patches%cmp_elem_patch()

    write(*,"(a,*(1x,g0))") "number of patches      :", patches%get_npatch()
    write(*,"(a,*(1x,g0))") "number of connections  :", patches%get_nconnection()
    write(*,"(a,*(1x,g0))") "connections are valid  :", patches%is_valid()
    write(*,"(a,*(1x,g0))") "local control dofs     :", size(dof_map)
    write(*,"(a,*(1x,g0))") "shared global dofs     :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements               :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element      :", size(elem, 2)
    write(*,"(a,*(1x,g0))") "element patch ids      :", elem_patch
    call patches%create(41, 41)
    call patches%export_Xc("vtk/surface_multipatch")
    call patches%export_Xg("vtk/surface_multipatch")
    call patches%export_Xth_in_Xg("vtk/surface_multipatch", 41)
    call patches%show(&
        vtkfile_Xc        = "vtk/surface_multipatch_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_multipatch_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_multipatch_Xth_*.vtk")

end program surface_multipatch
