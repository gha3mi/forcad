!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates mixed continuous and discontinuous interfaces in a surface strip.
program surface_multipatch_mixed

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface, multipatch_connection

    implicit none

    type(nurbs_surface) :: patch
    type(nurbs_multipatch_surface) :: patches
    type(multipatch_connection) :: conn
    integer, allocatable :: dof_map(:), elem(:,:), elem_patch(:)
    integer :: id(3), i

    do i = 1, 3
        call patch%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call patch%translate_Xc([real(i - 1, rk), 0.0_rk, 0.0_rk])
        call patches%add_patch(patch, id(i))
    end do

    call patches%connect(id(1), "u_max", id(2), "u_min", 1)
    call patches%connect(id(2), "u_max", id(3), "u_min", -1)

    dof_map = patches%cmp_dof_map()
    elem = patches%cmp_elem(shared=.true.)
    elem_patch = patches%cmp_elem_patch()

    write(*,"(a,*(1x,g0))") "surface strip patches      :", patches%get_npatch()
    write(*,"(a,*(1x,g0))") "surface strip connections  :", patches%get_nconnection()
    write(*,"(a,*(1x,g0))") "surface strip valid        :", patches%is_valid()
    write(*,"(a,*(1x,g0))") "local control dofs         :", size(dof_map)
    write(*,"(a,*(1x,g0))") "shared global dofs         :", maxval(dof_map)
    conn = patches%get_connection(1)
    write(*,"(a,1x,i0)") "First connection continuity :", conn%get_continuity()
    conn = patches%get_connection(2)
    write(*,"(a,1x,i0)") "Second connection continuity:", conn%get_continuity()
    write(*,"(a,*(1x,g0))") "elements                   :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element          :", size(elem, 2)
    write(*,"(a,*(1x,g0))") "element patch ids          :", elem_patch
    call patches%create(41, 41)
    call patches%export_Xc("vtk/surface_multipatch_mixed")
    call patches%export_Xg("vtk/surface_multipatch_mixed")
    call patches%export_Xth_in_Xg("vtk/surface_multipatch_mixed", 41)
    call patches%show(&
        vtkfile_Xc        = "vtk/surface_multipatch_mixed_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_multipatch_mixed_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_multipatch_mixed_Xth_*.vtk")

end program surface_multipatch_mixed
