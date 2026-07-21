!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates a conforming 2-by-2 multi-patch NURBS surface grid.
program surface_multipatch_grid

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface

    implicit none

    type(nurbs_surface) :: patch
    type(nurbs_multipatch_surface) :: patches
    integer, allocatable :: dof_map(:), elem(:,:), elem_patch(:), elem_local(:)
    integer :: id(2,2), i, j

    do j = 1, 2
        do i = 1, 2
            call patch%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
            call patch%translate_Xc([real(i - 1, rk), real(j - 1, rk), 0.0_rk])
            call patches%add_patch(patch, id(i,j))
        end do
    end do

    do j = 1, 2
        call patches%connect(id(1,j), "u_max", id(2,j), "u_min", 0)
    end do
    do i = 1, 2
        call patches%connect(id(i,1), "v_max", id(i,2), "v_min", 0)
    end do

    dof_map = patches%cmp_dof_map()
    elem = patches%cmp_elem(shared=.true.)
    elem_patch = patches%cmp_elem_patch()
    elem_local = patches%cmp_elem_local()

    write(*,"(a,*(1x,g0))") "surface grid patches      :", patches%get_npatch()
    write(*,"(a,*(1x,g0))") "surface grid connections  :", patches%get_nconnection()
    write(*,"(a,*(1x,g0))") "surface grid valid        :", patches%is_valid()
    write(*,"(a,*(1x,g0))") "local control dofs        :", size(dof_map)
    write(*,"(a,*(1x,g0))") "shared global dofs        :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "expected global grid dofs :", 9
    write(*,"(a,*(1x,g0))") "elements                  :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element         :", size(elem, 2)
    write(*,"(a,*(1x,g0))") "element patch ids         :", elem_patch
    write(*,"(a,*(1x,g0))") "element local ids         :", elem_local
    call patches%create(41, 41)
    call patches%export_Xc("vtk/surface_multipatch_grid")
    call patches%export_Xg("vtk/surface_multipatch_grid")
    call patches%export_Xth_in_Xg("vtk/surface_multipatch_grid", 41)
    call patches%show(&
        vtkfile_Xc        = "vtk/surface_multipatch_grid_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_multipatch_grid_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_multipatch_grid_Xth_*.vtk")

end program surface_multipatch_grid
