!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> L-bracket assembled from connected NURBS volume patches.
program volume_multipatch_bracket

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: bracket
    integer, allocatable :: dof_map(:), elem(:,:)
    integer, parameter :: nx = 5, ny = 5, nz = 2
    integer :: id(nx, ny, nz), i, j, k
    logical :: active(nx, ny, nz)
    real(rk), parameter :: dx = 1.0_rk, dy = 1.0_rk, dz = 0.55_rk

    id = 0
    active = .false.
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                active(i, j, k) = (j <= 2) .or. (i <= 2)
            end do
        end do
    end do

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                if (.not. active(i, j, k)) cycle
                call patch%set_hexahedron([dx, dy, dz], [2, 2, 2])
                call patch%translate_Xc([real(i - 1, rk)*dx, real(j - 1, rk)*dy, real(k - 1, rk)*dz])
                call bracket%add_patch(patch, id(i, j, k))
            end do
        end do
    end do

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx - 1
                if (id(i, j, k) > 0 .and. id(i + 1, j, k) > 0) then
                    call bracket%connect(id(i, j, k), "u_max", id(i + 1, j, k), "u_min", 0)
                end if
            end do
        end do
    end do
    do k = 1, nz
        do j = 1, ny - 1
            do i = 1, nx
                if (id(i, j, k) > 0 .and. id(i, j + 1, k) > 0) then
                    call bracket%connect(id(i, j, k), "v_max", id(i, j + 1, k), "v_min", 0)
                end if
            end do
        end do
    end do
    do k = 1, nz - 1
        do j = 1, ny
            do i = 1, nx
                if (id(i, j, k) > 0 .and. id(i, j, k + 1) > 0) then
                    call bracket%connect(id(i, j, k), "w_max", id(i, j, k + 1), "w_min", 0)
                end if
            end do
        end do
    end do

    dof_map = bracket%cmp_dof_map()
    elem = bracket%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "bracket patches        :", bracket%get_npatch()
    write(*,"(a,*(1x,g0))") "bracket connections    :", bracket%get_nconnection()
    write(*,"(a,*(1x,g0))") "bracket valid          :", bracket%is_valid()
    write(*,"(a,*(1x,g0))") "shared control nodes   :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "shared elements        :", size(elem, 1)
    call bracket%create(5, 5, 5)
    call bracket%export_Xc("vtk/volume_multipatch_bracket_geometry")
    call bracket%export_Xg("vtk/volume_multipatch_bracket_geometry")
    call bracket%export_Xth_in_Xg("vtk/volume_multipatch_bracket_geometry", 5)
    call bracket%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_bracket_geometry_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_bracket_geometry_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_bracket_geometry_Xth_*.vtk")
    call bracket%finalize()

end program volume_multipatch_bracket
