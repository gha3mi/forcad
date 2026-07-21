!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Geometry-only multi-patch clevis mount with base, rear flange, side rails, ribs, and lugs.
program volume_multipatch_mount

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: nx = 9, ny = 5, nz = 5
    real(rk), parameter :: dx = 0.55_rk, dy = 0.42_rk, dz = 0.32_rk

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: mount
    logical :: active(nx, ny, nz)
    integer :: id(nx, ny, nz), i, j, k
    integer, allocatable :: dof_map(:), elem(:,:)

    id = 0
    active = .false.
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                active(i, j, k) = &
                    k == 1 .or. &
                    (i <= 2 .and. k <= 4) .or. &
                    (i >= 3 .and. i <= 8 .and. (j == 1 .or. j == ny) .and. k <= 2) .or. &
                    (i >= 7 .and. (j == 2 .or. j == 4) .and. k <= 5) .or. &
                    (i >= 3 .and. i <= 7 .and. (j == 2 .or. j == 4) .and. k == 3) .or. &
                    (i == 5 .and. j == 3 .and. k <= 3)
            end do
        end do
    end do

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                if (.not. active(i, j, k)) cycle
                call patch%set_hexahedron([dx, dy, dz], [2, 2, 2])
                call patch%translate_Xc([real(i - 1, rk)*dx, real(j - 1, rk)*dy, real(k - 1, rk)*dz])
                call mount%add_patch(patch, id(i, j, k))
            end do
        end do
    end do

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx - 1
                if (id(i, j, k) > 0 .and. id(i + 1, j, k) > 0) then
                  call mount%connect(id(i, j, k), "u_max", id(i + 1, j, k), "u_min", 0)
                end if
            end do
        end do
    end do
    do k = 1, nz
        do j = 1, ny - 1
            do i = 1, nx
                if (id(i, j, k) > 0 .and. id(i, j + 1, k) > 0) then
                  call mount%connect(id(i, j, k), "v_max", id(i, j + 1, k), "v_min", 0)
                end if
            end do
        end do
    end do
    do k = 1, nz - 1
        do j = 1, ny
            do i = 1, nx
                if (id(i, j, k) > 0 .and. id(i, j, k + 1) > 0) then
                  call mount%connect(id(i, j, k), "w_max", id(i, j, k + 1), "w_min", 0)
                end if
            end do
        end do
    end do

    dof_map = mount%cmp_dof_map()
    elem = mount%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "mount patches         :", mount%get_npatch()
    write(*,"(a,*(1x,g0))") "mount connections     :", mount%get_nconnection()
    write(*,"(a,*(1x,g0))") "mount valid           :", mount%is_valid()
    write(*,"(a,*(1x,g0))") "shared global dofs    :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements              :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element     :", size(elem, 2)

    call mount%create(5, 5, 5)
    call mount%export_Xc("vtk/volume_multipatch_mount")
    call mount%export_Xg("vtk/volume_multipatch_mount")
    call mount%export_Xth_in_Xg("vtk/volume_multipatch_mount", 5)
    call mount%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_mount_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_mount_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_mount_Xth_*.vtk")
    call mount%finalize()

end program volume_multipatch_mount
