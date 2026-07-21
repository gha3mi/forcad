!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Geometry-only multi-patch C-frame press body with throat, base, crown, and rear ribs.
program volume_multipatch_frame

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: nx = 8, ny = 3, nz = 8
    real(rk), parameter :: dx = 0.48_rk, dy = 0.50_rk, dz = 0.42_rk

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: frame
    logical :: active(nx, ny, nz)
    integer :: id(nx, ny, nz), i, j, k
    integer, allocatable :: dof_map(:), elem(:,:)

    id = 0
    active = .false.
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                active(i, j, k) = &
                    k <= 2 .or. &
                    k >= 7 .or. &
                    i <= 2 .or. &
                    (j == 1 .and. i <= 5 .and. k >= 3 .and. k <= 6) .or. &
                    ((j == 1 .or. j == ny) .and. i >= 3 .and. i <= 4 .and. k >= 3 .and. k <= 6) .or. &
                    (i >= 7 .and. k >= 5)
            end do
        end do
    end do

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                if (.not. active(i, j, k)) cycle
                call patch%set_hexahedron([dx, dy, dz], [2, 2, 2])
                call patch%translate_Xc([real(i - 1, rk)*dx, real(j - 1, rk)*dy, real(k - 1, rk)*dz])
                call frame%add_patch(patch, id(i, j, k))
            end do
        end do
    end do

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx - 1
                if (id(i, j, k) > 0 .and. id(i + 1, j, k) > 0) then
                  call frame%connect(id(i, j, k), "u_max", id(i + 1, j, k), "u_min", 0)
                end if
            end do
        end do
    end do
    do k = 1, nz
        do j = 1, ny - 1
            do i = 1, nx
                if (id(i, j, k) > 0 .and. id(i, j + 1, k) > 0) then
                  call frame%connect(id(i, j, k), "v_max", id(i, j + 1, k), "v_min", 0)
                end if
            end do
        end do
    end do
    do k = 1, nz - 1
        do j = 1, ny
            do i = 1, nx
                if (id(i, j, k) > 0 .and. id(i, j, k + 1) > 0) then
                  call frame%connect(id(i, j, k), "w_max", id(i, j, k + 1), "w_min", 0)
                end if
            end do
        end do
    end do

    dof_map = frame%cmp_dof_map()
    elem = frame%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "frame patches         :", frame%get_npatch()
    write(*,"(a,*(1x,g0))") "frame connections     :", frame%get_nconnection()
    write(*,"(a,*(1x,g0))") "frame valid           :", frame%is_valid()
    write(*,"(a,*(1x,g0))") "shared global dofs    :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements              :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element     :", size(elem, 2)

    call frame%create(5, 5, 5)
    call frame%export_Xc("vtk/volume_multipatch_frame")
    call frame%export_Xg("vtk/volume_multipatch_frame")
    call frame%export_Xth_in_Xg("vtk/volume_multipatch_frame", 5)
    call frame%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_frame_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_frame_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_frame_Xth_*.vtk")
    call frame%finalize()

end program volume_multipatch_frame
