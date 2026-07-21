!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Geometry-only multi-patch gearbox housing with perimeter shell, central opening, bosses, and ribs.
program volume_multipatch_housing

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: nx = 8, ny = 6, nz = 4
    real(rk), parameter :: dx = 0.42_rk, dy = 0.42_rk, dz = 0.36_rk

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: housing
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
                    (i == 1 .or. i == nx .or. j == 1 .or. j == ny) .or. &
                    (k == nz .and. (i <= 2 .or. i >= nx - 1 .or. j <= 2 .or. j >= ny - 1)) .or. &
                    ((i <= 2 .or. i >= nx - 1) .and. (j <= 2 .or. j >= ny - 1)) .or. &
                    ((i == 4 .or. i == 5) .and. (j == 2 .or. j == ny - 1) .and. k <= 3) .or. &
                    ((j == 3 .or. j == 4) .and. (i == 2 .or. i == nx - 1) .and. k <= 3)
            end do
        end do
    end do

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                if (.not. active(i, j, k)) cycle
                call patch%set_hexahedron([dx, dy, dz], [2, 2, 2])
                call patch%translate_Xc([real(i - 1, rk)*dx, real(j - 1, rk)*dy, real(k - 1, rk)*dz])
                call housing%add_patch(patch, id(i, j, k))
            end do
        end do
    end do

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx - 1
                if (id(i, j, k) > 0 .and. id(i + 1, j, k) > 0) then
                  call housing%connect(id(i, j, k), "u_max", id(i + 1, j, k), "u_min", 0)
                end if
            end do
        end do
    end do
    do k = 1, nz
        do j = 1, ny - 1
            do i = 1, nx
                if (id(i, j, k) > 0 .and. id(i, j + 1, k) > 0) then
                  call housing%connect(id(i, j, k), "v_max", id(i, j + 1, k), "v_min", 0)
                end if
            end do
        end do
    end do
    do k = 1, nz - 1
        do j = 1, ny
            do i = 1, nx
                if (id(i, j, k) > 0 .and. id(i, j, k + 1) > 0) then
                  call housing%connect(id(i, j, k), "w_max", id(i, j, k + 1), "w_min", 0)
                end if
            end do
        end do
    end do

    dof_map = housing%cmp_dof_map()
    elem = housing%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "housing patches       :", housing%get_npatch()
    write(*,"(a,*(1x,g0))") "housing connections   :", housing%get_nconnection()
    write(*,"(a,*(1x,g0))") "housing valid         :", housing%is_valid()
    write(*,"(a,*(1x,g0))") "shared global dofs    :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements              :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element     :", size(elem, 2)

    call housing%create(5, 5, 5)
    call housing%export_Xc("vtk/volume_multipatch_housing")
    call housing%export_Xg("vtk/volume_multipatch_housing")
    call housing%export_Xth_in_Xg("vtk/volume_multipatch_housing", 5)
    call housing%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_housing_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_housing_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_housing_Xth_*.vtk")
    call housing%finalize()

end program volume_multipatch_housing
