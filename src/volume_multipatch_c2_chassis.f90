!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Exact C2 multi-patch volume chassis with a sparse high-continuity topology.
program volume_multipatch_c2_chassis

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: nx = 10, ny = 6, nz = 5, ncp = 3
    real(rk), parameter :: dx = 0.46_rk, dy = 0.42_rk, dz = 0.34_rk

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: chassis
    logical :: active(nx, ny, nz)
    real(rk), allocatable :: Xc(:,:), cval(:)
    integer, allocatable :: dof_map(:), elem(:,:), rowptr(:), col(:)
    integer :: id(nx, ny, nz), i, j, k, iu, iv, iw, p
    real(rk) :: u, v, w, gx, gy, gz

    id = 0
    active = .false.
    allocate(Xc(ncp*ncp*ncp, 3))

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                active(i,j,k) = &
                    k == 1 .or. &
                    ((j == 1 .or. j == ny) .and. k <= 4) .or. &
                    ((i == 1 .or. i == nx) .and. k <= 4) .or. &
                    (k == 3 .and. (j == 2 .or. j == ny - 1) .and. i >= 3 .and. i <= nx - 2) .or. &
                    (k == 4 .and. i >= 3 .and. i <= nx - 2 .and. j >= 2 .and. j <= ny - 1) .or. &
                    (i >= 4 .and. i <= 7 .and. j >= 3 .and. j <= 4 .and. k <= nz)
            end do
        end do
    end do

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                if (.not. active(i,j,k)) cycle
                do iw = 0, ncp - 1
                    w = real(iw, rk)/real(ncp - 1, rk)
                    do iv = 0, ncp - 1
                        v = real(iv, rk)/real(ncp - 1, rk)
                        do iu = 0, ncp - 1
                            u = real(iu, rk)/real(ncp - 1, rk)
                            gx = real(i - 1, rk) + u
                            gy = real(j - 1, rk) + v
                            gz = real(k - 1, rk) + w
                            p = iu + ncp*(iv + ncp*iw) + 1
                            Xc(p,1) = dx*(gx + 0.18_rk*gy)
                            Xc(p,2) = dy*(gy + 0.10_rk*gz)
                            Xc(p,3) = dz*(gz + 0.08_rk*gx)
                        end do
                    end do
                end do
                call patch%set([ncp, ncp, ncp], Xc)
                call chassis%add_patch(patch, id(i,j,k))
            end do
        end do
    end do

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx - 1
                if (id(i,j,k) > 0 .and. id(i + 1,j,k) > 0) then
                  call chassis%connect(id(i,j,k), "u_max", id(i + 1,j,k), "u_min", 2)
                end if
            end do
        end do
    end do
    do k = 1, nz
        do j = 1, ny - 1
            do i = 1, nx
                if (id(i,j,k) > 0 .and. id(i,j + 1,k) > 0) then
                  call chassis%connect(id(i,j,k), "v_max", id(i,j + 1,k), "v_min", 2)
                end if
            end do
        end do
    end do
    do k = 1, nz - 1
        do j = 1, ny
            do i = 1, nx
                if (id(i,j,k) > 0 .and. id(i,j,k + 1) > 0) then
                  call chassis%connect(id(i,j,k), "w_max", id(i,j,k + 1), "w_min", 2)
                end if
            end do
        end do
    end do

    dof_map = chassis%cmp_dof_map()
    elem = chassis%cmp_elem(shared=.true.)
    call chassis%cmp_dof_constraint(rowptr, col, cval)

    write(*,"(a,*(1x,g0))") "C2 chassis patches       :", chassis%get_npatch()
    write(*,"(a,*(1x,g0))") "C2 chassis connections   :", chassis%get_nconnection()
    write(*,"(a,*(1x,g0))") "C2 chassis valid         :", chassis%is_valid()
    write(*,"(a,*(1x,g0))") "local control dofs       :", size(dof_map)
    write(*,"(a,*(1x,g0))") "shared global dofs       :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "Cn constraint rows       :", size(rowptr) - 1
    write(*,"(a,*(1x,g0))") "Cn constraint nnz        :", size(cval)
    write(*,"(a,*(1x,g0))") "elements                 :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element        :", size(elem, 2)

    call chassis%create(5, 5, 5)
    call chassis%export_Xc("vtk/volume_multipatch_c2_chassis")
    call chassis%export_Xg("vtk/volume_multipatch_c2_chassis")
    call chassis%export_Xth_in_Xg("vtk/volume_multipatch_c2_chassis", 5)
    call chassis%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_c2_chassis_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_c2_chassis_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_c2_chassis_Xth_*.vtk")
    call chassis%finalize()

end program volume_multipatch_c2_chassis
