!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Geometry-only warped perforated volume plate with two through openings.
program volume_multipatch_perforated

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: nx = 10, ny = 6, ncp = 3
    real(rk), parameter :: pi = acos(-1.0_rk)

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: plate
    integer :: id(nx, ny), ix, iy
    integer, allocatable :: dof_map(:), elem(:,:)

    id = 0
    do iy = 1, ny
        do ix = 1, nx
            if (.not. is_opening(ix, iy)) then
                call set_plate_patch(patch, ix, iy)
                call plate%add_patch(patch, id(ix, iy))
            end if
        end do
    end do

    do iy = 1, ny
        do ix = 1, nx - 1
            if (id(ix, iy) > 0 .and. id(ix + 1, iy) > 0) then
              call plate%connect(id(ix, iy), "u_max", id(ix + 1, iy), "u_min", 0)
            end if
        end do
    end do
    do iy = 1, ny - 1
        do ix = 1, nx
            if (id(ix, iy) > 0 .and. id(ix, iy + 1) > 0) then
              call plate%connect(id(ix, iy), "v_max", id(ix, iy + 1), "v_min", 0)
            end if
        end do
    end do

    dof_map = plate%cmp_dof_map()
    elem = plate%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "perforated patches    :", plate%get_npatch()
    write(*,"(a,*(1x,g0))") "perforated connections:", plate%get_nconnection()
    write(*,"(a,*(1x,g0))") "perforated valid      :", plate%is_valid()
    write(*,"(a,*(1x,g0))") "shared global dofs    :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements              :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element     :", size(elem, 2)

    call plate%create(5, 5, 5)
    call plate%export_Xc("vtk/volume_multipatch_perforated")
    call plate%export_Xg("vtk/volume_multipatch_perforated")
    call plate%export_Xth_in_Xg("vtk/volume_multipatch_perforated", 5)
    call plate%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_perforated_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_perforated_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_perforated_Xth_*.vtk")
    call plate%finalize()

contains

    pure logical function is_opening(ix, iy) result(opening)
        integer, intent(in) :: ix, iy
        real(rk) :: x, y, h1, h2, slot

        x = -5.0_rk + 10.0_rk*(real(ix, rk) - 0.5_rk)/real(nx, rk)
        y = -3.0_rk + 6.0_rk*(real(iy, rk) - 0.5_rk)/real(ny, rk)
        h1 = ((x + 2.0_rk)/0.95_rk)**2 + (y/1.15_rk)**2
        h2 = ((x - 2.0_rk)/0.95_rk)**2 + (y/1.15_rk)**2
        slot = (x/0.55_rk)**2 + ((abs(y) - 1.55_rk)/0.80_rk)**2
        opening = h1 < 1.0_rk .or. h2 < 1.0_rk .or. slot < 1.0_rk
    end function
    !===============================================================================

    pure subroutine set_plate_patch(p, ix, iy)
        type(nurbs_volume), intent(inout) :: p
        integer, intent(in) :: ix, iy
        real(rk) :: Xc(ncp*ncp*ncp, 3)
        real(rk) :: u, v, w, x, y, z, thickness, relief
        integer :: a, b, c, n

        do c = 1, ncp
            w = real(c - 1, rk)/real(ncp - 1, rk)
            do b = 1, ncp
                v = real(b - 1, rk)/real(ncp - 1, rk)
                y = -3.0_rk + 6.0_rk*(real(iy - 1, rk) + v)/real(ny, rk)
                do a = 1, ncp
                    u = real(a - 1, rk)/real(ncp - 1, rk)
                    x = -5.0_rk + 10.0_rk*(real(ix - 1, rk) + u)/real(nx, rk)
                    thickness = 0.34_rk*(1.0_rk + 0.08_rk*cos(pi*x/5.0_rk)*cos(pi*y/3.0_rk))
                    relief = 0.10_rk*cos(2.0_rk*pi*x/5.0_rk)*cos(pi*y/1.5_rk)
                    z = thickness*(w - 0.5_rk) + relief*(2.0_rk*w - 1.0_rk)
                    n = a + (b - 1)*ncp + (c - 1)*ncp*ncp
                    Xc(n, 1) = x + 0.08_rk*sin(pi*y/3.0_rk)*sin(pi*x/5.0_rk)
                    Xc(n, 2) = y + 0.06_rk*sin(2.0_rk*pi*x/5.0_rk)
                    Xc(n, 3) = z
                end do
            end do
        end do
        call p%set([ncp, ncp, ncp], Xc)
    end subroutine
    !===============================================================================

end program volume_multipatch_perforated
