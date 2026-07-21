!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates a multi-patch arched bridge with deck, side arches, and ribs.
program surface_multipatch_bridge

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface

    implicit none

    integer, parameter :: nx = 8, ny = 3, nh = 2, nr = 7, ncp = 4
    real(rk), parameter :: pi = 3.1415926535897932384626433832795_rk
    real(rk), parameter :: length = 9.0_rk, width = 2.4_rk

    type(nurbs_surface) :: patch
    type(nurbs_multipatch_surface) :: patches
    real(rk), allocatable :: Xc(:,:)
    integer, allocatable :: dof_map(:), elem(:,:)
    integer :: deck_id(nx,ny), arch_id(2,nx,nh), rib_id(nr)
    integer :: i, j, k, n, side
    real(rk) :: u, v, q, t, w, x, y, z, side_sign, deck_z, arch_z, arch_shape

    do j = 1, ny
        do i = 1, nx
            call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
            Xc = patch%get_Xc()

            do n = 1, size(Xc, 1)
                u = Xc(n,1)
                v = Xc(n,2)
                t = (real(i - 1, rk) + u)/real(nx, rk)
                w = (real(j - 1, rk) + v)/real(ny, rk) - 0.5_rk
                x = length*(t - 0.5_rk)
                y = width*w
                z = 0.22_rk*sin(pi*t) - 0.13_rk*(2.0_rk*w)**2 + &
                    0.04_rk*sin(4.0_rk*pi*t)*cos(2.0_rk*pi*w)
                Xc(n,:) = [x, y, z]
            end do

            call patch%set([ncp, ncp], Xc)
            call patches%add_patch(patch, deck_id(i,j))
        end do
    end do

    do side = 1, 2
        side_sign = merge(-1.0_rk, 1.0_rk, side == 1)
        do k = 1, nh
            do i = 1, nx
                call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
                Xc = patch%get_Xc()

                do n = 1, size(Xc, 1)
                    u = Xc(n,1)
                    v = Xc(n,2)
                    t = (real(i - 1, rk) + u)/real(nx, rk)
                    q = (real(k - 1, rk) + v)/real(nh, rk)
                    w = 0.5_rk*side_sign
                    deck_z = 0.22_rk*sin(pi*t) - 0.13_rk*(2.0_rk*w)**2 + &
                        0.04_rk*sin(4.0_rk*pi*t)*cos(2.0_rk*pi*w)
                    arch_shape = sin(pi*t)
                    arch_z = deck_z + 1.05_rk + 2.15_rk*arch_shape**0.82_rk + &
                        0.12_rk*sin(5.0_rk*pi*t)**2
                    x = length*(t - 0.5_rk) + 0.05_rk*q*sin(2.0_rk*pi*t)
                    y = side_sign*(0.5_rk*width + 0.12_rk*q + 0.20_rk*q*arch_shape)
                    z = deck_z + q*(arch_z - deck_z)
                    Xc(n,:) = [x, y, z]
                end do

                call patch%set([ncp, ncp], Xc)
                call patches%add_patch(patch, arch_id(side,i,k))
            end do
        end do
    end do

    do i = 1, nr
        call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
        Xc = patch%get_Xc()

        do n = 1, size(Xc, 1)
            u = Xc(n,1)
            q = Xc(n,2)
            t = real(i, rk)/real(nr + 1, rk)
            w = u - 0.5_rk
            deck_z = 0.22_rk*sin(pi*t) - 0.13_rk*(2.0_rk*w)**2 + &
                0.04_rk*sin(4.0_rk*pi*t)*cos(2.0_rk*pi*w)
            arch_shape = max(0.0_rk, 1.0_rk - (2.0_rk*w)**2)
            x = length*(t - 0.5_rk) + 0.06_rk*(q - 0.5_rk)
            y = 1.22_rk*width*w
            z = deck_z + 0.12_rk + q*(0.55_rk + 1.10_rk*arch_shape)*(0.55_rk + 0.45_rk*sin(pi*t))
            Xc(n,:) = [x, y, z]
        end do

        call patch%set([ncp, ncp], Xc)
        call patches%add_patch(patch, rib_id(i))
    end do

    do j = 1, ny
        do i = 1, nx - 1
            call patches%connect(deck_id(i,j), "u_max", deck_id(i+1,j), "u_min", 0)
        end do
    end do
    do j = 1, ny - 1
        do i = 1, nx
            call patches%connect(deck_id(i,j), "v_max", deck_id(i,j+1), "v_min", 0)
        end do
    end do

    do side = 1, 2
        do k = 1, nh
            do i = 1, nx - 1
                call patches%connect(arch_id(side,i,k), "u_max", arch_id(side,i+1,k), "u_min", 0)
            end do
        end do
        do k = 1, nh - 1
            do i = 1, nx
                call patches%connect(arch_id(side,i,k), "v_max", arch_id(side,i,k+1), "v_min", 0)
            end do
        end do
        do i = 1, nx
            if (side == 1) then
                call patches%connect(deck_id(i,1), "v_min", arch_id(side,i,1), "v_min", 0)
            else
                call patches%connect(deck_id(i,ny), "v_max", arch_id(side,i,1), "v_min", 0)
            end if
        end do
    end do

    dof_map = patches%cmp_dof_map()
    elem = patches%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "bridge deck patches      :", nx*ny
    write(*,"(a,*(1x,g0))") "bridge arch patches      :", 2*nx*nh
    write(*,"(a,*(1x,g0))") "bridge rib patches       :", size(rib_id)
    write(*,"(a,*(1x,g0))") "bridge total patches     :", patches%get_npatch()
    write(*,"(a,*(1x,g0))") "bridge connections       :", patches%get_nconnection()
    write(*,"(a,*(1x,g0))") "bridge valid             :", patches%is_valid()
    write(*,"(a,*(1x,g0))") "local control dofs       :", size(dof_map)
    write(*,"(a,*(1x,g0))") "shared global dofs       :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements                 :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element        :", size(elem, 2)

    call patches%create(37, 31)
    call patches%export_Xc("vtk/surface_multipatch_bridge")
    call patches%export_Xg("vtk/surface_multipatch_bridge")
    call patches%export_Xth_in_Xg("vtk/surface_multipatch_bridge", 37)
    call patches%show(&
        vtkfile_Xc        = "vtk/surface_multipatch_bridge_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_multipatch_bridge_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_multipatch_bridge_Xth_*.vtk")

end program surface_multipatch_bridge
