!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Geometry-only serpentine hollow manifold built from curved connected NURBS volume patches.
program volume_multipatch_manifold

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: ns = 10, nr = 2, nt = 8, ncp = 3
    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: length = 6.50_rk
    real(rk), parameter :: radius_inner = 0.28_rk
    real(rk), parameter :: radius_outer = 0.52_rk

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: manifold
    integer :: id(ns, nr, nt), is, ir, it
    integer, allocatable :: dof_map(:), elem(:,:)

    do it = 1, nt
        do ir = 1, nr
            do is = 1, ns
                call set_manifold_patch(patch, is, ir, it)
                call manifold%add_patch(patch, id(is, ir, it))
            end do
        end do
    end do

    do it = 1, nt
        do ir = 1, nr
            do is = 1, ns - 1
                call manifold%connect(id(is, ir, it), "u_max", id(is + 1, ir, it), "u_min", 0)
            end do
        end do
    end do
    do it = 1, nt
        do ir = 1, nr - 1
            do is = 1, ns
                call manifold%connect(id(is, ir, it), "v_max", id(is, ir + 1, it), "v_min", 0)
            end do
        end do
    end do
    do it = 1, nt - 1
        do ir = 1, nr
            do is = 1, ns
                call manifold%connect(id(is, ir, it), "w_max", id(is, ir, it + 1), "w_min", 0)
            end do
        end do
    end do
    do ir = 1, nr
        do is = 1, ns
            call manifold%connect(id(is, ir, nt), "w_max", id(is, ir, 1), "w_min", 0)
        end do
    end do

    dof_map = manifold%cmp_dof_map()
    elem = manifold%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "manifold patches      :", manifold%get_npatch()
    write(*,"(a,*(1x,g0))") "manifold connections  :", manifold%get_nconnection()
    write(*,"(a,*(1x,g0))") "manifold valid        :", manifold%is_valid()
    write(*,"(a,*(1x,g0))") "shared global dofs    :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements              :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element     :", size(elem, 2)

    call manifold%create(6, 5, 8)
    call manifold%export_Xc("vtk/volume_multipatch_manifold")
    call manifold%export_Xg("vtk/volume_multipatch_manifold")
    call manifold%export_Xth_in_Xg("vtk/volume_multipatch_manifold", 8)
    call manifold%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_manifold_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_manifold_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_manifold_Xth_*.vtk")
    call manifold%finalize()

contains

    pure subroutine set_manifold_patch(p, is, ir, it)
        type(nurbs_volume), intent(inout) :: p
        integer, intent(in) :: is, ir, it
        real(rk) :: Xc(ncp*ncp*ncp, 3)
        real(rk) :: u, v, w, s, q, phi, radius, cx, cy, cz, tx, ty, mag, nx, ny, bx, by, bz
        real(rk) :: aellipse, bellipse, twist, cph, sph
        integer :: a, b, c, n

        do c = 1, ncp
            w = real(c - 1, rk)/real(ncp - 1, rk)
            phi = 2.0_rk*pi*(real(it - 1, rk) + w)/real(nt, rk)
            do b = 1, ncp
                v = real(b - 1, rk)/real(ncp - 1, rk)
                q = (real(ir - 1, rk) + v)/real(nr, rk)
                radius = radius_inner + (radius_outer - radius_inner)*q
                do a = 1, ncp
                    u = real(a - 1, rk)/real(ncp - 1, rk)
                    s = (real(is - 1, rk) + u)/real(ns, rk)
                    cx = length*(s - 0.5_rk)
                    cy = 0.80_rk*sin(2.0_rk*pi*s) + 0.22_rk*sin(6.0_rk*pi*s)
                    cz = 0.38_rk*sin(pi*s) + 0.18_rk*sin(4.0_rk*pi*s)
                    tx = length
                    ty = 1.60_rk*pi*cos(2.0_rk*pi*s) + 1.32_rk*pi*cos(6.0_rk*pi*s)
                    mag = sqrt(tx*tx + ty*ty)
                    nx = -ty/mag
                    ny = tx/mag
                    bx = 0.0_rk
                    by = 0.0_rk
                    bz = 1.0_rk
                    twist = phi + 0.90_rk*pi*s
                    cph = cos(twist)
                    sph = sin(twist)
                    aellipse = 1.0_rk + 0.18_rk*sin(2.0_rk*pi*s)
                    bellipse = 0.78_rk + 0.14_rk*cos(2.0_rk*pi*s)
                    n = a + (b - 1)*ncp + (c - 1)*ncp*ncp
                    Xc(n, 1) = cx + radius*aellipse*cph*nx + radius*bellipse*sph*bx
                    Xc(n, 2) = cy + radius*aellipse*cph*ny + radius*bellipse*sph*by
                    Xc(n, 3) = cz + radius*bellipse*sph*bz
                end do
            end do
        end do
        call p%set([ncp, ncp, ncp], Xc)
    end subroutine
    !===============================================================================

end program volume_multipatch_manifold
