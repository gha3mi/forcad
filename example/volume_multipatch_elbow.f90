!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Geometry-only curved hollow elbow built from connected non-hexahedral NURBS volume patches.
program volume_multipatch_elbow

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: ns = 8, nr = 2, nt = 8, ncp = 3
    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: bend_angle = 0.70_rk*pi
    real(rk), parameter :: center_radius = 2.20_rk
    real(rk), parameter :: radius_inner = 0.34_rk
    real(rk), parameter :: radius_outer = 0.58_rk

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: elbow
    integer :: id(ns, nr, nt), is, ir, it
    integer, allocatable :: dof_map(:), elem(:,:)

    do it = 1, nt
        do ir = 1, nr
            do is = 1, ns
                call set_elbow_patch(patch, is, ir, it)
                call elbow%add_patch(patch, id(is, ir, it))
            end do
        end do
    end do

    do it = 1, nt
        do ir = 1, nr
            do is = 1, ns - 1
                call elbow%connect(id(is, ir, it), "u_max", id(is + 1, ir, it), "u_min", 0)
            end do
        end do
    end do
    do it = 1, nt
        do ir = 1, nr - 1
            do is = 1, ns
                call elbow%connect(id(is, ir, it), "v_max", id(is, ir + 1, it), "v_min", 0)
            end do
        end do
    end do
    do it = 1, nt - 1
        do ir = 1, nr
            do is = 1, ns
                call elbow%connect(id(is, ir, it), "w_max", id(is, ir, it + 1), "w_min", 0)
            end do
        end do
    end do
    do ir = 1, nr
        do is = 1, ns
            call elbow%connect(id(is, ir, nt), "w_max", id(is, ir, 1), "w_min", 0)
        end do
    end do

    dof_map = elbow%cmp_dof_map()
    elem = elbow%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "elbow patches         :", elbow%get_npatch()
    write(*,"(a,*(1x,g0))") "elbow connections     :", elbow%get_nconnection()
    write(*,"(a,*(1x,g0))") "elbow valid           :", elbow%is_valid()
    write(*,"(a,*(1x,g0))") "shared global dofs    :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements              :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element     :", size(elem, 2)

    call elbow%create(6, 5, 8)
    call elbow%export_Xc("vtk/volume_multipatch_elbow")
    call elbow%export_Xg("vtk/volume_multipatch_elbow")
    call elbow%export_Xth_in_Xg("vtk/volume_multipatch_elbow", 8)
    call elbow%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_elbow_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_elbow_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_elbow_Xth_*.vtk")
    call elbow%finalize()

contains

    pure subroutine set_elbow_patch(p, is, ir, it)
        type(nurbs_volume), intent(inout) :: p
        integer, intent(in) :: is, ir, it
        real(rk) :: Xc(ncp*ncp*ncp, 3)
        real(rk) :: u, v, w, s, q, theta, phi, radius, sweep, cth, sth, cph, sph
        integer :: a, b, c, n

        do c = 1, ncp
            w = real(c - 1, rk)/real(ncp - 1, rk)
            phi = 2.0_rk*pi*(real(it - 1, rk) + w)/real(nt, rk)
            cph = cos(phi)
            sph = sin(phi)
            do b = 1, ncp
                v = real(b - 1, rk)/real(ncp - 1, rk)
                q = (real(ir - 1, rk) + v)/real(nr, rk)
                radius = radius_inner + (radius_outer - radius_inner)*q
                do a = 1, ncp
                    u = real(a - 1, rk)/real(ncp - 1, rk)
                    s = (real(is - 1, rk) + u)/real(ns, rk)
                    sweep = bend_angle*s
                    theta = sweep + 0.08_rk*sin(pi*s)*sph
                    cth = cos(theta)
                    sth = sin(theta)
                    n = a + (b - 1)*ncp + (c - 1)*ncp*ncp
                    Xc(n, 1) = (center_radius + radius*cph)*cth
                    Xc(n, 2) = (center_radius + radius*cph)*sth
                    Xc(n, 3) = radius*sph + 0.18_rk*sin(2.0_rk*pi*s)*q
                end do
            end do
        end do
        call p%set([ncp, ncp, ncp], Xc)
    end subroutine
    !===============================================================================

end program volume_multipatch_elbow
