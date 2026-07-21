!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Geometry-only spiral volute volume with an expanding tangential outlet.
program volume_multipatch_volute

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: ns = 14, nr = 2, nz = 2, nout = 4, ncp = 3
    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: theta0 = -0.35_rk*pi
    real(rk), parameter :: sweep = 1.72_rk*pi
    real(rk), parameter :: radius_inner = 0.62_rk
    real(rk), parameter :: height = 0.62_rk
    real(rk), parameter :: outlet_length = 2.35_rk

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: volute
    integer :: spiral_id(ns, nr, nz), outlet_id(nout, nr, nz)
    integer :: is, ir, iz, io
    integer, allocatable :: dof_map(:), elem(:,:)

    do iz = 1, nz
        do ir = 1, nr
            do is = 1, ns
                call set_spiral_patch(patch, is, ir, iz)
                call volute%add_patch(patch, spiral_id(is, ir, iz))
            end do
        end do
    end do

    do iz = 1, nz
        do ir = 1, nr
            do io = 1, nout
                call set_outlet_patch(patch, io, ir, iz)
                call volute%add_patch(patch, outlet_id(io, ir, iz))
            end do
        end do
    end do

    do iz = 1, nz
        do ir = 1, nr
            do is = 1, ns - 1
                call volute%connect(spiral_id(is, ir, iz), "u_max", spiral_id(is + 1, ir, iz), "u_min", 0)
            end do
            call volute%connect(spiral_id(ns, ir, iz), "u_max", outlet_id(1, ir, iz), "u_min", 0)
            do io = 1, nout - 1
                call volute%connect(outlet_id(io, ir, iz), "u_max", outlet_id(io + 1, ir, iz), "u_min", 0)
            end do
        end do
    end do
    do iz = 1, nz
        do is = 1, ns
            call volute%connect(spiral_id(is, 1, iz), "v_max", spiral_id(is, 2, iz), "v_min", 0)
        end do
        do io = 1, nout
            call volute%connect(outlet_id(io, 1, iz), "v_max", outlet_id(io, 2, iz), "v_min", 0)
        end do
    end do
    do ir = 1, nr
        do is = 1, ns
            call volute%connect(spiral_id(is, ir, 1), "w_max", spiral_id(is, ir, 2), "w_min", 0)
        end do
        do io = 1, nout
            call volute%connect(outlet_id(io, ir, 1), "w_max", outlet_id(io, ir, 2), "w_min", 0)
        end do
    end do

    dof_map = volute%cmp_dof_map()
    elem = volute%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "volute spiral patches :", ns*nr*nz
    write(*,"(a,*(1x,g0))") "volute outlet patches :", nout*nr*nz
    write(*,"(a,*(1x,g0))") "volute total patches  :", volute%get_npatch()
    write(*,"(a,*(1x,g0))") "volute connections    :", volute%get_nconnection()
    write(*,"(a,*(1x,g0))") "volute valid          :", volute%is_valid()
    write(*,"(a,*(1x,g0))") "shared global dofs    :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements              :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element     :", size(elem, 2)

    call volute%create(6, 5, 5)
    call volute%export_Xc("vtk/volume_multipatch_volute")
    call volute%export_Xg("vtk/volume_multipatch_volute")
    call volute%export_Xth_in_Xg("vtk/volume_multipatch_volute", 6)
    call volute%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_volute_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_volute_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_volute_Xth_*.vtk")
    call volute%finalize()

contains

    pure subroutine set_spiral_patch(p, is, ir, iz)
        type(nurbs_volume), intent(inout) :: p
        integer, intent(in) :: is, ir, iz
        real(rk) :: Xc(ncp*ncp*ncp, 3)
        real(rk) :: u, v, w, s, q, theta, radius, radius_outer, z
        integer :: a, b, c, n

        do c = 1, ncp
            w = real(c - 1, rk)/real(ncp - 1, rk)
            do b = 1, ncp
                v = real(b - 1, rk)/real(ncp - 1, rk)
                q = (real(ir - 1, rk) + v)/real(nr, rk)
                do a = 1, ncp
                    u = real(a - 1, rk)/real(ncp - 1, rk)
                    s = (real(is - 1, rk) + u)/real(ns, rk)
                    theta = theta0 + sweep*s
                    radius_outer = 1.18_rk + 0.78_rk*s + 0.08_rk*sin(pi*s)
                    radius = radius_inner + (radius_outer - radius_inner)*q
                    z = height*((real(iz - 1, rk) + w)/real(nz, rk) - 0.5_rk)*(1.0_rk + 0.08_rk*sin(theta))
                    n = a + (b - 1)*ncp + (c - 1)*ncp*ncp
                    Xc(n, 1) = radius*cos(theta)
                    Xc(n, 2) = radius*sin(theta)
                    Xc(n, 3) = z
                end do
            end do
        end do
        call p%set([ncp, ncp, ncp], Xc)
    end subroutine
    !===============================================================================

    pure subroutine set_outlet_patch(p, io, ir, iz)
        type(nurbs_volume), intent(inout) :: p
        integer, intent(in) :: io, ir, iz
        real(rk) :: Xc(ncp*ncp*ncp, 3)
        real(rk) :: u, v, w, s, q, theta, radius_outer, width, radius, z, radial(2), tangent(2)
        integer :: a, b, c, n

        theta = theta0 + sweep
        radial = [cos(theta), sin(theta)]
        tangent = [-sin(theta), cos(theta)]
        radius_outer = 1.18_rk + 0.78_rk + 0.08_rk*sin(pi)
        width = radius_outer - radius_inner

        do c = 1, ncp
            w = real(c - 1, rk)/real(ncp - 1, rk)
            do b = 1, ncp
                v = real(b - 1, rk)/real(ncp - 1, rk)
                q = (real(ir - 1, rk) + v)/real(nr, rk)
                do a = 1, ncp
                    u = real(a - 1, rk)/real(ncp - 1, rk)
                    s = (real(io - 1, rk) + u)/real(nout, rk)
                    radius = 0.5_rk*(radius_inner + radius_outer) + width*(q - 0.5_rk)*(1.0_rk + 0.22_rk*s)
                    radius = radius + 0.25_rk*s*s
                    z = height*((real(iz - 1, rk) + w)/real(nz, rk) - 0.5_rk)*(1.0_rk + 0.08_rk*sin(theta))
                    z = z*(1.0_rk + 0.16_rk*s)
                    n = a + (b - 1)*ncp + (c - 1)*ncp*ncp
                    Xc(n, 1) = radius*radial(1) + outlet_length*s*tangent(1)
                    Xc(n, 2) = radius*radial(2) + outlet_length*s*tangent(2)
                    Xc(n, 3) = z
                end do
            end do
        end do
        call p%set([ncp, ncp, ncp], Xc)
    end subroutine
    !===============================================================================

end program volume_multipatch_volute
