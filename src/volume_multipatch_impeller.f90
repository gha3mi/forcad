!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Geometry-only twisted impeller-style annular volume with curved hub and blade bands.
program volume_multipatch_impeller

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: ns = 12, nz = 2, nb = 3, ncp = 3
    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: hub_inner = 0.55_rk
    real(rk), parameter :: hub_outer = 1.05_rk
    real(rk), parameter :: blade_outer = 2.35_rk
    real(rk), parameter :: height = 0.56_rk

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: impeller
    integer :: hub_id(ns, nz), blade_id(ns, nb, nz), is, ir, iz
    integer, allocatable :: dof_map(:), elem(:,:)

    do iz = 1, nz
        do is = 1, ns
            call set_hub_patch(patch, is, iz)
            call impeller%add_patch(patch, hub_id(is, iz))
        end do
    end do

    do iz = 1, nz
        do ir = 1, nb
            do is = 1, ns
                call set_blade_patch(patch, is, ir, iz)
                call impeller%add_patch(patch, blade_id(is, ir, iz))
            end do
        end do
    end do

    do iz = 1, nz
        do is = 1, ns - 1
            call impeller%connect(hub_id(is, iz), "u_max", hub_id(is + 1, iz), "u_min", 0)
        end do
        call impeller%connect(hub_id(ns, iz), "u_max", hub_id(1, iz), "u_min", 0)
    end do
    do iz = 1, nz - 1
        do is = 1, ns
            call impeller%connect(hub_id(is, iz), "w_max", hub_id(is, iz + 1), "w_min", 0)
        end do
    end do

    do iz = 1, nz
        do ir = 1, nb
            do is = 1, ns - 1
                call impeller%connect(blade_id(is, ir, iz), "u_max", blade_id(is + 1, ir, iz), "u_min", 0)
            end do
            call impeller%connect(blade_id(ns, ir, iz), "u_max", blade_id(1, ir, iz), "u_min", 0)
        end do
    end do
    do iz = 1, nz
        do ir = 1, nb - 1
            do is = 1, ns
                call impeller%connect(blade_id(is, ir, iz), "v_max", blade_id(is, ir + 1, iz), "v_min", 0)
            end do
        end do
    end do
    do iz = 1, nz - 1
        do ir = 1, nb
            do is = 1, ns
                call impeller%connect(blade_id(is, ir, iz), "w_max", blade_id(is, ir, iz + 1), "w_min", 0)
            end do
        end do
    end do
    do iz = 1, nz
        do is = 1, ns
            call impeller%connect(hub_id(is, iz), "v_max", blade_id(is, 1, iz), "v_min", 0)
        end do
    end do

    dof_map = impeller%cmp_dof_map()
    elem = impeller%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "impeller hub patches  :", ns*nz
    write(*,"(a,*(1x,g0))") "impeller blade patches:", ns*nb*nz
    write(*,"(a,*(1x,g0))") "impeller total patches:", impeller%get_npatch()
    write(*,"(a,*(1x,g0))") "impeller connections  :", impeller%get_nconnection()
    write(*,"(a,*(1x,g0))") "impeller valid        :", impeller%is_valid()
    write(*,"(a,*(1x,g0))") "shared global dofs    :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements              :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element     :", size(elem, 2)

    call impeller%create(6, 6, 5)
    call impeller%export_Xc("vtk/volume_multipatch_impeller")
    call impeller%export_Xg("vtk/volume_multipatch_impeller")
    call impeller%export_Xth_in_Xg("vtk/volume_multipatch_impeller", 6)
    call impeller%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_impeller_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_impeller_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_impeller_Xth_*.vtk")
    call impeller%finalize()

contains

    pure subroutine set_hub_patch(p, is, iz)
        type(nurbs_volume), intent(inout) :: p
        integer, intent(in) :: is, iz
        real(rk) :: Xc(ncp*ncp*ncp, 3)
        real(rk) :: u, v, w, theta, radius, z
        integer :: a, b, c, n

        do c = 1, ncp
            w = real(c - 1, rk)/real(ncp - 1, rk)
            z = height*((real(iz - 1, rk) + w)/real(nz, rk) - 0.5_rk)
            do b = 1, ncp
                v = real(b - 1, rk)/real(ncp - 1, rk)
                radius = hub_inner + (hub_outer - hub_inner)*v
                do a = 1, ncp
                    u = real(a - 1, rk)/real(ncp - 1, rk)
                    theta = 2.0_rk*pi*(real(is - 1, rk) + u)/real(ns, rk)
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

    pure subroutine set_blade_patch(p, is, ir, iz)
        type(nurbs_volume), intent(inout) :: p
        integer, intent(in) :: is, ir, iz
        real(rk) :: Xc(ncp*ncp*ncp, 3)
        real(rk) :: u, v, w, q, phase, theta, radius, z, camber
        integer :: a, b, c, n

        do c = 1, ncp
            w = real(c - 1, rk)/real(ncp - 1, rk)
            z = height*((real(iz - 1, rk) + w)/real(nz, rk) - 0.5_rk)
            do b = 1, ncp
                v = real(b - 1, rk)/real(ncp - 1, rk)
                q = (real(ir - 1, rk) + v)/real(nb, rk)
                radius = hub_outer + (blade_outer - hub_outer)*q
                do a = 1, ncp
                    u = real(a - 1, rk)/real(ncp - 1, rk)
                    phase = 2.0_rk*pi*(real(is - 1, rk) + u)/real(ns, rk)
                    camber = 0.62_rk*q**1.25_rk + 0.08_rk*sin(pi*q)*sin(3.0_rk*phase)
                    theta = phase + camber + 0.08_rk*q*(z/height)
                    n = a + (b - 1)*ncp + (c - 1)*ncp*ncp
                    Xc(n, 1) = radius*cos(theta)
                    Xc(n, 2) = radius*sin(theta)
                    Xc(n, 3) = z + 0.10_rk*sin(pi*q)*cos(2.0_rk*phase)
                end do
            end do
        end do
        call p%set([ncp, ncp, ncp], Xc)
    end subroutine
    !===============================================================================

end program volume_multipatch_impeller
