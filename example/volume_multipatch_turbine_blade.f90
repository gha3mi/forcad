!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Twisted rational turbine blade volume with a dovetail root and platform.
program volume_multipatch_turbine_blade

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: nu = 4, nv = 2, nroot = 3, nplatform = 2, nblade = 9
    integer, parameter :: nw = nroot + nplatform + nblade, ncp = 3
    real(rk), parameter :: pi = acos(-1.0_rk), height = 7.0_rk
    real(rk), parameter :: root_top = height*real(nroot,rk)/real(nw,rk)
    real(rk), parameter :: blade_base = height*real(nroot+nplatform,rk)/real(nw,rk)

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: blade
    integer :: id(nu,nv,nw), i, j, k
    integer, allocatable :: dof_map(:), elem(:,:)

    do k = 1, nw
        do j = 1, nv
            do i = 1, nu
                call set_blade_patch(patch, i, j, k)
                call blade%add_patch(patch, id(i,j,k))
            end do
        end do
    end do

    do k = 1, nw
        do j = 1, nv
            do i = 1, nu - 1
                call blade%connect(id(i,j,k), "u_max", id(i+1,j,k), "u_min", 0)
            end do
        end do
        do j = 1, nv - 1
            do i = 1, nu
                call blade%connect(id(i,j,k), "v_max", id(i,j+1,k), "v_min", 0)
            end do
        end do
    end do
    do k = 1, nw - 1
        do j = 1, nv
            do i = 1, nu
                call blade%connect(id(i,j,k), "w_max", id(i,j,k+1), "w_min", 0)
            end do
        end do
    end do

    dof_map = blade%cmp_dof_map()
    elem = blade%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "root patches          :", nu*nv*nroot
    write(*,"(a,*(1x,g0))") "platform patches      :", nu*nv*nplatform
    write(*,"(a,*(1x,g0))") "airfoil patches       :", nu*nv*nblade
    write(*,"(a,*(1x,g0))") "total patches         :", blade%get_npatch()
    write(*,"(a,*(1x,g0))") "connections           :", blade%get_nconnection()
    write(*,"(a,*(1x,g0))") "valid connectivity    :", blade%is_valid()
    write(*,"(a,*(1x,g0))") "rational geometry     :", blade%is_rational()
    write(*,"(a,*(1x,g0))") "shared global dofs    :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements              :", size(elem,1)
    write(*,"(a,*(1x,g0))") "nodes per element     :", size(elem,2)

    call blade%create(7, 5, 7)
    call blade%export_Xc("vtk/volume_multipatch_turbine_blade")
    call blade%export_Xg("vtk/volume_multipatch_turbine_blade")
    call blade%export_Xth_in_Xg("vtk/volume_multipatch_turbine_blade", 7)
    call blade%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_turbine_blade_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_turbine_blade_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_turbine_blade_Xth_*.vtk")
    call blade%finalize()

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct one quadratic block from the global blade parameterization.
    pure subroutine set_blade_patch(p, iu, iv, iw)
        type(nurbs_volume), intent(inout) :: p
        integer, intent(in) :: iu, iv, iw
        real(rk) :: Xc(ncp*ncp*ncp,3), Wc(ncp*ncp*ncp)
        real(rk) :: u, v, w, z, xi, eta, q, flare, hx, hy, chord
        real(rk) :: thickness, camber, angle, sweep, lean, xlocal, ylocal
        real(rk) :: tip, ca, sa
        integer :: a, b, c, n

        do concurrent (c = 1:ncp, b = 1:ncp, a = 1:ncp) &
            local(u, v, w, z, xi, eta, q, flare, hx, hy, chord, thickness, &
                camber, angle, sweep, lean, xlocal, ylocal, tip, ca, sa, n)
            u = (real(iu-1,rk) + real(a-1,rk)/real(ncp-1,rk))/real(nu,rk)
            v = (real(iv-1,rk) + real(b-1,rk)/real(ncp-1,rk))/real(nv,rk)
            w = (real(iw-1,rk) + real(c-1,rk)/real(ncp-1,rk))/real(nw,rk)
            z = height*w
            xi = 2.0_rk*u - 1.0_rk
            eta = 2.0_rk*v - 1.0_rk

            if (z < root_top) then
                q = z/root_top
                flare = sin(3.0_rk*pi*q)
                hx = 0.70_rk + 0.20_rk*flare*flare + 0.08_rk*(1.0_rk-q)
                hy = 0.33_rk + 0.13_rk*flare*flare + 0.05_rk*(1.0_rk-q)
                xlocal = hx*xi*(1.0_rk-0.06_rk*eta*eta) + 0.07_rk*eta*flare
                ylocal = hy*eta*(1.0_rk-0.04_rk*xi*xi)
                sweep = 0.0_rk
                lean = 0.0_rk
                angle = 0.0_rk
            else if (z <= blade_base) then
                q = (z-root_top)/(blade_base-root_top)
                flare = sin(pi*q)**2
                hx = 0.70_rk + 0.12_rk*q + 0.90_rk*flare
                hy = 0.33_rk - 0.08_rk*q + 0.65_rk*flare
                xlocal = hx*xi*(1.0_rk-0.08_rk*eta*eta)
                ylocal = hy*eta*(1.0_rk-0.06_rk*xi*xi)
                sweep = 0.0_rk
                lean = 0.0_rk
                angle = 0.0_rk
            else
                q = (z-blade_base)/(height-blade_base)
                chord = 1.64_rk*(1.0_rk-0.42_rk*q+0.06_rk*sin(pi*q))
                thickness = chord*(0.045_rk + (0.17_rk-0.045_rk*q)* &
                    max(0.0_rk,sin(pi*u))**0.72_rk*(1.0_rk-0.25_rk*u))
                camber = chord*(0.08_rk*sin(pi*u)-0.025_rk*sin(2.0_rk*pi*u))*(1.0_rk-0.25_rk*q)
                xlocal = chord*(u-0.48_rk)
                ylocal = camber + eta*thickness
                angle = 1.05_rk*q + 0.12_rk*sin(pi*q)
                sweep = 0.35_rk*q**1.3_rk - 0.08_rk*sin(pi*q)
                lean = 0.42_rk*sin(0.70_rk*pi*q) + 0.10_rk*q
            end if

            ca = cos(angle)
            sa = sin(angle)
            tip = max(0.0_rk, min(1.0_rk, (z-blade_base)/(height-blade_base)))
            q = max(0.0_rk, (tip-0.82_rk)/0.18_rk)
            n = a + (b-1)*ncp + (c-1)*ncp*ncp
            Xc(n,1) = sweep + ca*xlocal - sa*ylocal
            Xc(n,2) = lean + sa*xlocal + ca*ylocal
            Xc(n,3) = z + 0.18_rk*q*q*(1.0_rk-xi*xi)*(1.0_rk-eta*eta)
            Wc(n) = 1.0_rk + sin(pi*u)**2*(0.04_rk*(1.0_rk-0.30_rk*eta*eta) + &
                0.06_rk*tip*(1.0_rk-eta*eta))
        end do

        call p%set([ncp,ncp,ncp], Xc, Wc)
    end subroutine
    !===============================================================================
end program volume_multipatch_turbine_blade
