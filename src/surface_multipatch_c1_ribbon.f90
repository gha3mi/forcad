!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Exact C1 multi-patch helical ribbon with derivative-compatible side rails.
program surface_multipatch_c1_ribbon

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface

    implicit none

    integer, parameter :: ns = 28, nb = 3, nu = 4, nv = 4
    real(rk), parameter :: pi = 3.1415926535897932384626433832795_rk
    real(rk), parameter :: turns = 1.55_rk, inner = 1.10_rk
    real(rk), parameter :: width = 1.75_rk, height = 5.20_rk, rail = 0.66_rk

    type(nurbs_surface) :: patch
    type(nurbs_multipatch_surface) :: ribbon
    real(rk), allocatable :: Xc(:,:), cval(:)
    integer, allocatable :: dof_map(:), elem(:,:), rowptr(:), col(:)
    integer :: deck_id(ns,nb), rail_id(ns,2)
    integer :: i, j, k
    real(rk) :: b0(3), b1(3), d0(3), d1(3), eta, q, t0, t1

    allocate(Xc(nu*nv, 3))

    do j = 1, nb
        do i = 1, ns
            t0 = real(i - 1, rk)/real(ns, rk)
            t1 = real(i, rk)/real(ns, rk)
            do k = 1, nv
                eta = real(k - 1, rk)/real(nv - 1, rk)
                q = (real(j - 1, rk) + eta)/real(nb, rk)
                call eval_deck(t0, q, b0, d0)
                call eval_deck(t1, q, b1, d1)
                Xc(1 + nu*(k - 1), :) = b0
                Xc(2 + nu*(k - 1), :) = b0 + d0/real(nu - 1, rk)
                Xc(3 + nu*(k - 1), :) = b1 - d1/real(nu - 1, rk)
                Xc(4 + nu*(k - 1), :) = b1
            end do
            call patch%set([nu, nv], Xc)
            call ribbon%add_patch(patch, deck_id(i,j))
        end do
    end do

    do j = 1, 2
        do i = 1, ns
            t0 = real(i - 1, rk)/real(ns, rk)
            t1 = real(i, rk)/real(ns, rk)
            do k = 1, nv
                eta = real(k - 1, rk)/real(nv - 1, rk)
                call eval_rail(t0, eta, j, b0, d0)
                call eval_rail(t1, eta, j, b1, d1)
                Xc(1 + nu*(k - 1), :) = b0
                Xc(2 + nu*(k - 1), :) = b0 + d0/real(nu - 1, rk)
                Xc(3 + nu*(k - 1), :) = b1 - d1/real(nu - 1, rk)
                Xc(4 + nu*(k - 1), :) = b1
            end do
            call patch%set([nu, nv], Xc)
            call ribbon%add_patch(patch, rail_id(i,j))
        end do
    end do

    do j = 1, nb
        do i = 1, ns - 1
            call ribbon%connect(deck_id(i,j), "u_max", deck_id(i + 1,j), "u_min", 1)
        end do
    end do
    do j = 1, nb - 1
        do i = 1, ns
            call ribbon%connect(deck_id(i,j), "v_max", deck_id(i,j + 1), "v_min", 0)
        end do
    end do
    do j = 1, 2
        do i = 1, ns - 1
            call ribbon%connect(rail_id(i,j), "u_max", rail_id(i + 1,j), "u_min", 1)
        end do
    end do
    do i = 1, ns
        call ribbon%connect(deck_id(i,1), "v_min", rail_id(i,1), "v_min", 0)
        call ribbon%connect(deck_id(i,nb), "v_max", rail_id(i,2), "v_min", 0)
    end do

    dof_map = ribbon%cmp_dof_map()
    elem = ribbon%cmp_elem(shared=.true.)
    call ribbon%cmp_dof_constraint(rowptr, col, cval)

    write(*,"(a,*(1x,g0))") "C1 ribbon deck patches   :", ns*nb
    write(*,"(a,*(1x,g0))") "C1 ribbon rail patches   :", 2*ns
    write(*,"(a,*(1x,g0))") "C1 ribbon total patches  :", ribbon%get_npatch()
    write(*,"(a,*(1x,g0))") "C1 ribbon connections    :", ribbon%get_nconnection()
    write(*,"(a,*(1x,g0))") "C1 ribbon valid          :", ribbon%is_valid()
    write(*,"(a,*(1x,g0))") "local control dofs       :", size(dof_map)
    write(*,"(a,*(1x,g0))") "shared global dofs       :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "Cn constraint rows       :", size(rowptr) - 1
    write(*,"(a,*(1x,g0))") "Cn constraint nnz        :", size(cval)
    write(*,"(a,*(1x,g0))") "elements                 :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element        :", size(elem, 2)

    call ribbon%create(31, 13)
    call ribbon%export_Xc("vtk/surface_multipatch_c1_ribbon")
    call ribbon%export_Xg("vtk/surface_multipatch_c1_ribbon")
    call ribbon%export_Xth_in_Xg("vtk/surface_multipatch_c1_ribbon", 31)
    call ribbon%show(&
        vtkfile_Xc        = "vtk/surface_multipatch_c1_ribbon_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_multipatch_c1_ribbon_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_multipatch_c1_ribbon_Xth_*.vtk")
    call ribbon%finalize()

contains

    pure subroutine eval_deck(t, q, b, d)
        real(rk), intent(in) :: t, q
        real(rk), intent(out) :: b(3), d(3)
        real(rk) :: theta, dtheta, radius, drdt, z, dzdt, shape

        dtheta = 2.0_rk*pi*turns
        theta = dtheta*t
        shape = q*(1.0_rk - q)
        radius = inner + width*q + 0.12_rk*shape*sin(6.0_rk*pi*t)
        drdt = 0.72_rk*pi*shape*cos(6.0_rk*pi*t)
        z = height*t + 0.18_rk*shape*cos(4.0_rk*pi*t)
        dzdt = height - 0.72_rk*pi*shape*sin(4.0_rk*pi*t)

        b = [radius*cos(theta), radius*sin(theta), z]
        d = [drdt*cos(theta) - radius*dtheta*sin(theta), &
            drdt*sin(theta) + radius*dtheta*cos(theta), dzdt]/real(ns, rk)
    end subroutine
    !===============================================================================

    pure subroutine eval_rail(t, eta, side, b, d)
        real(rk), intent(in) :: t, eta
        integer, intent(in) :: side
        real(rk), intent(out) :: b(3), d(3)
        real(rk) :: base(3), tangent(3), theta, dtheta, q, sgn, offset

        q = merge(0.0_rk, 1.0_rk, side == 1)
        sgn = merge(-1.0_rk, 1.0_rk, side == 1)
        offset = 0.16_rk*eta
        dtheta = 2.0_rk*pi*turns
        theta = dtheta*t

        call eval_deck(t, q, base, tangent)
        b = base + [sgn*offset*cos(theta), sgn*offset*sin(theta), rail*eta]
        d = tangent + sgn*offset*dtheta*[-sin(theta), cos(theta), 0.0_rk]/real(ns, rk)
    end subroutine
    !===============================================================================

end program surface_multipatch_c1_ribbon
