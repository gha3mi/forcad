!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Solve a uniformly loaded Euler-Bernoulli cantilever with two C1 IGA patches.
!! The dimensionless weak form is integral(u'' v'') dx = integral(v) dx on
!! [0,1], with u(0)=u'(0)=0 and natural free-end conditions at x=1.
program curve_iga_beam_c1

    use forcad, only: rk, nurbs_curve, nurbs_multipatch_curve, solve

    implicit none

    integer, parameter :: degree = 3, nc = 7, ngauss = 4
    real(rk), parameter :: knot(nc+degree+1) = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
        0.25_rk, 0.50_rk, 0.75_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    type(nurbs_curve) :: patch
    type(nurbs_multipatch_curve) :: beam, deflection
    integer :: patch_id(2), deflection_id(2), ip, is, ig, a, b, row, j, offset, ia, ib
    integer :: nraw, nconstraint, fixed(2)
    integer, allocatable :: offsets(:), rowptr(:), col(:), elem(:,:)
    real(rk), allocatable :: val(:), K(:,:), force(:), system(:,:), rhs(:,:), solution(:,:)
    real(rk), allocatable :: Xc(:,:), N(:), dN(:), dNdx(:,:), d2Ndx2(:)
    real(rk) :: dL, x, uh, exact
    real(rk) :: value_left, value_right, slope_left, slope_right
    real(rk) :: constraint_residual, row_residual, l2_error

    call set_line_patch(patch, 0.0_rk, 0.5_rk)
    call beam%add_patch(patch, patch_id(1))
    call set_line_patch(patch, 0.5_rk, 1.0_rk)
    call beam%add_patch(patch, patch_id(2))
    call beam%connect(patch_id(1), "right", patch_id(2), "left", 1)

    offsets = beam%cmp_dof_offsets()
    nraw = offsets(size(offsets))
    allocate(K(nraw,nraw), force(nraw), source=0.0_rk)

    do ip = 1, beam%get_npatch()
        patch = beam%get_patch(ip)
        Xc = patch%get_Xc()
        elem = patch%cmp_elem()
        offset = offsets(ip)
        do is = 1, size(elem,1)
            do ig = 1, ngauss
                call patch%ansatz(&
                    ie          = is,&
                    ig          = ig,&
                    Tgc         = N,&
                    dTgc_dXg    = dNdx,&
                    d2Tgc_dXg2  = d2Ndx2,&
                    dL          = dL,&
                    ngauss      = ngauss)
                do a = 1, size(elem,2)
                    ia = elem(is,a)
                    force(offset+ia) = force(offset+ia) + N(a)*dL
                end do
                do b = 1, size(elem,2)
                    ib = elem(is,b)
                    do a = 1, size(elem,2)
                        ia = elem(is,a)
                        K(offset+ia,offset+ib) = K(offset+ia,offset+ib) + d2Ndx2(a)*d2Ndx2(b)*dL
                    end do
                end do
            end do
        end do
    end do

    call beam%cmp_dof_constraint(rowptr, col, val)
    nconstraint = size(rowptr) - 1
    allocate(system(nraw+nconstraint,nraw+nconstraint), rhs(nraw+nconstraint,1), source=0.0_rk)
    system(1:nraw,1:nraw) = K
    rhs(1:nraw,1) = force
    do row = 1, nconstraint
        do j = rowptr(row), rowptr(row+1) - 1
            system(nraw+row,col(j)) = val(j)
            system(col(j),nraw+row) = val(j)
        end do
    end do

    fixed = [1, 2]
    do j = 1, size(fixed)
        system(fixed(j),:) = 0.0_rk
        system(:,fixed(j)) = 0.0_rk
        system(fixed(j),fixed(j)) = 1.0_rk
        rhs(fixed(j),1) = 0.0_rk
    end do
    solution = solve(system, rhs)

    constraint_residual = 0.0_rk
    do row = 1, nconstraint
        row_residual = 0.0_rk
        do j = rowptr(row), rowptr(row+1) - 1
            row_residual = row_residual + val(j)*solution(col(j),1)
        end do
        constraint_residual = max(constraint_residual, abs(row_residual))
    end do

    patch = beam%get_patch(1)
    Xc = patch%get_Xc()
    call patch%derivative(&
        Xt    = 1.0_rk,&
        dTgc = dN,&
        Tgc  = N)
    value_left = dot_product(N, solution(1:nc,1))
    slope_left = dot_product(dN, solution(1:nc,1))/dot_product(dN, Xc(:,1))
    patch = beam%get_patch(2)
    Xc = patch%get_Xc()
    call patch%derivative(&
        Xt    = 0.0_rk,&
        dTgc = dN,&
        Tgc  = N)
    value_right = dot_product(N, solution(offsets(2)+1:offsets(2)+nc,1))
    slope_right = dot_product(dN, solution(offsets(2)+1:offsets(2)+nc,1))/dot_product(dN, Xc(:,1))

    l2_error = 0.0_rk
    do ip = 1, beam%get_npatch()
        patch = beam%get_patch(ip)
        Xc = patch%get_Xc()
        elem = patch%cmp_elem()
        offset = offsets(ip)
        do is = 1, size(elem,1)
            do ig = 1, ngauss
                call patch%ansatz(&
                    ie         = is,&
                    ig         = ig,&
                    Tgc        = N,&
                    dTgc_dXg   = dNdx,&
                    dL         = dL,&
                    ngauss     = ngauss)
                x = 0.0_rk
                uh = 0.0_rk
                do a = 1, size(elem,2)
                    ia = elem(is,a)
                    x = x + N(a)*Xc(ia,1)
                    uh = uh + N(a)*solution(offset+ia,1)
                end do
                exact = x*x*(6.0_rk-4.0_rk*x+x*x)/24.0_rk
                l2_error = l2_error + (uh-exact)**2*dL
            end do
        end do
    end do

    write(*,"(a,l1)") "valid C1 geometry       : ", beam%is_valid()
    write(*,"(a,i0)") "raw patch dofs         : ", nraw
    write(*,"(a,i0)") "interface constraints  : ", nconstraint
    write(*,"(a,es12.4)") "constraint residual    : ", constraint_residual
    write(*,"(a,es12.4)") "interface value jump   : ", abs(value_left-value_right)
    write(*,"(a,es12.4)") "interface slope jump   : ", abs(slope_left-slope_right)
    write(*,"(a,es12.4)") "tip displacement       : ", solution(nraw,1)
    write(*,"(a,es12.4)") "exact tip displacement : ", 0.125_rk
    write(*,"(a,es12.4)") "L2 error               : ", sqrt(l2_error)

    do ip = 1, beam%get_npatch()
        patch = beam%get_patch(ip)
        Xc = patch%get_Xc()
        Xc(:,2) = solution(offsets(ip)+1:offsets(ip)+nc,1)
        call patch%set(knot, Xc)
        call deflection%add_patch(patch, deflection_id(ip))
    end do
    call deflection%connect(deflection_id(1), "right", deflection_id(2), "left", 1)
    write(*,"(a,l1)") "valid C1 solution       : ", deflection%is_valid()
    call deflection%create(81)
    call deflection%export_Xc("vtk/curve_iga_beam_c1")
    call deflection%export_Xg("vtk/curve_iga_beam_c1")
    call deflection%export_Xth_in_Xg("vtk/curve_iga_beam_c1", 17)
    call deflection%show(&
        vtkfile_Xc        = "vtk/curve_iga_beam_c1_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/curve_iga_beam_c1_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/curve_iga_beam_c1_Xth_*.vtk")

    call beam%finalize()
    call deflection%finalize()

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct an affine cubic B-spline patch from its Greville points.
    pure subroutine set_line_patch(curve, x0, x1)
        type(nurbs_curve), intent(inout) :: curve
        real(rk), intent(in) :: x0, x1
        real(rk) :: control_points(nc,2)
        integer :: i

        do concurrent (i = 1:nc)
            control_points(i,1) = x0 + (x1-x0)*sum(knot(i+1:i+degree))/real(degree,rk)
            control_points(i,2) = 0.0_rk
        end do
        call curve%set(knot, control_points)
    end subroutine
    !===============================================================================
end program curve_iga_beam_c1
