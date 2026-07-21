!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Rational twisted Cn sweep extracted exactly into volume patches.
!! Set continuity_order between zero and degree_u to select the continuity.
program volume_multipatch_cn_sweep

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: ns = 8, nu = 7, nv = 4, nw = 4
    integer, parameter :: degree_u = nu - 1, continuity_order = 3
    integer, parameter :: nlocal = nu*nv*nw
    real(rk), parameter :: pi = acos(-1.0_rk)

    type(nurbs_volume) :: master, patch
    type(nurbs_multipatch_volume) :: sweep
    integer :: patch_id(ns), is, a, b, c, i, j, k, d, n, source, offset, row
    integer :: ncu_refined
    integer, allocatable :: offsets(:), dof_map(:), elem(:,:), rowptr(:), col(:)
    real(rk) :: Xc(nlocal,3), Wc(nlocal)
    real(rk), allocatable :: refined_Xc(:,:), refined_Wc(:), cval(:), homogeneous(:,:)
    real(rk), allocatable :: master_point(:), patch_point(:)
    real(rk) :: s, eta, zeta, center_x, center_y, center_z, twist
    real(rk) :: radius_y, radius_z, local_y, local_z, xwarp
    real(rk) :: extraction_error, constraint_residual, residual

    do concurrent (c = 1:nw, b = 1:nv, a = 1:nu) &
        local(s, eta, zeta, center_x, center_y, center_z, twist, radius_y, &
            radius_z, local_y, local_z, xwarp, n)
        s = real(a-1,rk)/real(nu-1,rk)
        eta = 2.0_rk*real(b-1,rk)/real(nv-1,rk) - 1.0_rk
        zeta = 2.0_rk*real(c-1,rk)/real(nw-1,rk) - 1.0_rk
        center_x = 7.0_rk*(s-0.5_rk)
        center_y = 1.25_rk*sin(2.0_rk*pi*s) + 0.42_rk*sin(4.0_rk*pi*s)
        center_z = 0.90_rk*cos(pi*s) + 0.48_rk*sin(3.0_rk*pi*s)
        twist = 2.40_rk*pi*s + 0.30_rk*sin(2.0_rk*pi*s)
        radius_y = 0.72_rk*(1.0_rk-0.20_rk*s+0.16_rk*sin(pi*s))
        radius_z = 0.54_rk*(1.0_rk+0.20_rk*cos(2.0_rk*pi*s))
        local_y = radius_y*eta*(1.0_rk+0.12_rk*zeta)
        local_z = radius_z*zeta*(1.0_rk-0.10_rk*eta*eta)
        xwarp = 0.14_rk*eta*zeta*sin(pi*s) + 0.06_rk*(eta*eta-zeta*zeta)*sin(2.0_rk*pi*s)
        n = a + (b-1)*nu + (c-1)*nu*nv
        Xc(n,1) = center_x + xwarp
        Xc(n,2) = center_y + cos(twist)*local_y - sin(twist)*local_z
        Xc(n,3) = center_z + sin(twist)*local_y + cos(twist)*local_z
        Wc(n) = 1.0_rk + 0.10_rk*(1.0_rk-eta*eta)*(1.0_rk-zeta*zeta)*&
            (1.0_rk+0.20_rk*cos(2.0_rk*pi*s))
    end do

    call master%set(&
        nc = [nu, nv, nw],&
        Xc = Xc,&
        Wc = Wc)
    call master%insert_knots(&
        dir = 1,&
        Xth = [(real(i,rk)/real(ns,rk), i=1,ns-1)],&
        r   = [(degree_u, i=1,ns-1)])

    refined_Xc = master%get_Xc()
    refined_Wc = master%get_Wc()
    ncu_refined = master%get_nc(1)
    do is = 1, ns
        do c = 1, nw
            do b = 1, nv
                do a = 1, nu
                    n = a + (b-1)*nu + (c-1)*nu*nv
                    source = (is-1)*degree_u + a + (b-1)*ncu_refined +&
                        (c-1)*ncu_refined*nv
                    do d = 1, 3
                        Xc(n,d) = refined_Xc(source,d)
                    end do
                    Wc(n) = refined_Wc(source)
                end do
            end do
        end do
        call patch%set(&
            nc = [nu, nv, nw],&
            Xc = Xc,&
            Wc = Wc)
        call sweep%add_patch(&
            patch = patch,&
            id    = patch_id(is))
    end do
    do is = 1, ns - 1
        call sweep%connect(&
            patch_a    = patch_id(is),&
            side_a     = "u_max",&
            patch_b    = patch_id(is+1),&
            side_b     = "u_min",&
            continuity = continuity_order)
    end do

    extraction_error = 0.0_rk
    do is = 1, ns
        patch = sweep%get_patch(is)
        do i = 1, 3
            s = real(2*i-1,rk)/6.0_rk
            do j = 1, 2
                eta = real(j,rk)/3.0_rk
                do k = 1, 2
                    zeta = real(k,rk)/3.0_rk
                    master_point = master%cmp_Xg(&
                        [(real(is-1,rk)+s)/real(ns,rk), eta, zeta])
                    patch_point = patch%cmp_Xg([s, eta, zeta])
                    extraction_error = max(extraction_error, maxval(abs(master_point-patch_point)))
                end do
            end do
        end do
    end do

    offsets = sweep%cmp_dof_offsets()
    dof_map = sweep%cmp_dof_map()
    elem = sweep%cmp_elem(shared=.true.)
    call sweep%cmp_dof_constraint(&
        rowptr = rowptr,&
        col    = col,&
        val    = cval)
    allocate(homogeneous(offsets(size(offsets)),4), source=0.0_rk)
    do is = 1, ns
        patch = sweep%get_patch(is)
        Xc = patch%get_Xc()
        Wc = patch%get_Wc()
        offset = offsets(is)
        do i = 1, nlocal
            do d = 1, 3
                homogeneous(offset+i,d) = Xc(i,d)*Wc(i)
            end do
            homogeneous(offset+i,4) = Wc(i)
        end do
    end do
    constraint_residual = 0.0_rk
    do row = 1, size(rowptr) - 1
        do d = 1, 4
            residual = 0.0_rk
            do i = rowptr(row), rowptr(row+1) - 1
                residual = residual + cval(i)*homogeneous(col(i),d)
            end do
            constraint_residual = max(constraint_residual, abs(residual))
        end do
    end do

    write(*,"(a,i0)") "continuity order n      : ", continuity_order
    write(*,"(a,i0)") "maximum supported n     : ", degree_u
    write(*,"(a,i0)") "volume patches          : ", sweep%get_npatch()
    write(*,"(a,i0)") "connections             : ", sweep%get_nconnection()
    write(*,"(a,l1)") "rational geometry       : ", sweep%is_rational()
    write(*,"(a,l1)") "valid Cn geometry       : ", sweep%is_valid()
    write(*,"(a,es12.4)") "patch extraction error  : ", extraction_error
    write(*,"(a,es12.4)") "Cn constraint residual  : ", constraint_residual
    write(*,"(a,i0)") "raw patch dofs          : ", size(dof_map)
    write(*,"(a,i0)") "shared global dofs      : ", maxval(dof_map)
    write(*,"(a,i0)") "constraint rows         : ", size(rowptr)-1
    write(*,"(a,i0)") "constraint nonzeros     : ", size(cval)
    write(*,"(a,i0)") "elements                : ", size(elem,1)
    write(*,"(a,i0)") "nodes per element       : ", size(elem,2)

    call sweep%create(&
        res1 = 13,&
        res2 = 7,&
        res3 = 7)
    call sweep%export_Xc("vtk/volume_multipatch_cn_sweep")
    call sweep%export_Xg("vtk/volume_multipatch_cn_sweep")
    call sweep%export_Xth_in_Xg("vtk/volume_multipatch_cn_sweep", 9)
    call sweep%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_cn_sweep_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_cn_sweep_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_cn_sweep_Xth_*.vtk")

    call master%finalize()
    call sweep%finalize()

end program volume_multipatch_cn_sweep
