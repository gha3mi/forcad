!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Rational twisted Cn shell extracted exactly into a two-dimensional patch grid.
!! Set continuity_order between zero and min(degree_u, degree_v).
program surface_multipatch_cn_shell

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface

    implicit none

    integer, parameter :: ns = 6, nt = 4, nu = 7, nv = 6
    integer, parameter :: degree_u = nu - 1, degree_v = nv - 1
    integer, parameter :: continuity_order = 3, nlocal = nu*nv
    real(rk), parameter :: pi = acos(-1.0_rk)

    type(nurbs_surface) :: master, patch
    type(nurbs_multipatch_surface) :: shell
    integer :: patch_id(ns,nt), is, it, a, b, i, j, d, n, source, offset, row
    integer :: ncu_refined
    integer, allocatable :: offsets(:), dof_map(:), elem(:,:), rowptr(:), col(:)
    real(rk) :: Xc(nlocal,3), Wc(nlocal), master_point(3), patch_point(3)
    real(rk), allocatable :: refined_Xc(:,:), refined_Wc(:), cval(:), homogeneous(:,:)
    real(rk) :: s, eta, center_x, center_y, center_z, twist, half_width, camber
    real(rk) :: local_y, local_z, xwarp, extraction_error, constraint_residual, residual

    do concurrent (b = 1:nv, a = 1:nu) &
        local(s, eta, center_x, center_y, center_z, twist, half_width, camber, &
            local_y, local_z, xwarp, n)
        s = real(a-1,rk)/real(nu-1,rk)
        eta = 2.0_rk*real(b-1,rk)/real(nv-1,rk) - 1.0_rk
        center_x = 7.0_rk*(s-0.5_rk)
        center_y = 1.15_rk*sin(2.0_rk*pi*s) + 0.34_rk*sin(4.0_rk*pi*s)
        center_z = 0.72_rk*cos(pi*s) + 0.38_rk*sin(3.0_rk*pi*s)
        twist = 1.65_rk*pi*s + 0.28_rk*sin(2.0_rk*pi*s)
        half_width = 1.42_rk*(1.0_rk-0.16_rk*cos(2.0_rk*pi*s))
        camber = 0.46_rk*(1.0_rk-eta*eta)*(0.72_rk+0.28_rk*sin(pi*s))
        local_y = half_width*eta
        local_z = camber + 0.14_rk*eta*eta*sin(3.0_rk*pi*s)
        xwarp = 0.18_rk*eta*eta*sin(2.0_rk*pi*s) + 0.10_rk*eta*sin(pi*s)
        n = a + (b-1)*nu
        Xc(n,1) = center_x + xwarp
        Xc(n,2) = center_y + cos(twist)*local_y - sin(twist)*local_z
        Xc(n,3) = center_z + sin(twist)*local_y + cos(twist)*local_z
        Wc(n) = 1.0_rk + 0.14_rk*(1.0_rk-eta*eta)*&
            (1.0_rk+0.24_rk*cos(2.0_rk*pi*s))
    end do

    call master%set(&
        nc = [nu, nv],&
        Xc = Xc,&
        Wc = Wc)
    call master%insert_knots(&
        dir = 1,&
        Xth = [(real(i,rk)/real(ns,rk), i=1,ns-1)],&
        r   = [(degree_u, i=1,ns-1)])
    call master%insert_knots(&
        dir = 2,&
        Xth = [(real(i,rk)/real(nt,rk), i=1,nt-1)],&
        r   = [(degree_v, i=1,nt-1)])

    refined_Xc = master%get_Xc()
    refined_Wc = master%get_Wc()
    ncu_refined = master%get_nc(1)
    do it = 1, nt
        do is = 1, ns
            do b = 1, nv
                do a = 1, nu
                    n = a + (b-1)*nu
                    source = (is-1)*degree_u + a +&
                        ((it-1)*degree_v+b-1)*ncu_refined
                    do d = 1, 3
                        Xc(n,d) = refined_Xc(source,d)
                    end do
                    Wc(n) = refined_Wc(source)
                end do
            end do
            call patch%set(&
                nc = [nu, nv],&
                Xc = Xc,&
                Wc = Wc)
            call shell%add_patch(&
                patch = patch,&
                id    = patch_id(is,it))
        end do
    end do

    do it = 1, nt
        do is = 1, ns - 1
            call shell%connect(&
                patch_a    = patch_id(is,it),&
                side_a     = "u_max",&
                patch_b    = patch_id(is+1,it),&
                side_b     = "u_min",&
                continuity = continuity_order)
        end do
    end do
    do it = 1, nt - 1
        do is = 1, ns
            call shell%connect(&
                patch_a    = patch_id(is,it),&
                side_a     = "v_max",&
                patch_b    = patch_id(is,it+1),&
                side_b     = "v_min",&
                continuity = continuity_order)
        end do
    end do

    extraction_error = 0.0_rk
    do it = 1, nt
        do is = 1, ns
            patch = shell%get_patch(patch_id(is,it))
            do j = 1, 2
                eta = real(2*j-1,rk)/4.0_rk
                do i = 1, 3
                    s = real(2*i-1,rk)/6.0_rk
                    master_point = master%cmp_Xg([&
                        (real(is-1,rk)+s)/real(ns,rk),&
                        (real(it-1,rk)+eta)/real(nt,rk)])
                    patch_point = patch%cmp_Xg([s, eta])
                    extraction_error = max(&
                        extraction_error,&
                        maxval(abs(master_point-patch_point)))
                end do
            end do
        end do
    end do

    offsets = shell%cmp_dof_offsets()
    dof_map = shell%cmp_dof_map()
    elem = shell%cmp_elem(shared=.true.)
    call shell%cmp_dof_constraint(&
        rowptr = rowptr,&
        col    = col,&
        val    = cval)
    allocate(homogeneous(offsets(size(offsets)),4), source=0.0_rk)
    do it = 1, nt
        do is = 1, ns
            patch = shell%get_patch(patch_id(is,it))
            Xc = patch%get_Xc()
            Wc = patch%get_Wc()
            offset = offsets(patch_id(is,it))
            do i = 1, nlocal
                do d = 1, 3
                    homogeneous(offset+i,d) = Xc(i,d)*Wc(i)
                end do
                homogeneous(offset+i,4) = Wc(i)
            end do
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
    write(*,"(a,i0)") "maximum supported n     : ", min(degree_u,degree_v)
    write(*,"(a,i0)") "surface patches         : ", shell%get_npatch()
    write(*,"(a,i0)") "u-direction patches     : ", ns
    write(*,"(a,i0)") "v-direction patches     : ", nt
    write(*,"(a,i0)") "connections             : ", shell%get_nconnection()
    write(*,"(a,l1)") "rational geometry       : ", shell%is_rational()
    write(*,"(a,l1)") "valid Cn geometry       : ", shell%is_valid()
    write(*,"(a,es12.4)") "patch extraction error  : ", extraction_error
    write(*,"(a,es12.4)") "Cn constraint residual  : ", constraint_residual
    write(*,"(a,i0)") "raw patch dofs          : ", size(dof_map)
    write(*,"(a,i0)") "shared global dofs      : ", maxval(dof_map)
    write(*,"(a,i0)") "constraint rows         : ", size(rowptr)-1
    write(*,"(a,i0)") "constraint nonzeros     : ", size(cval)
    write(*,"(a,i0)") "elements                : ", size(elem,1)
    write(*,"(a,i0)") "nodes per element       : ", size(elem,2)

    call shell%create(&
        res1 = 13,&
        res2 = 11)
    call shell%export_Xc("vtk/surface_multipatch_cn_shell")
    call shell%export_Xg("vtk/surface_multipatch_cn_shell")
    call shell%export_Xth_in_Xg("vtk/surface_multipatch_cn_shell", 11)
    call shell%show(&
        vtkfile_Xc        = "vtk/surface_multipatch_cn_shell_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_multipatch_cn_shell_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_multipatch_cn_shell_Xth_*.vtk")

    call master%finalize()
    call shell%finalize()

end program surface_multipatch_cn_shell
