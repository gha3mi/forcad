!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Solve a clamped 3D biharmonic problem with two C1 IGA volume patches.
!! The weak form is integral(Hessian(u):Hessian(v)) dV = integral(f*v) dV.
!! Fixed DOFs and inactive constraints are eliminated before system allocation.
program volume_iga_c1

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume, solve

    implicit none

    integer, parameter :: degree = 4, ncu = 6, ncv = 5, ncw = 5
    integer, parameter :: ncp = ncu*ncv*ncw, ngauss = 5, npatch = 2
    real(rk), parameter :: visual_scale = 100.0_rk

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: domain, deformation
    integer :: patch_id(npatch), deformation_id(npatch)
    integer :: ip, is, ig, i, j, k, a, b, ia, ib, offset, base, row
    integer :: row_a, row_b, active_row, nraw, nfree, nconstraint, nactive, nelem
    integer, allocatable :: offsets(:), rowptr(:), col(:), elem(:,:), free_index(:)
    logical, allocatable :: fixed_dof(:)
    real(rk), allocatable :: val(:), system(:,:), rhs(:,:), reduced_solution(:,:), coefficient(:)
    real(rk), allocatable :: Xc(:,:), N(:), dN(:,:), d2N(:,:,:)
    real(rk) :: dV, x, y, z, uh, load_value, stiffness, row_norm
    real(rk) :: constraint_residual, row_residual, l2_error

    call patch%set_hexahedron(&
        L  = [0.5_rk, 1.0_rk, 1.0_rk],&
        nc = [2, 2, 2])
    do i = 1, 3
        call patch%elevate_degree(&
            dir = i,&
            t   = degree-1)
    end do
    call patch%insert_knots(&
        dir = 1,&
        Xth = [0.43_rk],&
        r   = [1])
    call domain%add_patch(&
        patch = patch,&
        id    = patch_id(1))
    call patch%translate_Xc([0.5_rk, 0.0_rk, 0.0_rk])
    call domain%add_patch(&
        patch = patch,&
        id    = patch_id(2))
    call domain%connect(&
        patch_a    = patch_id(1),&
        side_a     = "u_max",&
        patch_b    = patch_id(2),&
        side_b     = "u_min",&
        continuity = 1)

    offsets = domain%cmp_dof_offsets()
    nraw = offsets(size(offsets))
    call domain%cmp_dof_constraint(&
        rowptr = rowptr,&
        col    = col,&
        val    = val)
    nconstraint = size(rowptr) - 1

    allocate(fixed_dof(nraw), source=.false.)
    do k = 1, ncw
        do j = 1, ncv
            base = (k-1)*ncu*ncv + (j-1)*ncu
            fixed_dof(offsets(1)+base+1) = .true.
            fixed_dof(offsets(1)+base+2) = .true.
            fixed_dof(offsets(2)+base+ncu-1) = .true.
            fixed_dof(offsets(2)+base+ncu) = .true.
        end do
    end do
    do ip = 1, npatch
        offset = offsets(ip)
        do k = 1, ncw
            base = offset + (k-1)*ncu*ncv
            do i = 1, ncu
                fixed_dof(base+i) = .true.
                fixed_dof(base+ncu+i) = .true.
                fixed_dof(base+(ncv-2)*ncu+i) = .true.
                fixed_dof(base+(ncv-1)*ncu+i) = .true.
            end do
        end do
        do j = 1, ncv
            base = offset + (j-1)*ncu
            do i = 1, ncu
                fixed_dof(base+i) = .true.
                fixed_dof(base+ncu*ncv+i) = .true.
                fixed_dof(base+(ncw-2)*ncu*ncv+i) = .true.
                fixed_dof(base+(ncw-1)*ncu*ncv+i) = .true.
            end do
        end do
    end do

    allocate(free_index(nraw), source=0)
    nfree = 0
    do i = 1, nraw
        if (fixed_dof(i)) cycle
        nfree = nfree + 1
        free_index(i) = nfree
    end do
    nactive = 0
    do row = 1, nconstraint
        row_norm = 0.0_rk
        do j = rowptr(row), rowptr(row+1) - 1
            if (free_index(col(j)) > 0) row_norm = max(row_norm, abs(val(j)))
        end do
        if (row_norm > 0.0_rk) nactive = nactive + 1
    end do
    allocate(system(nfree+nactive,nfree+nactive), rhs(nfree+nactive,1), source=0.0_rk)

    nelem = 0
    do ip = 1, domain%get_npatch()
        patch = domain%get_patch(ip)
        Xc = patch%get_Xc()
        elem = patch%cmp_elem()
        offset = offsets(ip)
        nelem = nelem + size(elem,1)
        do is = 1, size(elem,1)
            do ig = 1, ngauss**3
                call patch%ansatz(&
                    ie          = is,&
                    ig          = ig,&
                    Tgc         = N,&
                    dTgc_dXg    = dN,&
                    d2Tgc_dXg2  = d2N,&
                    dV          = dV,&
                    ngauss      = [ngauss, ngauss, ngauss])
                x = 0.0_rk
                y = 0.0_rk
                z = 0.0_rk
                do a = 1, size(elem,2)
                    ia = elem(is,a)
                    x = x + N(a)*Xc(ia,1)
                    y = y + N(a)*Xc(ia,2)
                    z = z + N(a)*Xc(ia,3)
                end do
                load_value = source_term(x, y, z)*dV
                do a = 1, size(elem,2)
                    ia = elem(is,a)
                    row_a = free_index(offset+ia)
                    if (row_a == 0) cycle
                    rhs(row_a,1) = rhs(row_a,1) + N(a)*load_value
                    do b = 1, a
                        ib = elem(is,b)
                        row_b = free_index(offset+ib)
                        if (row_b == 0) cycle
                        stiffness = (&
                            d2N(a,1,1)*d2N(b,1,1) + d2N(a,1,2)*d2N(b,1,2) + d2N(a,1,3)*d2N(b,1,3) +&
                            d2N(a,2,1)*d2N(b,2,1) + d2N(a,2,2)*d2N(b,2,2) + d2N(a,2,3)*d2N(b,2,3) +&
                            d2N(a,3,1)*d2N(b,3,1) + d2N(a,3,2)*d2N(b,3,2) + d2N(a,3,3)*d2N(b,3,3))*dV
                        system(row_a,row_b) = system(row_a,row_b) + stiffness
                        if (row_a /= row_b) system(row_b,row_a) = system(row_b,row_a) + stiffness
                    end do
                end do
            end do
        end do
    end do

    active_row = 0
    do row = 1, nconstraint
        row_norm = 0.0_rk
        do j = rowptr(row), rowptr(row+1) - 1
            if (free_index(col(j)) > 0) row_norm = max(row_norm, abs(val(j)))
        end do
        if (row_norm <= 0.0_rk) cycle
        active_row = active_row + 1
        do j = rowptr(row), rowptr(row+1) - 1
            row_a = free_index(col(j))
            if (row_a == 0) cycle
            system(nfree+active_row,row_a) = system(nfree+active_row,row_a) + val(j)
            system(row_a,nfree+active_row) = system(row_a,nfree+active_row) + val(j)
        end do
    end do

    reduced_solution = solve(system, rhs)
    allocate(coefficient(nraw), source=0.0_rk)
    do i = 1, nraw
        if (free_index(i) > 0) coefficient(i) = reduced_solution(free_index(i),1)
    end do

    constraint_residual = 0.0_rk
    do row = 1, nconstraint
        row_residual = 0.0_rk
        do j = rowptr(row), rowptr(row+1) - 1
            row_residual = row_residual + val(j)*coefficient(col(j))
        end do
        constraint_residual = max(constraint_residual, abs(row_residual))
    end do

    l2_error = 0.0_rk
    do ip = 1, domain%get_npatch()
        patch = domain%get_patch(ip)
        Xc = patch%get_Xc()
        elem = patch%cmp_elem()
        offset = offsets(ip)
        do is = 1, size(elem,1)
            do ig = 1, ngauss**3
                call patch%ansatz(&
                    ie        = is,&
                    ig        = ig,&
                    Tgc       = N,&
                    dTgc_dXg = dN,&
                    dV        = dV,&
                    ngauss    = [ngauss, ngauss, ngauss])
                x = 0.0_rk
                y = 0.0_rk
                z = 0.0_rk
                uh = 0.0_rk
                do a = 1, size(elem,2)
                    ia = elem(is,a)
                    x = x + N(a)*Xc(ia,1)
                    y = y + N(a)*Xc(ia,2)
                    z = z + N(a)*Xc(ia,3)
                    uh = uh + N(a)*coefficient(offset+ia)
                end do
                l2_error = l2_error + (uh-exact_solution(x,y,z))**2*dV
            end do
        end do
    end do

    patch = domain%get_patch(1)
    call patch%basis(&
        Xt  = [1.0_rk, 0.5_rk, 0.5_rk],&
        Tgc = N)
    uh = dot_product(N, coefficient(1:ncp))

    write(*,"(a,l1)") "valid C1 geometry       : ", domain%is_valid()
    write(*,"(a,i0)") "patches                 : ", domain%get_npatch()
    write(*,"(a,i0)") "elements                : ", nelem
    write(*,"(a,i0)") "raw patch dofs          : ", nraw
    write(*,"(a,i0)") "fixed boundary dofs     : ", count(fixed_dof)
    write(*,"(a,i0)") "free field dofs         : ", nfree
    write(*,"(a,i0)") "interface constraints   : ", nconstraint
    write(*,"(a,i0)") "active constraints      : ", nactive
    write(*,"(a,i0)") "reduced system size     : ", nfree+nactive
    write(*,"(a,es12.4)") "constraint residual     : ", constraint_residual
    write(*,"(a,es12.4)") "center value            : ", uh
    write(*,"(a,es12.4)") "exact center value      : ", exact_solution(0.5_rk,0.5_rk,0.5_rk)
    write(*,"(a,es12.4)") "L2 error                : ", sqrt(l2_error)

    do ip = 1, npatch
        patch = domain%get_patch(ip)
        Xc = patch%get_Xc()
        offset = offsets(ip)
        do i = 1, ncp
            call patch%modify_Xc(&
                X   = Xc(i,3) + visual_scale*coefficient(offset+i),&
                num = i,&
                dir = 3)
        end do
        call deformation%add_patch(&
            patch = patch,&
            id    = deformation_id(ip))
    end do
    call deformation%connect(&
        patch_a    = deformation_id(1),&
        side_a     = "u_max",&
        patch_b    = deformation_id(2),&
        side_b     = "u_min",&
        continuity = 1)
    write(*,"(a,l1)") "valid C1 field          : ", deformation%is_valid()
    write(*,"(a,f6.1)") "visualization scale     : ", visual_scale

    call deformation%create(&
        res1 = 17,&
        res2 = 17,&
        res3 = 17)
    call deformation%export_Xc("vtk/volume_iga_c1")
    call deformation%export_Xg("vtk/volume_iga_c1")
    call deformation%export_Xth_in_Xg("vtk/volume_iga_c1", 9)
    call deformation%show(&
        vtkfile_Xc        = "vtk/volume_iga_c1_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_iga_c1_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_iga_c1_Xth_*.vtk")

    call domain%finalize()
    call deformation%finalize()

contains
    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure elemental real(rk) function exact_solution(x, y, z) result(u)
        real(rk), intent(in) :: x, y, z

        u = x*x*(1.0_rk-x)**2*y*y*(1.0_rk-y)**2*z*z*(1.0_rk-z)**2
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure elemental real(rk) function source_term(x, y, z) result(f)
        real(rk), intent(in) :: x, y, z
        real(rk) :: ax, ay, az, d2x, d2y, d2z

        ax = x*x*(1.0_rk-x)**2
        ay = y*y*(1.0_rk-y)**2
        az = z*z*(1.0_rk-z)**2
        d2x = 2.0_rk - 12.0_rk*x + 12.0_rk*x*x
        d2y = 2.0_rk - 12.0_rk*y + 12.0_rk*y*y
        d2z = 2.0_rk - 12.0_rk*z + 12.0_rk*z*z
        f = 24.0_rk*(ay*az + ax*az + ax*ay) +&
            2.0_rk*(d2x*d2y*az + d2x*ay*d2z + ax*d2y*d2z)
    end function
    !===============================================================================
end program volume_iga_c1
