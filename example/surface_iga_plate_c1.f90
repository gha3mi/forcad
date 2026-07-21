!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Solve a clamped biharmonic plate with two C1 IGA surface patches.
!! The weak form is integral(Hessian(u):Hessian(v)) dA = integral(f*v) dA.
!! A quartic manufactured solution belongs exactly to the degree-four space.
program surface_iga_plate_c1

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface, solve

    implicit none

    integer, parameter :: degree = 4, nc = 7, ngauss = 5, npatch = 2
    real(rk), parameter :: knot(nc+degree+1) = [&
        0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.35_rk, 0.72_rk,&
        1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
    real(rk), parameter :: visual_scale = 20.0_rk

    type(nurbs_surface) :: patch
    type(nurbs_multipatch_surface) :: plate, deflection
    integer :: patch_id(npatch), deflection_id(npatch)
    integer :: ip, is, ig, i, j, a, b, ia, ib, offset, row
    integer :: nraw, nconstraint, nactive, nelem
    integer, allocatable :: offsets(:), rowptr(:), col(:), elem(:,:)
    logical, allocatable :: fixed_dof(:)
    real(rk), allocatable :: val(:), system(:,:), rhs(:,:), solution(:,:)
    real(rk), allocatable :: Xc(:,:), N(:), dN(:,:), d2N(:,:,:)
    real(rk) :: Xplot(nc*nc,3)
    real(rk) :: dA, x, y, uh, load_value, stiffness, row_norm
    real(rk) :: constraint_residual, row_residual, l2_error

    call set_plane_patch(&
        surface = patch,&
        x0      = 0.0_rk,&
        x1      = 0.5_rk)
    call plate%add_patch(&
        patch = patch,&
        id    = patch_id(1))
    call set_plane_patch(&
        surface = patch,&
        x0      = 0.5_rk,&
        x1      = 1.0_rk)
    call plate%add_patch(&
        patch = patch,&
        id    = patch_id(2))
    call plate%connect(&
        patch_a    = patch_id(1),&
        side_a     = "u_max",&
        patch_b    = patch_id(2),&
        side_b     = "u_min",&
        continuity = 1)

    offsets = plate%cmp_dof_offsets()
    nraw = offsets(size(offsets))
    call plate%cmp_dof_constraint(&
        rowptr = rowptr,&
        col    = col,&
        val    = val)
    nconstraint = size(rowptr) - 1
    allocate(system(nraw+nconstraint,nraw+nconstraint), rhs(nraw+nconstraint,1), source=0.0_rk)

    nelem = 0
    do ip = 1, plate%get_npatch()
        patch = plate%get_patch(ip)
        Xc = patch%get_Xc()
        elem = patch%cmp_elem()
        offset = offsets(ip)
        nelem = nelem + size(elem,1)
        do is = 1, size(elem,1)
            do ig = 1, ngauss*ngauss
                call patch%ansatz(&
                    ie          = is,&
                    ig          = ig,&
                    Tgc         = N,&
                    dTgc_dXg    = dN,&
                    d2Tgc_dXg2  = d2N,&
                    dA          = dA,&
                    ngauss      = [ngauss, ngauss])
                x = 0.0_rk
                y = 0.0_rk
                do a = 1, size(elem,2)
                    ia = elem(is,a)
                    x = x + N(a)*Xc(ia,1)
                    y = y + N(a)*Xc(ia,2)
                end do
                load_value = source_term(x, y)*dA
                do a = 1, size(elem,2)
                    ia = elem(is,a)
                    rhs(offset+ia,1) = rhs(offset+ia,1) + N(a)*load_value
                    do b = 1, a
                        ib = elem(is,b)
                        stiffness = (&
                            d2N(a,1,1)*d2N(b,1,1) + d2N(a,1,2)*d2N(b,1,2) +&
                            d2N(a,2,1)*d2N(b,2,1) + d2N(a,2,2)*d2N(b,2,2))*dA
                        system(offset+ia,offset+ib) = system(offset+ia,offset+ib) + stiffness
                        if (a /= b) system(offset+ib,offset+ia) = system(offset+ib,offset+ia) + stiffness
                    end do
                end do
            end do
        end do
    end do

    do row = 1, nconstraint
        do j = rowptr(row), rowptr(row+1) - 1
            system(nraw+row,col(j)) = val(j)
            system(col(j),nraw+row) = val(j)
        end do
    end do

    allocate(fixed_dof(nraw), source=.false.)
    do j = 1, nc
        fixed_dof(offsets(1)+(j-1)*nc+1) = .true.
        fixed_dof(offsets(1)+(j-1)*nc+2) = .true.
        fixed_dof(offsets(2)+(j-1)*nc+nc-1) = .true.
        fixed_dof(offsets(2)+(j-1)*nc+nc) = .true.
    end do
    do ip = 1, npatch
        offset = offsets(ip)
        do i = 1, nc
            fixed_dof(offset+i) = .true.
            fixed_dof(offset+nc+i) = .true.
            fixed_dof(offset+(nc-2)*nc+i) = .true.
            fixed_dof(offset+(nc-1)*nc+i) = .true.
        end do
    end do
    do i = 1, nraw
        if (.not. fixed_dof(i)) cycle
        system(i,:) = 0.0_rk
        system(:,i) = 0.0_rk
        system(i,i) = 1.0_rk
        rhs(i,1) = 0.0_rk
    end do

    nactive = 0
    do row = 1, nconstraint
        row_norm = 0.0_rk
        do j = rowptr(row), rowptr(row+1) - 1
            if (.not. fixed_dof(col(j))) row_norm = max(row_norm, abs(val(j)))
        end do
        if (row_norm > 0.0_rk) then
            nactive = nactive + 1
        else
            system(nraw+row,nraw+row) = 1.0_rk
        end if
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

    l2_error = 0.0_rk
    do ip = 1, plate%get_npatch()
        patch = plate%get_patch(ip)
        Xc = patch%get_Xc()
        elem = patch%cmp_elem()
        offset = offsets(ip)
        do is = 1, size(elem,1)
            do ig = 1, ngauss*ngauss
                call patch%ansatz(&
                    ie        = is,&
                    ig        = ig,&
                    Tgc       = N,&
                    dTgc_dXg = dN,&
                    dA        = dA,&
                    ngauss    = [ngauss, ngauss])
                x = 0.0_rk
                y = 0.0_rk
                uh = 0.0_rk
                do a = 1, size(elem,2)
                    ia = elem(is,a)
                    x = x + N(a)*Xc(ia,1)
                    y = y + N(a)*Xc(ia,2)
                    uh = uh + N(a)*solution(offset+ia,1)
                end do
                l2_error = l2_error + (uh-exact_solution(x,y))**2*dA
            end do
        end do
    end do

    patch = plate%get_patch(1)
    call patch%basis(&
        Xt  = [1.0_rk, 0.5_rk],&
        Tgc = N)
    uh = dot_product(N, solution(1:nc*nc,1))

    write(*,"(a,l1)") "valid C1 geometry       : ", plate%is_valid()
    write(*,"(a,i0)") "patches                 : ", plate%get_npatch()
    write(*,"(a,i0)") "elements                : ", nelem
    write(*,"(a,i0)") "raw patch dofs          : ", nraw
    write(*,"(a,i0)") "fixed boundary dofs     : ", count(fixed_dof)
    write(*,"(a,i0)") "interface constraints   : ", nconstraint
    write(*,"(a,i0)") "active constraints      : ", nactive
    write(*,"(a,es12.4)") "constraint residual     : ", constraint_residual
    write(*,"(a,es12.4)") "center displacement     : ", uh
    write(*,"(a,es12.4)") "exact center value      : ", exact_solution(0.5_rk,0.5_rk)
    write(*,"(a,es12.4)") "L2 error                : ", sqrt(l2_error)

    do ip = 1, npatch
        patch = plate%get_patch(ip)
        Xc = patch%get_Xc()
        offset = offsets(ip)
        do concurrent (i = 1:nc*nc)
            Xplot(i,1) = Xc(i,1)
            Xplot(i,2) = Xc(i,2)
            Xplot(i,3) = visual_scale*solution(offset+i,1)
        end do
        call patch%set(&
            knot1  = knot,&
            knot2  = knot,&
            Xc     = Xplot,&
            degree = [degree, degree])
        call deflection%add_patch(&
            patch = patch,&
            id    = deflection_id(ip))
    end do
    call deflection%connect(&
        patch_a    = deflection_id(1),&
        side_a     = "u_max",&
        patch_b    = deflection_id(2),&
        side_b     = "u_min",&
        continuity = 1)
    write(*,"(a,l1)") "valid C1 solution       : ", deflection%is_valid()
    write(*,"(a,f6.1)") "visualization scale     : ", visual_scale

    call deflection%create(&
        res1 = 41,&
        res2 = 41)
    call deflection%export_Xc("vtk/surface_iga_plate_c1")
    call deflection%export_Xg("vtk/surface_iga_plate_c1")
    call deflection%export_Xth_in_Xg("vtk/surface_iga_plate_c1", 17)
    call deflection%show(&
        vtkfile_Xc        = "vtk/surface_iga_plate_c1_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_iga_plate_c1_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_iga_plate_c1_Xth_*.vtk")

    call plate%finalize()
    call deflection%finalize()

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_plane_patch(surface, x0, x1)
        type(nurbs_surface), intent(inout) :: surface
        real(rk), intent(in) :: x0, x1
        real(rk) :: control_points(nc*nc,2)
        integer :: i, j

        do concurrent (j = 1:nc, i = 1:nc)
            control_points((j-1)*nc+i,1) = x0 + (x1-x0)*sum(knot(i+1:i+degree))/real(degree,rk)
            control_points((j-1)*nc+i,2) = sum(knot(j+1:j+degree))/real(degree,rk)
        end do
        call surface%set(&
            knot1  = knot,&
            knot2  = knot,&
            Xc     = control_points,&
            degree = [degree, degree])
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure elemental real(rk) function exact_solution(x, y) result(u)
        real(rk), intent(in) :: x, y

        u = x*x*(1.0_rk-x)**2*y*y*(1.0_rk-y)**2
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure elemental real(rk) function source_term(x, y) result(f)
        real(rk), intent(in) :: x, y

        f = 24.0_rk*y*y*(1.0_rk-y)**2 + 24.0_rk*x*x*(1.0_rk-x)**2 +&
            2.0_rk*(2.0_rk-12.0_rk*x+12.0_rk*x*x)*(2.0_rk-12.0_rk*y+12.0_rk*y*y)
    end function
    !===============================================================================
end program surface_iga_plate_c1
