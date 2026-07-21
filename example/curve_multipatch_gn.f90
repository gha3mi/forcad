!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Spatial quartic curve chain with exact G3 and non-C3 patch interfaces.
program curve_multipatch_gn

    use forcad, only: rk, nurbs_curve, nurbs_multipatch_curve

    implicit none

    integer, parameter :: continuity_order = 3, degree = 7, nc = degree + 1
    integer, parameter :: npatch = 6, ncomponent = 3
    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: transition(continuity_order,npatch-1) = reshape([&
        0.86_rk, -0.040_rk,  0.018_rk,&
        1.12_rk,  0.050_rk, -0.020_rk,&
        0.91_rk, -0.035_rk,  0.015_rk,&
        1.08_rk,  0.045_rk, -0.018_rk,&
        0.94_rk, -0.030_rk,  0.012_rk], [continuity_order,npatch-1])
    type(nurbs_curve) :: patch
    type(nurbs_multipatch_curve) :: curve, parametric_curve
    integer :: patch_id(npatch), parametric_id(npatch)
    integer :: ipatch, component, i, row, constraint_set
    integer :: constraint_rows(2), constraint_nonzeros(2)
    integer, allocatable :: offsets(:), dof_map(:), elem(:,:), rowptr(:), col(:)
    real(rk) :: Xc(nc,ncomponent,npatch)
    real(rk) :: derivative_a1, derivative_a2, derivative_a3
    real(rk) :: derivative_b1, derivative_b2, derivative_b3
    real(rk) :: local_coordinate, chain_coordinate, row_residual, constraint_residual(2)
    real(rk), allocatable :: cval(:), dof(:,:)

    ! Degree seven leaves four independent control layers at each patch end.
    Xc = 0.0_rk
    do i = 1, nc
        chain_coordinate = real(i-1,rk)/real(degree,rk)
        Xc(i,1,1) = -4.0_rk + 1.25_rk*chain_coordinate +&
            0.25_rk*sin(0.70_rk*pi*chain_coordinate)
        Xc(i,2,1) = 1.35_rk*sin(0.85_rk*pi*chain_coordinate) +&
            0.30_rk*sin(1.80_rk*pi*chain_coordinate)
        Xc(i,3,1) = 0.95_rk*cos(0.55_rk*pi*chain_coordinate) +&
            0.35_rk*sin(1.35_rk*pi*chain_coordinate)
    end do

    do ipatch = 2, npatch
        do component = 1, ncomponent
            derivative_a1 = real(degree,rk)*(Xc(nc,component,ipatch-1) -&
                Xc(nc-1,component,ipatch-1))
            derivative_a2 = real(degree*(degree-1),rk)*(Xc(nc,component,ipatch-1) -&
                2.0_rk*Xc(nc-1,component,ipatch-1) + Xc(nc-2,component,ipatch-1))
            derivative_a3 = real(degree*(degree-1)*(degree-2),rk)*(&
                Xc(nc,component,ipatch-1) - 3.0_rk*Xc(nc-1,component,ipatch-1) +&
                3.0_rk*Xc(nc-2,component,ipatch-1) - Xc(nc-3,component,ipatch-1))

            derivative_b1 = derivative_a1/transition(1,ipatch-1)
            derivative_b2 = (derivative_a2-transition(2,ipatch-1)*derivative_b1)/&
                transition(1,ipatch-1)**2
            derivative_b3 = (derivative_a3-transition(3,ipatch-1)*derivative_b1 -&
                3.0_rk*transition(1,ipatch-1)*transition(2,ipatch-1)*derivative_b2)/&
                transition(1,ipatch-1)**3

            Xc(1,component,ipatch) = Xc(nc,component,ipatch-1)
            Xc(2,component,ipatch) = Xc(1,component,ipatch) +&
                derivative_b1/real(degree,rk)
            Xc(3,component,ipatch) = derivative_b2/real(degree*(degree-1),rk) +&
                2.0_rk*Xc(2,component,ipatch) - Xc(1,component,ipatch)
            Xc(4,component,ipatch) = derivative_b3/real(degree*(degree-1)*(degree-2),rk) +&
                3.0_rk*Xc(3,component,ipatch) - 3.0_rk*Xc(2,component,ipatch) +&
                Xc(1,component,ipatch)
        end do
        do i = 5, nc
            local_coordinate = real(i-1,rk)/real(degree,rk)
            chain_coordinate = real(ipatch-1,rk) + local_coordinate
            Xc(i,1,ipatch) = -4.0_rk + 1.25_rk*chain_coordinate +&
                0.25_rk*sin(0.70_rk*pi*chain_coordinate)
            Xc(i,2,ipatch) = 1.35_rk*sin(0.85_rk*pi*chain_coordinate) +&
                0.30_rk*sin(1.80_rk*pi*chain_coordinate)
            Xc(i,3,ipatch) = 0.95_rk*cos(0.55_rk*pi*chain_coordinate) +&
                0.35_rk*sin(1.35_rk*pi*chain_coordinate)
        end do
    end do

    do ipatch = 1, npatch
        call patch%set(Xc(:,:,ipatch))
        call curve%add_patch(&
            patch = patch,&
            id    = patch_id(ipatch))
        call parametric_curve%add_patch(&
            patch = patch,&
            id    = parametric_id(ipatch))
    end do
    do ipatch = 1, npatch - 1
        call curve%connect(&
            patch_a            = patch_id(ipatch),&
            side_a             = "right",&
            patch_b            = patch_id(ipatch+1),&
            side_b             = "left",&
            continuity         = continuity_order,&
            geometric          = .true.,&
            reparameterization = transition(:,ipatch))
        call parametric_curve%connect(&
            patch_a         = parametric_id(ipatch),&
            side_a          = "right",&
            patch_b         = parametric_id(ipatch+1),&
            side_b          = "left",&
            continuity      = continuity_order,&
            geometric          = .false.)
    end do

    offsets = curve%cmp_dof_offsets()
    dof_map = curve%cmp_dof_map()
    elem = curve%cmp_elem(shared=.true.)
    allocate(dof(offsets(size(offsets)),ncomponent), source=0.0_rk)
    do ipatch = 1, npatch
        do i = 1, nc
            do component = 1, ncomponent
                dof(offsets(ipatch)+i,component) = Xc(i,component,ipatch)
            end do
        end do
    end do

    do constraint_set = 1, 2
        if (constraint_set == 1) then
            call curve%cmp_dof_constraint(&
                rowptr          = rowptr,&
                col             = col,&
                val             = cval,&
                geometric          = .true.)
        else
            call parametric_curve%cmp_dof_constraint(&
                rowptr          = rowptr,&
                col             = col,&
                val             = cval,&
                geometric          = .false.)
        end if
        constraint_rows(constraint_set) = size(rowptr) - 1
        constraint_nonzeros(constraint_set) = size(cval)
        constraint_residual(constraint_set) = 0.0_rk
        do row = 1, constraint_rows(constraint_set)
            do component = 1, ncomponent
                row_residual = 0.0_rk
                do i = rowptr(row), rowptr(row+1) - 1
                    row_residual = row_residual + cval(i)*dof(col(i),component)
                end do
                constraint_residual(constraint_set) = max(&
                    constraint_residual(constraint_set),&
                    abs(row_residual))
            end do
        end do
    end do

    write(*,"(a,i0)") "continuity order n      : ", continuity_order
    write(*,"(a,i0)") "curve patches           : ", curve%get_npatch()
    write(*,"(a,i0)") "connections             : ", curve%get_nconnection()
    write(*,"(a,5(f7.3,1x))") "transition speeds       : ", transition(1,:)
    write(*,"(a,l1)") "valid G3 geometry       : ", curve%is_valid()
    write(*,"(a,l1)") "same patches valid C3   : ", parametric_curve%is_valid()
    write(*,"(a,es12.4)") "G3 constraint residual  : ", constraint_residual(1)
    write(*,"(a,es12.4)") "C3 constraint residual  : ", constraint_residual(2)
    write(*,"(a,i0)") "G3 constraint rows      : ", constraint_rows(1)
    write(*,"(a,i0)") "G3 constraint nonzeros  : ", constraint_nonzeros(1)
    write(*,"(a,i0)") "raw patch dofs          : ", size(dof_map)
    write(*,"(a,i0)") "shared global dofs      : ", maxval(dof_map)
    write(*,"(a,i0)") "elements                : ", size(elem,1)
    write(*,"(a,i0)") "nodes per element       : ", size(elem,2)

    call curve%create(res=81)
    call curve%export_Xc("vtk/curve_multipatch_gn")
    call curve%export_Xg("vtk/curve_multipatch_gn")
    call curve%export_Xth_in_Xg("vtk/curve_multipatch_gn", 41)
    call curve%show(&
        vtkfile_Xc        = "vtk/curve_multipatch_gn_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/curve_multipatch_gn_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/curve_multipatch_gn_Xth_*.vtk")

    call curve%finalize()
    call parametric_curve%finalize()

end program curve_multipatch_gn
