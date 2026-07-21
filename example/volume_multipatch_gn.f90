!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Rational twisted solid with exact G3 and non-C3 patch interfaces.
program volume_multipatch_gn

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: continuity_order = 3, nu = 8, nv = 4, nw = 4
    integer, parameter :: degree_u = nu - 1, npatch = 4, nlocal = nu*nv*nw
    integer, parameter :: nhomogeneous = 4
    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: transition(continuity_order,npatch-1) = reshape([&
        0.88_rk, -0.050_rk,  0.022_rk,&
        1.10_rk,  0.045_rk, -0.018_rk,&
        0.93_rk, -0.035_rk,  0.015_rk], [continuity_order,npatch-1])

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: solid, parametric_solid
    integer :: patch_id(npatch), parametric_id(npatch)
    integer :: ipatch, i, j, k, d, n, base, first_layer, row, constraint_set
    integer :: constraint_rows(2), constraint_nonzeros(2)
    integer, allocatable :: offsets(:), dof_map(:), elem(:,:), rowptr(:), col(:)
    real(rk) :: Hc(nlocal,nhomogeneous,npatch), Xc(nlocal,3), Wc(nlocal)
    real(rk) :: local_coordinate, chain_coordinate, eta, zeta, x, y, z, weight
    real(rk) :: center_x, center_y, center_z, twist, radius_y, radius_z, local_y, local_z
    real(rk) :: derivative_a1, derivative_a2, derivative_a3
    real(rk) :: derivative_b1, derivative_b2, derivative_b3
    real(rk) :: row_residual, constraint_residual(2)
    real(rk), allocatable :: cval(:), dof(:,:)

    ! Degree seven leaves four independent control layers at each patch end.
    Hc = 0.0_rk
    do ipatch = 1, npatch
        first_layer = 5
        if (ipatch == 1) first_layer = 1
        do k = 1, nw
            zeta = 2.0_rk*real(k-1,rk)/real(nw-1,rk) - 1.0_rk
            do j = 1, nv
                eta = 2.0_rk*real(j-1,rk)/real(nv-1,rk) - 1.0_rk
                do i = first_layer, nu
                    local_coordinate = real(i-1,rk)/real(degree_u,rk)
                    chain_coordinate = real(ipatch-1,rk) + local_coordinate
                    n = i + (j-1)*nu + (k-1)*nu*nv
                    center_x = -4.0_rk + 1.80_rk*chain_coordinate
                    center_y = 0.62_rk*sin(0.55_rk*pi*chain_coordinate) +&
                        0.18_rk*sin(1.25_rk*pi*chain_coordinate)
                    center_z = 0.45_rk*cos(0.45_rk*pi*chain_coordinate) +&
                        0.14_rk*sin(0.95_rk*pi*chain_coordinate)
                    twist = 0.62_rk*pi*chain_coordinate +&
                        0.16_rk*sin(0.80_rk*pi*chain_coordinate)
                    radius_y = 0.90_rk*(1.0_rk-0.025_rk*chain_coordinate +&
                        0.08_rk*sin(0.50_rk*pi*chain_coordinate))
                    radius_z = 0.62_rk*(1.0_rk+0.10_rk*cos(0.55_rk*pi*chain_coordinate))
                    local_y = radius_y*eta*(1.0_rk+0.10_rk*zeta)
                    local_z = radius_z*zeta*(1.0_rk-0.08_rk*eta*eta)
                    x = center_x + 0.12_rk*eta*zeta*sin(0.60_rk*pi*chain_coordinate) +&
                        0.05_rk*(eta*eta-zeta*zeta)*sin(1.10_rk*pi*chain_coordinate)
                    y = center_y + cos(twist)*local_y - sin(twist)*local_z
                    z = center_z + sin(twist)*local_y + cos(twist)*local_z +&
                        0.06_rk*eta*zeta*cos(0.45_rk*pi*chain_coordinate)
                    weight = 1.0_rk + 0.10_rk*(1.0_rk-eta*eta)*(1.0_rk-zeta*zeta)*(&
                        0.82_rk+0.18_rk*cos(0.50_rk*pi*chain_coordinate))
                    Hc(n,1,ipatch) = weight*x
                    Hc(n,2,ipatch) = weight*y
                    Hc(n,3,ipatch) = weight*z
                    Hc(n,4,ipatch) = weight
                end do
            end do
        end do
    end do

    do ipatch = 2, npatch
        do k = 1, nw
            do j = 1, nv
                base = (j-1)*nu + (k-1)*nu*nv
                do d = 1, nhomogeneous
                    derivative_a1 = real(degree_u,rk)*(Hc(base+nu,d,ipatch-1) -&
                        Hc(base+nu-1,d,ipatch-1))
                    derivative_a2 = real(degree_u*(degree_u-1),rk)*(&
                        Hc(base+nu,d,ipatch-1) - 2.0_rk*Hc(base+nu-1,d,ipatch-1) +&
                        Hc(base+nu-2,d,ipatch-1))
                    derivative_a3 = real(degree_u*(degree_u-1)*(degree_u-2),rk)*(&
                        Hc(base+nu,d,ipatch-1) - 3.0_rk*Hc(base+nu-1,d,ipatch-1) +&
                        3.0_rk*Hc(base+nu-2,d,ipatch-1) - Hc(base+nu-3,d,ipatch-1))

                    derivative_b1 = derivative_a1/transition(1,ipatch-1)
                    derivative_b2 = (derivative_a2-transition(2,ipatch-1)*derivative_b1)/&
                        transition(1,ipatch-1)**2
                    derivative_b3 = (derivative_a3-transition(3,ipatch-1)*derivative_b1 -&
                        3.0_rk*transition(1,ipatch-1)*transition(2,ipatch-1)*derivative_b2)/&
                        transition(1,ipatch-1)**3

                    Hc(base+1,d,ipatch) = Hc(base+nu,d,ipatch-1)
                    Hc(base+2,d,ipatch) = Hc(base+1,d,ipatch) +&
                        derivative_b1/real(degree_u,rk)
                    Hc(base+3,d,ipatch) = derivative_b2/real(degree_u*(degree_u-1),rk) +&
                        2.0_rk*Hc(base+2,d,ipatch) - Hc(base+1,d,ipatch)
                    Hc(base+4,d,ipatch) = derivative_b3/&
                        real(degree_u*(degree_u-1)*(degree_u-2),rk) +&
                        3.0_rk*Hc(base+3,d,ipatch) - 3.0_rk*Hc(base+2,d,ipatch) +&
                        Hc(base+1,d,ipatch)
                end do
            end do
        end do
    end do

    do ipatch = 1, npatch
        do n = 1, nlocal
            Wc(n) = Hc(n,4,ipatch)
            do d = 1, 3
                Xc(n,d) = Hc(n,d,ipatch)/Wc(n)
            end do
        end do
        call patch%set(&
            nc = [nu, nv, nw],&
            Xc = Xc,&
            Wc = Wc)
        call solid%add_patch(&
            patch = patch,&
            id    = patch_id(ipatch))
        call parametric_solid%add_patch(&
            patch = patch,&
            id    = parametric_id(ipatch))
    end do
    do ipatch = 1, npatch - 1
        call solid%connect(&
            patch_a            = patch_id(ipatch),&
            side_a             = "u_max",&
            patch_b            = patch_id(ipatch+1),&
            side_b             = "u_min",&
            continuity         = continuity_order,&
            geometric          = .true.,&
            reparameterization = transition(:,ipatch))
        call parametric_solid%connect(&
            patch_a         = parametric_id(ipatch),&
            side_a          = "u_max",&
            patch_b         = parametric_id(ipatch+1),&
            side_b          = "u_min",&
            continuity      = continuity_order,&
            geometric          = .false.)
    end do

    offsets = solid%cmp_dof_offsets()
    dof_map = solid%cmp_dof_map()
    elem = solid%cmp_elem(shared=.true.)
    allocate(dof(offsets(size(offsets)),nhomogeneous), source=0.0_rk)
    do ipatch = 1, npatch
        do n = 1, nlocal
            do d = 1, nhomogeneous
                dof(offsets(ipatch)+n,d) = Hc(n,d,ipatch)
            end do
        end do
    end do

    do constraint_set = 1, 2
        if (constraint_set == 1) then
            call solid%cmp_dof_constraint(&
                rowptr          = rowptr,&
                col             = col,&
                val             = cval,&
                geometric          = .true.)
        else
            call parametric_solid%cmp_dof_constraint(&
                rowptr          = rowptr,&
                col             = col,&
                val             = cval,&
                geometric          = .false.)
        end if
        constraint_rows(constraint_set) = size(rowptr) - 1
        constraint_nonzeros(constraint_set) = size(cval)
        constraint_residual(constraint_set) = 0.0_rk
        do row = 1, constraint_rows(constraint_set)
            do d = 1, nhomogeneous
                row_residual = 0.0_rk
                do i = rowptr(row), rowptr(row+1) - 1
                    row_residual = row_residual + cval(i)*dof(col(i),d)
                end do
                constraint_residual(constraint_set) = max(&
                    constraint_residual(constraint_set),&
                    abs(row_residual))
            end do
        end do
    end do

    write(*,"(a,i0)") "continuity order n      : ", continuity_order
    write(*,"(a,i0)") "volume patches          : ", solid%get_npatch()
    write(*,"(a,i0)") "connections             : ", solid%get_nconnection()
    write(*,"(a,3(f7.3,1x))") "transition speeds       : ", transition(1,:)
    write(*,"(a,l1)") "rational geometry       : ", solid%is_rational()
    write(*,"(a,2(f8.4,1x))") "control weight range    : ", minval(Hc(:,4,:)), maxval(Hc(:,4,:))
    write(*,"(a,l1)") "valid G3 geometry       : ", solid%is_valid()
    write(*,"(a,l1)") "same patches valid C3   : ", parametric_solid%is_valid()
    write(*,"(a,es12.4)") "G3 constraint residual  : ", constraint_residual(1)
    write(*,"(a,es12.4)") "C3 constraint residual  : ", constraint_residual(2)
    write(*,"(a,i0)") "G3 constraint rows      : ", constraint_rows(1)
    write(*,"(a,i0)") "G3 constraint nonzeros  : ", constraint_nonzeros(1)
    write(*,"(a,i0)") "raw patch dofs          : ", size(dof_map)
    write(*,"(a,i0)") "shared global dofs      : ", maxval(dof_map)
    write(*,"(a,i0)") "elements                : ", size(elem,1)
    write(*,"(a,i0)") "nodes per element       : ", size(elem,2)

    call solid%create(&
        res1 = 21,&
        res2 = 9,&
        res3 = 9)
    call solid%export_Xc("vtk/volume_multipatch_gn")
    call solid%export_Xg("vtk/volume_multipatch_gn")
    call solid%export_Xth_in_Xg("vtk/volume_multipatch_gn", 11)
    call solid%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_gn_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_gn_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_gn_Xth_*.vtk")

    call solid%finalize()
    call parametric_solid%finalize()

end program volume_multipatch_gn
