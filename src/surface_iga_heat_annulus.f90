!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Solve steady heat conduction with uniform generation in a rational C-ring.
!! The radial Dirichlet data and insulated cut faces admit an analytical
!! solution, providing a geometry-exact IGA verification problem.
program surface_iga_heat_annulus

    use forcad, only: rk, nurbs_surface, solve

    implicit none

    integer, parameter :: quadrature(2) = [3, 4]
    real(rk), parameter :: radius_inner = 0.45_rk
    real(rk), parameter :: radius_outer = 1.20_rk
    real(rk), parameter :: conductivity = 18.0_rk
    real(rk), parameter :: heat_source = 2.5e4_rk
    real(rk), parameter :: temperature_inner = 420.0_rk
    real(rk), parameter :: temperature_outer = 300.0_rk

    type(nurbs_surface) :: shield
    integer :: nc(2), ncontrol, nelement, nfixed
    integer :: ie, ig, a, b, i
    integer, allocatable :: element(:,:), fixed(:)
    real(rk), allocatable :: Xc(:,:), Xg(:,:), basis(:), gradient(:,:)
    real(rk), allocatable :: stiffness(:,:), rhs(:,:), temperature(:,:)
    real(rk), allocatable :: fixed_value(:), basis_plot(:,:), plot_field(:,:)
    real(rk) :: Xq(3), dA, radius, exact, numerical
    real(rk) :: coefficient_log, coefficient_constant
    real(rk) :: area, error_l2, temperature_scale

    call shield%set_C(&
        center  = [0.0_rk, 0.0_rk, 0.0_rk],&
        radius1 = radius_inner,&
        radius2 = radius_outer)
    call shield%elevate_degree(&
        dir = 2,&
        t   = 2)
    call shield%insert_knots(&
        dir = 1,&
        Xth = [0.25_rk, 0.75_rk],&
        r   = [1, 1])
    call shield%insert_knots(&
        dir = 2,&
        Xth = [0.20_rk, 0.40_rk, 0.60_rk, 0.80_rk],&
        r   = [1, 1, 1, 1])

    nc = shield%get_nc()
    Xc = shield%get_Xc()
    element = shield%cmp_elem()
    ncontrol = product(nc)
    nelement = size(element, 1)
    allocate(stiffness(ncontrol,ncontrol), rhs(ncontrol,1), source=0.0_rk)

    area = 0.0_rk
    do ie = 1, nelement
        do ig = 1, product(quadrature)
            call shield%ansatz(&
                ie       = ie,&
                ig       = ig,&
                Tgc      = basis,&
                dTgc_dXg = gradient,&
                dA       = dA,&
                ngauss   = quadrature)
            area = area + dA
            do a = 1, size(element, 2)
                rhs(element(ie,a),1) = rhs(element(ie,a),1) + heat_source*basis(a)*dA
            end do
            do b = 1, size(element, 2)
                do a = 1, size(element, 2)
                    stiffness(element(ie,a),element(ie,b)) = stiffness(element(ie,a),element(ie,b)) + &
                        conductivity*(gradient(a,1)*gradient(b,1) + gradient(a,2)*gradient(b,2) + &
                        gradient(a,3)*gradient(b,3))*dA
                end do
            end do
        end do
    end do

    nfixed = 2*nc(1)
    allocate(fixed(nfixed), fixed_value(nfixed))
    do i = 1, nc(1)
        fixed(i) = i
        fixed_value(i) = temperature_inner
        fixed(nc(1)+i) = (nc(2)-1)*nc(1) + i
        fixed_value(nc(1)+i) = temperature_outer
    end do

    do i = 1, nfixed
        do a = 1, ncontrol
            rhs(a,1) = rhs(a,1) - stiffness(a,fixed(i))*fixed_value(i)
        end do
    end do
    do i = 1, nfixed
        stiffness(fixed(i),:) = 0.0_rk
        stiffness(:,fixed(i)) = 0.0_rk
        stiffness(fixed(i),fixed(i)) = 1.0_rk
        rhs(fixed(i),1) = fixed_value(i)
    end do
    temperature = solve(stiffness, rhs)

    temperature_scale = log(radius_outer/radius_inner)
    coefficient_log = (temperature_outer-temperature_inner + &
        heat_source*(radius_outer**2-radius_inner**2)/(4.0_rk*conductivity))/temperature_scale
    coefficient_constant = temperature_inner + heat_source*radius_inner**2/(4.0_rk*conductivity) - &
        coefficient_log*log(radius_inner)

    error_l2 = 0.0_rk
    do ie = 1, nelement
        do ig = 1, product(quadrature)
            call shield%ansatz(&
                ie       = ie,&
                ig       = ig,&
                Tgc      = basis,&
                dTgc_dXg = gradient,&
                dA       = dA,&
                ngauss   = quadrature)
            Xq = 0.0_rk
            numerical = 0.0_rk
            do a = 1, size(element, 2)
                Xq(1) = Xq(1) + basis(a)*Xc(element(ie,a),1)
                Xq(2) = Xq(2) + basis(a)*Xc(element(ie,a),2)
                Xq(3) = Xq(3) + basis(a)*Xc(element(ie,a),3)
                numerical = numerical + basis(a)*temperature(element(ie,a),1)
            end do
            radius = sqrt(Xq(1)**2 + Xq(2)**2)
            exact = -heat_source*radius**2/(4.0_rk*conductivity) + coefficient_log*log(radius) + &
                coefficient_constant
            error_l2 = error_l2 + (numerical-exact)**2*dA
        end do
    end do
    error_l2 = sqrt(error_l2/area)

    call shield%create(&
        res1 = 81,&
        res2 = 31)
    call shield%basis(&
        Tgc = basis_plot)
    Xg = shield%get_Xg()
    allocate(plot_field(size(Xg,1),3))
    plot_field(:,1) = matmul(basis_plot,temperature(:,1))
    do i = 1, size(Xg,1)
        radius = sqrt(Xg(i,1)**2 + Xg(i,2)**2)
        plot_field(i,2) = -heat_source*radius**2/(4.0_rk*conductivity) + coefficient_log*log(radius) + &
            coefficient_constant
        plot_field(i,3) = abs(plot_field(i,1)-plot_field(i,2))
    end do

    write(*,"(a,l1)") "Valid rational geometry : ", shield%err%ok
    write(*,"(a,2(i0,1x))") "Control net             : ", nc
    write(*,"(a,i0)") "Elements                : ", nelement
    write(*,"(a,es12.4)") "Integrated area [m2]    : ", area
    write(*,"(a,f10.4)") "Minimum temperature [K] : ", minval(plot_field(:,1))
    write(*,"(a,f10.4)") "Maximum temperature [K] : ", maxval(plot_field(:,1))
    write(*,"(a,es12.4)") "RMS error [K]           : ", error_l2

    call shield%export_Xc(&
        filename    = "vtk/surface_iga_heat_annulus_Xc.vtk",&
        point_data  = temperature,&
        field_names = [character(len=24) :: "temperature_coefficient"])
    call shield%export_Xg(&
        filename    = "vtk/surface_iga_heat_annulus_Xg.vtk",&
        point_data  = plot_field,&
        field_names = [character(len=24) :: "temperature", "exact_temperature", "absolute_error"])
    call shield%export_Xth_in_Xg(&
        filename = "vtk/surface_iga_heat_annulus_Xth.vtk",&
        res      = 17)
    call shield%show(&
        vtkfile_Xc        = "vtk/surface_iga_heat_annulus_Xc.vtk",&
        vtkfile_Xg        = "vtk/surface_iga_heat_annulus_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_iga_heat_annulus_Xth.vtk")
    call shield%finalize()

end program surface_iga_heat_annulus
