!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Simulate transient heat diffusion in a twisted tapered cooling fin.
!! Backward Euler time integration uses a consistent heat-capacity matrix,
!! internal generation, an isothermal base, and insulated remaining faces.
program volume_iga_heat_transient

    use forcad, only: rk, nurbs_volume, solve

    implicit none

    integer, parameter :: nc(3) = [4, 4, 6]
    integer, parameter :: degree(3) = [2, 2, 3]
    integer, parameter :: quadrature(3) = [3, 3, 4]
    integer, parameter :: nstep = 12
    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: fin_length = 0.30_rk
    real(rk), parameter :: width_root = 0.10_rk
    real(rk), parameter :: depth_root = 0.045_rk
    real(rk), parameter :: twist_angle = 0.75_rk*pi
    real(rk), parameter :: conductivity = 16.0_rk
    real(rk), parameter :: density = 7800.0_rk
    real(rk), parameter :: heat_capacity = 500.0_rk
    real(rk), parameter :: source_scale = 5.0e4_rk
    real(rk), parameter :: coolant_temperature = 293.15_rk
    real(rk), parameter :: time_step = 600.0_rk

    type(nurbs_volume) :: fin
    integer :: i, j, k, n, ie, ig, a, b, id, ncontrol, nfixed, step
    integer, allocatable :: element(:,:), fixed(:)
    real(rk), allocatable :: Xc(:,:), Xg(:,:), basis(:), gradient(:,:)
    real(rk), allocatable :: mass(:,:), stiffness(:,:), system(:,:)
    real(rk), allocatable :: source(:), rhs(:,:), solution(:,:), temperature(:)
    real(rk), allocatable :: boundary_correction(:), capacity_weight(:)
    real(rk), allocatable :: basis_plot(:,:), plot_field(:,:)
    real(rk) :: s, xi, eta, angle, width, depth, local_x, local_y
    real(rk) :: center_x, center_y, dV, source_density, Xq(3)
    real(rk) :: elapsed_time, stored_energy, generated_energy, removed_energy

    ncontrol = product(nc)
    allocate(Xc(ncontrol,3))
    do k = 1, nc(3)
        s = real(k-1,rk)/real(nc(3)-1,rk)
        angle = twist_angle*s
        width = width_root*(1.0_rk-0.35_rk*s)
        depth = depth_root*(1.0_rk-0.20_rk*s)
        center_x = 0.018_rk*sin(1.2_rk*pi*s)
        center_y = 0.025_rk*s*s
        do j = 1, nc(2)
            eta = real(j-1,rk)/real(nc(2)-1,rk) - 0.5_rk
            local_y = depth*eta
            do i = 1, nc(1)
                xi = real(i-1,rk)/real(nc(1)-1,rk) - 0.5_rk
                local_x = width*xi
                n = i + (j-1)*nc(1) + (k-1)*nc(1)*nc(2)
                Xc(n,1) = center_x + cos(angle)*local_x - sin(angle)*local_y
                Xc(n,2) = center_y + sin(angle)*local_x + cos(angle)*local_y
                Xc(n,3) = fin_length*s
            end do
        end do
    end do

    call fin%set(&
        degree = degree,&
        nc     = nc,&
        Xc     = Xc)
    element = fin%cmp_elem()
    allocate(&
        mass(ncontrol,ncontrol),&
        stiffness(ncontrol,ncontrol),&
        source(ncontrol),&
        source=0.0_rk)

    do ie = 1, size(element,1)
        do ig = 1, product(quadrature)
            call fin%ansatz(&
                ie       = ie,&
                ig       = ig,&
                Tgc      = basis,&
                dTgc_dXg = gradient,&
                dV       = dV,&
                ngauss   = quadrature)
            Xq = 0.0_rk
            do a = 1, size(element,2)
                Xq(1) = Xq(1) + basis(a)*Xc(element(ie,a),1)
                Xq(2) = Xq(2) + basis(a)*Xc(element(ie,a),2)
                Xq(3) = Xq(3) + basis(a)*Xc(element(ie,a),3)
            end do
            source_density = source_scale*(0.65_rk+0.70_rk*(Xq(3)/fin_length)**2)
            do a = 1, size(element,2)
                id = element(ie,a)
                source(id) = source(id) + source_density*basis(a)*dV
            end do
            do b = 1, size(element,2)
                do a = 1, size(element,2)
                    id = element(ie,a)
                    mass(id,element(ie,b)) = mass(id,element(ie,b)) + &
                        density*heat_capacity*basis(a)*basis(b)*dV
                    stiffness(id,element(ie,b)) = stiffness(id,element(ie,b)) + &
                        conductivity*(gradient(a,1)*gradient(b,1) + gradient(a,2)*gradient(b,2) + &
                        gradient(a,3)*gradient(b,3))*dV
                end do
            end do
        end do
    end do

    system = stiffness + mass/time_step
    nfixed = nc(1)*nc(2)
    allocate(&
        fixed(nfixed),&
        boundary_correction(ncontrol),&
        capacity_weight(ncontrol),&
        temperature(ncontrol))
    do i = 1, nfixed
        fixed(i) = i
    end do
    boundary_correction = 0.0_rk
    capacity_weight = 0.0_rk
    do b = 1, ncontrol
        do a = 1, ncontrol
            capacity_weight(a) = capacity_weight(a) + mass(a,b)
        end do
    end do
    do i = 1, nfixed
        do a = 1, ncontrol
            boundary_correction(a) = boundary_correction(a) + system(a,fixed(i))*coolant_temperature
        end do
    end do
    do i = 1, nfixed
        system(fixed(i),:) = 0.0_rk
        system(:,fixed(i)) = 0.0_rk
        system(fixed(i),fixed(i)) = 1.0_rk
    end do

    allocate(rhs(ncontrol,1))
    temperature = coolant_temperature
    write(*,"(a)") "Transient cooling-fin solution"
    write(*,"(a)") "  step       time [s]       minimum [K]       maximum [K]"
    do step = 1, nstep
        do a = 1, ncontrol
            rhs(a,1) = source(a) - boundary_correction(a)
        end do
        do b = 1, ncontrol
            do a = 1, ncontrol
                rhs(a,1) = rhs(a,1) + mass(a,b)*temperature(b)/time_step
            end do
        end do
        do i = 1, nfixed
            rhs(fixed(i),1) = coolant_temperature
        end do
        solution = solve(system, rhs)
        temperature = solution(:,1)
        elapsed_time = real(step,rk)*time_step
        write(*,"(2x,i4,3(3x,es14.6))") step, elapsed_time, minval(temperature), maxval(temperature)
    end do

    stored_energy = 0.0_rk
    do i = 1, ncontrol
        stored_energy = stored_energy + capacity_weight(i)*(temperature(i)-coolant_temperature)
    end do
    generated_energy = elapsed_time*sum(source)
    removed_energy = generated_energy-stored_energy

    call fin%create(&
        res1 = 13,&
        res2 = 13,&
        res3 = 25)
    call fin%basis(&
        Tgc = basis_plot)
    Xg = fin%get_Xg()
    allocate(plot_field(size(Xg,1),1))
    plot_field(:,1) = matmul(basis_plot,temperature)

    write(*,"(a,l1)") "Valid fin geometry     : ", fin%err%ok
    write(*,"(a,3(i0,1x))") "Control lattice        : ", nc
    write(*,"(a,i0)") "Elements               : ", size(element,1)
    write(*,"(a,i0)") "Thermal dofs           : ", ncontrol
    write(*,"(a,es12.4)") "Stored energy [J]      : ", stored_energy
    write(*,"(a,es12.4)") "Generated heat [J]     : ", generated_energy
    write(*,"(a,es12.4)") "Heat removed [J]       : ", removed_energy

    call fin%export_Xc(&
        filename    = "vtk/volume_iga_heat_transient_Xc.vtk",&
        point_data  = solution,&
        field_names = [character(len=16) :: "temperature"])
    call fin%export_Xg(&
        filename    = "vtk/volume_iga_heat_transient_Xg.vtk",&
        point_data  = plot_field,&
        field_names = [character(len=16) :: "temperature"])
    call fin%export_Xth_in_Xg(&
        filename = "vtk/volume_iga_heat_transient_Xth.vtk",&
        res      = 11)
    call fin%show(&
        vtkfile_Xc        = "vtk/volume_iga_heat_transient_Xc.vtk",&
        vtkfile_Xg        = "vtk/volume_iga_heat_transient_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_iga_heat_transient_Xth.vtk")
    call fin%finalize()

end program volume_iga_heat_transient
