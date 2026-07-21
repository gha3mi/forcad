!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Solve plane-stress elasticity for a curved tapered beam under self-weight.
!! The example assembles interleaved vector degrees of freedom directly from
!! physical NURBS gradients and exports an exaggerated deformed configuration.
program surface_iga_elasticity_beam

    use forcad, only: rk, nurbs_surface, solve

    implicit none

    integer, parameter :: ncu = 8, ncv = 4
    integer, parameter :: degree(2) = [3, 2]
    integer, parameter :: quadrature(2) = [4, 3]
    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: length = 2.4_rk
    real(rk), parameter :: height_root = 0.34_rk
    real(rk), parameter :: height_tip = 0.19_rk
    real(rk), parameter :: camber = 0.24_rk
    real(rk), parameter :: thickness = 0.025_rk
    real(rk), parameter :: young_modulus = 70.0e9_rk
    real(rk), parameter :: poisson_ratio = 0.33_rk
    real(rk), parameter :: density = 2700.0_rk
    real(rk), parameter :: gravity = 9.81_rk

    type(nurbs_surface) :: beam, deformed
    integer :: ie, ig, i, j, a, b, id, jd, n, ncontrol, ndof
    integer :: ua, va, ub, vb, nfixed
    integer, allocatable :: element(:,:), fixed(:)
    real(rk), allocatable :: Xc(:,:), Xdef(:,:), basis(:), gradient(:,:)
    real(rk), allocatable :: stiffness(:,:), rhs(:,:), solution(:,:), displacement(:,:)
    real(rk), allocatable :: basis_plot(:,:), plot_field(:,:), control_field(:,:)
    real(rk) :: s, eta, center_y, slope, normal_scale, section_height
    real(rk) :: dA, integration_scale, c11, c12, c33
    real(rk) :: strain(3), stress(3), von_mises, max_von_mises
    real(rk) :: displacement_norm, max_displacement, tip_displacement
    real(rk) :: strain_energy, visual_scale

    ncontrol = ncu*ncv
    ndof = 2*ncontrol
    allocate(Xc(ncontrol,2))
    do j = 1, ncv
        eta = real(j-1,rk)/real(ncv-1,rk) - 0.5_rk
        do i = 1, ncu
            s = real(i-1,rk)/real(ncu-1,rk)
            center_y = camber*sin(pi*s)
            slope = camber*pi*cos(pi*s)/length
            normal_scale = sqrt(1.0_rk+slope*slope)
            section_height = height_root + (height_tip-height_root)*s
            n = i + (j-1)*ncu
            Xc(n,1) = length*s - eta*section_height*slope/normal_scale
            Xc(n,2) = center_y + eta*section_height/normal_scale
        end do
    end do

    call beam%set(&
        degree = degree,&
        nc     = [ncu, ncv],&
        Xc     = Xc)
    element = beam%cmp_elem()
    allocate(stiffness(ndof,ndof), rhs(ndof,1), source=0.0_rk)

    c11 = young_modulus/(1.0_rk-poisson_ratio**2)
    c12 = poisson_ratio*c11
    c33 = 0.5_rk*young_modulus/(1.0_rk+poisson_ratio)
    do ie = 1, size(element,1)
        do ig = 1, product(quadrature)
            call beam%ansatz(&
                ie       = ie,&
                ig       = ig,&
                Tgc      = basis,&
                dTgc_dXg = gradient,&
                dA       = dA,&
                ngauss   = quadrature)
            integration_scale = thickness*dA
            do a = 1, size(element,2)
                id = element(ie,a)
                va = 2*id
                rhs(va,1) = rhs(va,1) - density*gravity*basis(a)*integration_scale
            end do
            do b = 1, size(element,2)
                jd = element(ie,b)
                ub = 2*jd - 1
                vb = ub + 1
                do a = 1, size(element,2)
                    id = element(ie,a)
                    ua = 2*id - 1
                    va = ua + 1
                    stiffness(ua,ub) = stiffness(ua,ub) + integration_scale*(&
                        c11*gradient(a,1)*gradient(b,1) + c33*gradient(a,2)*gradient(b,2))
                    stiffness(ua,vb) = stiffness(ua,vb) + integration_scale*(&
                        c12*gradient(a,1)*gradient(b,2) + c33*gradient(a,2)*gradient(b,1))
                    stiffness(va,ub) = stiffness(va,ub) + integration_scale*(&
                        c12*gradient(a,2)*gradient(b,1) + c33*gradient(a,1)*gradient(b,2))
                    stiffness(va,vb) = stiffness(va,vb) + integration_scale*(&
                        c11*gradient(a,2)*gradient(b,2) + c33*gradient(a,1)*gradient(b,1))
                end do
            end do
        end do
    end do

    nfixed = 2*ncv
    allocate(fixed(nfixed))
    do j = 1, ncv
        id = 1 + (j-1)*ncu
        fixed(2*j-1) = 2*id - 1
        fixed(2*j) = 2*id
    end do
    do i = 1, nfixed
        stiffness(fixed(i),:) = 0.0_rk
        stiffness(:,fixed(i)) = 0.0_rk
        stiffness(fixed(i),fixed(i)) = 1.0_rk
        rhs(fixed(i),1) = 0.0_rk
    end do
    solution = solve(stiffness, rhs)

    allocate(displacement(ncontrol,2))
    do i = 1, ncontrol
        displacement(i,1) = solution(2*i-1,1)
        displacement(i,2) = solution(2*i,1)
    end do
    max_displacement = 0.0_rk
    do i = 1, ncontrol
        displacement_norm = sqrt(displacement(i,1)**2+displacement(i,2)**2)
        max_displacement = max(max_displacement,displacement_norm)
    end do
    tip_displacement = 0.0_rk
    do j = 1, ncv
        tip_displacement = tip_displacement + displacement(ncu+(j-1)*ncu,2)
    end do
    tip_displacement = tip_displacement/real(ncv,rk)
    strain_energy = 0.5_rk*dot_product(solution(:,1),rhs(:,1))

    max_von_mises = 0.0_rk
    do ie = 1, size(element,1)
        do ig = 1, product(quadrature)
            call beam%ansatz(&
                ie       = ie,&
                ig       = ig,&
                Tgc      = basis,&
                dTgc_dXg = gradient,&
                dA       = dA,&
                ngauss   = quadrature)
            strain = 0.0_rk
            do a = 1, size(element,2)
                id = element(ie,a)
                strain(1) = strain(1) + gradient(a,1)*displacement(id,1)
                strain(2) = strain(2) + gradient(a,2)*displacement(id,2)
                strain(3) = strain(3) + gradient(a,2)*displacement(id,1) + gradient(a,1)*displacement(id,2)
            end do
            stress(1) = c11*strain(1) + c12*strain(2)
            stress(2) = c12*strain(1) + c11*strain(2)
            stress(3) = c33*strain(3)
            von_mises = sqrt(stress(1)**2-stress(1)*stress(2)+stress(2)**2+3.0_rk*stress(3)**2)
            max_von_mises = max(max_von_mises,von_mises)
        end do
    end do

    visual_scale = 0.18_rk*length/max(max_displacement,tiny(1.0_rk))
    allocate(Xdef(ncontrol,2))
    Xdef = Xc + visual_scale*displacement
    call deformed%set(&
        degree = degree,&
        nc     = [ncu, ncv],&
        Xc     = Xdef)
    call deformed%create(&
        res1 = 81,&
        res2 = 25)
    call deformed%basis(&
        Tgc = basis_plot)
    allocate(plot_field(size(basis_plot,1),3), control_field(ncontrol,3))
    plot_field(:,1:2) = matmul(basis_plot,displacement)
    control_field(:,1:2) = displacement
    do i = 1, size(plot_field,1)
        plot_field(i,3) = sqrt(plot_field(i,1)**2+plot_field(i,2)**2)
    end do
    do i = 1, ncontrol
        control_field(i,3) = sqrt(displacement(i,1)**2+displacement(i,2)**2)
    end do

    write(*,"(a,l1)") "Valid beam geometry       : ", beam%err%ok
    write(*,"(a,2(i0,1x))") "Control net               : ", ncu, ncv
    write(*,"(a,i0)") "Elements                  : ", size(element,1)
    write(*,"(a,i0)") "Displacement dofs         : ", ndof
    write(*,"(a,es12.4)") "Tip displacement [m]      : ", tip_displacement
    write(*,"(a,es12.4)") "Maximum displacement [m]  : ", max_displacement
    write(*,"(a,es12.4)") "Maximum von Mises [Pa]    : ", max_von_mises
    write(*,"(a,es12.4)") "Strain energy [J]         : ", strain_energy
    write(*,"(a,es12.4)") "Visualization scale       : ", visual_scale

    call deformed%export_Xc(&
        filename    = "vtk/surface_iga_elasticity_beam_Xc.vtk",&
        point_data  = control_field,&
        field_names = [character(len=16) :: "ux", "uy", "u_magnitude"])
    call deformed%export_Xg(&
        filename    = "vtk/surface_iga_elasticity_beam_Xg.vtk",&
        point_data  = plot_field,&
        field_names = [character(len=16) :: "ux", "uy", "u_magnitude"])
    call deformed%export_Xth_in_Xg(&
        filename = "vtk/surface_iga_elasticity_beam_Xth.vtk",&
        res      = 17)
    call deformed%show(&
        vtkfile_Xc        = "vtk/surface_iga_elasticity_beam_Xc.vtk",&
        vtkfile_Xg        = "vtk/surface_iga_elasticity_beam_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_iga_elasticity_beam_Xth.vtk")
    call beam%finalize()
    call deformed%finalize()

end program surface_iga_elasticity_beam
