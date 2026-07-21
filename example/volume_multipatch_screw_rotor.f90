!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> C1-connected variable-pitch three-lobe screw rotor assembled from cubic volume patches.
program volume_multipatch_screw_rotor

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: nsegment = 18, nlobe = 3
    integer, parameter :: nu = 4, nv = 5, nw = 5, ncontrol = nu*nv*nw
    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: rotor_length = 12.0_rk, rotor_radius = 1.55_rk
    real(rk), parameter :: rotor_turns = 2.30_rk, pitch_variation = 0.20_rk
    real(rk), parameter :: lobe_amplitude = 0.20_rk, corner_rounding = 0.16_rk
    real(rk), parameter :: end_scale = 0.86_rk
    real(rk), parameter :: cross_parameter(5) = [-1.0_rk,-2.0_rk/3.0_rk,0.0_rk,2.0_rk/3.0_rk,1.0_rk]
    real(rk), parameter :: knot_axial(8) = [0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk]
    real(rk), parameter :: knot_cross(9) = [0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk]

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: rotor
    integer :: id(nsegment), segment
    integer, allocatable :: dof_map(:), elem(:,:), rowptr(:), column(:)
    real(rk), allocatable :: constraint(:)

    do segment = 1, nsegment
        call set_rotor_patch(patch, segment)
        call rotor%add_patch(&
            patch = patch,&
            id    = id(segment))
    end do

    do segment = 1, nsegment - 1
        call rotor%connect(&
            patch_a    = id(segment),&
            side_a     = "u_max",&
            patch_b    = id(segment+1),&
            side_b     = "u_min",&
            continuity = 1)
    end do

    dof_map = rotor%cmp_dof_map()
    elem = rotor%cmp_elem(shared=.true.)
    call rotor%cmp_dof_constraint(&
        rowptr = rowptr,&
        col    = column,&
        val    = constraint)

    write(*,"(a,i0)") "Rotor lobes            : ", nlobe
    write(*,"(a,f8.3)") "Rotor turns            : ", rotor_turns
    write(*,"(a,f8.3)") "Rotor length           : ", rotor_length
    write(*,"(a)") "Patch continuity       : C1"
    write(*,"(a,i0)") "Patches                : ", rotor%get_npatch()
    write(*,"(a,i0)") "Connections            : ", rotor%get_nconnection()
    write(*,"(a,l1)") "Valid connectivity     : ", rotor%is_valid()
    write(*,"(a,i0)") "Local control dofs     : ", size(dof_map)
    write(*,"(a,i0)") "Shared global dofs     : ", maxval(dof_map)
    write(*,"(a,i0)") "Constraint equations   : ", size(rowptr)-1
    write(*,"(a,i0)") "Constraint nonzeros    : ", size(constraint)
    write(*,"(a,i0)") "Elements               : ", size(elem,1)
    write(*,"(a,i0)") "Nodes per element      : ", size(elem,2)

    call rotor%create(&
        res1 = 9,&
        res2 = 21,&
        res3 = 21)
    call rotor%export_Xc("vtk/volume_multipatch_screw_rotor")
    call rotor%export_Xg("vtk/volume_multipatch_screw_rotor")
    call rotor%export_Xth_in_Xg("vtk/volume_multipatch_screw_rotor", 9)
    call rotor%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_screw_rotor_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_screw_rotor_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_screw_rotor_Xth_*.vtk")
    call rotor%finalize()

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct one cubic rotor segment from endpoint positions and axial tangents.
    pure subroutine set_rotor_patch(p, segment)
        type(nurbs_volume), intent(inout) :: p
        integer, intent(in) :: segment
        real(rk) :: Xc(ncontrol,3)
        real(rk) :: s0, s1, ds, s, xi, eta, xbase, ybase, xlocal, ylocal
        real(rk) :: theta, dtheta, scale, dscale, xrot, yrot, dxds, dyds
        integer :: a, b, c, n

        s0 = real(segment-1,rk)/real(nsegment,rk)
        s1 = real(segment,rk)/real(nsegment,rk)
        ds = s1 - s0

        do concurrent (c = 1:nw, b = 1:nv, a = 1:nu) &
            local(s, xi, eta, xbase, ybase, xlocal, ylocal, theta, dtheta, &
                scale, dscale, xrot, yrot, dxds, dyds, n)
            xi = cross_parameter(b)
            eta = cross_parameter(c)
            xbase = xi*(1.0_rk-corner_rounding*eta*eta)
            ybase = eta*(1.0_rk-corner_rounding*xi*xi)
            xlocal = rotor_radius*(xbase+lobe_amplitude*(xbase*xbase-ybase*ybase))
            ylocal = rotor_radius*(ybase-2.0_rk*lobe_amplitude*xbase*ybase)

            if (a <= 2) then
                s = s0
            else
                s = s1
            end if
            theta = 2.0_rk*pi*rotor_turns*s + pitch_variation*sin(2.0_rk*pi*s)
            dtheta = 2.0_rk*pi*(rotor_turns+pitch_variation*cos(2.0_rk*pi*s))
            scale = end_scale + (1.0_rk-end_scale)*sin(pi*s)**2
            dscale = (1.0_rk-end_scale)*pi*sin(2.0_rk*pi*s)
            xrot = cos(theta)*xlocal - sin(theta)*ylocal
            yrot = sin(theta)*xlocal + cos(theta)*ylocal
            dxds = dscale*xrot - scale*dtheta*yrot
            dyds = dscale*yrot + scale*dtheta*xrot

            n = a + (b-1)*nu + (c-1)*nu*nv
            Xc(n,1) = scale*xrot
            Xc(n,2) = scale*yrot
            Xc(n,3) = rotor_length*(s-0.5_rk)
            if (a == 2) then
                Xc(n,1) = Xc(n,1) + ds*dxds/3.0_rk
                Xc(n,2) = Xc(n,2) + ds*dyds/3.0_rk
                Xc(n,3) = Xc(n,3) + ds*rotor_length/3.0_rk
            else if (a == 3) then
                Xc(n,1) = Xc(n,1) - ds*dxds/3.0_rk
                Xc(n,2) = Xc(n,2) - ds*dyds/3.0_rk
                Xc(n,3) = Xc(n,3) - ds*rotor_length/3.0_rk
            end if
        end do

        call p%set(&
            knot1  = knot_axial,&
            knot2  = knot_cross,&
            knot3  = knot_cross,&
            Xc     = Xc,&
            degree = [3,3,3])
    end subroutine
    !===============================================================================
end program volume_multipatch_screw_rotor
