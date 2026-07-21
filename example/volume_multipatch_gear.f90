!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> C1-connected quadratic multipatch approximation of a complete involute helical gear.
program volume_multipatch_gear

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: gear_teeth = 20, tooth_elements = 10
    integer, parameter :: nr = 3, ns = 2*tooth_elements + 1, nz = 3
    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: pressure_angle = 20.0_rk*pi/180.0_rk
    real(rk), parameter :: helix_angle = 24.0_rk*pi/180.0_rk
    real(rk), parameter :: gear_module = 1.0_rk
    real(rk), parameter :: pitch_radius = 0.5_rk*gear_module*real(gear_teeth,rk)
    real(rk), parameter :: base_radius = pitch_radius*cos(pressure_angle)
    real(rk), parameter :: tip_radius = pitch_radius + gear_module
    real(rk), parameter :: root_radius = pitch_radius - 1.25_rk*gear_module
    real(rk), parameter :: bore_radius = 0.58_rk*pitch_radius
    real(rk), parameter :: face_width = 4.0_rk*gear_module
    real(rk), parameter :: half_pitch = pi/real(gear_teeth,rk)
    real(rk), parameter :: root_handle_angle = half_pitch/real(tooth_elements,rk)
    real(rk), parameter :: helix_twist = face_width*tan(helix_angle)/pitch_radius
    real(rk), parameter :: pitch_tangent = tan(pressure_angle)
    real(rk), parameter :: pitch_involute = pitch_tangent - atan(pitch_tangent)
    real(rk), parameter :: tip_tangent = sqrt((tip_radius/base_radius)**2 - 1.0_rk)
    real(rk), parameter :: tip_involute = tip_tangent - atan(tip_tangent)
    real(rk), parameter :: tip_half = 0.5_rk*half_pitch + pitch_involute - tip_involute
    real(rk), parameter :: base_half = 0.5_rk*half_pitch + pitch_involute
    real(rk), parameter :: root_half = 0.8_rk*half_pitch
    real(rk), parameter :: knot_bezier(6) = [0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk]
    real(rk), parameter :: smooth_knots(5) = [0.1_rk,0.3_rk,0.5_rk,0.7_rk,0.9_rk]
    integer, parameter :: smooth_removal(5) = 1

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: gear
    real(rk) :: knot_tooth(2*tooth_elements + 4)
    integer :: id(gear_teeth), i, tooth
    integer, allocatable :: dof_map(:), elem(:,:), rowptr(:), column(:)
    real(rk), allocatable :: constraint(:)

    knot_tooth(1:3) = 0.0_rk
    do i = 1, tooth_elements - 1
        knot_tooth(2*i+2:2*i+3) = real(i,rk)/real(tooth_elements,rk)
    end do
    knot_tooth(size(knot_tooth)-2:) = 1.0_rk

    call set_gear_patch(patch, knot_tooth)
    do tooth = 1, gear_teeth
        call gear%add_patch(patch, id(tooth))
        call patch%rotate_Xc(&
            alpha = 0.0_rk,&
            beta  = 0.0_rk,&
            theta = 360.0_rk/real(gear_teeth,rk))
    end do

    do tooth = 1, gear_teeth - 1
        call gear%connect(&
            patch_a    = id(tooth),&
            side_a     = "v_max",&
            patch_b    = id(tooth+1),&
            side_b     = "v_min",&
            continuity = 1)
    end do
    call gear%connect(&
        patch_a    = id(gear_teeth),&
        side_a     = "v_max",&
        patch_b    = id(1),&
        side_b     = "v_min",&
        continuity = 1)

    dof_map = gear%cmp_dof_map()
    elem = gear%cmp_elem(shared=.true.)
    call gear%cmp_dof_constraint(rowptr, column, constraint)

    write(*,"(a,i0)") "Gear teeth             : ", gear_teeth
    write(*,"(a,f8.3)") "Angular coverage       : ", 360.0_rk
    write(*,"(a,f8.3)") "Helix angle [degree]   : ", helix_angle*180.0_rk/pi
    write(*,"(a)") "Patch continuity       : C1"
    write(*,"(a,i0)") "Elements per tooth     : ", tooth_elements
    write(*,"(a,i0)") "Patches                : ", gear%get_npatch()
    write(*,"(a,i0)") "Connections            : ", gear%get_nconnection()
    write(*,"(a,l1)") "Valid connectivity     : ", gear%is_valid()
    write(*,"(a,i0)") "Shared global dofs     : ", maxval(dof_map)
    write(*,"(a,i0)") "Constraint equations   : ", size(rowptr)-1
    write(*,"(a,i0)") "Constraint nonzeros    : ", size(constraint)
    write(*,"(a,i0)") "Elements               : ", size(elem,1)
    write(*,"(a,i0)") "Nodes per element      : ", size(elem,2)

    call gear%create(7, 41, 7)
    call gear%export_Xc("vtk/volume_multipatch_gear_c1")
    call gear%export_Xg("vtk/volume_multipatch_gear_c1")
    call gear%export_Xth_in_Xg("vtk/volume_multipatch_gear_c1", 7)
    call gear%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_gear_c1_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_gear_c1_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_gear_c1_Xth_*.vtk")
    call gear%finalize()

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct one tooth patch of the bore-to-tip helical parameterization.
    pure subroutine set_gear_patch(p, knot_tooth)
        type(nurbs_volume), intent(inout) :: p
        real(rk), intent(in), contiguous :: knot_tooth(:)
        real(rk) :: Xc(nr*ns*nz,3)
        real(rk) :: u, v, w, phase, theta, radius, outside, z, delta
        integer :: a, b, c, d, n, left, middle, right

        do b = 1, ns
            v = real(b-1,rk)/real(ns-1,rk)
            phase = (2.0_rk*v - 1.0_rk)*half_pitch
            outside = tooth_radius(phase)
            do concurrent (c = 1:nz, a = 1:nr) local(u, w, theta, radius, z, n)
                u = real(a-1,rk)/real(nr-1,rk)
                w = real(c-1,rk)/real(nz-1,rk)
                radius = bore_radius + u*(outside-bore_radius)
                if (b == 2 .or. b == ns-1) radius = radius/cos(root_handle_angle)
                theta = phase + helix_twist*(w-0.5_rk)
                z = face_width*(w-0.5_rk)
                n = a + (b-1)*nr + (c-1)*nr*ns
                Xc(n,1) = radius*cos(theta)
                Xc(n,2) = radius*sin(theta)
                Xc(n,3) = z
            end do
        end do

        do concurrent (c = 1:nz, a = 1:nr) local(b, d, left, middle, right, delta)
            left = a + (2-1)*nr + (c-1)*nr*ns
            middle = a + (3-1)*nr + (c-1)*nr*ns
            right = a + (4-1)*nr + (c-1)*nr*ns
            do d = 1, 3
                Xc(right,d) = 2.0_rk*Xc(middle,d) - Xc(left,d)
            end do

            do b = 7, 15, 4
                left = a + (b-2)*nr + (c-1)*nr*ns
                middle = a + (b-1)*nr + (c-1)*nr*ns
                right = a + b*nr + (c-1)*nr*ns
                do d = 1, 3
                    delta = 0.5_rk*(Xc(right,d)-Xc(left,d))
                    Xc(left,d) = Xc(middle,d) - delta
                    Xc(right,d) = Xc(middle,d) + delta
                end do
            end do

            left = a + (18-1)*nr + (c-1)*nr*ns
            middle = a + (19-1)*nr + (c-1)*nr*ns
            right = a + (20-1)*nr + (c-1)*nr*ns
            do d = 1, 3
                Xc(left,d) = 2.0_rk*Xc(middle,d) - Xc(right,d)
            end do
        end do

        call p%set(knot_bezier, knot_tooth, knot_bezier, Xc, degree=[2,2,2])
        call p%remove_knots(&
            dir = 2,&
            Xth = smooth_knots,&
            r = smooth_removal)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the external radius from the standard involute polar relation.
    pure elemental real(rk) function tooth_radius(angle) result(radius)
        real(rk), intent(in) :: angle
        real(rk) :: target, t, q
        integer :: iteration

        if (abs(angle) <= tip_half) then
            radius = tip_radius
        else if (abs(angle) < base_half) then
            target = base_half - abs(angle)
            t = (3.0_rk*target)**(1.0_rk/3.0_rk)
            do iteration = 1, 6
                t = max(0.0_rk, t - (t-atan(t)-target)*(1.0_rk+t*t)/max(t*t,epsilon(1.0_rk)))
            end do
            radius = base_radius*sqrt(1.0_rk+t*t)
        else if (abs(angle) < root_half) then
            q = (abs(angle)-base_half)/(root_half-base_half)
            radius = base_radius + (root_radius-base_radius)*q*q*(3.0_rk-2.0_rk*q)
        else
            radius = root_radius
        end if
    end function
    !===============================================================================
end program volume_multipatch_gear
