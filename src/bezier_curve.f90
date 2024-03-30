!> This program demonstrates the usage of a Bezier curve object to create, and finalize a Bezier curve.
!> It sets up control points and weights, generates the curve, and exports the control points
!> and the curve to VTK files at various stages.

program example_bezier_curve

    use forcad, only: rk, bezier_curve

    implicit none
    type(bezier_curve) :: bezier              !! Declare a bezier curve object
    real(rk), allocatable :: Xc(:,:), Wc(:)   !! Arrays for control points and weights

    !-----------------------------------------------------------------------------
    ! Setting up the bezier curve
    !-----------------------------------------------------------------------------

    !> Define control points for the Bezier curve
    Xc = generate_Xc(5, 1.0_rk, 2.0_rk, 20)

    !> Define weights for the control points
    allocate(Wc(size(Xc, 1)), source = 1.0_rk)

    !> Set control points and weights for the bezier curve object
    call bezier%set(Xc=Xc, Wc=Wc)

    !> Export initial control points to a VTK file
    call bezier%export_Xc('vtk/bezier_curve_Xc.vtk')

    !-----------------------------------------------------------------------------
    ! Creating the bezier curve
    !-----------------------------------------------------------------------------

    !> Generate the Bezier curve with a resolution of 1000
    call bezier%create(res=500)

    !> Export the generated curve to a VTK file
    call bezier%export_Xg('vtk/bezier_curve_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the Bezier curve object
    call bezier%finalize()

contains

    !-----------------------------------------------------------------------------
    function generate_Xc(num_coils, radius, height, num_points_per_coil) result(control_points)
        integer, intent(in) :: num_coils, num_points_per_coil
        real(rk), intent(in) :: radius, height
        real(rk), dimension(:,:), allocatable :: control_points
        integer :: coil, i
        real(rk) :: theta, coil_height
        allocate(control_points(num_coils * num_points_per_coil, 3))
        do coil = 1, num_coils
            coil_height = height * (coil-1) / real(num_coils-1, rk)
            theta = 0.0_rk
            do i = 1, num_points_per_coil
                theta = theta + 2.0_rk * acos(-1.0_rk) / real(num_points_per_coil, rk)
                control_points((coil - 1) * num_points_per_coil + i, 1) = radius * cos(theta)
                control_points((coil - 1) * num_points_per_coil + i, 2) = radius * sin(theta)
                control_points((coil - 1) * num_points_per_coil + i, 3) = coil_height
            end do
        end do
    end function
    !-----------------------------------------------------------------------------

end program example_bezier_curve
