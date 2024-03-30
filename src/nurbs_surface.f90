!> This program demonstrates the usage of a NURBS (Non-Uniform Rational B-Spline) surface object to create  and finalize a NURBS surface.
!> It sets up control points, weights, and knot vectors for all three dimensions, generates the surface, and exports the control points and the surface to VTK files.

program example_nurbs_surface

    use forcad

    implicit none

    type(nurbs_surface) :: nurbs            !! Declare a NURBS surface object
    real(rk), allocatable :: Xc(:,:), Wc(:) !! Arrays for control points and weights
    real(rk) :: knot1(6), knot2(6)          !! Arrays for knot vectors in both dimensions

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS surface
    !-----------------------------------------------------------------------------

    !> Define control points for the NURBS surface
    Xc = generate_Xc(3, 3, 1.0_rk)

    !> Define weights for the control points
    allocate(Wc(size(Xc, 1)))
    Wc = 1.0_rk

    !> Define knot vectors for both dimensions
    knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
    knot2 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    !> Set knot vectors, control points, and weights for the NURBS surface object
    call nurbs%set(knot1, knot2, Xc, Wc)

    !> Export the control points to a VTK file
    call nurbs%export_Xc('vtk/nurbs_surface_Xc.vtk')

    !-----------------------------------------------------------------------------
    ! Creating the NURBS surface
    !-----------------------------------------------------------------------------

    !> Generate the NURBS surface with resolutions of 30 in both dimensions
    call nurbs%create(30, 30)

    !> Export the generated surface to a VTK file
    call nurbs%export_Xg('vtk/nurbs_surface_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS surface object
    call nurbs%finalize()

contains

    !-----------------------------------------------------------------------------
    function generate_Xc(num_rows, num_cols, peak_height) result(control_points)
        integer, intent(in) :: num_rows, num_cols
        real(rk), intent(in) :: peak_height
        real(rk), dimension(:,:), allocatable :: control_points
        integer :: i, j
        real(rk) :: x_spacing, y_spacing, x_offset, y_offset
        x_spacing = 1.0_rk / real(num_cols - 1)
        y_spacing = 1.0_rk / real(num_rows - 1)
        x_offset = -0.5_rk
        y_offset = -0.5_rk
        allocate(control_points(num_rows * num_cols, 3))
        do i = 1, num_rows
            do j = 1, num_cols
                control_points((i - 1) * num_cols + j, 1) = x_offset + real(j - 1) * x_spacing
                control_points((i - 1) * num_cols + j, 2) = y_offset + real(i - 1) * y_spacing
                control_points((i - 1) * num_cols + j, 3) = &
                    peak_height * exp(-((control_points((i - 1) * num_cols + j, 1) ** 2) &
                    + (control_points((i - 1) * num_cols + j, 2) ** 2))) + 0.5_rk * peak_height * 0.2_rk
            end do
        end do
    end function
    !-----------------------------------------------------------------------------

end program example_nurbs_surface
