!> This program demonstrates the usage of a NURBS (Non-Uniform Rational B-Spline) surface object to create  and finalize a NURBS surface.
!> It sets up control points, weights, and knot vectors for all three dimensions, generates the surface, and exports the control points and the surface to VTK files.

program example3_surface

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
    Wc(2) = 2.0_rk

    !> Define knot vectors for both dimensions
    knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
    knot2 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    !> Set knot vectors, control points, and weights for the NURBS surface object
    call nurbs%set(knot1, knot2, Xc, Wc)

    !> Deallocate local arrays
    deallocate(Xc, Wc)

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
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call nurbs%show('vtk/nurbs_surface_Xc.vtk','vtk/nurbs_surface_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Refinements
    !-----------------------------------------------------------------------------

    !> Print size of the knot vectors
    print*, size(nurbs%get_knot(1))
    print*, size(nurbs%get_knot(2))

    !> Insert knots 0.25, twice and 0.75, once in both directions
    call nurbs%insert_knots(1, [0.25_rk, 0.75_rk], [2,1]) ! direction 1
    call nurbs%insert_knots(2, [0.25_rk, 0.75_rk], [2,1]) ! direction 2

    !> Print size of the knot vectors after inserting knots
    print*, size(nurbs%get_knot(1))
    print*, size(nurbs%get_knot(2))

    !> Print the degrees
    print*, nurbs%get_degree()

    !> Elevate degree by 2 in both directions
    call nurbs%elevate_degree(1, 2) ! direction 1
    call nurbs%elevate_degree(2, 2) ! direction 2

    !> Print the degrees after elevating
    print*, nurbs%get_degree()

    !> Print size of the knot vectors
    print*, size(nurbs%get_knot(1))
    print*, size(nurbs%get_knot(2))

    !> Remove knots 0.25, twice and 0.75, once in both directions
    call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1]) ! direction 1
    call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1]) ! direction 2

    !> Print size of the knot vectors after removing knots
    print*, size(nurbs%get_knot(1))
    print*, size(nurbs%get_knot(2))

    !> Generate the refined NURBS surface with resolutions of 30 in both dimensions
    call nurbs%create()

    !> Export updated control points to a VTK file
    call nurbs%export_Xc('vtk/nurbs_surface_Xc2.vtk')

    !> Export the refined generated surface to a VTK file
    call nurbs%export_Xg('vtk/nurbs_surface_Xg2.vtk')

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call nurbs%show('vtk/nurbs_surface_Xc2.vtk','vtk/nurbs_surface_Xg2.vtk')

    !-----------------------------------------------------------------------------
    ! Transformations
    !-----------------------------------------------------------------------------

    !> Rotate the control points
    call nurbs%rotate_Xc(alpha=45.0_rk, beta=0.0_rk, theta=90.0_rk)

    !> Rotate the generated curve
    call nurbs%rotate_Xg(alpha=-45.0_rk, beta=0.0_rk, theta=-90.0_rk)

    !> Translate the control points
    call nurbs%translate_Xc([1.0_rk, 2.0_rk, -3.0_rk])

    !> Translate the generated curve
    call nurbs%translate_Xg([-1.0_rk, -2.0_rk, 3.0_rk])

    !> Export the transformed control points to a VTK file
    call nurbs%export_Xc('vtk/nurbs_surface_Xc3.vtk')

    !> Export the transformed generated volume to a VTK file
    call nurbs%export_Xg('vtk/nurbs_surface_Xg3.vtk')

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call nurbs%show('vtk/nurbs_surface_Xc3.vtk','vtk/nurbs_surface_Xg3.vtk')

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
        real(rk), allocatable :: control_points(:,:)
        integer :: i, j
        real(rk) :: x_spacing, y_spacing, x_offset, y_offset
        x_spacing = 1.0_rk / real(num_cols - 1, rk)
        y_spacing = 1.0_rk / real(num_rows - 1, rk)
        x_offset = -0.5_rk
        y_offset = -0.5_rk
        allocate(control_points(num_rows * num_cols, 3))
        do i = 1, num_rows
            do j = 1, num_cols
                control_points((i - 1) * num_cols + j, 1) = x_offset + real(j - 1, rk) * x_spacing
                control_points((i - 1) * num_cols + j, 2) = y_offset + real(i - 1, rk) * y_spacing
                control_points((i - 1) * num_cols + j, 3) = &
                    peak_height * exp(-((control_points((i - 1) * num_cols + j, 1) ** 2) &
                    + (control_points((i - 1) * num_cols + j, 2) ** 2))) + 0.5_rk * peak_height * 0.2_rk
            end do
        end do
    end function
    !-----------------------------------------------------------------------------

end program example3_surface
