!> This program demonstrates the usage of a NURBS curve object to create, and finalize a NURBS curve.
!> It sets up control points and weights, generates the curve, and exports the control points
!> and the curve to VTK files at various stages.

program example_nurbs_curve

    use forcad, only: rk, nurbs_curve

    implicit none
    type(nurbs_curve) :: nurbs                !! Declare a NURBS curve object
    real(rk), allocatable :: Xc(:,:), Wc(:)   !! Arrays for control points and weights

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS curve
    !-----------------------------------------------------------------------------

    !> Define control points for the NURBS curve
    Xc = generate_Xc(5, 1.0_rk, 2.0_rk, 20)

    !> Define weights for the control points
    allocate(Wc(size(Xc, 1)), source = 1.0_rk)

    !> Set control points and weights for the NURBS curve object
    call nurbs%set(Xc, Wc)

    !> Deallocate local arrays
    deallocate(Xc, Wc)

    !> Export initial control points to a VTK file
    call nurbs%export_Xc('vtk/demo_curve_Xc.vtk')

    !-----------------------------------------------------------------------------
    ! Creating the NURBS curve
    !-----------------------------------------------------------------------------

    !> Generate the NURBS curve with a resolution of 500
    call nurbs%create(res=500)

    !> Export the generated curve to a VTK file
    call nurbs%export_Xg('vtk/demo_curve_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call nurbs%show('vtk/demo_curve_Xc.vtk','vtk/demo_curve_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS curve object
    call nurbs%finalize()

contains

    !-----------------------------------------------------------------------------
    function generate_Xc(num_coils, radius, height, num_points_per_coil) result(control_points)
        integer, intent(in) :: num_coils, num_points_per_coil
        real(rk), intent(in) :: radius, height
        real(rk), allocatable :: control_points(:,:)
        integer :: coil, i
        real(rk) :: theta, coil_height
        allocate(control_points(num_coils * num_points_per_coil, 3))
        do coil = 1, num_coils
            coil_height = height * real(coil-1, rk) / real(num_coils-1, rk)
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

end program example_nurbs_curve
