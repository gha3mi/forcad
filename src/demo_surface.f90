!> This program demonstrates the usage of a NURBS surface object to create, and finalize a NURBS surface.
!> It sets up control points and weights, generates the surface, and exports the control points
!> and the surface to VTK files at various stages.

program example_nurbs_surface

    use forcad, only: rk, nurbs_surface

    implicit none
    type(nurbs_surface) :: nurbs             !! Declare a NURBS surface object
    real(rk), allocatable :: Xc(:,:), Wc(:)  !! Arrays for control points and weights

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS surface
    !-----------------------------------------------------------------------------

    !> Define control points for the NURBS surface
    Xc = generate_Xc(10, 10, 1.5_rk)

    !> Define weights for the control points
    allocate(Wc(size(Xc,1)), source=1.0_rk)

    !> Set control points and weights for the NURBS surface object
    call nurbs%set([10,10],Xc,Wc)

    !> Deallocate local arrays
    deallocate(Xc, Wc)

    !-----------------------------------------------------------------------------
    ! Refinement
    !-----------------------------------------------------------------------------

    !> Elevate the degree of the NURBS surface in the first and second directions
    call nurbs%elevate_degree(1,3)
    call nurbs%elevate_degree(2,3)

    !> Insert knots into the NURBS surface in the first and second directions
    call nurbs%insert_knots(1,[0.5_rk], [1])
    call nurbs%insert_knots(2,[0.5_rk], [1])


    !> Export control points to a VTK file
    call nurbs%export_Xc('vtk/demo_surface_Xc.vtk')

    !-----------------------------------------------------------------------------
    ! Creating the NURBS surface
    !-----------------------------------------------------------------------------

    !> Generate the NURBS surface with a resolution of 30x30
    call nurbs%create(res1=30, res2=30)


    !> Export the generated surface to a VTK file
    call nurbs%export_Xg('vtk/demo_surface_Xg.vtk')

    !> Export the parameter space to a VTK file
    call nurbs%export_Xth_in_Xg('vtk/demo_surface_Xth_in_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry, geometry and parameters using PyVista
    call nurbs%show('vtk/demo_surface_Xc.vtk', 'vtk/demo_surface_Xg.vtk', 'vtk/demo_surface_Xth_in_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS surface object
    call nurbs%finalize()

contains

    !-----------------------------------------------------------------------------
    pure function generate_Xc(num_rows, num_cols, peak_height) result(control_points)
        integer,  intent(in) :: num_rows, num_cols
        real(rk), intent(in) :: peak_height
        real(rk), allocatable :: control_points(:,:)

        integer  :: i, j, idx
        real(rk) :: x_spacing, y_spacing, x_offset, y_offset
        real(rk) :: x, y, z

        x_spacing = merge(0.0_rk, 1.0_rk/real(num_cols-1, rk), num_cols == 1)
        y_spacing = merge(0.0_rk, 1.0_rk/real(num_rows-1, rk), num_rows == 1)

        x_offset = -0.5_rk
        y_offset = -0.5_rk

        allocate(control_points(num_rows*num_cols, 3))

        do concurrent (i = 1:num_rows, j = 1:num_cols) local(idx, x, y, z)
            idx = (i-1)*num_cols+j

            x = x_offset+real(j-1, rk)*x_spacing
            y = y_offset+real(i-1, rk)*y_spacing
            z = peak_height*exp(-(x**2+y**2))+0.1_rk*peak_height

            control_points(idx, 1) = x
            control_points(idx, 2) = y
            control_points(idx, 3) = z
        end do
    end function
    !-----------------------------------------------------------------------------

end program example_nurbs_surface
