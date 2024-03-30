!> This program demonstrates the usage of a Bezier surface object to create, and finalize a Bezier surface.
!> It sets up control points and weights, generates the surface, and exports the control points
!> and the surface to VTK files at various stages.

program example_bezier_surface

    use forcad, only: rk, bezier_surface

    implicit none
    type(bezier_surface) :: bezier           !! Declare a bezier surface object
    real(rk), allocatable :: Xc(:,:), Wc(:)  !! Arrays for control points and weights

    !-----------------------------------------------------------------------------
    ! Setting up the bezier surface
    !-----------------------------------------------------------------------------

    !> Define control points for the Bezier surface
    Xc = generate_Xc(10, 10, 1.5_rk)

    !> Define weights for the control points
    allocate(Wc(size(Xc,1)), source=1.0_rk)

    !> Set control points and weights for the bezier surface object
    call bezier%set([10,10],Xc,Wc)

    !> Export initial control points to a VTK file
    call bezier%export_Xc('vtk/bezier_surface_Xc.vtk')

    !-----------------------------------------------------------------------------
    ! Creating the bezier surface
    !-----------------------------------------------------------------------------

    !> Generate the Bezier surface with a resolution of 100x100
    call bezier%create(res1=30, res2=30)

    !> Export the generated surface to a VTK file
    call bezier%export_Xg('vtk/bezier_surface_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the Bezier surface object
    call bezier%finalize()

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

end program example_bezier_surface
