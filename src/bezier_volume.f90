!> This program demonstrates the usage of a Bezier volume object to create, and finalize a Bezier volume.
!> It sets up control points and weights, generates the volume, and exports the control points
!> and the volume to VTK files at various stages.

program example_bezier_volume

    use forcad, only: rk, bezier_volume

    implicit none
    type(bezier_volume) :: bezier            !! Declare a bezier volume object
    real(rk), allocatable :: Xc(:,:), Wc(:)  !! Arrays for control points and weights

    !-----------------------------------------------------------------------------
    ! Setting up the bezier volume
    !-----------------------------------------------------------------------------

    !> Define control points for the Bezier volume
    Xc = generate_Xc(1.0_rk)

    !> Define weights for the control points
    allocate(Wc(size(Xc, 1)), source=1.0_rk)

    !> Set control points and weights for the bezier volume object
    call bezier%set([2,2,2], Xc, Wc)

    !> Export initial control points to a VTK file
    call bezier%export_Xc('vtk/bezier_volume_Xc.vtk')

    !-----------------------------------------------------------------------------
    ! Creating the bezier volume
    !-----------------------------------------------------------------------------

    !> Generate the Bezier volume with a resolution of 10x10x10
    call bezier%create(15, 15, 15)

    !> Export the generated volume to a VTK file
    call bezier%export_Xg('vtk/bezier_volume_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the Bezier volume object
    call bezier%finalize()

contains

    !-----------------------------------------------------------------------------
    function generate_Xc(L) result(control_points)
        implicit none
        real(rk), intent(in) :: L
        real(rk), dimension(:,:), allocatable :: control_points
        real(rk) :: L2
        L2 = L / 2.0_rk
        allocate(control_points(8, 3))
        control_points(1,:) = [ L2, -L2,  L2]
        control_points(2,:) = [ L2, -L2, -L2]
        control_points(3,:) = [-L2, -L2,  L2]
        control_points(4,:) = [-L2, -L2, -L2]
        control_points(5,:) = [ L2,  L2,  L2]
        control_points(6,:) = [ L2,  L2, -L2]
        control_points(7,:) = [-L2,  L2,  L2]
        control_points(8,:) = [-L2,  L2, -L2]
    end function
    !-----------------------------------------------------------------------------

end program example_bezier_volume
