!> This program demonstrates the usage of a NURBS (Non-Uniform Rational B-Spline) volume object to create  and finalize a NURBS volume.
!> It sets up control points, weights, and knot vectors for all three dimensions, generates the volume, and exports the control points and the volume to VTK files.

program example_nurbs_volume

    use forcad, only: rk, nurbs_volume

    implicit none
    type(nurbs_volume) :: nurbs              !! Declare a NURBS volume object
    real(rk), allocatable :: Xc(:,:), Wc(:)  !! Arrays for control points and weights
    real(rk) :: knot1(4), knot2(4), knot3(4) !! Arrays for knot vectors in all three dimensions

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS volume
    !-----------------------------------------------------------------------------

    !> Define the control points for the NURBS volume
    Xc = generate_Xc(5.0_rk)

    !> Define weights for the control points (optional)
    allocate(Wc(size(Xc,1)), source=1.0_rk)
    Wc(2) = 5.0_rk

    !> Define knot vectors for all three dimensions
    knot1 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
    knot2 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
    knot3 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

    !> Set knot vectors, control points, and weights for the NURBS volume object
    !> Wc is optional.
    call nurbs%set(knot1, knot2, knot3, Xc, Wc)

    !> Export the control points to a VTK file
    call nurbs%export_Xc('vtk/nurbs_volume_Xc.vtk')

    !-----------------------------------------------------------------------------
    ! Creating the NURBS volume
    !-----------------------------------------------------------------------------

    !> Generate the NURBS volume with resolutions of 20, 20, and 20 in the three dimensions
    call nurbs%create(20, 20, 20)

    !> Export the generated volume to a VTK file
    call nurbs%export_Xg('vtk/nurbs_volume_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Refinements
    !-----------------------------------------------------------------------------

    ! Insert knots 0.25 and 0.75 in all three directions
    call nurbs%insert_knots(1, [0.25_rk, 0.75_rk], [1,1]) ! direction 1
    call nurbs%insert_knots(2, [0.25_rk, 0.75_rk], [1,1]) ! direction 2
    call nurbs%insert_knots(3, [0.25_rk, 0.75_rk], [1,1]) ! direction 3

    ! Elevate degree by 2 in all three directions
    call nurbs%elevate_degree(1, 2) ! direction 1
    call nurbs%elevate_degree(2, 2) ! direction 2
    call nurbs%elevate_degree(3, 2) ! direction 3

    ! Export updated control points to a VTK file
    call nurbs%export_Xc('vtk/nurbs_volume_Xc2.vtk')

    ! Export the refined generated volume to a VTK file
    call nurbs%export_Xg('vtk/nurbs_volume_Xg2.vtk')

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS volume object
    call nurbs%finalize()

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

end program example_nurbs_volume
