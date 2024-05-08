!> This program demonstrates the usage of a NURBS volume object to create, and finalize a NURBS volume.
!> It sets up control points and weights, generates the volume, and exports the control points
!> and the volume to VTK files at various stages.

program example_nurbs_volume

    use forcad, only: rk, nurbs_volume

    implicit none
    type(nurbs_volume) :: nurbs              !! Declare a NURBS volume object
    real(rk), allocatable :: Xc(:,:), Wc(:)  !! Arrays for control points and weights

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS volume
    !-----------------------------------------------------------------------------

    !> Define control points for the NURBS volume
    Xc = generate_Xc(1.0_rk)

    !> Define weights for the control points
    allocate(Wc(size(Xc, 1)), source=1.0_rk)

    !> Set control points and weights for the NURBS volume object
    call nurbs%set([2,2,2], Xc, Wc)

    !> Deallocate local arrays
    deallocate(Xc, Wc)

    !> Export initial control points to a VTK file
    call nurbs%export_Xc('vtk/demo_volume_Xc.vtk')

    !-----------------------------------------------------------------------------
    ! Creating the NURBS volume
    !-----------------------------------------------------------------------------

    !> Generate the NURBS volume with a resolution of 15X15X15
    call nurbs%create(15, 15, 15)

    !> Export the generated volume to a VTK file
    call nurbs%export_Xg('vtk/demo_volume_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call nurbs%show('vtk/demo_volume_Xc.vtk','vtk/demo_volume_Xg.vtk')

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
        real(rk), allocatable :: control_points(:,:)
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
