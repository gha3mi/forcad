!> This program demonstrates the usage of a NURBS (Non-Uniform Rational B-Spline) volume object to create  and finalize a NURBS volume.
!> It sets up control points, weights, and knot vectors for all three dimensions, generates the volume, and exports the control points and the volume to VTK files.

program example3_volume

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

    !> Deallocate local arrays
    deallocate(Xc, Wc)

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
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call nurbs%show('vtk/nurbs_volume_Xc.vtk','vtk/nurbs_volume_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Refinements
    !-----------------------------------------------------------------------------

    !> Print size of knot vectors
    print*, size(nurbs%get_knot(1))
    print*, size(nurbs%get_knot(2))
    print*, size(nurbs%get_knot(3))

    !> Insert knots 0.25 and 0.75 in all three directions
    call nurbs%insert_knots(1, [0.25_rk, 0.75_rk], [1,1]) ! direction 1
    call nurbs%insert_knots(2, [0.25_rk, 0.75_rk], [1,1]) ! direction 2
    call nurbs%insert_knots(3, [0.25_rk, 0.75_rk], [1,1]) ! direction 3

    !> Print size of knot vectors after inserting knots
    print*, size(nurbs%get_knot(1))
    print*, size(nurbs%get_knot(2))
    print*, size(nurbs%get_knot(3))

    !> Print degrees
    print*, nurbs%get_degree()

    !> Elevate degree by 2 in all three directions
    call nurbs%elevate_degree(1, 2) ! direction 1
    call nurbs%elevate_degree(2, 2) ! direction 2
    call nurbs%elevate_degree(3, 2) ! direction 3

    !> Print degrees after elevating
    print*, nurbs%get_degree()

    !> Print size of knot vectors
    print*, size(nurbs%get_knot(1))
    print*, size(nurbs%get_knot(2))
    print*, size(nurbs%get_knot(3))

    call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [1,1]) ! direction 1
    call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [1,1]) ! direction 2
    call nurbs%remove_knots(3, [0.25_rk, 0.75_rk], [1,1]) ! direction 3

    !> Print size of knot vectors after removing knots
    print*, size(nurbs%get_knot(1))
    print*, size(nurbs%get_knot(2))
    print*, size(nurbs%get_knot(3))

    !> Generate the refined NURBS volume with resolutions of 40, 40, and 40 in the three dimensions
    call nurbs%create()

    !> Export updated control points to a VTK file
    call nurbs%export_Xc('vtk/nurbs_volume_Xc2.vtk')

    !> Export the refined generated volume to a VTK file
    call nurbs%export_Xg('vtk/nurbs_volume_Xg2.vtk')

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call nurbs%show('vtk/nurbs_volume_Xc2.vtk','vtk/nurbs_volume_Xg2.vtk')

    !-----------------------------------------------------------------------------
    ! Transformations
    !-----------------------------------------------------------------------------

    !> Rotate the control points
    call nurbs%rotate_Xc(alpha=-45.0_rk, beta=0.0_rk, theta=90.0_rk)

    !> Rotate the generated curve
    call nurbs%rotate_Xg(alpha=-45.0_rk, beta=0.0_rk, theta=90.0_rk)

    !> Translate the control points
    call nurbs%translate_Xc([1.0_rk, 2.0_rk, -3.0_rk])

    !> Translate the generated curve
    call nurbs%translate_Xg([1.0_rk, 2.0_rk, -3.0_rk])

    !> Export the transformed control points to a VTK file
    call nurbs%export_Xc('vtk/nurbs_volume_Xc3.vtk')

    !> Export the transformed generated volume to a VTK file
    call nurbs%export_Xg('vtk/nurbs_volume_Xg3.vtk')

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call nurbs%show('vtk/nurbs_volume_Xc3.vtk','vtk/nurbs_volume_Xg3.vtk')

    !-----------------------------------------------------------------------------
    ! Extract faces
    !-----------------------------------------------------------------------------

    !> first compute and set the connectivities of volume elements
    call nurbs%set_elem(nurbs%cmp_elem())

    !> get the connectivity of the face1 of the first element
    print*, 'Face 1 of element 1:', nurbs%cmp_elemFace(elem=1, face=1)
    print*, 'Face 2 of element 1:', nurbs%cmp_elemFace(elem=1, face=2)
    print*, 'Face 3 of element 1:', nurbs%cmp_elemFace(elem=1, face=3)
    print*, 'Face 4 of element 1:', nurbs%cmp_elemFace(elem=1, face=4)
    print*, 'Face 5 of element 1:', nurbs%cmp_elemFace(elem=1, face=5)
    print*, 'Face 6 of element 1:', nurbs%cmp_elemFace(elem=1, face=6)

    !> get the degree of the faces
    print*, 'Degree of face 1:', nurbs%cmp_degreeFace(face=1)
    print*, 'Degree of face 2:', nurbs%cmp_degreeFace(face=2)
    print*, 'Degree of face 3:', nurbs%cmp_degreeFace(face=3)
    print*, 'Degree of face 4:', nurbs%cmp_degreeFace(face=4)
    print*, 'Degree of face 5:', nurbs%cmp_degreeFace(face=5)
    print*, 'Degree of face 6:', nurbs%cmp_degreeFace(face=6)
    
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

end program example3_volume
