!> This program demonstrates the usage of a NURBS (Non-Uniform Rational B-Spline) curve object to create  and finalize a NURBS curve.
!> It sets up control points, weights, and knot vectors for all three dimensions, generates the curve, and exports the control points and the curve to VTK files.

program example1_curve

    use forcad, only: rk, nurbs_curve

    implicit none
    type(nurbs_curve) :: nurbs              !! Declare a NURBS curve object
    real(rk), allocatable :: Xc(:,:), Wc(:) !! Arrays for control points and weights
    real(rk) :: knot(6)                     !! Array for knot vector

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS curve
    !-----------------------------------------------------------------------------

    !> Define control points for the NURBS curve
    allocate(Xc(3, 3))
    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [0.0_rk, 5.0_rk, 0.0_rk]
    Xc(3,:) = [5.0_rk, 5.0_rk, 0.0_rk]

    !> Define weights for the control points (optional)
    allocate(Wc(3))
    Wc = [1.0_rk, 2.0_rk, 0.3_rk]

    !> Define knot vector
    knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    !> Set knot vector, control points, and weights for the NURBS curve object.
    !> Wc is optional
    call nurbs%set(knot, Xc, Wc)

    !> Deallocate local arrays
    deallocate(Xc, Wc)

    !> Export control points to a VTK file
    call nurbs%export_Xc('vtk/nurbs_curve_Xc.vtk')

    !-----------------------------------------------------------------------------
    ! Creating the NURBS curve
    !-----------------------------------------------------------------------------

    !> Generate the NURBS curve with a resolution of 20
    call nurbs%create(res = 20)

    !> Export the generated curve to a VTK file
    call nurbs%export_Xg('vtk/nurbs_curve_Xg.vtk')


    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call nurbs%show('vtk/nurbs_curve_Xc.vtk','vtk/nurbs_curve_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Refinements
    !-----------------------------------------------------------------------------

    !> Print size of the knot vector
    print*, size(nurbs%get_knot())

    !> Insert knots 0.25, twice and 0.75, once
    call nurbs%insert_knots([0.25_rk, 0.75_rk], [2,1])

    !> Print size of the updated knot vector
    print*, size(nurbs%get_knot())

    !> Print the degree of the curve
    print*, nurbs%get_degree()

    !> Elevate the degree of the curve (2 times)
    call nurbs%elevate_degree(2)

    !> Print the updated degree of the curve
    print*, nurbs%get_degree()

    !> Print size of the knot vector
    print*, size(nurbs%get_knot())

    !> Remove knots 0.25, twice and 0.75, once
    call nurbs%remove_knots([0.25_rk, 0.75_rk], [2,1])

    !> Print size of the updated knot vector
    print*, size(nurbs%get_knot())

    !> Generate the refined curve with a resolution of 20
    call nurbs%create()

    !> Export updated control points to a VTK file
    call nurbs%export_Xc('vtk/nurbs_curve_Xc2.vtk')

    !> Export the refined generated curve to a VTK file
    call nurbs%export_Xg('vtk/nurbs_curve_Xg2.vtk')

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call nurbs%show('vtk/nurbs_curve_Xc2.vtk','vtk/nurbs_curve_Xg2.vtk')

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
    call nurbs%export_Xc('vtk/nurbs_curve_Xc3.vtk')

    !> Export the transformed generated volume to a VTK file
    call nurbs%export_Xg('vtk/nurbs_curve_Xg3.vtk')

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call nurbs%show('vtk/nurbs_curve_Xc3.vtk','vtk/nurbs_curve_Xg3.vtk')

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS curve object
    call nurbs%finalize()

end program example1_curve
