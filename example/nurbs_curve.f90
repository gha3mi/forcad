!> This program demonstrates the usage of a NURBS (Non-Uniform Rational B-Spline) curve object to create  and finalize a NURBS curve.
!> It sets up control points, weights, and knot vectors for all three dimensions, generates the curve, and exports the control points and the curve to VTK files.

program example_nurbs_curve

    use forcad, only: rk, nurbs_curve

    implicit none
    type(nurbs_curve) :: nurbs              !! Declare a NURBS curve object
    real(rk), allocatable :: Xc(:,:), Wc(:) !! Arrays for control points and weights
    real(rk) :: knot(6)                     !! Array for knot vector

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS curve
    !-----------------------------------------------------------------------------

    ! Define control points for the NURBS curve
    allocate(Xc(3, 3))
    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [0.0_rk, 5.0_rk, 0.0_rk]
    Xc(3,:) = [5.0_rk, 5.0_rk, 0.0_rk]

    ! Define weights for the control points
    allocate(Wc(3))
    Wc = [1.0_rk, 2.0_rk, 0.3_rk]

    ! Define knot vector
    knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    ! Set knot vector, control points, and weights for the NURBS curve object
    call nurbs%set(knot, Xc, Wc)

    ! Export control points to a VTK file
    call nurbs%export_Xc('vtk/nurbs_curve_Xc.vtk')

    !-----------------------------------------------------------------------------
    ! Creating the NURBS curve
    !-----------------------------------------------------------------------------

    ! Generate the NURBS curve with a resolution of 20
    call nurbs%create(res = 20)

    ! Export the generated curve to a VTK file
    call nurbs%export_Xg('vtk/nurbs_curve_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Refinements
    !-----------------------------------------------------------------------------

    ! Insert knots
    call nurbs%insert_knot([0.25_rk, 0.75_rk])

    ! Elevate the degree of the curve (2 times)
    call nurbs%elevate_degree(2)

    ! Export updated control points to a VTK file
    call nurbs%export_Xc('vtk/nurbs_curve_Xc2.vtk')

    ! Export the refined generated curve to a VTK file
    call nurbs%export_Xg('vtk/nurbs_curve_Xg2.vtk')

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    ! Finalize the NURBS curve object
    call nurbs%finalize()

end program example_nurbs_curve
