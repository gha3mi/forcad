program nearest_point_1d

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: shape              !! Declare a NURBS curve object
    real(rk), allocatable :: Xc(:,:), Wc(:) !! Arrays for control points and weights
    real(rk) :: knot(6)                     !! Array for knot vector
    real(rk), allocatable :: nearest_Xg(:)  !! Array for the nearest point on the curve
    real(rk) :: nearest_Xt                  !! Array for the parametric coordinates of the nearest point
    integer :: id                           !! Variable for the id of the nearest point

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
    Wc = [1.0_rk, 1.1_rk, 1.0_rk]

    !> Define knot vector
    knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    !> Set knot vector, control points, and weights for the NURBS curve object.
    !> Wc is optional
    call shape%set(knot, Xc, Wc)

    !-----------------------------------------------------------------------------
    ! Creating the NURBS curve
    !-----------------------------------------------------------------------------

    !> Generate the NURBS curve with a resolution of 20
    call shape%create(20)

    !-----------------------------------------------------------------------------
    ! Nearest point on the curve
    !-----------------------------------------------------------------------------

    !> Find the nearest point on the curve to a given point
    ! nearest_Xg: Coordinates of the nearest point on the curve (optional)
    ! nearest_Xt: Corresponding parametric coordinates of the nearest point (optional)
    ! id: id of the nearest point (optional)
    call shape%nearest_point([4.5_rk, 4.5_rk, 5.0_rk], nearest_Xg, nearest_Xt, id)
    print '(a,1x,g0,2x,g0,2x,g0,a,2x,g0,2x,a,1x,g0)',&
        'Nearest point on the curve:', nearest_Xg, ' with parametric coordinates:', nearest_Xt, ' and id:', id

    !-----------------------------------------------------------------------------
    ! Nearest point on the curve (Optimization)
    !-----------------------------------------------------------------------------

    !> Find the nearest point on the curve to a given point
    !> The optimization method is used to find the nearest point
    !> The optimization method is based on the Newton-Raphson method
    ! nearest_Xt: Corresponding parametric coordinates of the nearest point
    ! nearest_Xg: Coordinates of the nearest point on the curve (optional)
    call shape%nearest_point2([4.5_rk, 4.5_rk, 5.0_rk], 1.0e-11_rk, 30, nearest_Xt, nearest_Xg)
    print '(a,1x,g0,2x,g0,a,2x,g0,2x,g0)',&
        'Nearest point on the curve:', nearest_Xg, ' with parametric coordinates:', nearest_Xt

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS curve object
    call shape%finalize()
    deallocate(nearest_Xg, Xc, Wc)

end program
