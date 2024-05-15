program nearest_point_1d

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: shape             !! Declare a NURBS curve object
    real(rk), allocatable :: nearest_Xg(:) !! Coordinates of the nearest point on the curve
    real(rk) :: nearest_Xt                 !! Corresponding parametric coordinates of the nearest point
    integer :: id                          !! id of the nearest point

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS circle
    !-----------------------------------------------------------------------------

    !> Set a circle with radius 2.0 and center at [0.0, 0.0, 0.0]
    call shape%set_circle(center = [0.0_rk, 0.0_rk, 0.0_rk], radius = 2.0_rk)

    !-----------------------------------------------------------------------------
    ! Creating circle
    !-----------------------------------------------------------------------------

    !> Generate the NURBS circle with a resolution of 100
    call shape%create(res = 100)

    !-----------------------------------------------------------------------------
    ! Nearest point on the curve
    !-----------------------------------------------------------------------------

    !> Find the nearest point on the curve to a given point
    ! nearest_Xg: Coordinates of the nearest point on the curve (optional)
    ! nearest_Xt: Corresponding parametric coordinates of the nearest point (optional)
    ! id: id of the nearest point (optional)
    call shape%nearest_point([2.0_rk, 3.0_rk, 5.0_rk], nearest_Xg, nearest_Xt, id)
    print *, 'Nearest point on the curve:', nearest_Xg, 'with parametric coordinates:', nearest_Xt, 'and id:', id

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS curve object
    call shape%finalize()

end program
