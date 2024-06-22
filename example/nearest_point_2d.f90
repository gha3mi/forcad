program nearest_point_2d

    use forcad, only: rk, nurbs_surface

    implicit none

    type(nurbs_surface) :: shape           !! Declare a NURBS surface object
    real(rk), allocatable :: nearest_Xg(:) !! Coordinates of the nearest point on the surface
    real(rk), allocatable :: nearest_Xt(:) !! Corresponding parametric coordinates of the nearest point
    integer :: id                          !! id of the nearest point
    real(rk) :: Xc(4,3)                    !! Control points
    real(rk) :: Wc(4)                      !! Weights of the control points

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS tetrangon
    !-----------------------------------------------------------------------------

    !> Set a surface with 4 control points
    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [2.0_rk, 0.0_rk, 0.0_rk]
    Xc(3,:) = [0.0_rk, 2.0_rk, 0.0_rk]
    Xc(4,:) = [2.0_rk, 2.0_rk, 0.0_rk]

    !> The weights of the control points (Wc) are optional.
    Wc = [1.0_rk, 1.1_rk, 0.7_rk, 1.0_rk]

    call shape%set(knot1=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], knot2=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xc=Xc, Wc=Wc)

    !-----------------------------------------------------------------------------
    ! Creating the NURBS surface
    !-----------------------------------------------------------------------------

    !> Generate the NURBS surface with resolutions of 30 in both dimensions
    call shape%create(30, 30)

    !-----------------------------------------------------------------------------
    ! Nearest point on the surface (Approximation)
    !-----------------------------------------------------------------------------

    !> Find the nearest point on the surface to a given point
    ! nearest_Xg: Coordinates of the nearest point on the surface (optional)
    ! nearest_Xt: Corresponding parametric coordinates of the nearest point (optional)
    ! id: id of the nearest point (optional)
    call shape%nearest_point([1.3_rk, 1.0_rk, 1.999999999_rk], nearest_Xg, nearest_Xt, id)
    print '(a,1x,g0,2x,g0,2x,g0,a,2x,g0,2x,g0,2x,a,1x,g0)',&
        'Nearest point on the surface:', nearest_Xg, ' with parametric coordinates:', nearest_Xt, ' and id:', id

    !-----------------------------------------------------------------------------
    ! Nearest point on the surface (Optimization)
    !-----------------------------------------------------------------------------

    !> Find the nearest point on the surface to a given point
    !> The optimization method is used to find the nearest point
    !> The optimization method is based on the Newton-Raphson method
    ! nearest_Xt: Corresponding parametric coordinates of the nearest point
    ! nearest_Xg: Coordinates of the nearest point on the surface (optional)
    call shape%nearest_point2([1.3_rk, 1.0_rk, 1.999999999_rk], 1.0e-11_rk, 30, nearest_Xt, nearest_Xg)
    print '(a,1x,g0,2x,g0,2x,g0,a,2x,g0,2x,g0)',&
        'Nearest point on the surface:', nearest_Xg, ' with parametric coordinates:', nearest_Xt

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS surface object
    call shape%finalize()
    ! deallocate(nearest_Xg, nearest_Xt)

end program
