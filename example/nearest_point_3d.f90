program nearest_point_3d

    use forcad, only: rk, nurbs_volume

    implicit none

    type(nurbs_volume) :: shape            !! Declare a NURBS volume object
    real(rk), allocatable :: nearest_Xg(:) !! Coordinates of the nearest point on the volume
    real(rk), allocatable :: nearest_Xt(:) !! Corresponding parametric coordinates of the nearest point
    integer :: id                          !! id of the nearest point
    real(rk) :: Xc(8,3)                    !! Control points
    real(rk) :: Wc(8)                      !! Weights of the control points

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS hexahedron
    !-----------------------------------------------------------------------------

    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [2.0_rk, 0.0_rk, 0.0_rk]
    Xc(3,:) = [0.0_rk, 4.0_rk, 0.0_rk]
    Xc(4,:) = [2.0_rk, 4.0_rk, 0.0_rk]
    Xc(5,:) = [0.0_rk, 0.0_rk, 2.0_rk]
    Xc(6,:) = [2.0_rk, 0.0_rk, 2.0_rk]
    Xc(7,:) = [0.0_rk, 4.0_rk, 2.0_rk]
    Xc(8,:) = [2.0_rk, 4.0_rk, 2.0_rk]

    !> The weights of the control points (Wc) are optional.
    Wc = [1.0_rk, 1.1_rk, 1.11_rk, 1.0_rk, 0.5_rk, 0.5_rk, 1.2_rk, 1.0_rk]

    call shape%set(&
        knot1=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
        knot2=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
        knot3=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
        Xc=Xc, Wc=Wc)

    !-----------------------------------------------------------------------------
    ! Creating the NURBS volume
    !-----------------------------------------------------------------------------

    !> Generate the NURBS volume with resolutions of 20, 20, 20
    call shape%create(30, 30, 30)

    !-----------------------------------------------------------------------------
    ! Nearest point on the volume
    !-----------------------------------------------------------------------------

    !> Find the nearest point on the volume to a given point
    ! nearest_Xg: Coordinates of the nearest point on the volume (optional)
    ! nearest_Xt: Corresponding parametric coordinates of the nearest point (optional)
    ! id: id of the nearest point (optional)
    call shape%nearest_point([1.5_rk, 3.5_rk, 1.1_rk], nearest_Xg, nearest_Xt, id)
    print '(a,1x,g0,2x,g0,2x,g0,a,2x,g0,2x,g0,2x,g0,2x,a,1x,g0)',&
        'Nearest point on the volume:', nearest_Xg, ' with parametric coordinates:', nearest_Xt, ' and id:', id

    !-----------------------------------------------------------------------------
    ! Nearest point on the volume (Optimization)
    !-----------------------------------------------------------------------------

    !> Find the nearest point on the volume to a given point
    !> The optimization method is used to find the nearest point
    !> The optimization method is based on the Newton-Raphson method
    ! nearest_Xt: Corresponding parametric coordinates of the nearest point
    ! nearest_Xg: Coordinates of the nearest point on the volume (optional)
    call shape%nearest_point2([1.5_rk, 3.5_rk, 1.1_rk], 1.0e-11_rk, 500, nearest_Xt, nearest_Xg)
    print '(a,1x,g0,2x,g0,2x,g0,a,2x,g0,2x,g0,2x,g0)',&
        'Nearest point on the volume:', nearest_Xg, ' with parametric coordinates:', nearest_Xt

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS volume object
    call shape%finalize()
    ! deallocate(nearest_Xg, nearest_Xt)

end program
