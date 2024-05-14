program nearest_point_2d

    use forcad

    implicit none

    type(nurbs_surface) :: shape             !! Declare a NURBS surface object
    real(rk), allocatable :: nearest(:)   !! Array for nearest point on the surface
    integer :: id                         !! id of the nearest point

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS tetrangon
    !-----------------------------------------------------------------------------

    !> Set a tetragon with lengths of 2.0 and 3.0 and 3 and 4 control points in each direction
    !> The weights of the control points (Wc) are optional.
    call shape%set_tetragon(L=[2.0_rk, 3.0_rk], nc=[3,4])

    !-----------------------------------------------------------------------------
    ! Creating the NURBS surface
    !-----------------------------------------------------------------------------

    !> Generate the NURBS surface with resolutions of 30 in both dimensions
    call shape%create(30, 30)

    !-----------------------------------------------------------------------------
    ! Nearest point on the surface
    !-----------------------------------------------------------------------------

    !> Find the nearest point on the surface to a given point
    call shape%nearest_point([2.0_rk, 3.0_rk, 5.0_rk], nearest, id)
    print *, 'Nearest point on the surface:', nearest, 'with id:', id

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS surface object
    call shape%finalize()

end program
