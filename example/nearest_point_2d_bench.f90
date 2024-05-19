program nearest_point_2d_bench

    use forcad, only: rk, nurbs_surface
    use fortime

    implicit none

    type(nurbs_surface) :: shape           !! Declare a NURBS surface object
    real(rk), allocatable :: nearest_Xg(:) !! Coordinates of the nearest point on the surface
    real(rk), allocatable :: nearest_Xt(:) !! Corresponding parametric coordinates of the nearest point
    integer :: id                          !! id of the nearest point
    real(rk), allocatable :: points(:,:)
    integer :: i, j
    type(timer) :: t

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
    call shape%create(100, 100)

    !-----------------------------------------------------------------------------
    ! Nearest point on the surface
    !-----------------------------------------------------------------------------

    !> Find the nearest point on the surface to a given point
    ! nearest_Xg: Coordinates of the nearest point on the surface (optional)
    ! nearest_Xt: Corresponding parametric coordinates of the nearest point (optional)
    ! id: id of the nearest point (optional)
    do j = 1, 40
        allocate(points(j*1000, 3))
        print*, j*1000
        call random_number(points)
        call t%timer_start()
        do concurrent (i = 1: size(points,1))
            call shape%nearest_point(points(i,:), nearest_Xg, nearest_Xt, id)
        end do
        call t%timer_stop()
        deallocate(points)
    end do
    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS surface object
    call shape%finalize()
    deallocate(nearest_Xg, nearest_Xt)

end program
