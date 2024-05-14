program nearest_point_3d

    use forcad, only: rk, nurbs_volume

    implicit none

    type(nurbs_volume) :: shape            !! Declare a NURBS volume object
    real(rk), allocatable :: nearest_Xg(:) !! Coordinates of the nearest point on the volume
    real(rk), allocatable :: nearest_Xt(:) !! Corresponding parametric coordinates of the nearest point
    integer :: id                          !! id of the nearest point

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS hexahedron
    !-----------------------------------------------------------------------------

    !> Set up a hexahedron shape with dimensions L = [2.0, 4.0, 8.0] and a specified number of control points nc = [4, 6, 8].
    !> The weights of the control points (Wc) are optional.
    call shape%set_hexahedron(L=[2.0_rk, 4.0_rk, 8.0_rk], nc=[4,6,8])

    !-----------------------------------------------------------------------------
    ! Creating the NURBS volume
    !-----------------------------------------------------------------------------

    !> Generate the NURBS volume with resolutions of 8, 16 and 32
    call shape%create(8, 16, 32)

    !-----------------------------------------------------------------------------
    ! Nearest point on the volume
    !-----------------------------------------------------------------------------

    !> Find the nearest point on the volume to a given point
    ! nearest_Xg: Coordinates of the nearest point on the volume (optional)
    ! nearest_Xt: Corresponding parametric coordinates of the nearest point (optional)
    ! id: id of the nearest point (optional)
    call shape%nearest_point([2.0_rk, 3.0_rk, 5.0_rk], nearest_Xg, nearest_Xt, id)
    print *, 'Nearest point on the volume:', nearest_Xg, 'with parametric coordinates:', nearest_Xt, 'and id:', id

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS volume object
    call shape%finalize()
    deallocate(nearest_Xg, nearest_Xt)

end program
