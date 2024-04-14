program shape_tetragon

    use forcad

    implicit none

    type(nurbs_surface) :: shape                !! Declare a NURBS surface object
    real(rk), allocatable :: Xc(:,:), Wc(:)     !! Arrays for control points and weights
    real(rk), allocatable :: knot1(:), knot2(:) !! Arrays for knot vectors in both dimensions

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS tetrangon
    !-----------------------------------------------------------------------------

    !> Set a tetragon with lengths of 2.0 and 3.0 and 3 and 4 control points in each direction
    !> The weights of the control points (Wc) are optional.
    call shape%set_tetragon(L=[2.0_rk, 3.0_rk], nc=[3,4])

    !> Export the control points to a VTK file
    call shape%export_Xc('vtk/shape_tetragon_Xc.vtk')

    !-----------------------------------------------------------------------------
    ! Creating the NURBS surface
    !-----------------------------------------------------------------------------

    !> Generate the NURBS surface with resolutions of 30 in both dimensions
    call shape%create(30, 30)

    !> Export the generated surface to a VTK file
    call shape%export_Xg('vtk/shape_tetragon_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS surface object
    call shape%finalize()

end program
