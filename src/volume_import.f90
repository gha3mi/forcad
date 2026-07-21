!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program volume_import

    use forcad, only: rk, nurbs_volume, hexahedron_Xc

    implicit none

    type(nurbs_volume) :: control_shape
    real(rk), allocatable :: X(:,:)
    integer, allocatable :: elem(:,:)
    integer:: i, nunit

    !> You can create your shape or use a predefined one
    !> Read coordinates from file
    allocate(X(23200,3))
    open(newunit=nunit, file="example/points.txt")
    do i = 1,23200
        read(nunit,*) X(i,1), X(i,2), X(i,3)
    end do
    close(nunit)

    !> Read element connectivities from file
    allocate(elem(20577,8))
    open(newunit=nunit, file="example/elements.txt")
    do i = 1,20577
        read(nunit,*) elem(i,1), elem(i,2), elem(i,4), elem(i,3), elem(i,5), elem(i,6), elem(i,8), elem(i,7)
    end do
    close(nunit)

    !> Set a control shape that will be used to put the shape into
    !> The contol shape is a hexahedron with 100x40x10 with 10x5x3 number of control points
    !> By modifying the control shape you can modify the shape
    call control_shape%set(nc=[10,5,3], Xc=hexahedron_Xc(L=[100.0_rk, 40.0_rk, 10.0_rk], nc=[10,5,3]))

    !> Map the shape into the shape
    call control_shape%put_to_nurbs(X, elem)

    !> Deallocate local variables
    deallocate(X, elem)

    !> Export the shape and the control shape to vtk files
    call control_shape%export_Xc("vtk/volume_import_Xc.vtk")
    call control_shape%export_Xg("vtk/volume_import_Xg.vtk")

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call control_shape%export_Xth_in_Xg("vtk/volume_import_Xth.vtk")
    call control_shape%show("vtk/volume_import_Xc.vtk", "vtk/volume_import_Xg.vtk", "vtk/volume_import_Xth.vtk")

    !> Finalize the control shape
    call control_shape%finalize()

end program volume_import
