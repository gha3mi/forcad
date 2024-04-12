program example_morph

    use forcad
    use forcad_utils, only: hexahedron_Xc

    implicit none

    type(nurbs_volume) :: shape, control_shape
    real(rk), allocatable :: X(:,:)
    integer, allocatable :: elem(:,:)
    integer:: i

    !> You can create your shape or use a predefined one
    !> Read coordinates from file
    allocate(X(23200,3))
    open(unit=10, file='example/put_to_nurbs_X.txt')
    do i = 1,23200
        read(10,*) X(i,1), X(i,2), X(i,3)
    end do
    close(10)

    !> Read element connectivities from file
    allocate(elem(20577,8))
    open(unit=10, file='example/put_to_nurbs_elem.txt')
    do i = 1,20577
        read(10,*) elem(i,1), elem(i,2), elem(i,4), elem(i,3), elem(i,5), elem(i,6), elem(i,8), elem(i,7)
    end do

    !> Set a control shape that will be used to put the shape into
    !> The contol shape is a hexahedron with 100x40x10 with 10x5x3 number of control points
    !> By modifying the control shape you can modify the shape
    call control_shape%set(nc=[10,5,3], Xc=hexahedron_Xc(L=[100.0_rk, 40.0_rk, 10.0_rk], nc=[10,5,3]))

    !> Map the shape into the shape
    call control_shape%put_to_nurbs(X, elem)

    !> Export the shape and the control shape to vtk files
    call control_shape%export_Xc('vtk/control_shape.vtk')
    call control_shape%export_Xg('vtk/shape.vtk')

end program