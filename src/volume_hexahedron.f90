!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program volume_hexahedron

    use forcad, only: rk, nurbs_volume

    implicit none
    type(nurbs_volume) :: shape


    !> Set up a hexahedron shape with dimensions L = [2.0, 4.0, 8.0] and a specified number of control points nc = [4, 6, 8].
    !> The weights of the control points (Wc) are optional.
    call shape%set_hexahedron(L=[2.0_rk, 4.0_rk, 8.0_rk], nc=[4,6,8])

    ! Additional modifications can be made to control points and weights, or the NURBS can be refined using knot insertion or degree elevation.
    ! call shape%insert_knots(...)
    ! call shape%elevate_degree(...)
    ! ...

    !> Export the control points to a VTK file for visualization.
    call shape%export_Xc("vtk/volume_hexahedron_Xc.vtk")

    !> Create the shape using the specified number of elements in each direction.
    call shape%create(8, 16, 32)

    !> Export the geometry to a VTK file for visualization.
    call shape%export_Xg("vtk/volume_hexahedron_Xg.vtk")

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call shape%export_Xth_in_Xg("vtk/volume_hexahedron_Xth.vtk")
    call shape%show("vtk/volume_hexahedron_Xc.vtk", "vtk/volume_hexahedron_Xg.vtk", "vtk/volume_hexahedron_Xth.vtk")

    !> Finalize and clean up the shape object.
    call shape%finalize()

end program volume_hexahedron
