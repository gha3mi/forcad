program shape_ring_2d

    use forcad, only: rk, nurbs_surface

    implicit none
    type(nurbs_surface) :: shape


    !> Set up a ring shape with inner radius 1.0 and outer radius 2.0.
    call shape%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)

    !> Export the control points to a VTK file for visualization.
    call shape%export_Xc('vtk/shape_ring_2d_Xc.vtk')

    !> Create the shape using the specified number of elements in each direction.
    call shape%create(60, 15)

    !> Export the geometry to a VTK file for visualization.
    call shape%export_Xg('vtk/shape_ring_2d_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call shape%show('vtk/shape_ring_2d_Xc.vtk','vtk/shape_ring_2d_Xg.vtk')

    !> Finalize and clean up the shape object.
    call shape%finalize()

end program
