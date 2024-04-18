program shape_half_circle

    use forcad, only: rk, nurbs_curve

    implicit none
    type(nurbs_curve) :: shape


    !> Set up a half circle shape centered at the 0,0,0 with a radius of 1
    call shape%set_half_circle([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk)

    !> Export the control points to a VTK file for visualization.
    call shape%export_Xc('vtk/shape_half_circle_Xc.vtk')

    !> Create the shape using the specified number of elements in each direction.
    call shape%create(60)

    !> Export the geometry to a VTK file for visualization.
    call shape%export_Xg('vtk/shape_half_circle_Xg.vtk')

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call shape%show('vtk/shape_half_circle_Xc.vtk','vtk/shape_half_circle_Xg.vtk')

    !> Finalize and clean up the shape object.
    call shape%finalize()

end program
