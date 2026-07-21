!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program surface_half_ring

    use forcad, only: rk, nurbs_surface

    implicit none
    type(nurbs_surface) :: shape


    !> Set up a half ring shape centered at 0,0,0 with inner radius 1 and outer radius 2.
    call shape%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)

    !> Export the control points to a VTK file for visualization.
    call shape%export_Xc("vtk/surface_half_ring_Xc.vtk")

    !> Create the shape using the specified number of elements in each direction.
    call shape%create(60, 15)

    !> Export the geometry to a VTK file for visualization.
    call shape%export_Xg("vtk/surface_half_ring_Xg.vtk")

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call shape%export_Xth_in_Xg("vtk/surface_half_ring_Xth.vtk")
    call shape%show("vtk/surface_half_ring_Xc.vtk", "vtk/surface_half_ring_Xg.vtk", "vtk/surface_half_ring_Xth.vtk")

    !> Finalize and clean up the shape object.
    call shape%finalize()

end program surface_half_ring
