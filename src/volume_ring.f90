!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program volume_ring

    use forcad, only: rk, nurbs_volume

    implicit none
    type(nurbs_volume) :: shape


    !> Set up a ring shape centered at 0,0,0 with inner radius 1, outer radius 2, and length 1.
    call shape%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 1.0_rk)

    !> Export the control points to a VTK file for visualization.
    call shape%export_Xc("vtk/volume_ring_Xc.vtk")

    !> Create the shape using the specified number of elements in each direction.
    call shape%create(60, 15, 10)

    !> Export the geometry to a VTK file for visualization.
    call shape%export_Xg("vtk/volume_ring_Xg.vtk")

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call shape%export_Xth_in_Xg("vtk/volume_ring_Xth.vtk")
    call shape%show("vtk/volume_ring_Xc.vtk", "vtk/volume_ring_Xg.vtk", "vtk/volume_ring_Xth.vtk")

    !> Finalize and clean up the shape object.
    call shape%finalize()

end program volume_ring
