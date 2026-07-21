!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program curve_shape_c

    use forcad, only: rk, nurbs_curve

    implicit none
    type(nurbs_curve) :: shape

    !-----------------------------------------------------------------------------
    ! Setting up NURBS C-shape
    !-----------------------------------------------------------------------------

    !> Set a C-shape with radius 2.0 and center at [0.0, 0.0, 0.0]
    call shape%set_C(center = [0.0_rk, 0.0_rk, 0.0_rk], radius = 2.0_rk)

    !> Export control points to a VTK file
    call shape%export_Xc("vtk/curve_shape_c_Xc.vtk")

    !-----------------------------------------------------------------------------
    ! Creating C-shape
    !-----------------------------------------------------------------------------

    !> Generate the NURBS C-shape with a resolution of 100
    call shape%create(res = 100)

    !> Export the generated cirlce to a VTK file
    call shape%export_Xg("vtk/curve_shape_c_Xg.vtk")

    !-----------------------------------------------------------------------------
    ! Visualization using PyVista
    ! Note: PyVista is required for visualization. Install it using `pip install pyvista`
    !-----------------------------------------------------------------------------

    !> Show the control geometry and geometry using PyVista
    call shape%export_Xth_in_Xg("vtk/curve_shape_c_Xth.vtk")
    call shape%show("vtk/curve_shape_c_Xc.vtk", "vtk/curve_shape_c_Xg.vtk", "vtk/curve_shape_c_Xth.vtk")

    !-----------------------------------------------------------------------------
    ! Finalizing
    !-----------------------------------------------------------------------------

    !> Finalize the NURBS curve object
    call shape%finalize()

end program curve_shape_c
