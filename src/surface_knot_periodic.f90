!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Construct a rational surface periodic in its first parameter direction.
program surface_knot_periodic

    use forcad, only: rk, nurbs_surface

    implicit none

    type(nurbs_surface) :: surface
    real(rk) :: Xc(12,3)
    real(rk), parameter :: knot1(9) = [&
        0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
    real(rk), parameter :: knot2(4) = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
    real(rk), parameter :: Wc(12) = [&
        1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk,&
        1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]

    Xc(1,:)  = [ 1.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:)  = [ 0.0_rk, 1.0_rk, 0.0_rk]
    Xc(3,:)  = [-1.0_rk, 0.0_rk, 0.0_rk]
    Xc(4,:)  = [ 0.0_rk,-1.0_rk, 0.0_rk]
    Xc(5,:)  = Xc(1,:)
    Xc(6,:)  = Xc(2,:)
    Xc(7,:)  = [ 1.3_rk, 0.0_rk, 1.0_rk]
    Xc(8,:)  = [ 0.0_rk, 1.3_rk, 1.0_rk]
    Xc(9,:)  = [-1.3_rk, 0.0_rk, 1.0_rk]
    Xc(10,:) = [ 0.0_rk,-1.3_rk, 1.0_rk]
    Xc(11,:) = Xc(7,:)
    Xc(12,:) = Xc(8,:)

    call surface%set(&
        knot1           = knot1,&
        knot2           = knot2,&
        Xc              = Xc,&
        Wc              = Wc,&
        degree          = [2, 1],&
        wrap_parameters = [.true., .false.])
    call surface%create(&
        res1 = 121,&
        res2 = 31)

    write(*,"(a,2(i0,1x))") "Degrees           : ", surface%get_degree()
    write(*,"(a,2(i0,1x))") "Control points    : ", surface%get_nc()
    write(*,"(a,2(l1,1x))") "Parameter wrapping: ", [&
        surface%get_parameter_wrapping(1), surface%get_parameter_wrapping(2)]

    call surface%export_Xc("vtk/surface_knot_periodic_Xc.vtk")
    call surface%export_Xg("vtk/surface_knot_periodic_Xg.vtk")
    call surface%export_Xth_in_Xg("vtk/surface_knot_periodic_Xth.vtk")
    call surface%show(&
        vtkfile_Xc        = "vtk/surface_knot_periodic_Xc.vtk",&
        vtkfile_Xg        = "vtk/surface_knot_periodic_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_knot_periodic_Xth.vtk")

end program surface_knot_periodic
