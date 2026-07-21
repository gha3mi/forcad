!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Construct a rational annular volume periodic in its first parameter direction.
program volume_knot_periodic

    use forcad, only: rk, nurbs_volume

    implicit none

    type(nurbs_volume) :: volume
    real(rk) :: Xc(24,3), Wc(24), ring(6,2), radius
    real(rk), parameter :: knot1(9) = [&
        0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
    real(rk), parameter :: knot2(4) = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
    real(rk), parameter :: knot3(4) = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
    real(rk), parameter :: ring_Wc(6) = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
    integer :: i1, i2, i3, n

    ring(1,:) = [ 1.0_rk, 0.0_rk]
    ring(2,:) = [ 0.0_rk, 1.0_rk]
    ring(3,:) = [-1.0_rk, 0.0_rk]
    ring(4,:) = [ 0.0_rk,-1.0_rk]
    ring(5,:) = ring(1,:)
    ring(6,:) = ring(2,:)

    do i3 = 1, 2
        do i2 = 1, 2
            radius = 1.0_rk + 0.6_rk*real(i2 - 1, rk)
            do i1 = 1, 6
                n = i1 + 6*(i2 - 1) + 12*(i3 - 1)
                Xc(n,:) = [radius*ring(i1,1), radius*ring(i1,2), real(i3 - 1, rk)]
                Wc(n) = ring_Wc(i1)
            end do
        end do
    end do

    call volume%set(&
        knot1           = knot1,&
        knot2           = knot2,&
        knot3           = knot3,&
        Xc              = Xc,&
        Wc              = Wc,&
        degree          = [2, 1, 1],&
        wrap_parameters = [.true., .false., .false.])
    call volume%create(&
        res1 = 81,&
        res2 = 13,&
        res3 = 13)

    write(*,"(a,3(i0,1x))") "Degrees           : ", volume%get_degree()
    write(*,"(a,3(i0,1x))") "Control points    : ", volume%get_nc()
    write(*,"(a,3(l1,1x))") "Parameter wrapping: ", [&
        volume%get_parameter_wrapping(1), volume%get_parameter_wrapping(2),&
        volume%get_parameter_wrapping(3)]

    call volume%export_Xc("vtk/volume_knot_periodic_Xc.vtk")
    call volume%export_Xg("vtk/volume_knot_periodic_Xg.vtk")
    call volume%export_Xth_in_Xg("vtk/volume_knot_periodic_Xth.vtk")
    call volume%show(&
        vtkfile_Xc        = "vtk/volume_knot_periodic_Xc.vtk",&
        vtkfile_Xg        = "vtk/volume_knot_periodic_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_knot_periodic_Xth.vtk")

end program volume_knot_periodic
