!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program volume_knot_clamped

    use forcad, only: rk, nurbs_volume

    implicit none

    type(nurbs_volume) :: volume
    real(rk) :: knot1(7), knot2(7), knot3(7)
    real(rk), allocatable :: Xc(:,:), Wc(:)
    integer :: i, j, k, id
    real(rk) :: x, y, z

    knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.35_rk, 1.0_rk, 1.0_rk, 1.0_rk]
    knot2 = [0.0_rk, 0.0_rk, 0.0_rk, 0.55_rk, 1.0_rk, 1.0_rk, 1.0_rk]
    knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.70_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    allocate(Xc(64, 3), Wc(64))
    do concurrent (k = 1:4, j = 1:4, i = 1:4) local(id, x, y, z)
        id = (k - 1)*16 + (j - 1)*4 + i
        x = real(i - 1, rk)/3.0_rk
        y = real(j - 1, rk)/3.0_rk
        z = real(k - 1, rk)/3.0_rk
        Xc(id,:) = [x, y, z + 0.10_rk*sin(3.141592653589793_rk*x)*sin(3.141592653589793_rk*y)]
        Wc(id) = 1.0_rk
    end do

    call volume%set(knot1, knot2, knot3, Xc, Wc, degree=[2, 2, 2])
    call volume%create(24, 20, 16)
    call volume%export_Xc("vtk/volume_knot_clamped_Xc.vtk")
    call volume%export_Xg("vtk/volume_knot_clamped_Xg.vtk")
    call volume%export_Xth_in_Xg("vtk/volume_knot_clamped_Xthg.vtk")
    call volume%show("vtk/volume_knot_clamped_Xc.vtk", "vtk/volume_knot_clamped_Xg.vtk", "vtk/volume_knot_clamped_Xthg.vtk")

    write(*,"(a,3(i0,1x))") "degree:", volume%get_degree()

    call volume%finalize()

end program volume_knot_clamped
