!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program surface_knot_clamped

    use forcad, only: rk, nurbs_surface

    implicit none

    type(nurbs_surface) :: surface
    real(rk) :: knot1(7), knot2(7)
    real(rk), allocatable :: Xc(:,:), Wc(:)
    integer :: i, j, id
    real(rk) :: x, y

    knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.35_rk, 1.0_rk, 1.0_rk, 1.0_rk]
    knot2 = [0.0_rk, 0.0_rk, 0.0_rk, 0.65_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    allocate(Xc(16, 3), Wc(16))
    do concurrent (j = 1:4, i = 1:4) local(id, x, y)
        id = (j - 1)*4 + i
        x = real(i - 1, rk)/3.0_rk
        y = real(j - 1, rk)/3.0_rk
        Xc(id,:) = [x, y, 0.20_rk*sin(3.141592653589793_rk*x)*cos(3.141592653589793_rk*y)]
        Wc(id) = 1.0_rk
    end do

    call surface%set(knot1, knot2, Xc, Wc, degree=[2, 2])
    call surface%create(64, 48)
    call surface%export_Xc("vtk/surface_knot_clamped_Xc.vtk")
    call surface%export_Xg("vtk/surface_knot_clamped_Xg.vtk")
    call surface%export_Xth_in_Xg("vtk/surface_knot_clamped_Xthg.vtk")
    call surface%show("vtk/surface_knot_clamped_Xc.vtk", "vtk/surface_knot_clamped_Xg.vtk", "vtk/surface_knot_clamped_Xthg.vtk")

    write(*,"(a,2(i0,1x))") "degree:", surface%get_degree()

    call surface%finalize()

end program surface_knot_clamped
