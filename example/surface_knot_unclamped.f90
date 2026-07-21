!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program surface_knot_unclamped

    use forcad, only: rk, nurbs_surface

    implicit none

    type(nurbs_surface) :: surface
    real(rk) :: knot1(7), knot2(7)
    real(rk), allocatable :: Xc(:,:), Wc(:)
    integer :: i, j, id
    real(rk) :: x, y

    knot1 = [-0.4_rk, -0.1_rk, 0.2_rk, 0.7_rk, 1.4_rk, 1.9_rk, 2.3_rk]
    knot2 = [-0.2_rk, 0.1_rk, 0.5_rk, 0.8_rk, 1.5_rk, 2.0_rk, 2.6_rk]

    allocate(Xc(16, 3), Wc(16))
    do concurrent (j = 1:4, i = 1:4) local(id, x, y)
        id = (j - 1)*4 + i
        x = real(i - 1, rk)/3.0_rk
        y = real(j - 1, rk)/3.0_rk
        Xc(id,:) = [x, y, 0.15_rk*sin(6.283185307179586_rk*x*y)]
        Wc(id) = 1.0_rk + 0.05_rk*real(i + j - 2, rk)/6.0_rk
    end do

    call surface%set(knot1, knot2, Xc, Wc, degree=[2, 2])
    call surface%create(64, 48)
    call surface%export_Xc("vtk/surface_knot_unclamped_Xc.vtk")
    call surface%export_Xg("vtk/surface_knot_unclamped_Xg.vtk")
    call surface%export_Xth_in_Xg("vtk/surface_knot_unclamped_Xthg.vtk")
    call surface%show(&
        "vtk/surface_knot_unclamped_Xc.vtk", &
        "vtk/surface_knot_unclamped_Xg.vtk", &
        "vtk/surface_knot_unclamped_Xthg.vtk")

    write(*,"(a,2(i0,1x))") "degree:", surface%get_degree()

    call surface%finalize()

end program surface_knot_unclamped
