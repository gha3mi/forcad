program test_surface

    use forcad, only: rk, nurbs_surface
    use forunittest, only: unit_test

    implicit none
    type(nurbs_surface) :: nurbs, bsp
    real(rk) :: knot1(6), knot2(6)
    real(rk), allocatable :: Xc(:,:), Wc(:)
    real(rk), allocatable :: Xg(:,:), Xgb(:,:)
    type(unit_test) :: ut

    Xc = generate_Xc(3, 3, 1.0_rk)

    allocate(Wc(size(Xc, 1)))
    Wc = 1.0_rk
    Wc(2) = 2.0_rk

    knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
    knot2 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    call nurbs%set(knot1, knot2, Xc, Wc)
    call bsp%set(knot1, knot2, Xc)

    call nurbs%create(30, 30)
    call bsp%create(30, 30)

    Xg = nurbs%get_Xg()
    Xgb = bsp%get_Xg()

    call nurbs%insert_knots(1, [0.25_rk, 0.75_rk], [2,1])
    call nurbs%insert_knots(2, [0.25_rk, 0.75_rk], [2,1])

    call bsp%insert_knots(1, [0.25_rk, 0.75_rk], [2,1])
    call bsp%insert_knots(2, [0.25_rk, 0.75_rk], [2,1])

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test 07")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test 08")

    call nurbs%elevate_degree(1, 2)
    call nurbs%elevate_degree(2, 2)

    call bsp%elevate_degree(1, 2)
    call bsp%elevate_degree(2, 2)

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test 09")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test 10")

    call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
    call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

    call bsp%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
    call bsp%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test 11")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test 12")

    call nurbs%finalize()
    call bsp%finalize()

contains

    !-----------------------------------------------------------------------------
    function generate_Xc(num_rows, num_cols, peak_height) result(control_points)
        integer, intent(in) :: num_rows, num_cols
        real(rk), intent(in) :: peak_height
        real(rk), allocatable :: control_points(:,:)
        integer :: i, j
        real(rk) :: x_spacing, y_spacing, x_offset, y_offset
        x_spacing = 1.0_rk / real(num_cols - 1)
        y_spacing = 1.0_rk / real(num_rows - 1)
        x_offset = -0.5_rk
        y_offset = -0.5_rk
        allocate(control_points(num_rows * num_cols, 3))
        do i = 1, num_rows
            do j = 1, num_cols
                control_points((i - 1) * num_cols + j, 1) = x_offset + real(j - 1) * x_spacing
                control_points((i - 1) * num_cols + j, 2) = y_offset + real(i - 1) * y_spacing
                control_points((i - 1) * num_cols + j, 3) = &
                    peak_height * exp(-((control_points((i - 1) * num_cols + j, 1) ** 2) &
                    + (control_points((i - 1) * num_cols + j, 2) ** 2))) + 0.5_rk * peak_height * 0.2_rk
            end do
        end do
    end function
    !-----------------------------------------------------------------------------

end program
