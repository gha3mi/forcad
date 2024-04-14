program test_volume

    use forcad, only: rk, nurbs_volume

    implicit none
    type(nurbs_volume) :: nurbs, bsp
    real(rk), allocatable :: Xc(:,:), Wc(:)
    real(rk), allocatable :: Xg(:,:), Xgb(:,:)
    real(rk) :: knot1(4), knot2(4), knot3(4)

    Xc = generate_Xc(5.0_rk)

    allocate(Wc(size(Xc,1)), source=1.0_rk)
    Wc(2) = 5.0_rk

    knot1 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
    knot2 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
    knot3 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

    call nurbs%set(knot1, knot2, knot3, Xc, Wc)
    call bsp%set(knot1, knot2, knot3, Xc)

    call nurbs%create(20, 20, 20)
    call bsp%create(20, 20, 20)

    Xg = nurbs%get_Xg()
    Xgb = bsp%get_Xg()

    call nurbs%insert_knots(1, [0.25_rk, 0.75_rk], [1,1])
    call nurbs%insert_knots(2, [0.25_rk, 0.75_rk], [1,1])
    call nurbs%insert_knots(3, [0.25_rk, 0.75_rk], [1,1])

    call bsp%insert_knots(1, [0.25_rk, 0.75_rk], [1,1])
    call bsp%insert_knots(2, [0.25_rk, 0.75_rk], [1,1])
    call bsp%insert_knots(3, [0.25_rk, 0.75_rk], [1,1])

    call nurbs%create()
    call bsp%create()

    print*,'test 13: ', norm2(Xg - nurbs%get_Xg())
    print*,'test 14: ', norm2(Xgb - bsp%get_Xg())

    call nurbs%elevate_degree(1, 2)
    call nurbs%elevate_degree(2, 2)
    call nurbs%elevate_degree(3, 2)

    call bsp%elevate_degree(1, 2)
    call bsp%elevate_degree(2, 2)
    call bsp%elevate_degree(3, 2)

    call nurbs%create()
    call bsp%create()

    print*,'test 15: ', norm2(Xg - nurbs%get_Xg())
    print*,'test 16: ', norm2(Xgb - bsp%get_Xg())

    call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [1,1])
    call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [1,1])
    call nurbs%remove_knots(3, [0.25_rk, 0.75_rk], [1,1])

    call bsp%remove_knots(1, [0.25_rk, 0.75_rk], [1,1])
    call bsp%remove_knots(2, [0.25_rk, 0.75_rk], [1,1])
    call bsp%remove_knots(3, [0.25_rk, 0.75_rk], [1,1])

    call nurbs%create()
    call bsp%create()

    print*,'test 17: ', norm2(Xg - nurbs%get_Xg())
    print*,'test 18: ', norm2(Xgb - bsp%get_Xg())

    call nurbs%finalize()
    call bsp%finalize()

contains

    !-----------------------------------------------------------------------------
    function generate_Xc(L) result(control_points)
        implicit none
        real(rk), intent(in) :: L
        real(rk), allocatable :: control_points(:,:)
        real(rk) :: L2
        L2 = L / 2.0_rk
        allocate(control_points(8, 3))
        control_points(1,:) = [ L2, -L2,  L2]
        control_points(2,:) = [ L2, -L2, -L2]
        control_points(3,:) = [-L2, -L2,  L2]
        control_points(4,:) = [-L2, -L2, -L2]
        control_points(5,:) = [ L2,  L2,  L2]
        control_points(6,:) = [ L2,  L2, -L2]
        control_points(7,:) = [-L2,  L2,  L2]
        control_points(8,:) = [-L2,  L2, -L2]
    end function
    !-----------------------------------------------------------------------------

end program
