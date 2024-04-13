program test_curve

    use forcad, only: rk, nurbs_curve

    implicit none
    type(nurbs_curve) :: nurbs, bsp
    real(rk), allocatable :: Xc(:,:), Wc(:)
    real(rk), allocatable :: Xg(:,:), Xgb(:,:)
    real(rk) :: knot(6)

    allocate(Xc(3, 3))
    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [0.0_rk, 5.0_rk, 0.0_rk]
    Xc(3,:) = [5.0_rk, 5.0_rk, 0.0_rk]

    allocate(Wc(3))
    Wc = [1.0_rk, 2.0_rk, 0.3_rk]

    knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    call nurbs%set(knot, Xc, Wc)
    call bsp%set(knot, Xc)

    call nurbs%create(res = 20)
    call bsp%create(res = 20)

    Xg = nurbs%get_Xg()
    Xgb = bsp%get_Xg()

    call nurbs%insert_knots([0.25_rk, 0.75_rk], [2,1])
    call bsp%insert_knots([0.25_rk, 0.75_rk], [2,1])

    call nurbs%create()
    call bsp%create()

    print*,'test 01: ', norm2(Xg - nurbs%get_Xg())
    print*,'test 02: ', norm2(Xgb - bsp%get_Xg())

    call nurbs%elevate_degree(2)
    call bsp%elevate_degree(2)

    call nurbs%create()
    call bsp%create()

    print*,'test 03: ', norm2(Xg - nurbs%get_Xg())
    print*,'test 04: ', norm2(Xgb - bsp%get_Xg())

    call nurbs%remove_knots([0.25_rk, 0.75_rk], [2,1])
    call bsp%remove_knots([0.25_rk, 0.75_rk], [2,1])

    call nurbs%create()
    call bsp%create()

    print*,'test 05: ', norm2(Xg - nurbs%get_Xg())
    print*,'test 06: ', norm2(Xgb - bsp%get_Xg())

    call nurbs%finalize()
    call bsp%finalize()

end program
