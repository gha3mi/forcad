program test_curve

    use forcad, only: rk, nurbs_curve
    use forunittest, only: unit_test

    implicit none
    type(nurbs_curve) :: nurbs, bsp
    real(rk), allocatable :: Xc(:,:), Wc(:)
    real(rk), allocatable :: Xg(:,:), Xgb(:,:)
    real(rk) :: knot(6)
    type(unit_test) :: ut

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

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test 01")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test 02")

    call nurbs%elevate_degree(2)
    call bsp%elevate_degree(2)

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test 03")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test 04")

    call nurbs%remove_knots([0.25_rk, 0.75_rk], [2,1])
    call bsp%remove_knots([0.25_rk, 0.75_rk], [2,1])

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test 05")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test 06")

    call nurbs%finalize()
    call bsp%finalize()

end program
