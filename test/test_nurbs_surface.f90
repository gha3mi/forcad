program test_nurbs_surface

    use forcad, only: rk, nurbs_surface
    use forunittest, only: unit_test

    implicit none
    type(nurbs_surface) :: nurbs, bsp
    real(rk) :: knot1(4), knot2(4), area, areab
    real(rk), allocatable :: Xc(:,:), Wc(:)
    real(rk), allocatable :: Xg(:,:), Xgb(:,:)
    integer :: id
    real(rk), allocatable :: nearest_Xg(:), nearest_Xt(:)
    type(unit_test) :: ut

    allocate(Xc(4, 3))
    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [5.0_rk, 0.0_rk, 0.0_rk]
    Xc(3,:) = [0.0_rk, 5.0_rk, 0.0_rk]
    Xc(4,:) = [5.0_rk, 5.0_rk, 0.0_rk]

    allocate(Wc(size(Xc, 1)))
    Wc = 1.0_rk
    Wc(2) = 2.0_rk

    knot1 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
    knot2 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

    call nurbs%set(knot1, knot2, Xc, Wc)
    call bsp%set(knot1, knot2, Xc)

    call nurbs%create(30, 30)
    call bsp%create(30, 30)

    call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 01")
    call ut%check(res=nearest_Xt, expected=[0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 02")
    call ut%check(res=id, expected=1, msg="test_nurbs_surface: 03")

    call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 04")
    call ut%check(res=nearest_Xt, expected=[0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 05")
    call ut%check(res=id, expected=1, msg="test_nurbs_surface: 06")

    call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 07")
    call ut%check(res=nearest_Xt, expected=[0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 08")

    call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 09")
    call ut%check(res=nearest_Xt, expected=[0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 10")

    Xg = nurbs%get_Xg()
    Xgb = bsp%get_Xg()

    call nurbs%insert_knots(1, [0.25_rk, 0.75_rk], [1,1])
    call nurbs%insert_knots(2, [0.25_rk, 0.75_rk], [1,1])

    call bsp%insert_knots(1, [0.25_rk, 0.75_rk], [1,1])
    call bsp%insert_knots(2, [0.25_rk, 0.75_rk], [1,1])

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 11")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 12")

    call nurbs%elevate_degree(1, 2)
    call nurbs%elevate_degree(2, 2)

    call bsp%elevate_degree(1, 2)
    call bsp%elevate_degree(2, 2)

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 13")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 14")

    call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
    call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

    call bsp%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
    call bsp%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 15")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 16")

    call nurbs%cmp_area(area)
    call bsp%cmp_area(areab)

    call ut%check(res=area,  expected=25.0_rk, tol=1e-5_rk, msg="test_nurbs_surface: 17")
    call ut%check(res=areab, expected=25.0_rk, tol=1e-5_rk, msg="test_nurbs_surface: 18")

    call nurbs%finalize()
    call bsp%finalize()
    deallocate(Xc, Wc, Xg, Xgb)

end program
