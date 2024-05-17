program test_curve

    use forcad, only: rk, nurbs_curve
    use forunittest, only: unit_test

    implicit none
    type(nurbs_curve) :: nurbs, bsp
    real(rk), allocatable :: Xc(:,:), Wc(:)
    real(rk), allocatable :: Xg(:,:), Xgb(:,:)
    real(rk) :: knot(6)
    integer, allocatable :: elemConn(:,:)
    real(rk), allocatable :: Tgc(:,:), dTgc(:,:)
    integer :: i
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

    call nurbs%create(res = 23)
    call bsp%create(res = 23)

    Xg = nurbs%get_Xg()
    Xgb = bsp%get_Xg()

    call nurbs%set([0.0_rk, 1.0_rk], 2, [-1, -1], Xc, Wc)
    call bsp%set([0.0_rk, 1.0_rk], 2, [-1, -1], Xc)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_curve: 01")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_curve: 02")

    call nurbs%set(Xc, Wc)
    call bsp%set(Xc)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_curve: 03")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_curve: 04")

    call nurbs%create(Xt = nurbs%get_Xt())
    call bsp%create(Xt = bsp%get_Xt())

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_curve: 05")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_curve: 06")

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_curve: 07")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_curve: 08")

    call ut%check(res=nurbs%get_Xc(1), expected=Xc(1,:),  tol=1e-5_rk, msg="test_curve: 09")
    call ut%check(res=bsp%get_Xc(1),   expected=Xc(1,:), tol=1e-5_rk, msg="test_curve: 10")

    call ut%check(res=nurbs%get_Xc(1,1), expected=Xc(1,1),  tol=1e-5_rk, msg="test_curve: 11")
    call ut%check(res=bsp%get_Xc(1,1),   expected=Xc(1,1), tol=1e-5_rk, msg="test_curve: 12")

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_curve: 13")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_curve: 14")

    call ut%check(res=nurbs%get_Xg(1), expected=Xg(1,:),  tol=1e-5_rk, msg="test_curve: 15")
    call ut%check(res=bsp%get_Xg(1),   expected=Xgb(1,:), tol=1e-5_rk, msg="test_curve: 16")

    call ut%check(res=nurbs%get_Xg(1,1), expected=Xg(1,1),  tol=1e-5_rk, msg="test_curve: 17")
    call ut%check(res=bsp%get_Xg(1,1),   expected=Xgb(1,1), tol=1e-5_rk, msg="test_curve: 18")

    call ut%check(res=nurbs%get_Wc(), expected=Wc,  tol=1e-5_rk, msg="test_curve: 19")

    call ut%check(res=nurbs%get_Wc(1), expected=Wc(1),  tol=1e-5_rk, msg="test_curve: 20")

    call ut%check(res=nurbs%get_knot(), expected=knot,  tol=1e-5_rk, msg="test_curve: 21")
    call ut%check(res=bsp%get_knot(), expected=knot,  tol=1e-5_rk, msg="test_curve: 22")

    call ut%check(res=nurbs%get_knot(1), expected=knot(1),  tol=1e-5_rk, msg="test_curve: 23")
    call ut%check(res=bsp%get_knot(1), expected=knot(1),  tol=1e-5_rk, msg="test_curve: 24")

    call ut%check(res=nurbs%get_ng(), expected=size(Xg,1), msg="test_curve: 25")
    call ut%check(res=bsp%get_ng(), expected=size(Xgb,1), msg="test_curve: 26")

    call ut%check(res=nurbs%get_degree(), expected=2, msg="test_curve: 27")
    call ut%check(res=bsp%get_degree(), expected=2, msg="test_curve: 28")

    call ut%check(res=nurbs%get_multiplicity(), expected=[3,3], msg="test_curve: 29")
    call ut%check(res=bsp%get_multiplicity(), expected=[3,3], msg="test_curve: 30")

    call ut%check(res=nurbs%get_continuity(), expected=[-1,-1], msg="test_curve: 31")
    call ut%check(res=bsp%get_continuity(), expected=[-1,-1], msg="test_curve: 32")

    call ut%check(res=nurbs%get_nc(), expected=size(Xc,1), msg="test_curve: 33")
    call ut%check(res=bsp%get_nc(), expected=size(Xc,1), msg="test_curve: 34")

    call nurbs%cmp_nc()
    call bsp%cmp_nc()

    call ut%check(res=nurbs%get_nc(), expected=size(Xc,1), msg="test_curve: 35")
    call ut%check(res=bsp%get_nc(), expected=size(Xc,1), msg="test_curve: 36")


    elemConn = nurbs%cmp_elem_Xc_vis(2)
    call nurbs%set_elem_Xc_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xc_vis(), expected=elemConn, msg="test_curve: 37")
    deallocate(elemConn)
    elemConn = nurbs%cmp_elem_Xc_vis()
    call nurbs%set_elem_Xc_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xc_vis(), expected=elemConn, msg="test_curve: 38")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xc_vis(2)
    call bsp%set_elem_Xc_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xc_vis(), expected=elemConn, msg="test_curve: 39")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xc_vis()
    call bsp%set_elem_Xc_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xc_vis(), expected=elemConn, msg="test_curve: 40")
    deallocate(elemConn)

    elemConn = nurbs%cmp_elem_Xg_vis(2)
    call nurbs%set_elem_Xg_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xg_vis(), expected=elemConn, msg="test_curve: 41")
    deallocate(elemConn)
    elemConn = nurbs%cmp_elem_Xg_vis()
    call nurbs%set_elem_Xg_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xg_vis(), expected=elemConn, msg="test_curve: 42")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xg_vis(2)
    call bsp%set_elem_Xg_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xg_vis(), expected=elemConn, msg="test_curve: 43")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xg_vis()
    call bsp%set_elem_Xg_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xg_vis(), expected=elemConn, msg="test_curve: 44")
    deallocate(elemConn)

    elemConn = nurbs%cmp_elem()
    call nurbs%set_elem(elemConn)
    call ut%check(res=nurbs%get_elem(), expected=elemConn, msg="test_curve: 45")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem()
    call bsp%set_elem(elemConn)
    call ut%check(res=bsp%get_elem(), expected=elemConn, msg="test_curve: 46")
    deallocate(elemConn)

    call nurbs%modify_Xc(Xc(1,1), 1,1)
    call bsp%modify_Xc(Xc(1,1), 1,1)

    call nurbs%modify_Wc(Wc(1),1)

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_curve: 47")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_curve: 48")

    call nurbs%basis(res=23, Tgc=Tgc)
    call bsp%basis(res=23, Tgc=Tgc)

    call nurbs%basis(Xt=[(real(i-1, rk) / real(23-1, rk), i=1, 23)], Tgc=Tgc)
    call bsp%basis(Xt=[(real(i-1, rk) / real(23-1, rk), i=1, 23)], Tgc=Tgc)

    call nurbs%derivative(res=23, dTgc=dTgc)
    call bsp%derivative(res=23, dTgc=dTgc)

    call nurbs%derivative(Xt=[(real(i-1, rk) / real(23-1, rk), i=1, 23)], dTgc=dTgc)
    call bsp%derivative(Xt=[(real(i-1, rk) / real(23-1, rk), i=1, 23)], dTgc=dTgc)

    call nurbs%rotate_Xc(45.0_rk, 0.0_rk, 0.0_rk)
    call nurbs%rotate_Xc(-45.0_rk, 0.0_rk, 0.0_rk)

    call bsp%rotate_Xc(45.0_rk, 0.0_rk, 0.0_rk)
    call bsp%rotate_Xc(-45.0_rk, 0.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_curve: 49")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_curve: 50")

    call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
    call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

    call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
    call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_curve: 51")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_curve: 52")

    call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
    call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

    call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
    call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_curve: 53")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_curve: 54")

    call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
    call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

    call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
    call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_curve: 55")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_curve: 56")

    call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
    call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

    call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
    call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_curve: 57")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_curve: 58")

    call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
    call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

    call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
    call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_curve: 59")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_curve: 60")

    call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
    call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

    call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
    call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_curve: 61")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_curve: 62")

    call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
    call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

    call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
    call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_curve: 63")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_curve: 64")

    call nurbs%export_Xc("vtk/test_curve_Xc.vtk")
    call nurbs%export_Xg("vtk/test_curve_Xg.vtk")

    call bsp%export_Xc("vtk/test_curve_Xc.vtk")
    call bsp%export_Xg("vtk/test_curve_Xg.vtk")

    call nurbs%insert_knots([0.25_rk, 0.75_rk], [2,1])
    call bsp%insert_knots([0.25_rk, 0.75_rk], [2,1])

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_curve: 65")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_curve: 66")

    call nurbs%elevate_degree(2)
    call bsp%elevate_degree(2)

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_curve: 67")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_curve: 68")

    call nurbs%remove_knots([0.25_rk, 0.75_rk], [2,1])
    call bsp%remove_knots([0.25_rk, 0.75_rk], [2,1])

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_curve: 69")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_curve: 70")

    call nurbs%finalize()
    call bsp%finalize()
    deallocate(Xc, Wc, Xg, Xgb)

end program
