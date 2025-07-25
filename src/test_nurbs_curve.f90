program test_nurbs_curve

    use forcad, only: rk, nurbs_curve
    use forunittest, only: unit_test

    implicit none
    type(nurbs_curve) :: nurbs, bsp
    real(rk), allocatable :: Xc(:,:), Wc(:)
    real(rk), allocatable :: Xg(:,:), Xgb(:,:)
    real(rk) :: knot(6)
    integer, allocatable :: elemConn(:,:)
    real(rk), allocatable :: Tgc(:,:), dTgc(:,:), Tgcb(:,:), dTgcb(:,:), d2Tgc(:,:), d2Tgcb(:,:)
    real(rk), allocatable :: Tgc1(:), dTgc1(:), Tgc1b(:), dTgc1b(:), d2Tgc1(:), d2Tgc1b(:)
    integer :: i, id
    real(rk), allocatable :: nearest_Xg(:)
    real(rk) :: nearest_Xt, length, lengthb
    type(unit_test) :: ut

    allocate(Xc(3, 3))
    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
    Xc(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]

    allocate(Wc(3))
    Wc = [1.0_rk, 0.9_rk, 0.8_rk]

    knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    call nurbs%set(knot, Xc, Wc)
    call bsp%set(knot, Xc)

    call nurbs%set(degree=nurbs%get_degree(), nc=nurbs%get_nc(), Xc=nurbs%get_Xc(), Wc=nurbs%get_Wc())
    call bsp%set(degree=nurbs%get_degree(), nc=nurbs%get_nc(), Xc=nurbs%get_Xc())

    call nurbs%create(res = 23)
    call bsp%create(res = 23)

    call nurbs%export_Xc("vtk/test_nurbs_curve_Xc.vtk")
    call bsp%export_Xc("vtk/test_bsp_curve_Xc.vtk")

    call nurbs%export_Xg("vtk/test_nurbs_curve_Xg.vtk")
    call bsp%export_Xg("vtk/test_bsp_curve_Xg.vtk")

    call nurbs%export_Xth("vtk/test_nurbs_curve_Xth.vtk")
    call bsp%export_Xth("vtk/test_bsp_curve_Xth.vtk")

    call nurbs%export_iges('iges/test_nurbs_curve.iges')
    call bsp%export_iges('iges/test_bsp_curve.iges')

    call nurbs%cmp_length(length)
    call bsp%cmp_length(lengthb)

    call ut%check(res=length, expected=2.0_rk,  tol=1e-5_rk, msg="test_nurbs_curve: 01")
    call ut%check(res=lengthb, expected=2.0_rk,  tol=1e-5_rk, msg="test_nurbs_curve: 02")

    call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_curve: 03")
    call ut%check(res=nearest_Xt, expected=0.0_rk,  tol=1e-5_rk, msg="test_nurbs_curve: 04")
    call ut%check(res=id, expected=1, msg="test_nurbs_curve: 05")

    call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_curve: 06")
    call ut%check(res=nearest_Xt, expected=0.0_rk,  tol=1e-5_rk, msg="test_nurbs_curve: 07")
    call ut%check(res=id, expected=1, msg="test_nurbs_curve: 08")

    call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_curve: 09")
    call ut%check(res=nearest_Xt, expected=0.0_rk,  tol=1e-5_rk, msg="test_nurbs_curve: 10")

    call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_curve: 11")
    call ut%check(res=nearest_Xt, expected=0.0_rk,  tol=1e-5_rk, msg="test_nurbs_curve: 12")

    Xg = nurbs%get_Xg()
    Xgb = bsp%get_Xg()

    call nurbs%set([0.0_rk, 1.0_rk], 2, [-1, -1], Xc, Wc)
    call bsp%set([0.0_rk, 1.0_rk], 2, [-1, -1], Xc)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_curve: 13")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_curve: 14")

    call nurbs%set(Xc, Wc)
    call bsp%set(Xc)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_curve: 15")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_curve: 16")

    call nurbs%create(Xt = nurbs%get_Xt())
    call bsp%create(Xt = bsp%get_Xt())

    call nurbs%export_Xc("vtk/test_nurbs_curve_Xc.vtk")
    call bsp%export_Xc("vtk/test_bsp_curve_Xc.vtk")

    call nurbs%export_Xg("vtk/test_nurbs_curve_Xg.vtk")
    call bsp%export_Xg("vtk/test_bsp_curve_Xg.vtk")

    call nurbs%export_Xth("vtk/test_nurbs_curve_Xth.vtk")
    call bsp%export_Xth("vtk/test_bsp_curve_Xth.vtk")

    call nurbs%export_iges('iges/test_nurbs_curve.iges')
    call bsp%export_iges('iges/test_bsp_curve.iges')

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_curve: 17")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_curve: 18")

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_curve: 19")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_curve: 20")

    call ut%check(res=nurbs%get_Xc(1), expected=Xc(1,:),  tol=1e-5_rk, msg="test_nurbs_curve: 21")
    call ut%check(res=bsp%get_Xc(1),   expected=Xc(1,:), tol=1e-5_rk, msg="test_nurbs_curve: 22")

    call ut%check(res=nurbs%get_Xc(1,1), expected=Xc(1,1),  tol=1e-5_rk, msg="test_nurbs_curve: 23")
    call ut%check(res=bsp%get_Xc(1,1),   expected=Xc(1,1), tol=1e-5_rk, msg="test_nurbs_curve: 24")

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_curve: 25")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_curve: 26")

    call ut%check(res=nurbs%get_Xg(1), expected=Xg(1,:),  tol=1e-5_rk, msg="test_nurbs_curve: 27")
    call ut%check(res=bsp%get_Xg(1),   expected=Xgb(1,:), tol=1e-5_rk, msg="test_nurbs_curve: 28")

    call ut%check(res=nurbs%get_Xg(1,1), expected=Xg(1,1),  tol=1e-5_rk, msg="test_nurbs_curve: 29")
    call ut%check(res=bsp%get_Xg(1,1),   expected=Xgb(1,1), tol=1e-5_rk, msg="test_nurbs_curve: 30")

    call ut%check(res=nurbs%get_Wc(), expected=Wc,  tol=1e-5_rk, msg="test_nurbs_curve: 31")

    call ut%check(res=nurbs%get_Wc(1), expected=Wc(1),  tol=1e-5_rk, msg="test_nurbs_curve: 32")

    call ut%check(res=nurbs%get_knot(), expected=knot,  tol=1e-5_rk, msg="test_nurbs_curve: 33")
    call ut%check(res=bsp%get_knot(), expected=knot,  tol=1e-5_rk, msg="test_nurbs_curve: 34")

    call ut%check(res=nurbs%get_knot(1), expected=knot(1),  tol=1e-5_rk, msg="test_nurbs_curve: 35")
    call ut%check(res=bsp%get_knot(1), expected=knot(1),  tol=1e-5_rk, msg="test_nurbs_curve: 36")

    call ut%check(res=nurbs%get_ng(), expected=size(Xg,1), msg="test_nurbs_curve: 37")
    call ut%check(res=bsp%get_ng(), expected=size(Xgb,1), msg="test_nurbs_curve: 38")

    call ut%check(res=nurbs%get_degree(), expected=2, msg="test_nurbs_curve: 39")
    call ut%check(res=bsp%get_degree(), expected=2, msg="test_nurbs_curve: 40")

    call ut%check(res=nurbs%get_multiplicity(), expected=[3,3], msg="test_nurbs_curve: 41")
    call ut%check(res=bsp%get_multiplicity(), expected=[3,3], msg="test_nurbs_curve: 42")

    call ut%check(res=nurbs%get_continuity(), expected=[-1,-1], msg="test_nurbs_curve: 43")
    call ut%check(res=bsp%get_continuity(), expected=[-1,-1], msg="test_nurbs_curve: 44")

    call ut%check(res=nurbs%get_nc(), expected=size(Xc,1), msg="test_nurbs_curve: 45")
    call ut%check(res=bsp%get_nc(), expected=size(Xc,1), msg="test_nurbs_curve: 46")

    call nurbs%cmp_nc()
    call bsp%cmp_nc()

    call ut%check(res=nurbs%get_nc(), expected=size(Xc,1), msg="test_nurbs_curve: 47")
    call ut%check(res=bsp%get_nc(), expected=size(Xc,1), msg="test_nurbs_curve: 48")

    elemConn = nurbs%cmp_elem_Xc_vis(2)
    call nurbs%set_elem_Xc_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xc_vis(), expected=elemConn, msg="test_nurbs_curve: 49")
    deallocate(elemConn)
    elemConn = nurbs%cmp_elem_Xc_vis()
    call nurbs%set_elem_Xc_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xc_vis(), expected=elemConn, msg="test_nurbs_curve: 50")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xc_vis(2)
    call bsp%set_elem_Xc_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xc_vis(), expected=elemConn, msg="test_nurbs_curve: 51")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xc_vis()
    call bsp%set_elem_Xc_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xc_vis(), expected=elemConn, msg="test_nurbs_curve: 52")
    deallocate(elemConn)

    elemConn = nurbs%cmp_elem_Xg_vis(2)
    call nurbs%set_elem_Xg_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xg_vis(), expected=elemConn, msg="test_nurbs_curve: 53")
    deallocate(elemConn)
    elemConn = nurbs%cmp_elem_Xg_vis()
    call nurbs%set_elem_Xg_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xg_vis(), expected=elemConn, msg="test_nurbs_curve: 54")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xg_vis(2)
    call bsp%set_elem_Xg_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xg_vis(), expected=elemConn, msg="test_nurbs_curve: 55")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xg_vis()
    call bsp%set_elem_Xg_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xg_vis(), expected=elemConn, msg="test_nurbs_curve: 56")
    deallocate(elemConn)

    elemConn = nurbs%cmp_elem()
    call nurbs%set_elem(elemConn)
    call ut%check(res=nurbs%get_elem(), expected=elemConn, msg="test_nurbs_curve: 57")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem()
    call bsp%set_elem(elemConn)
    call ut%check(res=bsp%get_elem(), expected=elemConn, msg="test_nurbs_curve: 58")
    deallocate(elemConn)

    call nurbs%modify_Xc(Xc(1,1), 1,1)
    call bsp%modify_Xc(Xc(1,1), 1,1)

    call nurbs%modify_Wc(Wc(1),1)

    call nurbs%create()
    call bsp%create()

    call nurbs%export_Xc("vtk/test_nurbs_curve_Xc.vtk")
    call bsp%export_Xc("vtk/test_bsp_curve_Xc.vtk")

    call nurbs%export_Xg("vtk/test_nurbs_curve_Xg.vtk")
    call bsp%export_Xg("vtk/test_bsp_curve_Xg.vtk")

    call nurbs%export_Xth("vtk/test_nurbs_curve_Xth.vtk")
    call bsp%export_Xth("vtk/test_bsp_curve_Xth.vtk")

    call nurbs%export_iges('iges/test_nurbs_curve.iges')
    call bsp%export_iges('iges/test_bsp_curve.iges')

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_curve: 59")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_curve: 60")

    call nurbs%basis(res=23, Tgc=Tgc)
    call bsp%basis(res=23, Tgc=Tgc)

    call nurbs%basis(Xt=0.0_rk, Tgc=Tgc1)
    call bsp%basis(Xt=0.0_rk, Tgc=Tgc1b)

    call nurbs%basis(Xt=[(real(i-1, rk) / real(23-1, rk), i=1, 23)], Tgc=Tgc)
    call bsp%basis(Xt=[(real(i-1, rk) / real(23-1, rk), i=1, 23)], Tgc=Tgc)

    call nurbs%derivative(res=23, dTgc=dTgc, Tgc=Tgc)
    call bsp%derivative(res=23, dTgc=dTgcb, Tgc=Tgcb)

    call nurbs%derivative(Xt=[(real(i-1, rk) / real(23-1, rk), i=1, 23)], dTgc=dTgc, Tgc=Tgc)
    call bsp%derivative(Xt=[(real(i-1, rk) / real(23-1, rk), i=1, 23)], dTgc=dTgcb, Tgc=Tgcb)

    call nurbs%derivative(Xt=0.0_rk, dTgc=dTgc1, Tgc=Tgc1)
    call bsp%derivative(Xt=0.0_rk, dTgc=dTgc1b, Tgc=Tgc1b)

    call nurbs%derivative(Xt=0.0_rk, dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
    call bsp%derivative(Xt=0.0_rk, dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

    call nurbs%derivative2(res=23, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
    call bsp%derivative2(res=23, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

    call nurbs%derivative2(Xt=[(real(i-1, rk) / real(23-1, rk), i=1, 23)], d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
    call bsp%derivative2(Xt=[(real(i-1, rk) / real(23-1, rk), i=1, 23)], d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

    call nurbs%derivative2(Xt=0.0_rk, d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
    call bsp%derivative2(Xt=0.0_rk, d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

    call nurbs%derivative2(Xt=0.0_rk, d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
    call bsp%derivative2(Xt=0.0_rk, d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

    call nurbs%rotate_Xc(45.0_rk, 0.0_rk, 0.0_rk)
    call nurbs%rotate_Xc(-45.0_rk, 0.0_rk, 0.0_rk)

    call bsp%rotate_Xc(45.0_rk, 0.0_rk, 0.0_rk)
    call bsp%rotate_Xc(-45.0_rk, 0.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_curve: 61")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_curve: 62")

    call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
    call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

    call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
    call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_curve: 63")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_curve: 64")

    call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
    call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

    call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
    call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_curve: 65")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_curve: 66")

    call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
    call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

    call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
    call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_curve: 67")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_curve: 68")

    call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
    call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

    call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
    call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_curve: 69")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_curve: 70")

    call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
    call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

    call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
    call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_curve: 71")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_curve: 72")

    call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
    call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

    call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
    call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_curve: 73")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_curve: 74")

    call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
    call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

    call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
    call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_curve: 75")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_curve: 76")

    call nurbs%export_Xc("vtk/test_nurbs_curve_Xc.vtk")
    call nurbs%export_Xg("vtk/test_nurbs_curve_Xg.vtk")

    call bsp%export_Xc("vtk/test_nurbs_curve_Xc.vtk")
    call bsp%export_Xg("vtk/test_nurbs_curve_Xg.vtk")

    call nurbs%insert_knots([0.25_rk, 0.75_rk], [2,1])
    call bsp%insert_knots([0.25_rk, 0.75_rk], [2,1])

    call nurbs%create()
    call bsp%create()

    call nurbs%export_Xc("vtk/test_nurbs_curve_Xc.vtk")
    call bsp%export_Xc("vtk/test_bsp_curve_Xc.vtk")

    call nurbs%export_Xg("vtk/test_nurbs_curve_Xg.vtk")
    call bsp%export_Xg("vtk/test_bsp_curve_Xg.vtk")

    call nurbs%export_Xth("vtk/test_nurbs_curve_Xth.vtk")
    call bsp%export_Xth("vtk/test_bsp_curve_Xth.vtk")

    call nurbs%export_iges('iges/test_nurbs_curve.iges')
    call bsp%export_iges('iges/test_bsp_curve.iges')

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_curve: 77")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_curve: 78")

    call nurbs%elevate_degree(2)
    call bsp%elevate_degree(2)

    call nurbs%create()
    call bsp%create()

    call nurbs%export_Xc("vtk/test_nurbs_curve_Xc.vtk")
    call bsp%export_Xc("vtk/test_bsp_curve_Xc.vtk")

    call nurbs%export_Xg("vtk/test_nurbs_curve_Xg.vtk")
    call bsp%export_Xg("vtk/test_bsp_curve_Xg.vtk")

    call nurbs%export_Xth("vtk/test_nurbs_curve_Xth.vtk")
    call bsp%export_Xth("vtk/test_bsp_curve_Xth.vtk")

    call nurbs%export_iges('iges/test_nurbs_curve.iges')
    call bsp%export_iges('iges/test_bsp_curve.iges')

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_curve: 79")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_curve: 80")

    call nurbs%remove_knots([0.25_rk, 0.75_rk], [2,1])
    call bsp%remove_knots([0.25_rk, 0.75_rk], [2,1])

    call nurbs%create()
    call bsp%create()

    call nurbs%export_Xc("vtk/test_nurbs_curve_Xc.vtk")
    call bsp%export_Xc("vtk/test_bsp_curve_Xc.vtk")

    call nurbs%export_Xg("vtk/test_nurbs_curve_Xg.vtk")
    call bsp%export_Xg("vtk/test_bsp_curve_Xg.vtk")

    call nurbs%export_Xth("vtk/test_nurbs_curve_Xth.vtk")
    call bsp%export_Xth("vtk/test_bsp_curve_Xth.vtk")

    call nurbs%export_iges('iges/test_nurbs_curve.iges')
    call bsp%export_iges('iges/test_bsp_curve.iges')

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_curve: 81")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_curve: 82")


    call nurbs%set_circle([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk)
    call nurbs%set_half_circle([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk)
    call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk)

    call nurbs%finalize()
    call bsp%finalize()
    deallocate(Xc, Wc, Xg, Xgb)

    call nurbs%set(knot=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xc=[0.0_rk,2.0_rk], Wc=[1.0_rk, 0.9_rk])
    call bsp%set(knot=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xc=[0.0_rk,2.0_rk])

    call nurbs%finalize()
    call bsp%finalize()

    !============================================================================
    ! Test: Least-squares B-spline curve fitting
    !============================================================================
    block
       use forcad, only: rk, nurbs_curve
       use forunittest, only: unit_test

       type(nurbs_curve) :: bsp
       integer :: n, i
       real(rk), parameter :: pi = acos(-1.0_rk)
       real(rk), allocatable :: Xt(:), Xdata(:,:), Xg_eval(:,:)
       real(rk) :: err1, err2, err3, rms
       type(unit_test) :: ut

       n = 42
       allocate(Xt(n), Xdata(n, 3))

       do i = 1, n
          Xt(i) = real(i - 1, rk) / real(n - 1, rk)
          Xdata(i,1) = Xt(i)
          Xdata(i,2) = 0.3_rk * sin(4.0_rk * pi * Xt(i))
          Xdata(i,3) = 0.3_rk * cos(4.0_rk * pi * Xt(i))
       end do

       call bsp%set(&
          degree     = 5,&
          Xth_dir    = [ (real(i - 1, rk)/10.0_rk, i=1,11) ],&
          continuity = [ -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1 ])

       call bsp%lsq_fit_bspline(Xt, Xdata, n)
       call bsp%create(res = n)
       Xg_eval = bsp%get_Xg()

       err1 = norm2(Xg_eval(:,1) - Xdata(:,1)) / norm2(Xdata(:,1))
       err2 = norm2(Xg_eval(:,2) - Xdata(:,2)) / norm2(Xdata(:,2))
       err3 = norm2(Xg_eval(:,3) - Xdata(:,3)) / norm2(Xdata(:,3))
       rms  = sqrt((err1**2 + err2**2 + err3**2) / 3.0_rk)

       call ut%check(res=rms, expected=0.0_rk, tol=1e-6_rk, msg="test_nurbs_curve: 83")
       call bsp%finalize()
       deallocate(Xt, Xdata, Xg_eval)
    end block

end program
