program test_nurbs_surface

    use forcad, only: rk, nurbs_surface
    use forunittest, only: unit_test

    implicit none
    type(nurbs_surface) :: nurbs, bsp
    real(rk) :: knot1(4), knot2(4), area, areab
    real(rk), allocatable :: Xc(:,:), Wc(:)
    real(rk), allocatable :: Xg(:,:), Xgb(:,:)
    integer, allocatable :: elemConn(:,:)
    real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
    real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
    integer :: i, id
    real(rk), allocatable :: nearest_Xg(:), nearest_Xt(:)
    type(unit_test) :: ut

    allocate(Xc(4, 3))
    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [5.0_rk, 0.0_rk, 0.0_rk]
    Xc(3,:) = [0.0_rk, 5.0_rk, 0.0_rk]
    Xc(4,:) = [5.0_rk, 5.0_rk, 0.0_rk]

    allocate(Wc(size(Xc, 1)))
    Wc = 1.0_rk
    Wc(2) = 0.9_rk

    knot1 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
    knot2 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

    call nurbs%set(knot1, knot2, Xc, Wc)
    call bsp%set(knot1, knot2, Xc)

    call nurbs%set(degree=nurbs%get_degree(), nc=nurbs%get_nc(), Xc=nurbs%get_Xc(), Wc=nurbs%get_Wc())
    call bsp%set(degree=nurbs%get_degree(), nc=nurbs%get_nc(), Xc=nurbs%get_Xc())

    call nurbs%create(30, 30)
    call bsp%create(30, 30)

    call nurbs%cmp_area(area)
    call bsp%cmp_area(areab)

    call ut%check(res=area, expected=25.0_rk,  tol=1e-5_rk, msg="test_nurbs_surface: 01")
    call ut%check(res=areab, expected=25.0_rk,  tol=1e-5_rk, msg="test_nurbs_surface: 02")

    call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 03")
    call ut%check(res=nearest_Xt, expected=[0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 04")
    call ut%check(res=id, expected=1, msg="test_nurbs_surface: 05")

    call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 06")
    call ut%check(res=nearest_Xt, expected=[0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 07")
    call ut%check(res=id, expected=1, msg="test_nurbs_surface: 08")

    call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 09")
    call ut%check(res=nearest_Xt, expected=[0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 10")

    call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 11")
    call ut%check(res=nearest_Xt, expected=[0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_surface: 12")

    Xg = nurbs%get_Xg()
    Xgb = bsp%get_Xg()

    call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
    call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 13")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 14")

    call nurbs%set([2,2], Xc, Wc)
    call bsp%set([2,2], Xc)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 15")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 16")


    call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
    call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = nurbs%get_Xt(2))

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 17")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 18")

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_surface: 19")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_surface: 20")

    call ut%check(res=nurbs%get_Xc(1), expected=Xc(1,:),  tol=1e-5_rk, msg="test_nurbs_surface: 21")
    call ut%check(res=bsp%get_Xc(1),   expected=Xc(1,:), tol=1e-5_rk, msg="test_nurbs_surface: 22")

    call ut%check(res=nurbs%get_Xc(1,1), expected=Xc(1,1),  tol=1e-5_rk, msg="test_nurbs_surface: 23")
    call ut%check(res=bsp%get_Xc(1,1),   expected=Xc(1,1), tol=1e-5_rk, msg="test_nurbs_surface: 24")

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 25")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 26")

    call ut%check(res=nurbs%get_Xg(1), expected=Xg(1,:),  tol=1e-5_rk, msg="test_nurbs_surface: 27")
    call ut%check(res=bsp%get_Xg(1),   expected=Xgb(1,:), tol=1e-5_rk, msg="test_nurbs_surface: 28")

    call ut%check(res=nurbs%get_Xg(1,1), expected=Xg(1,1),  tol=1e-5_rk, msg="test_nurbs_surface: 29")
    call ut%check(res=bsp%get_Xg(1,1),   expected=Xgb(1,1), tol=1e-5_rk, msg="test_nurbs_surface: 30")

    call ut%check(res=nurbs%get_Wc(), expected=Wc,  tol=1e-5_rk, msg="test_nurbs_surface: 31")

    call ut%check(res=nurbs%get_Wc(1), expected=Wc(1),  tol=1e-5_rk, msg="test_nurbs_surface: 32")

    call ut%check(res=nurbs%get_knot(1), expected=knot1,  tol=1e-5_rk, msg="test_nurbs_surface: 33")
    call ut%check(res=bsp%get_knot(1), expected=knot1,  tol=1e-5_rk, msg="test_nurbs_surface: 34")

    call ut%check(res=nurbs%get_knot(1,1), expected=knot1(1),  tol=1e-5_rk, msg="test_nurbs_surface: 35")
    call ut%check(res=bsp%get_knot(1,1), expected=knot1(1),  tol=1e-5_rk, msg="test_nurbs_surface: 36")

    ! call ut%check(res=nurbs%get_ng(), expected=size(Xg,1), msg="test_nurbs_surface: 37")
    ! call ut%check(res=bsp%get_ng(), expected=size(Xgb,1), msg="test_nurbs_surface: 38")

    call ut%check(res=nurbs%get_degree(1), expected=1, msg="test_nurbs_surface: 39")
    call ut%check(res=bsp%get_degree(1), expected=1, msg="test_nurbs_surface: 40")

    call ut%check(res=nurbs%get_multiplicity(1), expected=[2,2], msg="test_nurbs_surface: 41")
    call ut%check(res=bsp%get_multiplicity(1), expected=[2,2], msg="test_nurbs_surface: 42")

    call ut%check(res=nurbs%get_continuity(1), expected=[-1,-1], msg="test_nurbs_surface: 43")
    call ut%check(res=bsp%get_continuity(1), expected=[-1,-1], msg="test_nurbs_surface: 44")

    ! call ut%check(res=nurbs%get_nc(), expected=size(Xc,1), msg="test_nurbs_surface: 45")
    ! call ut%check(res=bsp%get_nc(), expected=size(Xc,1), msg="test_nurbs_surface: 46")

    call nurbs%cmp_nc()
    call bsp%cmp_nc()

    ! call ut%check(res=nurbs%get_nc(), expected=size(Xc,1), msg="test_nurbs_surface: 47")
    ! call ut%check(res=bsp%get_nc(), expected=size(Xc,1), msg="test_nurbs_surface: 48")

    elemConn = nurbs%cmp_elem_Xc_vis([1,1])
    call nurbs%set_elem_Xc_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xc_vis(), expected=elemConn, msg="test_nurbs_surface: 49")
    deallocate(elemConn)
    elemConn = nurbs%cmp_elem_Xc_vis()
    call nurbs%set_elem_Xc_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xc_vis(), expected=elemConn, msg="test_nurbs_surface: 50")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xc_vis([1,1])
    call bsp%set_elem_Xc_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xc_vis(), expected=elemConn, msg="test_nurbs_surface: 51")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xc_vis()
    call bsp%set_elem_Xc_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xc_vis(), expected=elemConn, msg="test_nurbs_surface: 52")
    deallocate(elemConn)

    elemConn = nurbs%cmp_elem_Xg_vis([1,1])
    call nurbs%set_elem_Xg_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xg_vis(), expected=elemConn, msg="test_nurbs_surface: 53")
    deallocate(elemConn)
    elemConn = nurbs%cmp_elem_Xg_vis()
    call nurbs%set_elem_Xg_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xg_vis(), expected=elemConn, msg="test_nurbs_surface: 54")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xg_vis([1,1])
    call bsp%set_elem_Xg_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xg_vis(), expected=elemConn, msg="test_nurbs_surface: 55")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xg_vis()
    call bsp%set_elem_Xg_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xg_vis(), expected=elemConn, msg="test_nurbs_surface: 56")
    deallocate(elemConn)

    elemConn = nurbs%cmp_elem()
    call nurbs%set_elem(elemConn)
    call ut%check(res=nurbs%get_elem(), expected=elemConn, msg="test_nurbs_surface: 57")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem()
    call bsp%set_elem(elemConn)
    call ut%check(res=bsp%get_elem(), expected=elemConn, msg="test_nurbs_surface: 58")
    deallocate(elemConn)

    call nurbs%modify_Xc(Xc(1,1), 1,1)
    call bsp%modify_Xc(Xc(1,1), 1,1)

    call nurbs%modify_Wc(Wc(1),1)

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 59")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 60")

    call nurbs%basis(res1=30, res2=30, Tgc=Tgc)
    call bsp%basis(res1=30, res2=30, Tgc=Tgc)

    call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
    call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)

    call nurbs%basis(Xt1=[(real(i-1, rk) / real(30-1, rk), i=1, 30)], Xt2=[(real(i-1, rk) / real(30-1, rk), i=1, 30)], Tgc=Tgc)
    call bsp%basis(Xt1=[(real(i-1, rk) / real(30-1, rk), i=1, 30)], Xt2=[(real(i-1, rk) / real(30-1, rk), i=1, 30)], Tgc=Tgc)

    call nurbs%derivative(res1=30, res2=30, dTgc=dTgc, Tgc=Tgc)
    call bsp%derivative(res1=30, res2=30, dTgc=dTgcb, Tgc=Tgcb)

    call nurbs%derivative(&
        Xt1=[(real(i-1, rk) / real(30-1, rk), i=1, 30)], Xt2=[(real(i-1, rk) / real(30-1, rk), i=1, 30)], dTgc=dTgc, Tgc=Tgc)
    call bsp%derivative(&
        Xt1=[(real(i-1, rk) / real(30-1, rk), i=1, 30)], Xt2=[(real(i-1, rk) / real(30-1, rk), i=1, 30)], dTgc=dTgcb, Tgc=Tgcb)

    call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)
    call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

    call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])
    call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

    call nurbs%derivative2(res1=30, res2=30, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
    call bsp%derivative2(res1=30, res2=30, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

    call nurbs%derivative2(&
        Xt1=[(real(i-1, rk) / real(30-1, rk), i=1, 30)], Xt2=[(real(i-1, rk) / real(30-1, rk), i=1, 30)],&
        d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
    call bsp%derivative2(&
        Xt1=[(real(i-1, rk) / real(30-1, rk), i=1, 30)], Xt2=[(real(i-1, rk) / real(30-1, rk), i=1, 30)],&
        d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

    call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1b, Tgc=Tgc1b)
    call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

    call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1b, Tgc=Tgc1b)
    call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

    call nurbs%rotate_Xc(45.0_rk, 0.0_rk, 0.0_rk)
    call nurbs%rotate_Xc(-45.0_rk, 0.0_rk, 0.0_rk)

    call bsp%rotate_Xc(45.0_rk, 0.0_rk, 0.0_rk)
    call bsp%rotate_Xc(-45.0_rk, 0.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_surface: 61")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_surface: 62")

    call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
    call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

    call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
    call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_surface: 63")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_surface: 64")

    call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
    call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

    call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
    call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_surface: 65")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_surface: 66")

    call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
    call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

    call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
    call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 67")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 68")

    call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
    call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

    call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
    call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 69")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 70")

    call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
    call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

    call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
    call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 71")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 72")

    call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
    call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

    call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
    call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_surface: 73")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_surface: 74")

    call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
    call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

    call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
    call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 75")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 76")

    call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
    call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")

    call bsp%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
    call bsp%export_Xg("vtk/test_nurbs_surface_Xg.vtk")

    call nurbs%insert_knots(1, [0.25_rk, 0.75_rk], [1,1])
    call nurbs%insert_knots(2, [0.25_rk, 0.75_rk], [1,1])

    call bsp%insert_knots(1, [0.25_rk, 0.75_rk], [1,1])
    call bsp%insert_knots(2, [0.25_rk, 0.75_rk], [1,1])

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 77")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 78")

    call nurbs%elevate_degree(1, 2)
    call nurbs%elevate_degree(2, 2)

    call bsp%elevate_degree(1, 2)
    call bsp%elevate_degree(2, 2)

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 79")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 80")

    call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
    call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

    call bsp%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
    call bsp%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

    call nurbs%create()
    call bsp%create()

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_surface: 81")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_surface: 82")


    call nurbs%set_tetragon([2.0_rk, 2.0_rk], [2,2])
    call bsp%set_tetragon([2.0_rk, 2.0_rk], [2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk])
    call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
    call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
    call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)

    call nurbs%finalize()
    call bsp%finalize()
    deallocate(Xc, Wc, Xg, Xgb)

end program
