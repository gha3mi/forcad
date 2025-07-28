program test_nurbs_volume

    use forcad, only: rk, nurbs_volume
    use forunittest, only: unit_test

    implicit none
    type(nurbs_volume) :: nurbs, bsp
    real(rk), allocatable :: Xc(:,:), Wc(:)
    real(rk), allocatable :: Xg(:,:), Xgb(:,:)
    integer, allocatable :: elemConn(:,:)
    real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
    real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
    real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
    integer :: i, id
    real(rk), allocatable :: nearest_Xg(:), nearest_Xt(:)
    type(unit_test) :: ut

    allocate(Xc(8, 3))
    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [5.0_rk, 0.0_rk, 0.0_rk]
    Xc(3,:) = [0.0_rk, 5.0_rk, 0.0_rk]
    Xc(4,:) = [5.0_rk, 5.0_rk, 0.0_rk]
    Xc(5,:) = [0.0_rk, 0.0_rk, 5.0_rk]
    Xc(6,:) = [5.0_rk, 0.0_rk, 5.0_rk]
    Xc(7,:) = [0.0_rk, 5.0_rk, 5.0_rk]
    Xc(8,:) = [5.0_rk, 5.0_rk, 5.0_rk]

    allocate(Wc(size(Xc, 1)))
    Wc = 1.0_rk
    Wc(2) = 0.9_rk

    knot1 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
    knot2 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
    knot3 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

    call nurbs%set(knot1, knot2, knot3, Xc, Wc)
    call bsp%set(knot1, knot2, knot3, Xc)

    call nurbs%set(degree=nurbs%get_degree(), nc=nurbs%get_nc(), Xc=nurbs%get_Xc(), Wc=nurbs%get_Wc())
    call bsp%set(degree=nurbs%get_degree(), nc=nurbs%get_nc(), Xc=nurbs%get_Xc())

    call nurbs%create(20, 20, 20)
    call bsp%create(20, 20, 20)

    call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
    call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

    call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
    call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

    call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
    call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

    call nurbs%cmp_volume(volume)
    call bsp%cmp_volume(volumeb)

    call ut%check(res=volume, expected=125.0_rk,  tol=1e-5_rk, msg="test_nurbs_volume: 01")
    call ut%check(res=volumeb, expected=125.0_rk,  tol=1e-5_rk, msg="test_nurbs_volume: 02")

    call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_volume: 03")
    call ut%check(res=nearest_Xt, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_volume: 04")
    call ut%check(res=id, expected=1, msg="test_nurbs_volume: 05")

    call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_volume: 06")
    call ut%check(res=nearest_Xt, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_volume: 07")
    call ut%check(res=id, expected=1, msg="test_nurbs_volume: 08")

    call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_volume: 09")
    call ut%check(res=nearest_Xt, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_volume: 10")

    call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
    call ut%check(res=nearest_Xg, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_volume: 11")
    call ut%check(res=nearest_Xt, expected=[0.0_rk, 0.0_rk, 0.0_rk],  tol=1e-5_rk, msg="test_nurbs_volume: 12")

    Xg = nurbs%get_Xg()
    Xgb = bsp%get_Xg()

    call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
    call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_volume: 13")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_volume: 14")

    call nurbs%set([2,2,2], Xc, Wc)
    call bsp%set([2,2,2], Xc)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_volume: 15")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_volume: 16")


    call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
    call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = nurbs%get_Xt(2))

    call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
    call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

    call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
    call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

    call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
    call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_volume: 17")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_volume: 18")

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_volume: 19")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_volume: 20")

    call ut%check(res=nurbs%get_Xc(1), expected=Xc(1,:),  tol=1e-5_rk, msg="test_nurbs_volume: 21")
    call ut%check(res=bsp%get_Xc(1),   expected=Xc(1,:), tol=1e-5_rk, msg="test_nurbs_volume: 22")

    call ut%check(res=nurbs%get_Xc(1,1), expected=Xc(1,1),  tol=1e-5_rk, msg="test_nurbs_volume: 23")
    call ut%check(res=bsp%get_Xc(1,1),   expected=Xc(1,1), tol=1e-5_rk, msg="test_nurbs_volume: 24")

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_volume: 25")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_volume: 26")

    call ut%check(res=nurbs%get_Xg(1), expected=Xg(1,:),  tol=1e-5_rk, msg="test_nurbs_volume: 27")
    call ut%check(res=bsp%get_Xg(1),   expected=Xgb(1,:), tol=1e-5_rk, msg="test_nurbs_volume: 28")

    call ut%check(res=nurbs%get_Xg(1,1), expected=Xg(1,1),  tol=1e-5_rk, msg="test_nurbs_volume: 29")
    call ut%check(res=bsp%get_Xg(1,1),   expected=Xgb(1,1), tol=1e-5_rk, msg="test_nurbs_volume: 30")

    call ut%check(res=nurbs%get_Wc(), expected=Wc,  tol=1e-5_rk, msg="test_nurbs_volume: 31")

    call ut%check(res=nurbs%get_Wc(1), expected=Wc(1),  tol=1e-5_rk, msg="test_nurbs_volume: 32")

    call ut%check(res=nurbs%get_knot(1), expected=knot1,  tol=1e-5_rk, msg="test_nurbs_volume: 33")
    call ut%check(res=bsp%get_knot(1), expected=knot1,  tol=1e-5_rk, msg="test_nurbs_volume: 34")

    call ut%check(res=nurbs%get_knot(1,1), expected=knot1(1),  tol=1e-5_rk, msg="test_nurbs_volume: 35")
    call ut%check(res=bsp%get_knot(1,1), expected=knot1(1),  tol=1e-5_rk, msg="test_nurbs_volume: 36")

    ! call ut%check(res=nurbs%get_ng(), expected=size(Xg,1), msg="test_nurbs_volume: 37")
    ! call ut%check(res=bsp%get_ng(), expected=size(Xgb,1), msg="test_nurbs_volume: 38")

    call ut%check(res=nurbs%get_degree(1), expected=1, msg="test_nurbs_volume: 39")
    call ut%check(res=bsp%get_degree(1), expected=1, msg="test_nurbs_volume: 40")

    call ut%check(res=nurbs%get_multiplicity(1), expected=[2,2], msg="test_nurbs_volume: 41")
    call ut%check(res=bsp%get_multiplicity(1), expected=[2,2], msg="test_nurbs_volume: 42")

    call ut%check(res=nurbs%get_continuity(1), expected=[-1,-1], msg="test_nurbs_volume: 43")
    call ut%check(res=bsp%get_continuity(1), expected=[-1,-1], msg="test_nurbs_volume: 44")

    ! call ut%check(res=nurbs%get_nc(), expected=size(Xc,1), msg="test_nurbs_volume: 45")
    ! call ut%check(res=bsp%get_nc(), expected=size(Xc,1), msg="test_nurbs_volume: 46")

    call nurbs%cmp_nc()
    call bsp%cmp_nc()

    ! call ut%check(res=nurbs%get_nc(), expected=size(Xc,1), msg="test_nurbs_volume: 47")
    ! call ut%check(res=bsp%get_nc(), expected=size(Xc,1), msg="test_nurbs_volume: 48")

    elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
    call nurbs%set_elem_Xc_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xc_vis(), expected=elemConn, msg="test_nurbs_volume: 49")
    deallocate(elemConn)
    elemConn = nurbs%cmp_elem_Xc_vis()
    call nurbs%set_elem_Xc_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xc_vis(), expected=elemConn, msg="test_nurbs_volume: 50")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
    call bsp%set_elem_Xc_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xc_vis(), expected=elemConn, msg="test_nurbs_volume: 51")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xc_vis()
    call bsp%set_elem_Xc_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xc_vis(), expected=elemConn, msg="test_nurbs_volume: 52")
    deallocate(elemConn)

    elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
    call nurbs%set_elem_Xg_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xg_vis(), expected=elemConn, msg="test_nurbs_volume: 53")
    deallocate(elemConn)
    elemConn = nurbs%cmp_elem_Xg_vis()
    call nurbs%set_elem_Xg_vis(elemConn)
    call ut%check(res=nurbs%get_elem_Xg_vis(), expected=elemConn, msg="test_nurbs_volume: 54")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
    call bsp%set_elem_Xg_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xg_vis(), expected=elemConn, msg="test_nurbs_volume: 55")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem_Xg_vis()
    call bsp%set_elem_Xg_vis(elemConn)
    call ut%check(res=bsp%get_elem_Xg_vis(), expected=elemConn, msg="test_nurbs_volume: 56")
    deallocate(elemConn)

    elemConn = nurbs%cmp_elem()
    call nurbs%set_elem(elemConn)
    call ut%check(res=nurbs%get_elem(), expected=elemConn, msg="test_nurbs_volume: 57")
    deallocate(elemConn)
    elemConn = bsp%cmp_elem()
    call bsp%set_elem(elemConn)
    call ut%check(res=bsp%get_elem(), expected=elemConn, msg="test_nurbs_volume: 58")
    deallocate(elemConn)

    call nurbs%modify_Xc(Xc(1,1), 1,1)
    call bsp%modify_Xc(Xc(1,1), 1,1)

    call nurbs%modify_Wc(Wc(1),1)

    call nurbs%create()
    call bsp%create()

    call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
    call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

    call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
    call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

    call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
    call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_volume: 59")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_volume: 60")

    call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
    call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

    call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
    call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)

    call nurbs%basis(&
        Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        Xt3=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
    call bsp%basis(&
        Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        Xt3=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

    call nurbs%derivative(res1=20, res2=20, res3=20, dTgc=dTgc, Tgc=Tgc)
    call bsp%derivative(res1=20, res2=20, res3=20, dTgc=dTgcb, Tgc=Tgcb)

    call nurbs%derivative(&
        Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        Xt3=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
    call bsp%derivative(&
        Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        Xt3=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

    call nurbs%derivative(Xt=[0.0_rk,0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
    call bsp%derivative(Xt=[0.0_rk,0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

    call nurbs%derivative(Xt=[0.0_rk,0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
    call bsp%derivative(Xt=[0.0_rk,0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

    call nurbs%derivative2(res1=20, res2=20, res3=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
    call bsp%derivative2(res1=20, res2=20, res3=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

    call nurbs%derivative2(&
        Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        Xt3=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
    call bsp%derivative2(&
        Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        Xt3=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
        d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

    call nurbs%derivative2(Xt=[0.0_rk,0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1b, Tgc=Tgc1b)
    call bsp%derivative2(Xt=[0.0_rk,0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

    call nurbs%derivative2(Xt=[0.0_rk,0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1b, Tgc=Tgc1b)
    call bsp%derivative2(Xt=[0.0_rk,0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

    call nurbs%rotate_Xc(45.0_rk, 0.0_rk, 0.0_rk)
    call nurbs%rotate_Xc(-45.0_rk, 0.0_rk, 0.0_rk)

    call bsp%rotate_Xc(45.0_rk, 0.0_rk, 0.0_rk)
    call bsp%rotate_Xc(-45.0_rk, 0.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_volume: 61")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_volume: 62")

    call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
    call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

    call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
    call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_volume: 63")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_volume: 64")

    call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
    call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

    call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
    call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_volume: 65")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_volume: 66")

    call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
    call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

    call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
    call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_volume: 67")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_volume: 68")

    call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
    call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

    call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
    call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_volume: 69")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_volume: 70")

    call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
    call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

    call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
    call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_volume: 71")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_volume: 72")

    call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
    call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

    call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
    call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

    call ut%check(res=nurbs%get_Xc(), expected=Xc,  tol=1e-5_rk, msg="test_nurbs_volume: 73")
    call ut%check(res=bsp%get_Xc(),   expected=Xc, tol=1e-5_rk, msg="test_nurbs_volume: 74")

    call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
    call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

    call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
    call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_volume: 75")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_volume: 76")

    call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
    call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")

    call bsp%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
    call bsp%export_Xg("vtk/test_nurbs_volume_Xg.vtk")

    call nurbs%insert_knots(1, [0.25_rk, 0.75_rk], [1,1])
    call nurbs%insert_knots(2, [0.25_rk, 0.75_rk], [1,1])
    call nurbs%insert_knots(3, [0.25_rk, 0.75_rk], [1,1])

    call bsp%insert_knots(1, [0.25_rk, 0.75_rk], [1,1])
    call bsp%insert_knots(2, [0.25_rk, 0.75_rk], [1,1])
    call bsp%insert_knots(3, [0.25_rk, 0.75_rk], [1,1])

    call nurbs%create()
    call bsp%create()

    call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
    call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

    call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
    call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

    call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
    call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_volume: 77")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_volume: 78")

    call nurbs%elevate_degree(1, 2)
    call nurbs%elevate_degree(2, 2)
    call nurbs%elevate_degree(3, 2)

    call bsp%elevate_degree(1, 2)
    call bsp%elevate_degree(2, 2)
    call bsp%elevate_degree(3, 2)

    call nurbs%create()
    call bsp%create()

    call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
    call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

    call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
    call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

    call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
    call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_volume: 79")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_volume: 80")

    call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
    call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])
    call nurbs%remove_knots(3, [0.25_rk, 0.75_rk], [2,1])

    ! call bsp%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
    ! call bsp%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])
    ! call bsp%remove_knots(3, [0.25_rk, 0.75_rk], [2,1])

    call nurbs%create()
    call bsp%create()

    call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
    call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

    call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
    call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

    call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
    call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

    call ut%check(res=nurbs%get_Xg(), expected=Xg,  tol=1e-5_rk, msg="test_nurbs_volume: 81")
    call ut%check(res=bsp%get_Xg(),   expected=Xgb, tol=1e-5_rk, msg="test_nurbs_volume: 82")


    call nurbs%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2])
    call bsp%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk,1.0_rk,1.0_rk,1.0_rk,0.9_rk])
    call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
    call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
    call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)

    call nurbs%finalize()
    call bsp%finalize()
    deallocate(Xc, Wc, Xg, Xgb)

    !============================================================================
    ! Least-squares B-spline volume fitting
    !============================================================================
    block
       use forcad, only: rk, nurbs_volume
       use forcad_utils, only: ndgrid
       use forunittest, only: unit_test

       type(nurbs_volume) :: bsp
       integer :: n(3), ndata, i
       real(rk), parameter :: pi = acos(-1.0_rk)
       real(rk), allocatable :: Xdata(:,:)
       real(rk), allocatable :: Xt1(:), Xt2(:), Xt3(:), Xt(:,:)
       real(rk), allocatable :: Xg_eval(:,:)
       real(rk) :: err1, err2, err3, rms

       n = [6,6,6]

       allocate(Xt1(n(1)), Xt2(n(2)), Xt3(n(3)))
       do concurrent (i = 1:n(1))
          Xt1(i) = real(i - 1, rk) / real(n(1) - 1, rk)
       end do
       do concurrent (i = 1:n(2))
          Xt2(i) = real(i - 1, rk) / real(n(2) - 1, rk)
       end do
       do concurrent (i = 1:n(3))
          Xt3(i) = real(i - 1, rk) / real(n(3) - 1, rk)
       end do
       call ndgrid(Xt1, Xt2, Xt3, Xt)

       ndata = n(1) * n(2) * n(3)
       allocate(Xdata(ndata, 3))
       do i = 1, ndata
          Xdata(i,1) = Xt(i,1) + 0.1_rk * sin(2.0_rk * pi * Xt(i,2))
          Xdata(i,2) = Xt(i,2) + 0.1_rk * sin(2.0_rk * pi * Xt(i,3))
          Xdata(i,3) = Xt(i,3) + 0.1_rk * sin(2.0_rk * pi * Xt(i,1))
       end do

       call bsp%set(&
          degree      = [2, 2, 2],&
          Xth_dir1    = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk],&
          Xth_dir2    = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk],&
          Xth_dir3    = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk],&
          continuity1 = [ -1   ,   1    ,   1   ,   1    ,  -1   ],&
          continuity2 = [ -1   ,   1    ,   1   ,   1    ,  -1   ],&
          continuity3 = [ -1   ,   1    ,   1   ,   1    ,  -1   ])


       call bsp%lsq_fit_bspline(Xt, Xdata, n)
       call bsp%create(Xt1=Xt1, Xt2=Xt2, Xt3=Xt3)
       Xg_eval = bsp%get_Xg()

       err1 = norm2(Xg_eval(:,1) - Xdata(:,1)) / norm2(Xdata(:,1))
       err2 = norm2(Xg_eval(:,2) - Xdata(:,2)) / norm2(Xdata(:,2))
       err3 = norm2(Xg_eval(:,3) - Xdata(:,3)) / norm2(Xdata(:,3))
       rms  = sqrt((err1**2 + err2**2 + err3**2) / 3.0_rk)

       call ut%check(res=rms, expected=0.0_rk, tol=1e-6_rk, msg="test_nurbs_volume: 83")

       call bsp%finalize()
       deallocate(Xt1, Xt2, Xt3, Xt, Xdata, Xg_eval)
    end block

end program
