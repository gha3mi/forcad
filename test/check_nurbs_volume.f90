module forcad_test_nurbs_volume

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use forcad_kinds, only: rk
    use forcad_nurbs_volume, only: nurbs_volume
    use forcad_utils, only: ndgrid, rotation, linspace, kron_eye, active_knots, active_span_count
    use forunittest, only: unit_tests

    implicit none

    private
    public :: run_nurbs_volume_tests

contains

    subroutine forcad_nurbs_volume_0001(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb

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
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS volume measure",&
            res      = volume,&
            expected = 125.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Volume measure is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0001


    subroutine forcad_nurbs_volume_0002(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb

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
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline volume measure",&
            res      = volumeb,&
            expected = 125.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Volume measure is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0002


    subroutine forcad_nurbs_volume_0003(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS sampled nearest-point geometry",&
            res      = nearest_Xg,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Sampled nearest-point geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0003


    subroutine forcad_nurbs_volume_0004(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS sampled nearest-point parameters",&
            res      = nearest_Xt,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Sampled nearest-point parameters are incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0004


    subroutine forcad_nurbs_volume_0005(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS sampled nearest-point index",&
            res      = id,&
            expected = 1,&
            msg      = "Sampled nearest-point index is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0005


    subroutine forcad_nurbs_volume_0006(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline sampled nearest-point geometry",&
            res      = nearest_Xg,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Sampled nearest-point geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0006


    subroutine forcad_nurbs_volume_0007(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline sampled nearest-point parameters",&
            res      = nearest_Xt,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Sampled nearest-point parameters are incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0007


    subroutine forcad_nurbs_volume_0008(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline sampled nearest-point index",&
            res      = id,&
            expected = 1,&
            msg      = "Sampled nearest-point index is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0008


    subroutine forcad_nurbs_volume_0009(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS iterative nearest-point geometry",&
            res      = nearest_Xg,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Iterative nearest-point geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0009


    subroutine forcad_nurbs_volume_0010(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS iterative nearest-point parameters",&
            res      = nearest_Xt,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Iterative nearest-point parameters are incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0010


    subroutine forcad_nurbs_volume_0011(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline iterative nearest-point geometry",&
            res      = nearest_Xg,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Iterative nearest-point geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0011


    subroutine forcad_nurbs_volume_0012(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline iterative nearest-point parameters",&
            res      = nearest_Xt,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Iterative nearest-point parameters are incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0012


    subroutine forcad_nurbs_volume_0013(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS explicit-knot setter geometry",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Explicit-knot setter geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0013


    subroutine forcad_nurbs_volume_0014(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline explicit-knot setter geometry",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Explicit-knot setter geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0014


    subroutine forcad_nurbs_volume_0015(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS degree-control setter geometry",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Degree-control setter geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0015


    subroutine forcad_nurbs_volume_0016(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline degree-control setter geometry",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Degree-control setter geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0016


    subroutine forcad_nurbs_volume_0017(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS explicit-grid geometry recreation",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Explicit-grid geometry recreation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0017


    subroutine forcad_nurbs_volume_0018(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline explicit-grid geometry recreation",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Explicit-grid geometry recreation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0018


    subroutine forcad_nurbs_volume_0019(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS control-net getter",&
            res      = nurbs%get_Xc(),&
            expected = Xc,&
            tol      = 1e-5_rk,&
            msg      = "Control-net getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0019


    subroutine forcad_nurbs_volume_0020(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline control-net getter",&
            res      = bsp%get_Xc(),&
            expected = Xc,&
            tol      = 1e-5_rk,&
            msg      = "Control-net getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0020


    subroutine forcad_nurbs_volume_0021(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS control-point getter",&
            res      = nurbs%get_Xc(1),&
            expected = Xc(1,:),&
            tol      = 1e-5_rk,&
            msg      = "Control-point getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0021


    subroutine forcad_nurbs_volume_0022(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline control-point getter",&
            res      = bsp%get_Xc(1),&
            expected = Xc(1,:),&
            tol      = 1e-5_rk,&
            msg      = "Control-point getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0022


    subroutine forcad_nurbs_volume_0023(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS control-coordinate getter",&
            res      = nurbs%get_Xc(1,1),&
            expected = Xc(1,1),&
            tol      = 1e-5_rk,&
            msg      = "Control-coordinate getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0023


    subroutine forcad_nurbs_volume_0024(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline control-coordinate getter",&
            res      = bsp%get_Xc(1,1),&
            expected = Xc(1,1),&
            tol      = 1e-5_rk,&
            msg      = "Control-coordinate getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0024


    subroutine forcad_nurbs_volume_0025(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS geometry getter",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Geometry getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0025


    subroutine forcad_nurbs_volume_0026(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline geometry getter",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Geometry getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0026


    subroutine forcad_nurbs_volume_0027(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS geometry-point getter",&
            res      = nurbs%get_Xg(1),&
            expected = Xg(1,:),&
            tol      = 1e-5_rk,&
            msg      = "Geometry-point getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0027


    subroutine forcad_nurbs_volume_0028(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline geometry-point getter",&
            res      = bsp%get_Xg(1),&
            expected = Xgb(1,:),&
            tol      = 1e-5_rk,&
            msg      = "Geometry-point getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0028


    subroutine forcad_nurbs_volume_0029(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS geometry-coordinate getter",&
            res      = nurbs%get_Xg(1,1),&
            expected = Xg(1,1),&
            tol      = 1e-5_rk,&
            msg      = "Geometry-coordinate getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0029


    subroutine forcad_nurbs_volume_0030(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline geometry-coordinate getter",&
            res      = bsp%get_Xg(1,1),&
            expected = Xgb(1,1),&
            tol      = 1e-5_rk,&
            msg      = "Geometry-coordinate getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0030


    subroutine forcad_nurbs_volume_0031(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS weight-vector getter",&
            res      = nurbs%get_Wc(),&
            expected = Wc,&
            tol      = 1e-5_rk,&
            msg      = "Weight-vector getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0031


    subroutine forcad_nurbs_volume_0032(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS weight getter",&
            res      = nurbs%get_Wc(1),&
            expected = Wc(1),&
            tol      = 1e-5_rk,&
            msg      = "Weight getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0032


    subroutine forcad_nurbs_volume_0033(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS directional knot-vector getter",&
            res      = nurbs%get_knot(1),&
            expected = knot1,&
            tol      = 1e-5_rk,&
            msg      = "Directional knot-vector getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0033


    subroutine forcad_nurbs_volume_0034(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline directional knot-vector getter",&
            res      = bsp%get_knot(1),&
            expected = knot1,&
            tol      = 1e-5_rk,&
            msg      = "Directional knot-vector getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0034


    subroutine forcad_nurbs_volume_0035(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS knot getter",&
            res      = nurbs%get_knot(1,1),&
            expected = knot1(1),&
            tol      = 1e-5_rk,&
            msg      = "The NURBS knot getter result does not match the expected value.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0035


    subroutine forcad_nurbs_volume_0036(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline knot getter",&
            res      = bsp%get_knot(1,1),&
            expected = knot1(1),&
            tol      = 1e-5_rk,&
            msg      = "Knot getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0036


    subroutine forcad_nurbs_volume_0037(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS directional degree getter",&
            res      = nurbs%get_degree(1),&
            expected = 1,&
            msg      = "Directional degree getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0037


    subroutine forcad_nurbs_volume_0038(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline directional degree getter",&
            res      = bsp%get_degree(1),&
            expected = 1,&
            msg      = "Directional degree getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0038


    subroutine forcad_nurbs_volume_0039(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS knot multiplicity getter",&
            res      = nurbs%get_multiplicity(1),&
            expected = [2,2],&
            msg      = "Knot multiplicity getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0039


    subroutine forcad_nurbs_volume_0040(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline knot multiplicity getter",&
            res      = bsp%get_multiplicity(1),&
            expected = [2,2],&
            msg      = "Knot multiplicity getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0040


    subroutine forcad_nurbs_volume_0041(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS knot continuity getter",&
            res      = nurbs%get_continuity(1),&
            expected = [-1,-1],&
            msg      = "Knot continuity getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0041


    subroutine forcad_nurbs_volume_0042(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline knot continuity getter",&
            res      = bsp%get_continuity(1),&
            expected = [-1,-1],&
            msg      = "Knot continuity getter is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0042


    subroutine forcad_nurbs_volume_0043(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS control-net shape",&
            res      = nurbs%get_nc(),&
            expected = [2, 2, 2],&
            msg      = "Control-net shape is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0043


    subroutine forcad_nurbs_volume_0044(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline control-net shape",&
            res      = bsp%get_nc(),&
            expected = [2, 2, 2],&
            msg      = "Control-net shape is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0044


    subroutine forcad_nurbs_volume_0045(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS recomputed control-net shape",&
            res      = nurbs%get_nc(),&
            expected = [2, 2, 2],&
            msg      = "Recomputed control-net shape is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0045


    subroutine forcad_nurbs_volume_0046(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline recomputed control-net shape",&
            res      = bsp%get_nc(),&
            expected = [2, 2, 2],&
            msg      = "Recomputed control-net shape is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0046


    subroutine forcad_nurbs_volume_0047(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS custom control-net visualization connectivity",&
            res      = nurbs%get_elem_Xc_vis(),&
            expected = elemConn,&
            msg      = "Custom control-net visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0047


    subroutine forcad_nurbs_volume_0048(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS default control-net visualization connectivity",&
            res      = nurbs%get_elem_Xc_vis(),&
            expected = elemConn,&
            msg      = "Default control-net visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0048


    subroutine forcad_nurbs_volume_0049(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline custom control-net visualization connectivity",&
            res      = bsp%get_elem_Xc_vis(),&
            expected = elemConn,&
            msg      = "Custom control-net visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0049


    subroutine forcad_nurbs_volume_0050(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline default control-net visualization connectivity",&
            res      = bsp%get_elem_Xc_vis(),&
            expected = elemConn,&
            msg      = "Default control-net visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0050


    subroutine forcad_nurbs_volume_0051(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS custom geometry visualization connectivity",&
            res      = nurbs%get_elem_Xg_vis(),&
            expected = elemConn,&
            msg      = "Custom geometry visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0051


    subroutine forcad_nurbs_volume_0052(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS default geometry visualization connectivity",&
            res      = nurbs%get_elem_Xg_vis(),&
            expected = elemConn,&
            msg      = "Default geometry visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0052


    subroutine forcad_nurbs_volume_0053(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline custom geometry visualization connectivity",&
            res      = bsp%get_elem_Xg_vis(),&
            expected = elemConn,&
            msg      = "Custom geometry visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0053


    subroutine forcad_nurbs_volume_0054(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline default geometry visualization connectivity",&
            res      = bsp%get_elem_Xg_vis(),&
            expected = elemConn,&
            msg      = "Default geometry visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0054


    subroutine forcad_nurbs_volume_0055(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS analysis element connectivity",&
            res      = nurbs%get_elem(),&
            expected = elemConn,&
            msg      = "Analysis element connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0055


    subroutine forcad_nurbs_volume_0056(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline analysis element connectivity",&
            res      = bsp%get_elem(),&
            expected = elemConn,&
            msg      = "Analysis element connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0056


    subroutine forcad_nurbs_volume_0057(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS modified-control geometry recreation",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Modified-control geometry recreation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0057


    subroutine forcad_nurbs_volume_0058(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline modified-control geometry recreation",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Modified-control geometry recreation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0058


    subroutine forcad_nurbs_volume_0059(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:)
        real(rk), allocatable :: Tgc1(:)
        real(rk), allocatable :: Tgc1b(:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS boundary active basis values",&
            res      = Tgc1e,&
            expected = Tgc1(1:3),&
            tol      = 1e-5_rk,&
            msg      = "Boundary active basis values are incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0059


    subroutine forcad_nurbs_volume_0060(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:)
        real(rk), allocatable :: Tgc1(:)
        real(rk), allocatable :: Tgc1b(:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline boundary active basis values",&
            res      = Tgc1be,&
            expected = Tgc1b(1:3),&
            tol      = 1e-5_rk,&
            msg      = "Boundary active basis values are incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0060


    subroutine forcad_nurbs_volume_0061(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:)
        real(rk), allocatable :: Tgc1(:)
        real(rk), allocatable :: Tgc1b(:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS boundary element basis selection",&
            res      = maxval(abs(Tgc1e - [Tgc1(1), Tgc1(8)])),&
            expected = 0.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Boundary element basis selection is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0061


    subroutine forcad_nurbs_volume_0062(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:)
        real(rk), allocatable :: Tgc1(:)
        real(rk), allocatable :: Tgc1b(:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline boundary element basis selection",&
            res      = maxval(abs(Tgc1be - [Tgc1b(1), Tgc1b(8)])),&
            expected = 0.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Boundary element basis selection is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0062


    subroutine forcad_nurbs_volume_0063(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS control-net x-rotation round trip",&
            res      = nurbs%get_Xc(),&
            expected = Xc,&
            tol      = 1e-5_rk,&
            msg      = "Control-net x-rotation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0063


    subroutine forcad_nurbs_volume_0064(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline control-net x-rotation round trip",&
            res      = bsp%get_Xc(),&
            expected = Xc,&
            tol      = 1e-5_rk,&
            msg      = "Control-net x-rotation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0064


    subroutine forcad_nurbs_volume_0065(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS control-net y-rotation round trip",&
            res      = nurbs%get_Xc(),&
            expected = Xc,&
            tol      = 1e-5_rk,&
            msg      = "Control-net y-rotation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0065


    subroutine forcad_nurbs_volume_0066(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline control-net y-rotation round trip",&
            res      = bsp%get_Xc(),&
            expected = Xc,&
            tol      = 1e-5_rk,&
            msg      = "Control-net y-rotation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0066


    subroutine forcad_nurbs_volume_0067(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS control-net z-rotation round trip",&
            res      = nurbs%get_Xc(),&
            expected = Xc,&
            tol      = 1e-5_rk,&
            msg      = "Control-net z-rotation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0067


    subroutine forcad_nurbs_volume_0068(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline control-net z-rotation round trip",&
            res      = bsp%get_Xc(),&
            expected = Xc,&
            tol      = 1e-5_rk,&
            msg      = "Control-net z-rotation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0068


    subroutine forcad_nurbs_volume_0069(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS geometry x-rotation round trip",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Geometry x-rotation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0069


    subroutine forcad_nurbs_volume_0070(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline geometry x-rotation round trip",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Geometry x-rotation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0070


    subroutine forcad_nurbs_volume_0071(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS geometry y-rotation round trip",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Geometry y-rotation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0071


    subroutine forcad_nurbs_volume_0072(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline geometry y-rotation round trip",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Geometry y-rotation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0072


    subroutine forcad_nurbs_volume_0073(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS geometry z-rotation round trip",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Geometry z-rotation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0073


    subroutine forcad_nurbs_volume_0074(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline geometry z-rotation round trip",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Geometry z-rotation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0074


    subroutine forcad_nurbs_volume_0075(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS control-net translation round trip",&
            res      = nurbs%get_Xc(),&
            expected = Xc,&
            tol      = 1e-5_rk,&
            msg      = "Control-net translation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0075


    subroutine forcad_nurbs_volume_0076(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline control-net translation round trip",&
            res      = bsp%get_Xc(),&
            expected = Xc,&
            tol      = 1e-5_rk,&
            msg      = "Control-net translation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0076


    subroutine forcad_nurbs_volume_0077(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS geometry translation round trip",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Geometry translation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0077


    subroutine forcad_nurbs_volume_0078(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline geometry translation round trip",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Geometry translation round trip is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0078


    subroutine forcad_nurbs_volume_0079(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

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
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS knot insertion geometry preservation",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Knot insertion geometry preservation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0079


    subroutine forcad_nurbs_volume_0080(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

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
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline knot insertion geometry preservation",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Knot insertion geometry preservation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0080


    subroutine forcad_nurbs_volume_0081(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

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
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS degree elevation geometry preservation",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Degree elevation geometry preservation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0081


    subroutine forcad_nurbs_volume_0082(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

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
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline degree elevation geometry preservation",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Degree elevation geometry preservation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0082


    subroutine forcad_nurbs_volume_0083(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

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

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(3, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS knot removal geometry preservation",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Knot removal geometry preservation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0083


    subroutine forcad_nurbs_volume_0084(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

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

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(3, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline knot removal geometry preservation",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Knot removal geometry preservation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0084


    subroutine forcad_nurbs_volume_0085(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:), faceConn(:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

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

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(3, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2])
        call bsp%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk,1.0_rk,1.0_rk,1.0_rk,0.9_rk])
        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
        call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
        call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)

        call nurbs%set_hexahedron([2.0_rk, 3.0_rk, 4.0_rk], [2,2,2])
        call nurbs%create(2, 2, 2)
        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        faceConn = nurbs%cmp_elemFace(1, 6)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "volume analysis face connectivity",&
            res      = faceConn,&
            expected = [2,4,6,8],&
            msg      = "Volume analysis face connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0085


    subroutine forcad_nurbs_volume_0086(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:), faceConn(:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

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

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(3, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2])
        call bsp%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk,1.0_rk,1.0_rk,1.0_rk,0.9_rk])
        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
        call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
        call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)

        call nurbs%set_hexahedron([2.0_rk, 3.0_rk, 4.0_rk], [2,2,2])
        call nurbs%create(2, 2, 2)
        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        faceConn = nurbs%cmp_elemFace(1, 6)
        deallocate(faceConn, elemConn)

        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        faceConn = nurbs%cmp_elemFace_Xc_vis(1, 6)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "volume control-net face connectivity",&
            res      = faceConn,&
            expected = [2,4,6,8],&
            msg      = "Volume control-net face connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0086


    subroutine forcad_nurbs_volume_0087(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:), faceConn(:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

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

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(3, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2])
        call bsp%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk,1.0_rk,1.0_rk,1.0_rk,0.9_rk])
        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
        call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
        call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)

        call nurbs%set_hexahedron([2.0_rk, 3.0_rk, 4.0_rk], [2,2,2])
        call nurbs%create(2, 2, 2)
        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        faceConn = nurbs%cmp_elemFace(1, 6)
        deallocate(faceConn, elemConn)

        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        faceConn = nurbs%cmp_elemFace_Xc_vis(1, 6)
        deallocate(faceConn, elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        faceConn = nurbs%cmp_elemFace_Xg_vis(1, 6)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "volume geometry face connectivity",&
            res      = faceConn,&
            expected = [2,4,6,8],&
            msg      = "Volume geometry face connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0087


    subroutine forcad_nurbs_volume_0088(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:), faceConn(:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

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

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(3, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2])
        call bsp%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk,1.0_rk,1.0_rk,1.0_rk,0.9_rk])
        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
        call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
        call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)

        call nurbs%set_hexahedron([2.0_rk, 3.0_rk, 4.0_rk], [2,2,2])
        call nurbs%create(2, 2, 2)
        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        faceConn = nurbs%cmp_elemFace(1, 6)
        deallocate(faceConn, elemConn)

        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        faceConn = nurbs%cmp_elemFace_Xc_vis(1, 6)
        deallocate(faceConn, elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        faceConn = nurbs%cmp_elemFace_Xg_vis(1, 6)
        deallocate(faceConn, elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "volume face degree",&
            res      = nurbs%cmp_degreeFace(5),&
            expected = [0,1,1],&
            msg      = "The Volume face degree result does not match the expected value.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0088


    subroutine forcad_nurbs_volume_0089(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:), faceConn(:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

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

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(3, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2])
        call bsp%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk,1.0_rk,1.0_rk,1.0_rk,0.9_rk])
        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
        call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
        call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)

        call nurbs%set_hexahedron([2.0_rk, 3.0_rk, 4.0_rk], [2,2,2])
        call nurbs%create(2, 2, 2)
        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        faceConn = nurbs%cmp_elemFace(1, 6)
        deallocate(faceConn, elemConn)

        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        faceConn = nurbs%cmp_elemFace_Xc_vis(1, 6)
        deallocate(faceConn, elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        faceConn = nurbs%cmp_elemFace_Xg_vis(1, 6)
        deallocate(faceConn, elemConn)

        call nurbs%put_to_nurbs(nurbs%get_Xc(), reshape([1,2,3,4,5,6,7,8], [1,8]))
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "trilinear volume geometry",&
            res      = maxval(abs(nurbs%get_Xg() - nurbs%get_Xc())),&
            expected = 0.0_rk,&
            tol      = 1e-10_rk,&
            msg      = "Trilinear volume geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0089


    subroutine forcad_nurbs_volume_0090(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_volume) :: nurbs, bsp
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:), faceConn(:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: knot1(4), knot2(4), knot3(4), volume, volumeb
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

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

        call nurbs%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, -0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, -0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1,1], [-1, -1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2,2], Xc, Wc)
        call bsp%set([2,2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2), Xt3 = nurbs%get_Xt(3))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2), Xt3 = bsp%get_Xt(3))

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis()
        call bsp%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem()
        call bsp%set_elem(elemConn)
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

        call nurbs%basis(res1=20, res2=20, res3=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, res3=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1e, elem=[1,8])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk, 0.75_rk], Tgc=Tgc1be, elem=[1,8])

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

        call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
        call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
        call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)

        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
        call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)

        call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])

        call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

        call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
        call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])

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

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(3, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_volume_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_volume_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_volume_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_volume_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_volume_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_volume_Xth.vtk")

        call nurbs%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2])
        call bsp%set_hexahedron([2.0_rk, 2.0_rk, 2.0_rk], [2,2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk,1.0_rk,1.0_rk,1.0_rk,0.9_rk])
        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
        call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)
        call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 2.0_rk)

        call nurbs%set_hexahedron([2.0_rk, 3.0_rk, 4.0_rk], [2,2,2])
        call nurbs%create(2, 2, 2)
        elemConn = nurbs%cmp_elem()
        call nurbs%set_elem(elemConn)
        faceConn = nurbs%cmp_elemFace(1, 6)
        deallocate(faceConn, elemConn)

        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        faceConn = nurbs%cmp_elemFace_Xc_vis(1, 6)
        deallocate(faceConn, elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        faceConn = nurbs%cmp_elemFace_Xg_vis(1, 6)
        deallocate(faceConn, elemConn)

        call nurbs%put_to_nurbs(nurbs%get_Xc(), reshape([1,2,3,4,5,6,7,8], [1,8]))
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "volume geometry visualization connectivity shape",&
            res      = shape(nurbs%get_elem_Xg_vis()),&
            expected = [1,8],&
            msg      = "Volume geometry visualization connectivity shape is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0090


    subroutine forcad_nurbs_volume_0091(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

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
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline volume least-squares fit",&
            res      = rms,&
            expected = 0.0_rk,&
            tol      = 1e-6_rk,&
            msg      = "Volume least-squares fit is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0091


    subroutine forcad_nurbs_volume_0092(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nonopen
        real(rk), allocatable :: Xc_nonopen(:,:)
        real(rk), allocatable :: Xt_nonopen(:)
        real(rk) :: knot_nonopen(7)
        integer :: ii, jj, kk, idx

        allocate(Xc_nonopen(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ((kk - 1)*4 + jj - 1)*4 + ii
                    Xc_nonopen(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), real(kk - 1, rk)]
                end do
            end do
        end do
        knot_nonopen = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call nonopen%set(knot_nonopen, knot_nonopen, knot_nonopen, Xc_nonopen, degree=[2, 2, 2])
        call nonopen%create(3, 3, 3)
        Xt_nonopen = nonopen%get_Xt(1)
        call nonopen%err%print()

        call ut%test(ti)%check(&
            name     = "unclamped volume degree",&
            res      = nonopen%get_degree(),&
            expected = [2, 2, 2],&
            msg      = "Unclamped volume degree is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0092


    subroutine forcad_nurbs_volume_0093(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nonopen
        real(rk), allocatable :: Xc_nonopen(:,:)
        real(rk), allocatable :: Xt_nonopen(:)
        real(rk) :: knot_nonopen(7)
        integer :: ii, jj, kk, idx

        allocate(Xc_nonopen(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ((kk - 1)*4 + jj - 1)*4 + ii
                    Xc_nonopen(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), real(kk - 1, rk)]
                end do
            end do
        end do
        knot_nonopen = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call nonopen%set(knot_nonopen, knot_nonopen, knot_nonopen, Xc_nonopen, degree=[2, 2, 2])
        call nonopen%create(3, 3, 3)
        Xt_nonopen = nonopen%get_Xt(1)
        call nonopen%err%print()

        call ut%test(ti)%check(&
            name     = "unclamped volume control-net shape",&
            res      = nonopen%get_nc(),&
            expected = [4, 4, 4],&
            msg      = "Unclamped volume control-net shape is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0093


    subroutine forcad_nurbs_volume_0094(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nonopen
        real(rk), allocatable :: Xc_nonopen(:,:)
        real(rk), allocatable :: Xt_nonopen(:)
        real(rk) :: knot_nonopen(7)
        integer :: ii, jj, kk, idx

        allocate(Xc_nonopen(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ((kk - 1)*4 + jj - 1)*4 + ii
                    Xc_nonopen(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), real(kk - 1, rk)]
                end do
            end do
        end do
        knot_nonopen = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call nonopen%set(knot_nonopen, knot_nonopen, knot_nonopen, Xc_nonopen, degree=[2, 2, 2])
        call nonopen%create(3, 3, 3)
        Xt_nonopen = nonopen%get_Xt(1)
        call nonopen%err%print()

        call ut%test(ti)%check(&
            name     = "unclamped volume active domain",&
            res      = [Xt_nonopen(1), Xt_nonopen(size(Xt_nonopen))],&
            expected = [2.0_rk, 4.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Unclamped volume active domain is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0094


    subroutine forcad_nurbs_volume_0095(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nonopen
        real(rk), allocatable :: Xc_nonopen(:,:)
        real(rk), allocatable :: Xt_nonopen(:)
        real(rk), allocatable :: Tgc_nonopen(:)
        real(rk) :: knot_nonopen(7), pt(3)
        integer :: ii, jj, kk, idx

        allocate(Xc_nonopen(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ((kk - 1)*4 + jj - 1)*4 + ii
                    Xc_nonopen(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), real(kk - 1, rk)]
                end do
            end do
        end do
        knot_nonopen = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call nonopen%set(knot_nonopen, knot_nonopen, knot_nonopen, Xc_nonopen, degree=[2, 2, 2])
        call nonopen%create(3, 3, 3)
        Xt_nonopen = nonopen%get_Xt(1)

        pt = [2.5_rk, 3.0_rk, 3.5_rk]
        call nonopen%basis(pt, Tgc_nonopen)
        call nonopen%err%print()

        call ut%test(ti)%check(&
            name     = "unclamped volume partition of unity",&
            res      = sum(Tgc_nonopen),&
            expected = 1.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Unclamped volume partition of unity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0095


    subroutine forcad_nurbs_volume_0096(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nonopen
        real(rk), allocatable :: Xc_nonopen(:,:)
        real(rk), allocatable :: Xt_nonopen(:)
        real(rk), allocatable :: Tgc_nonopen(:)
        real(rk) :: knot_nonopen(7), pt(3)
        integer :: ii, jj, kk, idx

        allocate(Xc_nonopen(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ((kk - 1)*4 + jj - 1)*4 + ii
                    Xc_nonopen(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), real(kk - 1, rk)]
                end do
            end do
        end do
        knot_nonopen = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call nonopen%set(knot_nonopen, knot_nonopen, knot_nonopen, Xc_nonopen, degree=[2, 2, 2])
        call nonopen%create(3, 3, 3)
        Xt_nonopen = nonopen%get_Xt(1)

        pt = [2.5_rk, 3.0_rk, 3.5_rk]
        call nonopen%basis(pt, Tgc_nonopen)

        pt = [2.5_rk, 4.5_rk, 3.0_rk]
        call nonopen%basis(pt, Tgc_nonopen)
        call nonopen%err%print()

        call ut%test(ti)%check(&
            name     = "inactive unclamped volume span basis",&
            res      = sum(Tgc_nonopen),&
            expected = 0.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Inactive unclamped volume span basis is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0096


    subroutine forcad_nurbs_volume_0097(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: nonopen
        real(rk), allocatable :: Xc_nonopen(:,:), Xt_nonopen(:), Tgc_nonopen(:), dTgc_nonopen(:,:)
        real(rk) :: knot_nonopen(7), pt(3)
        integer :: ii, jj, kk, idx

        allocate(Xc_nonopen(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ((kk - 1)*4 + jj - 1)*4 + ii
                    Xc_nonopen(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), real(kk - 1, rk)]
                end do
            end do
        end do
        knot_nonopen = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call nonopen%set(knot_nonopen, knot_nonopen, knot_nonopen, Xc_nonopen, degree=[2, 2, 2])
        call nonopen%create(3, 3, 3)
        Xt_nonopen = nonopen%get_Xt(1)

        pt = [2.5_rk, 3.0_rk, 3.5_rk]
        call nonopen%basis(pt, Tgc_nonopen)

        pt = [2.5_rk, 4.5_rk, 3.0_rk]
        call nonopen%basis(pt, Tgc_nonopen)

        pt = [2.5_rk, 3.0_rk, 3.5_rk]
        call nonopen%derivative(pt, dTgc_nonopen, Tgc_nonopen, elem=[0, 65])
        call nonopen%err%print()

        call ut%test(ti)%check(&
            name     = "out-of-domain unclamped volume derivatives",&
            res      = max(maxval(abs(Tgc_nonopen)), maxval(abs(dTgc_nonopen))),&
            expected = 0.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Out-of-domain unclamped volume derivatives are incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0097


    subroutine forcad_nurbs_volume_0098(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: rational
        real(rk), allocatable :: Xc_out(:,:), Wc_out(:), Xg_out(:), Tgc_out(:)
        real(rk) :: kout(4)
        integer :: ii, jj, kk, idx

        allocate(Xc_out(8, 3), Wc_out(8))
        do kk = 1, 2
            do jj = 1, 2
                do ii = 1, 2
                    idx = ((kk - 1)*2 + jj - 1)*2 + ii
                    Xc_out(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), real(kk - 1, rk)]
                end do
            end do
        end do
        Wc_out = 1.0_rk
        Wc_out(2) = 0.9_rk
        kout = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

        call rational%set(kout, kout, kout, Xc_out, Wc_out)
        Xg_out = rational%cmp_Xg([-0.25_rk, 0.0_rk, 0.0_rk])
        call rational%basis([-0.25_rk, 0.0_rk, 0.0_rk], Tgc_out)
        call rational%err%print()

        call ut%test(ti)%check(&
            name     = "out-of-domain volume geometry",&
            res      = all(abs(Xg_out) <= 1e-5_rk),&
            expected = .true.,&
            msg      = "Out-of-domain volume geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0098


    subroutine forcad_nurbs_volume_0099(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: rational
        real(rk), allocatable :: Xc_out(:,:), Wc_out(:), Xg_out(:), Tgc_out(:)
        real(rk) :: kout(4)
        integer :: ii, jj, kk, idx

        allocate(Xc_out(8, 3), Wc_out(8))
        do kk = 1, 2
            do jj = 1, 2
                do ii = 1, 2
                    idx = ((kk - 1)*2 + jj - 1)*2 + ii
                    Xc_out(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), real(kk - 1, rk)]
                end do
            end do
        end do
        Wc_out = 1.0_rk
        Wc_out(2) = 0.9_rk
        kout = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

        call rational%set(kout, kout, kout, Xc_out, Wc_out)
        Xg_out = rational%cmp_Xg([-0.25_rk, 0.0_rk, 0.0_rk])
        call rational%basis([-0.25_rk, 0.0_rk, 0.0_rk], Tgc_out)
        call rational%err%print()

        call ut%test(ti)%check(&
            name     = "out-of-domain volume basis",&
            res      = all(abs(Tgc_out) <= 1e-5_rk),&
            expected = .true.,&
            msg      = "Out-of-domain volume basis is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0099


    subroutine forcad_nurbs_volume_0100(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: repeated
        real(rk), allocatable :: Xc_repeated(:,:)
        real(rk), allocatable :: Xt_repeated(:)
        real(rk) :: knot_repeated(9)
        integer :: ii, jj, kk, idx

        allocate(Xc_repeated(216, 3), source=0.0_rk)
        do kk = 1, 6
            do jj = 1, 6
                do ii = 1, 6
                    idx = ii + 6*((jj - 1) + 6*(kk - 1))
                    Xc_repeated(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), real(kk - 1, rk)]
                end do
            end do
        end do
        knot_repeated = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 3.0_rk, 3.0_rk, 4.0_rk, 6.0_rk, 9.0_rk]

        call repeated%set(knot_repeated, knot_repeated, knot_repeated, Xc_repeated, degree=[2, 2, 2])
        call repeated%create(4, 4, 4)
        Xt_repeated = repeated%get_Xt(1)
        call repeated%err%print()

        call ut%test(ti)%check(&
            name     = "repeated-knot volume active domain",&
            res      = [Xt_repeated(1), Xt_repeated(size(Xt_repeated))],&
            expected = [2.0_rk, 4.0_rk],&
            tol      = 1.0e-10_rk,&
            msg      = "Repeated-knot volume active domain is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0100


    subroutine forcad_nurbs_volume_0101(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: repeated
        real(rk), allocatable :: Xc_repeated(:,:), Xt_repeated(:), Tgc_repeated(:)
        real(rk) :: knot_repeated(9), pt(3)
        integer :: ii, jj, kk, idx

        allocate(Xc_repeated(216, 3), source=0.0_rk)
        do kk = 1, 6
            do jj = 1, 6
                do ii = 1, 6
                    idx = ii + 6*((jj - 1) + 6*(kk - 1))
                    Xc_repeated(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), real(kk - 1, rk)]
                end do
            end do
        end do
        knot_repeated = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 3.0_rk, 3.0_rk, 4.0_rk, 6.0_rk, 9.0_rk]

        call repeated%set(knot_repeated, knot_repeated, knot_repeated, Xc_repeated, degree=[2, 2, 2])
        call repeated%create(4, 4, 4)
        Xt_repeated = repeated%get_Xt(1)

        pt = [2.5_rk, 3.25_rk, 3.75_rk]
        call repeated%basis(pt, Tgc_repeated)
        call repeated%err%print()

        call ut%test(ti)%check(&
            name     = "repeated-knot volume partition of unity",&
            res      = sum(Tgc_repeated),&
            expected = 1.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "Repeated-knot volume partition of unity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0101


    subroutine forcad_nurbs_volume_0102(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume_ref, volume_scaled
        real(rk) :: knot_scale(4), pt(3)
        real(rk), allocatable :: Xc_scale(:,:), W_scale(:)
        real(rk), allocatable :: Xg_ref(:), Xg_scaled(:), T_ref(:), T_scaled(:), dT_ref(:,:), dT_scaled(:,:)
        integer :: ii, jj, kk, idx

        allocate(Xc_scale(8, 3), W_scale(8))
        do kk = 1, 2
            do jj = 1, 2
                do ii = 1, 2
                    idx = ii + 2*((jj - 1) + 2*(kk - 1))
                    Xc_scale(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), real(kk - 1, rk)]
                end do
            end do
        end do
        W_scale = [1.0_rk, 0.7_rk, 1.3_rk, 0.9_rk, 1.1_rk, 0.8_rk, 1.2_rk, 0.6_rk]
        knot_scale = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        pt = [0.37_rk, 0.61_rk, 0.29_rk]

        call volume_ref%set(knot_scale, knot_scale, knot_scale, Xc_scale, W_scale)
        call volume_scaled%set(knot_scale, knot_scale, knot_scale, Xc_scale, 1.0e-20_rk*W_scale)
        Xg_ref = volume_ref%cmp_Xg(pt)
        Xg_scaled = volume_scaled%cmp_Xg(pt)
        call volume_ref%derivative(pt, dT_ref, T_ref)
        call volume_scaled%derivative(pt, dT_scaled, T_scaled)
        call volume_ref%err%print()
        call volume_scaled%err%print()

        call ut%test(ti)%check(&
            name     = "affine knot scaling volume invariance",&
            res      = maxval(abs(Xg_ref - Xg_scaled)) <= 1.0e-10_rk .and. maxval(abs(T_ref - T_scaled)) <= 1.0e-10_rk .and. &
                maxval(abs(dT_ref - dT_scaled)) <= 1.0e-10_rk,&
            expected = .true.,&
            msg      = "Affine knot scaling volume invariance is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0102


    subroutine forcad_nurbs_volume_0103(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume_fd
        real(rk), allocatable :: Wc_fd(:), T0(:), Tm(:), Tp(:), dT0(:,:), dTm(:,:), dTp(:,:)
        real(rk), allocatable :: d2T0(:,:), d2Tm(:,:), d2Tp(:,:), cfd(:,:), cfd2(:,:)
        real(rk) :: h, pt(3), pm(3), pp(3)
        integer :: dir

        allocate(Wc_fd(64), source=1.0_rk)
        Wc_fd(10) = 0.5_rk
        h = 1.0e-6_rk
        pt = [0.5_rk, 0.3_rk, 0.7_rk]

        call volume_fd%set_hexahedron(L=[2.0_rk, 4.0_rk, 8.0_rk], nc=[4,4,4], Wc=Wc_fd)
        call volume_fd%derivative2(Xt=pt, d2Tgc=d2T0, dTgc=dT0, Tgc=T0)
        allocate(cfd(size(T0),3), cfd2(size(d2T0,1),size(d2T0,2)))
        do dir = 1, 3
            pm = pt
            pp = pt
            pm(dir) = pm(dir) - h
            pp(dir) = pp(dir) + h
            call volume_fd%derivative2(Xt=pm, d2Tgc=d2Tm, dTgc=dTm, Tgc=Tm)
            call volume_fd%derivative2(Xt=pp, d2Tgc=d2Tp, dTgc=dTp, Tgc=Tp)
            cfd(:,dir) = (Tp - Tm)/(2.0_rk*h)
            cfd2(:,dir) = reshape((dTp - dTm)/(2.0_rk*h), [size(d2T0,1)])
        end do
        call volume_fd%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS volume first-derivative finite difference",&
            res      = norm2(cfd - dT0),&
            expected = 0.0_rk,&
            tol      = 1.0e-7_rk,&
            msg      = "Volume first-derivative finite difference is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0103


    subroutine forcad_nurbs_volume_0104(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume_fd
        real(rk), allocatable :: Wc_fd(:), T0(:), Tm(:), Tp(:), dT0(:,:), dTm(:,:), dTp(:,:)
        real(rk), allocatable :: d2T0(:,:), d2Tm(:,:), d2Tp(:,:), cfd(:,:), cfd2(:,:)
        real(rk) :: h, pt(3), pm(3), pp(3)
        integer :: dir

        allocate(Wc_fd(64), source=1.0_rk)
        Wc_fd(10) = 0.5_rk
        h = 1.0e-6_rk
        pt = [0.5_rk, 0.3_rk, 0.7_rk]

        call volume_fd%set_hexahedron(L=[2.0_rk, 4.0_rk, 8.0_rk], nc=[4,4,4], Wc=Wc_fd)
        call volume_fd%derivative2(Xt=pt, d2Tgc=d2T0, dTgc=dT0, Tgc=T0)
        allocate(cfd(size(T0),3), cfd2(size(d2T0,1),size(d2T0,2)))
        do dir = 1, 3
            pm = pt
            pp = pt
            pm(dir) = pm(dir) - h
            pp(dir) = pp(dir) + h
            call volume_fd%derivative2(Xt=pm, d2Tgc=d2Tm, dTgc=dTm, Tgc=Tm)
            call volume_fd%derivative2(Xt=pp, d2Tgc=d2Tp, dTgc=dTp, Tgc=Tp)
            cfd(:,dir) = (Tp - Tm)/(2.0_rk*h)
            cfd2(:,dir) = reshape((dTp - dTm)/(2.0_rk*h), [size(d2T0,1)])
        end do
        call volume_fd%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS volume second-derivative finite difference",&
            res      = norm2(cfd2 - d2T0),&
            expected = 0.0_rk,&
            tol      = 1.0e-6_rk,&
            msg      = "Volume second-derivative finite difference is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0104


    subroutine forcad_nurbs_volume_0105(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume_fd
        real(rk), allocatable :: Wc_fd(:), T0(:), Tm(:), Tp(:), dT0(:,:), dTm(:,:), dTp(:,:)
        real(rk), allocatable :: d2T0(:,:), d2Tm(:,:), d2Tp(:,:), cfd(:,:), cfd2(:,:)
        real(rk) :: h, pt(3), pm(3), pp(3)
        integer :: dir

        allocate(Wc_fd(64), source=1.0_rk)
        Wc_fd(10) = 0.5_rk
        h = 1.0e-6_rk
        pt = [0.5_rk, 0.3_rk, 0.7_rk]

        call volume_fd%set_hexahedron(L=[2.0_rk, 4.0_rk, 8.0_rk], nc=[4,4,4], Wc=Wc_fd)
        call volume_fd%derivative2(Xt=pt, d2Tgc=d2T0, dTgc=dT0, Tgc=T0)
        allocate(cfd(size(T0),3), cfd2(size(d2T0,1),size(d2T0,2)))
        do dir = 1, 3
            pm = pt
            pp = pt
            pm(dir) = pm(dir) - h
            pp(dir) = pp(dir) + h
            call volume_fd%derivative2(Xt=pm, d2Tgc=d2Tm, dTgc=dTm, Tgc=Tm)
            call volume_fd%derivative2(Xt=pp, d2Tgc=d2Tp, dTgc=dTp, Tgc=Tp)
            cfd(:,dir) = (Tp - Tm)/(2.0_rk*h)
            cfd2(:,dir) = reshape((dTp - dTm)/(2.0_rk*h), [size(d2T0,1)])
        end do

        call volume_fd%finalize()
        call volume_fd%set_hexahedron(L=[2.0_rk, 4.0_rk, 8.0_rk], nc=[4,4,4])
        call volume_fd%derivative2(Xt=pt, d2Tgc=d2T0, dTgc=dT0, Tgc=T0)
        do dir = 1, 3
            pm = pt
            pp = pt
            pm(dir) = pm(dir) - h
            pp(dir) = pp(dir) + h
            call volume_fd%derivative2(Xt=pm, d2Tgc=d2Tm, dTgc=dTm, Tgc=Tm)
            call volume_fd%derivative2(Xt=pp, d2Tgc=d2Tp, dTgc=dTp, Tgc=Tp)
            cfd(:,dir) = (Tp - Tm)/(2.0_rk*h)
            cfd2(:,dir) = reshape((dTp - dTm)/(2.0_rk*h), [size(d2T0,1)])
        end do
        call volume_fd%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline volume first-derivative finite difference",&
            res      = norm2(cfd - dT0),&
            expected = 0.0_rk,&
            tol      = 1.0e-7_rk,&
            msg      = "Volume first-derivative finite difference is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0105


    subroutine forcad_nurbs_volume_0106(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume_fd
        real(rk), allocatable :: Wc_fd(:), T0(:), Tm(:), Tp(:), dT0(:,:), dTm(:,:), dTp(:,:)
        real(rk), allocatable :: d2T0(:,:), d2Tm(:,:), d2Tp(:,:), cfd(:,:), cfd2(:,:)
        real(rk) :: h, pt(3), pm(3), pp(3)
        integer :: dir

        allocate(Wc_fd(64), source=1.0_rk)
        Wc_fd(10) = 0.5_rk
        h = 1.0e-6_rk
        pt = [0.5_rk, 0.3_rk, 0.7_rk]

        call volume_fd%set_hexahedron(L=[2.0_rk, 4.0_rk, 8.0_rk], nc=[4,4,4], Wc=Wc_fd)
        call volume_fd%derivative2(Xt=pt, d2Tgc=d2T0, dTgc=dT0, Tgc=T0)
        allocate(cfd(size(T0),3), cfd2(size(d2T0,1),size(d2T0,2)))
        do dir = 1, 3
            pm = pt
            pp = pt
            pm(dir) = pm(dir) - h
            pp(dir) = pp(dir) + h
            call volume_fd%derivative2(Xt=pm, d2Tgc=d2Tm, dTgc=dTm, Tgc=Tm)
            call volume_fd%derivative2(Xt=pp, d2Tgc=d2Tp, dTgc=dTp, Tgc=Tp)
            cfd(:,dir) = (Tp - Tm)/(2.0_rk*h)
            cfd2(:,dir) = reshape((dTp - dTm)/(2.0_rk*h), [size(d2T0,1)])
        end do

        call volume_fd%finalize()
        call volume_fd%set_hexahedron(L=[2.0_rk, 4.0_rk, 8.0_rk], nc=[4,4,4])
        call volume_fd%derivative2(Xt=pt, d2Tgc=d2T0, dTgc=dT0, Tgc=T0)
        do dir = 1, 3
            pm = pt
            pp = pt
            pm(dir) = pm(dir) - h
            pp(dir) = pp(dir) + h
            call volume_fd%derivative2(Xt=pm, d2Tgc=d2Tm, dTgc=dTm, Tgc=Tm)
            call volume_fd%derivative2(Xt=pp, d2Tgc=d2Tp, dTgc=dTp, Tgc=Tp)
            cfd(:,dir) = (Tp - Tm)/(2.0_rk*h)
            cfd2(:,dir) = reshape((dTp - dTm)/(2.0_rk*h), [size(d2T0,1)])
        end do
        call volume_fd%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline volume second-derivative finite difference",&
            res      = norm2(cfd2 - d2T0),&
            expected = 0.0_rk,&
            tol      = 1.0e-6_rk,&
            msg      = "Volume second-derivative finite difference is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0106


    subroutine forcad_nurbs_volume_0107(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume_ref
        real(rk), allocatable :: Xc_ref(:,:)
        real(rk), allocatable :: Xg_ref(:,:)
        real(rk) :: knot_ref(7), Xt_ref(3)
        integer :: ii, jj, kk, idx

        allocate(Xc_ref(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ii + 4*((jj - 1) + 4*(kk - 1))
                    Xc_ref(idx, 1) = real(ii - 1, rk)
                    Xc_ref(idx, 2) = real(jj - 1, rk)
                    Xc_ref(idx, 3) = real(kk - 1, rk) + 0.05_rk*real((ii - 1)*(jj - 1), rk)
                end do
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 3.0_rk, 4.0_rk]

        call volume_ref%set(knot_ref, knot_ref, knot_ref, Xc_ref, degree=[2, 2, 2])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)
        Xg_ref = volume_ref%get_Xg()

        call volume_ref%insert_knots(1, [2.5_rk], [1])
        call volume_ref%insert_knots(2, [3.0_rk], [1])
        call volume_ref%insert_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)
        call volume_ref%err%print()

        call ut%test(ti)%check(&
            name     = "volume knot insertion geometry preservation",&
            res      = maxval(abs(volume_ref%get_Xg() - Xg_ref)),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "Volume knot insertion geometry preservation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0107


    subroutine forcad_nurbs_volume_0108(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume_ref
        real(rk), allocatable :: Xc_ref(:,:)
        real(rk), allocatable :: Xg_ref(:,:)
        real(rk) :: knot_ref(7), Xt_ref(3)
        integer :: ii, jj, kk, idx

        allocate(Xc_ref(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ii + 4*((jj - 1) + 4*(kk - 1))
                    Xc_ref(idx, 1) = real(ii - 1, rk)
                    Xc_ref(idx, 2) = real(jj - 1, rk)
                    Xc_ref(idx, 3) = real(kk - 1, rk) + 0.05_rk*real((ii - 1)*(jj - 1), rk)
                end do
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 3.0_rk, 4.0_rk]

        call volume_ref%set(knot_ref, knot_ref, knot_ref, Xc_ref, degree=[2, 2, 2])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)
        Xg_ref = volume_ref%get_Xg()

        call volume_ref%insert_knots(1, [2.5_rk], [1])
        call volume_ref%insert_knots(2, [3.0_rk], [1])
        call volume_ref%insert_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)

        call volume_ref%remove_knots(1, [2.5_rk], [1])
        call volume_ref%remove_knots(2, [3.0_rk], [1])
        call volume_ref%remove_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)
        call volume_ref%err%print()

        call ut%test(ti)%check(&
            name     = "volume insertion-removal geometry preservation",&
            res      = maxval(abs(volume_ref%get_Xg() - Xg_ref)),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "Volume insertion-removal geometry preservation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0108


    subroutine forcad_nurbs_volume_0109(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume_ref
        real(rk), allocatable :: Xc_ref(:,:), Xg_ref(:,:), Xt1_after(:), Xt2_after(:), Xt3_after(:)
        real(rk) :: knot_ref(7), Xt_ref(3)
        integer :: ii, jj, kk, idx

        allocate(Xc_ref(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ii + 4*((jj - 1) + 4*(kk - 1))
                    Xc_ref(idx, 1) = real(ii - 1, rk)
                    Xc_ref(idx, 2) = real(jj - 1, rk)
                    Xc_ref(idx, 3) = real(kk - 1, rk) + 0.05_rk*real((ii - 1)*(jj - 1), rk)
                end do
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 3.0_rk, 4.0_rk]

        call volume_ref%set(knot_ref, knot_ref, knot_ref, Xc_ref, degree=[2, 2, 2])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)
        Xg_ref = volume_ref%get_Xg()

        call volume_ref%insert_knots(1, [2.5_rk], [1])
        call volume_ref%insert_knots(2, [3.0_rk], [1])
        call volume_ref%insert_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)

        call volume_ref%remove_knots(1, [2.5_rk], [1])
        call volume_ref%remove_knots(2, [3.0_rk], [1])
        call volume_ref%remove_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)

        call volume_ref%elevate_degree(1, 1)
        call volume_ref%elevate_degree(2, 1)
        call volume_ref%elevate_degree(3, 1)
        call volume_ref%create(res1=3, res2=3, res3=3)
        Xt1_after = volume_ref%get_Xt(1)
        Xt2_after = volume_ref%get_Xt(2)
        Xt3_after = volume_ref%get_Xt(3)
        call volume_ref%err%print()

        call ut%test(ti)%check(&
            name     = "elevated volume direction-one domain",&
            res      = [Xt1_after(1), Xt1_after(size(Xt1_after))],&
            expected = [2.0_rk, 4.0_rk],&
            tol      = 1.0e-10_rk,&
            msg      = "Elevated volume direction-one domain is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0109


    subroutine forcad_nurbs_volume_0110(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume_ref
        real(rk), allocatable :: Xc_ref(:,:), Xg_ref(:,:), Xt1_after(:), Xt2_after(:), Xt3_after(:)
        real(rk) :: knot_ref(7), Xt_ref(3)
        integer :: ii, jj, kk, idx

        allocate(Xc_ref(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ii + 4*((jj - 1) + 4*(kk - 1))
                    Xc_ref(idx, 1) = real(ii - 1, rk)
                    Xc_ref(idx, 2) = real(jj - 1, rk)
                    Xc_ref(idx, 3) = real(kk - 1, rk) + 0.05_rk*real((ii - 1)*(jj - 1), rk)
                end do
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 3.0_rk, 4.0_rk]

        call volume_ref%set(knot_ref, knot_ref, knot_ref, Xc_ref, degree=[2, 2, 2])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)
        Xg_ref = volume_ref%get_Xg()

        call volume_ref%insert_knots(1, [2.5_rk], [1])
        call volume_ref%insert_knots(2, [3.0_rk], [1])
        call volume_ref%insert_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)

        call volume_ref%remove_knots(1, [2.5_rk], [1])
        call volume_ref%remove_knots(2, [3.0_rk], [1])
        call volume_ref%remove_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)

        call volume_ref%elevate_degree(1, 1)
        call volume_ref%elevate_degree(2, 1)
        call volume_ref%elevate_degree(3, 1)
        call volume_ref%create(res1=3, res2=3, res3=3)
        Xt1_after = volume_ref%get_Xt(1)
        Xt2_after = volume_ref%get_Xt(2)
        Xt3_after = volume_ref%get_Xt(3)
        call volume_ref%err%print()

        call ut%test(ti)%check(&
            name     = "elevated volume direction-two domain",&
            res      = [Xt2_after(1), Xt2_after(size(Xt2_after))],&
            expected = [2.0_rk, 4.0_rk],&
            tol      = 1.0e-10_rk,&
            msg      = "Elevated volume direction-two domain is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0110


    subroutine forcad_nurbs_volume_0111(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume_ref
        real(rk), allocatable :: Xc_ref(:,:), Xg_ref(:,:), Xt1_after(:), Xt2_after(:), Xt3_after(:)
        real(rk) :: knot_ref(7), Xt_ref(3)
        integer :: ii, jj, kk, idx

        allocate(Xc_ref(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ii + 4*((jj - 1) + 4*(kk - 1))
                    Xc_ref(idx, 1) = real(ii - 1, rk)
                    Xc_ref(idx, 2) = real(jj - 1, rk)
                    Xc_ref(idx, 3) = real(kk - 1, rk) + 0.05_rk*real((ii - 1)*(jj - 1), rk)
                end do
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 3.0_rk, 4.0_rk]

        call volume_ref%set(knot_ref, knot_ref, knot_ref, Xc_ref, degree=[2, 2, 2])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)
        Xg_ref = volume_ref%get_Xg()

        call volume_ref%insert_knots(1, [2.5_rk], [1])
        call volume_ref%insert_knots(2, [3.0_rk], [1])
        call volume_ref%insert_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)

        call volume_ref%remove_knots(1, [2.5_rk], [1])
        call volume_ref%remove_knots(2, [3.0_rk], [1])
        call volume_ref%remove_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)

        call volume_ref%elevate_degree(1, 1)
        call volume_ref%elevate_degree(2, 1)
        call volume_ref%elevate_degree(3, 1)
        call volume_ref%create(res1=3, res2=3, res3=3)
        Xt1_after = volume_ref%get_Xt(1)
        Xt2_after = volume_ref%get_Xt(2)
        Xt3_after = volume_ref%get_Xt(3)
        call volume_ref%err%print()

        call ut%test(ti)%check(&
            name     = "elevated volume direction-three domain",&
            res      = [Xt3_after(1), Xt3_after(size(Xt3_after))],&
            expected = [2.0_rk, 4.0_rk],&
            tol      = 1.0e-10_rk,&
            msg      = "Elevated volume direction-three domain is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0111


    subroutine forcad_nurbs_volume_0112(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume_ref
        real(rk), allocatable :: Xc_ref(:,:), Xg_ref(:,:), Xt1_after(:), Xt2_after(:), Xt3_after(:)
        real(rk) :: knot_ref(7), Xt_ref(3)
        integer :: ii, jj, kk, idx

        allocate(Xc_ref(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ii + 4*((jj - 1) + 4*(kk - 1))
                    Xc_ref(idx, 1) = real(ii - 1, rk)
                    Xc_ref(idx, 2) = real(jj - 1, rk)
                    Xc_ref(idx, 3) = real(kk - 1, rk) + 0.05_rk*real((ii - 1)*(jj - 1), rk)
                end do
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 3.0_rk, 4.0_rk]

        call volume_ref%set(knot_ref, knot_ref, knot_ref, Xc_ref, degree=[2, 2, 2])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)
        Xg_ref = volume_ref%get_Xg()

        call volume_ref%insert_knots(1, [2.5_rk], [1])
        call volume_ref%insert_knots(2, [3.0_rk], [1])
        call volume_ref%insert_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)

        call volume_ref%remove_knots(1, [2.5_rk], [1])
        call volume_ref%remove_knots(2, [3.0_rk], [1])
        call volume_ref%remove_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)

        call volume_ref%elevate_degree(1, 1)
        call volume_ref%elevate_degree(2, 1)
        call volume_ref%elevate_degree(3, 1)
        call volume_ref%create(res1=3, res2=3, res3=3)
        Xt1_after = volume_ref%get_Xt(1)
        Xt2_after = volume_ref%get_Xt(2)
        Xt3_after = volume_ref%get_Xt(3)
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)
        call volume_ref%err%print()

        call ut%test(ti)%check(&
            name     = "volume degree elevation geometry preservation",&
            res      = maxval(abs(volume_ref%get_Xg() - Xg_ref)),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "Volume degree elevation geometry preservation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0112


    subroutine forcad_nurbs_volume_0113(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume_ref
        real(rk), allocatable :: Xc_ref(:,:), Xg_ref(:,:), Xt1_after(:), Xt2_after(:), Xt3_after(:)
        real(rk) :: knot_ref(7), Xt_ref(3)
        integer :: ii, jj, kk, idx

        allocate(Xc_ref(64, 3), source=0.0_rk)
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    idx = ii + 4*((jj - 1) + 4*(kk - 1))
                    Xc_ref(idx, 1) = real(ii - 1, rk)
                    Xc_ref(idx, 2) = real(jj - 1, rk)
                    Xc_ref(idx, 3) = real(kk - 1, rk) + 0.05_rk*real((ii - 1)*(jj - 1), rk)
                end do
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 3.0_rk, 4.0_rk]

        call volume_ref%set(knot_ref, knot_ref, knot_ref, Xc_ref, degree=[2, 2, 2])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)
        Xg_ref = volume_ref%get_Xg()

        call volume_ref%insert_knots(1, [2.5_rk], [1])
        call volume_ref%insert_knots(2, [3.0_rk], [1])
        call volume_ref%insert_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)

        call volume_ref%remove_knots(1, [2.5_rk], [1])
        call volume_ref%remove_knots(2, [3.0_rk], [1])
        call volume_ref%remove_knots(3, [3.5_rk], [1])
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)

        call volume_ref%elevate_degree(1, 1)
        call volume_ref%elevate_degree(2, 1)
        call volume_ref%elevate_degree(3, 1)
        call volume_ref%create(res1=3, res2=3, res3=3)
        Xt1_after = volume_ref%get_Xt(1)
        Xt2_after = volume_ref%get_Xt(2)
        Xt3_after = volume_ref%get_Xt(3)
        call volume_ref%create(Xt1=Xt_ref, Xt2=Xt_ref, Xt3=Xt_ref)
        call volume_ref%err%print()

        call ut%test(ti)%check(&
            name     = "elevated volume degree",&
            res      = volume_ref%get_degree(),&
            expected = [3, 3, 3],&
            msg      = "Elevated volume degree is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0113


    subroutine forcad_nurbs_volume_0114(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: sh0
        type(nurbs_volume) :: shr
        type(nurbs_volume) :: shb
        real(rk), allocatable :: Xc_s(:,:)
        real(rk), allocatable :: Wc_s(:)
        real(rk), allocatable :: Xc0(:,:)
        real(rk), allocatable :: knot1_s(:)
        real(rk), allocatable :: knot2_s(:)
        real(rk), allocatable :: knot3_s(:)
        real(rk), allocatable :: S1(:,:)
        real(rk), allocatable :: Bdirect(:,:)
        integer :: ncoord
        integer :: nc0
        integer :: ndof_old

        allocate(Xc_s(8,3), Wc_s(8), knot1_s(4), knot2_s(4), knot3_s(4))
        Xc_s(1,:) = [ 2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(2,:) = [ 2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(3,:) = [-2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(4,:) = [-2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(5,:) = [ 2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(6,:) = [ 2.5_rk,  2.5_rk, -2.5_rk]
        Xc_s(7,:) = [-2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(8,:) = [-2.5_rk,  2.5_rk, -2.5_rk]
        Wc_s = 1.0_rk
        Wc_s(2) = 0.5_rk
        knot1_s = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot2_s = knot1_s
        knot3_s = knot1_s

        call shr%set(knot1_s, knot2_s, knot3_s, Xc_s, Wc_s)
        sh0 = shr
        Xc0 = sh0%get_Xc()
        ncoord = size(Xc0,2)
        nc0 = size(Xc0,1)
        ndof_old = nc0*ncoord

        call shr%elevate_degree(1, 4, Bs=S1)
        call shb%set(sh0%get_knot(1), sh0%get_knot(2), sh0%get_knot(3), Xc0, sh0%get_Wc())
        call shb%elevate_degree(1, 4, B=Bdirect)
        call sh0%err%print()
        call shr%err%print()
        call shb%err%print()

        call ut%test(ti)%check(&
            name     = "volume direction-one degree-elevation map",&
            res      = maxval(abs(Bdirect - kron_eye(S1, ncoord))),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "volume elevate dir1 B-only map",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0114


    subroutine forcad_nurbs_volume_0115(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: sh0
        type(nurbs_volume) :: shr
        type(nurbs_volume) :: shb
        real(rk), allocatable :: Xc_s(:,:)
        real(rk), allocatable :: Wc_s(:)
        real(rk), allocatable :: Xc0(:,:)
        real(rk), allocatable :: knot1_s(:)
        real(rk), allocatable :: knot2_s(:)
        real(rk), allocatable :: knot3_s(:)
        real(rk), allocatable :: S1(:,:)
        real(rk), allocatable :: S2(:,:)
        real(rk), allocatable :: Bdirect(:,:)
        integer :: ncoord
        integer :: nc0
        integer :: ndof_old

        allocate(Xc_s(8,3), Wc_s(8), knot1_s(4), knot2_s(4), knot3_s(4))
        Xc_s(1,:) = [ 2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(2,:) = [ 2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(3,:) = [-2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(4,:) = [-2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(5,:) = [ 2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(6,:) = [ 2.5_rk,  2.5_rk, -2.5_rk]
        Xc_s(7,:) = [-2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(8,:) = [-2.5_rk,  2.5_rk, -2.5_rk]
        Wc_s = 1.0_rk
        Wc_s(2) = 0.5_rk
        knot1_s = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot2_s = knot1_s
        knot3_s = knot1_s

        call shr%set(knot1_s, knot2_s, knot3_s, Xc_s, Wc_s)
        sh0 = shr
        Xc0 = sh0%get_Xc()
        ncoord = size(Xc0,2)
        nc0 = size(Xc0,1)
        ndof_old = nc0*ncoord

        call shr%elevate_degree(1, 4, Bs=S1)
        call shb%set(sh0%get_knot(1), sh0%get_knot(2), sh0%get_knot(3), Xc0, sh0%get_Wc())
        call shb%elevate_degree(1, 4, B=Bdirect)
        call shr%elevate_degree(2, 3, Bs=S2)
        call shb%elevate_degree(2, 3, B=Bdirect)
        call sh0%err%print()
        call shr%err%print()
        call shb%err%print()

        call ut%test(ti)%check(&
            name     = "volume direction-two degree-elevation map",&
            res      = maxval(abs(Bdirect - kron_eye(S2, ncoord))),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "volume elevate dir2 B-only map",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0115


    subroutine forcad_nurbs_volume_0116(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: sh0
        type(nurbs_volume) :: shr
        type(nurbs_volume) :: shb
        real(rk), allocatable :: Xc_s(:,:)
        real(rk), allocatable :: Wc_s(:)
        real(rk), allocatable :: Xc0(:,:)
        real(rk), allocatable :: knot1_s(:)
        real(rk), allocatable :: knot2_s(:)
        real(rk), allocatable :: knot3_s(:)
        real(rk), allocatable :: S1(:,:)
        real(rk), allocatable :: S2(:,:)
        real(rk), allocatable :: S3(:,:)
        real(rk), allocatable :: Bdirect(:,:)
        integer :: ncoord
        integer :: nc0
        integer :: ndof_old

        allocate(Xc_s(8,3), Wc_s(8), knot1_s(4), knot2_s(4), knot3_s(4))
        Xc_s(1,:) = [ 2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(2,:) = [ 2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(3,:) = [-2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(4,:) = [-2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(5,:) = [ 2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(6,:) = [ 2.5_rk,  2.5_rk, -2.5_rk]
        Xc_s(7,:) = [-2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(8,:) = [-2.5_rk,  2.5_rk, -2.5_rk]
        Wc_s = 1.0_rk
        Wc_s(2) = 0.5_rk
        knot1_s = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot2_s = knot1_s
        knot3_s = knot1_s

        call shr%set(knot1_s, knot2_s, knot3_s, Xc_s, Wc_s)
        sh0 = shr
        Xc0 = sh0%get_Xc()
        ncoord = size(Xc0,2)
        nc0 = size(Xc0,1)
        ndof_old = nc0*ncoord

        call shr%elevate_degree(1, 4, Bs=S1)
        call shb%set(sh0%get_knot(1), sh0%get_knot(2), sh0%get_knot(3), Xc0, sh0%get_Wc())
        call shb%elevate_degree(1, 4, B=Bdirect)
        call shr%elevate_degree(2, 3, Bs=S2)
        call shb%elevate_degree(2, 3, B=Bdirect)
        call shr%elevate_degree(3, 2, Bs=S3)
        call shb%elevate_degree(3, 2, B=Bdirect)
        call sh0%err%print()
        call shr%err%print()
        call shb%err%print()

        call ut%test(ti)%check(&
            name     = "volume direction-three degree-elevation map",&
            res      = maxval(abs(Bdirect - kron_eye(S3, ncoord))),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "volume elevate dir3 B-only map",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0116


    subroutine forcad_nurbs_volume_0117(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: sh0
        type(nurbs_volume) :: shr
        type(nurbs_volume) :: shb
        real(rk), allocatable :: Xc_s(:,:)
        real(rk), allocatable :: Wc_s(:)
        real(rk), allocatable :: Xc0(:,:)
        real(rk), allocatable :: knot1_s(:)
        real(rk), allocatable :: knot2_s(:)
        real(rk), allocatable :: knot3_s(:)
        real(rk), allocatable :: S1(:,:)
        real(rk), allocatable :: S2(:,:)
        real(rk), allocatable :: S3(:,:)
        real(rk), allocatable :: S4(:,:)
        real(rk), allocatable :: Bdirect(:,:)
        real(rk), allocatable :: u1(:), u2(:), u3(:)
        integer, allocatable :: r1(:), r2(:), r3(:)
        integer :: ncoord
        integer :: nc0
        integer :: ndof_old

        allocate(Xc_s(8,3), Wc_s(8), knot1_s(4), knot2_s(4), knot3_s(4))
        Xc_s(1,:) = [ 2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(2,:) = [ 2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(3,:) = [-2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(4,:) = [-2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(5,:) = [ 2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(6,:) = [ 2.5_rk,  2.5_rk, -2.5_rk]
        Xc_s(7,:) = [-2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(8,:) = [-2.5_rk,  2.5_rk, -2.5_rk]
        Wc_s = 1.0_rk
        Wc_s(2) = 0.5_rk
        knot1_s = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot2_s = knot1_s
        knot3_s = knot1_s

        call shr%set(knot1_s, knot2_s, knot3_s, Xc_s, Wc_s)
        sh0 = shr
        Xc0 = sh0%get_Xc()
        ncoord = size(Xc0,2)
        nc0 = size(Xc0,1)
        ndof_old = nc0*ncoord

        call shr%elevate_degree(1, 4, Bs=S1)
        call shb%set(sh0%get_knot(1), sh0%get_knot(2), sh0%get_knot(3), Xc0, sh0%get_Wc())
        call shb%elevate_degree(1, 4, B=Bdirect)
        call shr%elevate_degree(2, 3, Bs=S2)
        call shb%elevate_degree(2, 3, B=Bdirect)
        call shr%elevate_degree(3, 2, Bs=S3)
        call shb%elevate_degree(3, 2, B=Bdirect)
        u1 = linspace(0.0_rk, 1.0_rk, 4)
        u1 = u1(2:3)
        u2 = linspace(0.0_rk, 1.0_rk, 5)
        u2 = u2(2:4)
        u3 = linspace(0.0_rk, 1.0_rk, 6)
        u3 = u3(2:5)
        allocate(r1(size(u1)), source=2)
        allocate(r2(size(u2)), source=1)
        allocate(r3(size(u3)), source=2)
        call shr%insert_knots(1, u1, r1, Bs=S4)
        call shb%insert_knots(1, u1, r1, B=Bdirect)
        call sh0%err%print()
        call shr%err%print()
        call shb%err%print()

        call ut%test(ti)%check(&
            name     = "volume direction-one knot-insertion map",&
            res      = maxval(abs(Bdirect - kron_eye(S4, ncoord))),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "volume insert dir1 B-only map",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0117


    subroutine forcad_nurbs_volume_0118(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: sh0
        type(nurbs_volume) :: shr
        type(nurbs_volume) :: shb
        real(rk), allocatable :: Xc_s(:,:)
        real(rk), allocatable :: Wc_s(:)
        real(rk), allocatable :: Xc0(:,:)
        real(rk), allocatable :: knot1_s(:)
        real(rk), allocatable :: knot2_s(:)
        real(rk), allocatable :: knot3_s(:)
        real(rk), allocatable :: S1(:,:)
        real(rk), allocatable :: S2(:,:)
        real(rk), allocatable :: S3(:,:)
        real(rk), allocatable :: S4(:,:)
        real(rk), allocatable :: S5(:,:)
        real(rk), allocatable :: Bdirect(:,:)
        real(rk), allocatable :: u1(:), u2(:), u3(:)
        integer, allocatable :: r1(:), r2(:), r3(:)
        integer :: ncoord
        integer :: nc0
        integer :: ndof_old

        allocate(Xc_s(8,3), Wc_s(8), knot1_s(4), knot2_s(4), knot3_s(4))
        Xc_s(1,:) = [ 2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(2,:) = [ 2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(3,:) = [-2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(4,:) = [-2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(5,:) = [ 2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(6,:) = [ 2.5_rk,  2.5_rk, -2.5_rk]
        Xc_s(7,:) = [-2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(8,:) = [-2.5_rk,  2.5_rk, -2.5_rk]
        Wc_s = 1.0_rk
        Wc_s(2) = 0.5_rk
        knot1_s = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot2_s = knot1_s
        knot3_s = knot1_s

        call shr%set(knot1_s, knot2_s, knot3_s, Xc_s, Wc_s)
        sh0 = shr
        Xc0 = sh0%get_Xc()
        ncoord = size(Xc0,2)
        nc0 = size(Xc0,1)
        ndof_old = nc0*ncoord

        call shr%elevate_degree(1, 4, Bs=S1)
        call shb%set(sh0%get_knot(1), sh0%get_knot(2), sh0%get_knot(3), Xc0, sh0%get_Wc())
        call shb%elevate_degree(1, 4, B=Bdirect)
        call shr%elevate_degree(2, 3, Bs=S2)
        call shb%elevate_degree(2, 3, B=Bdirect)
        call shr%elevate_degree(3, 2, Bs=S3)
        call shb%elevate_degree(3, 2, B=Bdirect)
        u1 = linspace(0.0_rk, 1.0_rk, 4)
        u1 = u1(2:3)
        u2 = linspace(0.0_rk, 1.0_rk, 5)
        u2 = u2(2:4)
        u3 = linspace(0.0_rk, 1.0_rk, 6)
        u3 = u3(2:5)
        allocate(r1(size(u1)), source=2)
        allocate(r2(size(u2)), source=1)
        allocate(r3(size(u3)), source=2)
        call shr%insert_knots(1, u1, r1, Bs=S4)
        call shb%insert_knots(1, u1, r1, B=Bdirect)
        call shr%insert_knots(2, u2, r2, Bs=S5)
        call shb%insert_knots(2, u2, r2, B=Bdirect)
        call sh0%err%print()
        call shr%err%print()
        call shb%err%print()

        call ut%test(ti)%check(&
            name     = "volume direction-two knot-insertion map",&
            res      = maxval(abs(Bdirect - kron_eye(S5, ncoord))),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "volume insert dir2 B-only map",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0118


    subroutine forcad_nurbs_volume_0119(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: sh0
        type(nurbs_volume) :: shr
        type(nurbs_volume) :: shb
        real(rk), allocatable :: Xc_s(:,:)
        real(rk), allocatable :: Wc_s(:)
        real(rk), allocatable :: Xc0(:,:)
        real(rk), allocatable :: knot1_s(:)
        real(rk), allocatable :: knot2_s(:)
        real(rk), allocatable :: knot3_s(:)
        real(rk), allocatable :: S1(:,:)
        real(rk), allocatable :: S2(:,:)
        real(rk), allocatable :: S3(:,:)
        real(rk), allocatable :: S4(:,:)
        real(rk), allocatable :: S5(:,:)
        real(rk), allocatable :: S6(:,:)
        real(rk), allocatable :: Bdirect(:,:)
        real(rk), allocatable :: u1(:), u2(:), u3(:)
        integer, allocatable :: r1(:), r2(:), r3(:)
        integer :: ncoord
        integer :: nc0
        integer :: ndof_old

        allocate(Xc_s(8,3), Wc_s(8), knot1_s(4), knot2_s(4), knot3_s(4))
        Xc_s(1,:) = [ 2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(2,:) = [ 2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(3,:) = [-2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(4,:) = [-2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(5,:) = [ 2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(6,:) = [ 2.5_rk,  2.5_rk, -2.5_rk]
        Xc_s(7,:) = [-2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(8,:) = [-2.5_rk,  2.5_rk, -2.5_rk]
        Wc_s = 1.0_rk
        Wc_s(2) = 0.5_rk
        knot1_s = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot2_s = knot1_s
        knot3_s = knot1_s

        call shr%set(knot1_s, knot2_s, knot3_s, Xc_s, Wc_s)
        sh0 = shr
        Xc0 = sh0%get_Xc()
        ncoord = size(Xc0,2)
        nc0 = size(Xc0,1)
        ndof_old = nc0*ncoord

        call shr%elevate_degree(1, 4, Bs=S1)
        call shb%set(sh0%get_knot(1), sh0%get_knot(2), sh0%get_knot(3), Xc0, sh0%get_Wc())
        call shb%elevate_degree(1, 4, B=Bdirect)
        call shr%elevate_degree(2, 3, Bs=S2)
        call shb%elevate_degree(2, 3, B=Bdirect)
        call shr%elevate_degree(3, 2, Bs=S3)
        call shb%elevate_degree(3, 2, B=Bdirect)
        u1 = linspace(0.0_rk, 1.0_rk, 4)
        u1 = u1(2:3)
        u2 = linspace(0.0_rk, 1.0_rk, 5)
        u2 = u2(2:4)
        u3 = linspace(0.0_rk, 1.0_rk, 6)
        u3 = u3(2:5)
        allocate(r1(size(u1)), source=2)
        allocate(r2(size(u2)), source=1)
        allocate(r3(size(u3)), source=2)
        call shr%insert_knots(1, u1, r1, Bs=S4)
        call shb%insert_knots(1, u1, r1, B=Bdirect)
        call shr%insert_knots(2, u2, r2, Bs=S5)
        call shb%insert_knots(2, u2, r2, B=Bdirect)
        call shr%insert_knots(3, u3, r3, Bs=S6)
        call shb%insert_knots(3, u3, r3, B=Bdirect)
        call sh0%err%print()
        call shr%err%print()
        call shb%err%print()

        call ut%test(ti)%check(&
            name     = "volume direction-three knot-insertion map",&
            res      = maxval(abs(Bdirect - kron_eye(S6, ncoord))),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "volume insert dir3 B-only map",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0119


    subroutine forcad_nurbs_volume_0120(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: sh0, shr, shfd, shb
        real(rk), allocatable :: Xc_s(:,:), Wc_s(:), Xc0(:,:), Xp(:,:), Xm(:,:)
        real(rk), allocatable :: Xcp_vec(:), Xcm_vec(:), knot1_s(:), knot2_s(:), knot3_s(:)
        real(rk), allocatable :: S1(:,:), S2(:,:), S3(:,:), S4(:,:), S5(:,:), S6(:,:), Bs(:,:), B(:,:), Bdirect(:,:), Bfd(:,:)
        real(rk), allocatable :: u1(:), u2(:), u3(:)
        integer, allocatable :: r1(:), r2(:), r3(:)
        real(rk) :: rel_err
        integer :: ncoord, nc0, ndof_old, ndof_new, idx, ic, dc

        allocate(Xc_s(8,3), Wc_s(8), knot1_s(4), knot2_s(4), knot3_s(4))
        Xc_s(1,:) = [ 2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(2,:) = [ 2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(3,:) = [-2.5_rk, -2.5_rk,  2.5_rk]
        Xc_s(4,:) = [-2.5_rk, -2.5_rk, -2.5_rk]
        Xc_s(5,:) = [ 2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(6,:) = [ 2.5_rk,  2.5_rk, -2.5_rk]
        Xc_s(7,:) = [-2.5_rk,  2.5_rk,  2.5_rk]
        Xc_s(8,:) = [-2.5_rk,  2.5_rk, -2.5_rk]
        Wc_s = 1.0_rk
        Wc_s(2) = 0.5_rk
        knot1_s = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot2_s = knot1_s
        knot3_s = knot1_s

        call shr%set(knot1_s, knot2_s, knot3_s, Xc_s, Wc_s)
        sh0 = shr
        Xc0 = sh0%get_Xc()
        ncoord = size(Xc0,2)
        nc0 = size(Xc0,1)
        ndof_old = nc0*ncoord

        call shr%elevate_degree(1, 4, Bs=S1)
        call shb%set(sh0%get_knot(1), sh0%get_knot(2), sh0%get_knot(3), Xc0, sh0%get_Wc())
        call shb%elevate_degree(1, 4, B=Bdirect)
        call shr%elevate_degree(2, 3, Bs=S2)
        call shb%elevate_degree(2, 3, B=Bdirect)
        call shr%elevate_degree(3, 2, Bs=S3)
        call shb%elevate_degree(3, 2, B=Bdirect)
        u1 = linspace(0.0_rk, 1.0_rk, 4)
        u1 = u1(2:3)
        u2 = linspace(0.0_rk, 1.0_rk, 5)
        u2 = u2(2:4)
        u3 = linspace(0.0_rk, 1.0_rk, 6)
        u3 = u3(2:5)
        allocate(r1(size(u1)), source=2)
        allocate(r2(size(u2)), source=1)
        allocate(r3(size(u3)), source=2)
        call shr%insert_knots(1, u1, r1, Bs=S4)
        call shb%insert_knots(1, u1, r1, B=Bdirect)
        call shr%insert_knots(2, u2, r2, Bs=S5)
        call shb%insert_knots(2, u2, r2, B=Bdirect)
        call shr%insert_knots(3, u3, r3, Bs=S6)
        call shb%insert_knots(3, u3, r3, B=Bdirect)
        Bs = matmul(S6, matmul(S5, matmul(S4, matmul(S3, matmul(S2, S1)))))
        B = kron_eye(Bs, ncoord)
        ndof_new = size(shr%get_Xc(),1)*ncoord

        allocate(Xp(nc0,ncoord), Xm(nc0,ncoord), Xcp_vec(ndof_new), Xcm_vec(ndof_new), Bfd(ndof_new,ndof_old))
        do idx = 1, ndof_old
            Xp = Xc0
            Xm = Xc0
            ic = (idx - 1)/ncoord + 1
            dc = mod(idx - 1, ncoord) + 1
            Xp(ic,dc) = Xp(ic,dc) + 1.0e-5_rk
            Xm(ic,dc) = Xm(ic,dc) - 1.0e-5_rk

            call shfd%set(sh0%get_knot(1), sh0%get_knot(2), sh0%get_knot(3), Xp, sh0%get_Wc())
            call shfd%elevate_degree(1, 4)
            call shfd%elevate_degree(2, 3)
            call shfd%elevate_degree(3, 2)
            call shfd%insert_knots(1, u1, r1)
            call shfd%insert_knots(2, u2, r2)
            call shfd%insert_knots(3, u3, r3)
            Xcp_vec = reshape(transpose(shfd%get_Xc()), [ndof_new])

            call shfd%set(sh0%get_knot(1), sh0%get_knot(2), sh0%get_knot(3), Xm, sh0%get_Wc())
            call shfd%elevate_degree(1, 4)
            call shfd%elevate_degree(2, 3)
            call shfd%elevate_degree(3, 2)
            call shfd%insert_knots(1, u1, r1)
            call shfd%insert_knots(2, u2, r2)
            call shfd%insert_knots(3, u3, r3)
            Xcm_vec = reshape(transpose(shfd%get_Xc()), [ndof_new])

            Bfd(:,idx) = (Xcp_vec - Xcm_vec) * 5.0e4_rk
        end do
        rel_err = norm2(Bfd - B)/norm2(Bfd)
        call sh0%err%print()
        call shr%err%print()
        call shfd%err%print()
        call shb%err%print()

        call ut%test(ti)%check(&
            name     = "volume refinement-map finite difference",&
            res      = rel_err,&
            expected = 0.0_rk,&
            tol      = 1.0e-7_rk,&
            msg      = "Volume refinement-map finite difference is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0120


    subroutine forcad_nurbs_volume_0121(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vproj
        real(rk) :: Xproj(8,3), Xt_out(3), Xg_out(3)

        Xproj(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xproj(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xproj(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xproj(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xproj(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xproj(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xproj(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xproj(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        call vproj%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xproj)
        call vproj%nearest_point2([0.2_rk, 0.3_rk, 0.4_rk], 1.0e-12_rk, 0, Xt_out, Xg_out)
        call vproj%err%print()

        call ut%test(ti)%check(&
            name     = "zero-iteration volume projection outputs",&
            res      = all(Xt_out >= 0.0_rk .and. Xt_out <= 1.0_rk) .and. all(ieee_is_finite(Xg_out)),&
            expected = .true.,&
            msg      = "Zero-iteration volume projection outputs are incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0121


    subroutine forcad_nurbs_volume_0122(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vbad
        real(rk) :: Xbad(8,3)
        real(rk) :: Xnew(12,3)
        real(rk) :: Wbad(8)
        real(rk) :: Wbad12(12)
        real(rk) :: Wgood(8)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xbad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xbad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xbad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xbad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Xnew(7,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xnew(8,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xnew(9,:) = [2.0_rk, 0.0_rk, 1.0_rk]
        Xnew(10,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xnew(11,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(12,:) = [2.0_rk, 1.0_rk, 1.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad12 = 1.0_rk
        Wbad12(2) = 0.0_rk
        Wgood = 1.0_rk

        call vbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)
        call vbad%err%print()

        call ut%test(ti)%check(&
            name     = "volume rejects nonpositive weights",&
            res      = vbad%err%ok,&
            expected = .false.,&
            msg      = "A volume must reject nonpositive weights.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0122


    subroutine forcad_nurbs_volume_0123(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vbad
        real(rk) :: Xbad(8,3)
        real(rk) :: Xnew(12,3)
        real(rk) :: Wbad(8)
        real(rk) :: Wbad12(12)
        real(rk) :: Wgood(8)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xbad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xbad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xbad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xbad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Xnew(7,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xnew(8,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xnew(9,:) = [2.0_rk, 0.0_rk, 1.0_rk]
        Xnew(10,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xnew(11,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(12,:) = [2.0_rk, 1.0_rk, 1.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad12 = 1.0_rk
        Wbad12(2) = 0.0_rk
        Wgood = 1.0_rk

        call vbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)
        call vbad%err%print()

        call ut%test(ti)%check(&
            name     = "invalid-weight volume rational state",&
            res      = vbad%is_rational(),&
            expected = .false.,&
            msg      = "Invalid-weight volume rational state is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0123


    subroutine forcad_nurbs_volume_0124(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vbad
        type(nurbs_volume) :: vmod
        real(rk) :: Xbad(8,3)
        real(rk) :: Xnew(12,3)
        real(rk) :: Wbad(8)
        real(rk) :: Wbad12(12)
        real(rk) :: Wgood(8)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xbad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xbad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xbad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xbad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Xnew(7,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xnew(8,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xnew(9,:) = [2.0_rk, 0.0_rk, 1.0_rk]
        Xnew(10,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xnew(11,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(12,:) = [2.0_rk, 1.0_rk, 1.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad12 = 1.0_rk
        Wbad12(2) = 0.0_rk
        Wgood = 1.0_rk

        call vbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call vmod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call vmod%modify_Wc(0.5_rk, 1)
        call vbad%err%print()
        call vmod%err%print()

        call ut%test(ti)%check(&
            name     = "modified volume rational state",&
            res      = vmod%is_rational(),&
            expected = .true.,&
            msg      = "volume rational cache update",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0124


    subroutine forcad_nurbs_volume_0125(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vbad
        type(nurbs_volume) :: vmod
        real(rk) :: Xbad(8,3)
        real(rk) :: Xnew(12,3)
        real(rk) :: Wbad(8)
        real(rk) :: Wbad12(12)
        real(rk) :: Wgood(8)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xbad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xbad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xbad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xbad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Xnew(7,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xnew(8,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xnew(9,:) = [2.0_rk, 0.0_rk, 1.0_rk]
        Xnew(10,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xnew(11,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(12,:) = [2.0_rk, 1.0_rk, 1.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad12 = 1.0_rk
        Wbad12(2) = 0.0_rk
        Wgood = 1.0_rk

        call vbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call vmod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call vmod%modify_Wc(0.5_rk, 1)
        call vmod%modify_Wc(1.0_rk, 1)
        call vbad%err%print()
        call vmod%err%print()

        call ut%test(ti)%check(&
            name     = "restored uniform volume rational state",&
            res      = vmod%is_rational(),&
            expected = .false.,&
            msg      = "volume uniform cache restore",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0125


    subroutine forcad_nurbs_volume_0126(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vbad
        type(nurbs_volume) :: vmod
        real(rk) :: Xbad(8,3)
        real(rk) :: Xnew(12,3)
        real(rk) :: Wbad(8)
        real(rk) :: Wbad12(12)
        real(rk) :: Wgood(8)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xbad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xbad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xbad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xbad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Xnew(7,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xnew(8,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xnew(9,:) = [2.0_rk, 0.0_rk, 1.0_rk]
        Xnew(10,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xnew(11,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(12,:) = [2.0_rk, 1.0_rk, 1.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad12 = 1.0_rk
        Wbad12(2) = 0.0_rk
        Wgood = 1.0_rk

        call vbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call vmod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call vmod%modify_Wc(0.5_rk, 1)
        call vmod%modify_Wc(1.0_rk, 1)
        call vmod%modify_Wc(-1.0_rk, 1)
        call vbad%err%print()
        call vmod%err%print()

        call ut%test(ti)%check(&
            name     = "volume rejects negative weight modification",&
            res      = vmod%err%ok,&
            expected = .false.,&
            msg      = "Volume rejects negative weight modification is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0126


    subroutine forcad_nurbs_volume_0127(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vbad
        type(nurbs_volume) :: vmod
        type(nurbs_volume) :: viges
        real(rk) :: Xbad(8,3)
        real(rk) :: Xnew(12,3)
        real(rk) :: Wbad(8)
        real(rk) :: Wbad12(12)
        real(rk) :: Wgood(8)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xbad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xbad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xbad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xbad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Xnew(7,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xnew(8,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xnew(9,:) = [2.0_rk, 0.0_rk, 1.0_rk]
        Xnew(10,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xnew(11,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(12,:) = [2.0_rk, 1.0_rk, 1.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad12 = 1.0_rk
        Wbad12(2) = 0.0_rk
        Wgood = 1.0_rk

        call vbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call vmod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call vmod%modify_Wc(0.5_rk, 1)
        call vmod%modify_Wc(1.0_rk, 1)
        call vmod%modify_Wc(-1.0_rk, 1)

        call viges%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)
        call viges%export_iges("vtk/forcad_volume_unsupported.iges")
        call vbad%err%print()
        call vmod%err%print()
        call viges%err%print()

        call ut%test(ti)%check(&
            name     = "unsupported volume IGES export rejection",&
            res      = viges%err%ok,&
            expected = .false.,&
            msg      = "Unsupported volume IGES export rejection is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0127


    subroutine forcad_nurbs_volume_0128(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vbad
        type(nurbs_volume) :: vmod
        type(nurbs_volume) :: viges
        type(nurbs_volume) :: vstate
        real(rk) :: Xbad(8,3), Xnew(12,3), Wbad(8), Wbad12(12), Wgood(8), Wstate(8)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xbad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xbad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xbad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xbad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Xnew(7,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xnew(8,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xnew(9,:) = [2.0_rk, 0.0_rk, 1.0_rk]
        Xnew(10,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xnew(11,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(12,:) = [2.0_rk, 1.0_rk, 1.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad12 = 1.0_rk
        Wbad12(2) = 0.0_rk
        Wgood = 1.0_rk

        call vbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call vmod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call vmod%modify_Wc(0.5_rk, 1)
        call vmod%modify_Wc(1.0_rk, 1)
        call vmod%modify_Wc(-1.0_rk, 1)

        call viges%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)
        call viges%export_iges("vtk/forcad_volume_unsupported.iges")

        Wstate = 1.0_rk
        Wstate(2) = 0.5_rk
        call vstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wstate)
        call vstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)
        call vbad%err%print()
        call vmod%err%print()
        call viges%err%print()
        call vstate%err%print()

        call ut%test(ti)%check(&
            name     = "unweighted volume reset clears rational state",&
            res      = vstate%is_rational(),&
            expected = .false.,&
            msg      = "Unweighted volume reset clears rational state is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0128


    subroutine forcad_nurbs_volume_0129(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vbad
        type(nurbs_volume) :: vmod
        type(nurbs_volume) :: viges
        type(nurbs_volume) :: vstate
        type(nurbs_volume) :: vdegree1
        real(rk) :: Xbad(8,3), Xnew(12,3), Wbad(8), Wbad12(12), Wgood(8), Wstate(8)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xbad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xbad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xbad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xbad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Xnew(7,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xnew(8,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xnew(9,:) = [2.0_rk, 0.0_rk, 1.0_rk]
        Xnew(10,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xnew(11,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(12,:) = [2.0_rk, 1.0_rk, 1.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad12 = 1.0_rk
        Wbad12(2) = 0.0_rk
        Wgood = 1.0_rk

        call vbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call vmod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call vmod%modify_Wc(0.5_rk, 1)
        call vmod%modify_Wc(1.0_rk, 1)
        call vmod%modify_Wc(-1.0_rk, 1)

        call viges%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)
        call viges%export_iges("vtk/forcad_volume_unsupported.iges")

        Wstate = 1.0_rk
        Wstate(2) = 0.5_rk
        call vstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wstate)
        call vstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)

        call vdegree1%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, degree=[1, 1])
        call vbad%err%print()
        call vmod%err%print()
        call viges%err%print()
        call vstate%err%print()
        call vdegree1%err%print()

        call ut%test(ti)%check(&
            name     = "volume rejects incomplete degree vector",&
            res      = vdegree1%err%ok,&
            expected = .false.,&
            msg      = "Volume rejects incomplete degree vector is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0129


    subroutine forcad_nurbs_volume_0130(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vbad
        type(nurbs_volume) :: vmod
        type(nurbs_volume) :: viges
        type(nurbs_volume) :: vstate
        type(nurbs_volume) :: vdegree1
        type(nurbs_volume) :: vdegree2
        real(rk) :: Xbad(8,3), Xnew(12,3), Wbad(8), Wbad12(12), Wgood(8), Wstate(8)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xbad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xbad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xbad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xbad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Xnew(7,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xnew(8,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xnew(9,:) = [2.0_rk, 0.0_rk, 1.0_rk]
        Xnew(10,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xnew(11,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(12,:) = [2.0_rk, 1.0_rk, 1.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad12 = 1.0_rk
        Wbad12(2) = 0.0_rk
        Wgood = 1.0_rk

        call vbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call vmod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call vmod%modify_Wc(0.5_rk, 1)
        call vmod%modify_Wc(1.0_rk, 1)
        call vmod%modify_Wc(-1.0_rk, 1)

        call viges%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)
        call viges%export_iges("vtk/forcad_volume_unsupported.iges")

        Wstate = 1.0_rk
        Wstate(2) = 0.5_rk
        call vstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wstate)
        call vstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)

        call vdegree1%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, degree=[1, 1])

        call vdegree2%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], &
            [1, 1], [-1, -1], [-1, -1], [-1, -1])
        call vbad%err%print()
        call vmod%err%print()
        call viges%err%print()
        call vstate%err%print()
        call vdegree1%err%print()
        call vdegree2%err%print()

        call ut%test(ti)%check(&
            name     = "volume rejects inconsistent parameter metadata",&
            res      = vdegree2%err%ok,&
            expected = .false.,&
            msg      = "Volume rejects inconsistent parameter metadata is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0130


    subroutine forcad_nurbs_volume_0131(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vbad, vmod, viges, vstate, vdegree1, vdegree2, vatomic
        real(rk) :: Xbad(8,3), Xnew(12,3), Wbad(8), Wbad12(12), Wgood(8), Wstate(8)
        real(rk), allocatable :: Xkeep(:,:)
        logical :: rejected
        integer :: preserved

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xbad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xbad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xbad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xbad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Xnew(7,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xnew(8,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xnew(9,:) = [2.0_rk, 0.0_rk, 1.0_rk]
        Xnew(10,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xnew(11,:) = [1.0_rk, 1.0_rk, 1.0_rk]
        Xnew(12,:) = [2.0_rk, 1.0_rk, 1.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad12 = 1.0_rk
        Wbad12(2) = 0.0_rk
        Wgood = 1.0_rk

        call vbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call vmod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call vmod%modify_Wc(0.5_rk, 1)
        call vmod%modify_Wc(1.0_rk, 1)
        call vmod%modify_Wc(-1.0_rk, 1)

        call viges%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)
        call viges%export_iges("vtk/forcad_volume_unsupported.iges")

        Wstate = 1.0_rk
        Wstate(2) = 0.5_rk
        call vstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wstate)
        call vstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)

        call vdegree1%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, degree=[1, 1])

        call vdegree2%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], &
            [1, 1], [-1, -1], [-1, -1], [-1, -1])

        call vatomic%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wstate)
        call vatomic%set([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xnew, Wbad12, degree=[2, 1, 1])
        rejected = .not. vatomic%err%ok
        call vatomic%err%print()
        call vatomic%err%reset()
        preserved = 0
        if (rejected .and. all(vatomic%get_nc() == [2, 2, 2]) .and. all(vatomic%get_degree() == [1, 1, 1])) then
            Xkeep = vatomic%get_Xc()
            if (size(Xkeep, 1) == size(Xbad, 1) .and. size(Xkeep, 2) == size(Xbad, 2)) then
                if (maxval(abs(Xkeep - Xbad)) <= 1.0e-12_rk .and. vatomic%is_rational()) preserved = 1
            end if
        end if
        call vbad%err%print()
        call vmod%err%print()
        call viges%err%print()
        call vstate%err%print()
        call vdegree1%err%print()
        call vdegree2%err%print()
        call vatomic%err%print()

        call ut%test(ti)%check(&
            name     = "failed volume reset preserves object state",&
            res      = preserved,&
            expected = 1,&
            msg      = "Failed volume reset preserves object state is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0131


    subroutine forcad_nurbs_volume_0132(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vinsert
        real(rk) :: knot_bad(4), Xc_bad(8,3)

        knot_bad = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        Xc_bad = 0.0_rk
        Xc_bad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xc_bad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xc_bad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xc_bad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xc_bad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xc_bad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xc_bad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xc_bad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]

        call vinsert%set(knot_bad, knot_bad, knot_bad, Xc_bad)
        call vinsert%insert_knots(1, [0.5_rk], [3])
        call vinsert%err%print()

        call ut%test(ti)%check(&
            name     = "volume rejects excessive knot insertion",&
            res      = vinsert%err%ok,&
            expected = .false.,&
            msg      = "volume bad insert mult",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0132


    subroutine forcad_nurbs_volume_0133(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vinsert
        type(nurbs_volume) :: vremove
        real(rk) :: knot_bad(4), Xc_bad(8,3)

        knot_bad = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        Xc_bad = 0.0_rk
        Xc_bad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xc_bad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xc_bad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xc_bad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xc_bad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xc_bad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xc_bad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xc_bad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]

        call vinsert%set(knot_bad, knot_bad, knot_bad, Xc_bad)
        call vinsert%insert_knots(1, [0.5_rk], [3])

        call vremove%set(knot_bad, knot_bad, knot_bad, Xc_bad)
        call vremove%remove_knots(1, [0.5_rk], [1])
        call vinsert%err%print()
        call vremove%err%print()

        call ut%test(ti)%check(&
            name     = "volume rejects absent knot removal",&
            res      = vremove%err%ok,&
            expected = .false.,&
            msg      = "volume absent remove knot",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0133


    subroutine forcad_nurbs_volume_0134(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vinsert, vremove, velev
        real(rk) :: knot_bad(4), Xc_bad(8,3)

        knot_bad = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        Xc_bad = 0.0_rk
        Xc_bad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xc_bad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xc_bad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xc_bad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xc_bad(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        Xc_bad(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xc_bad(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        Xc_bad(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]

        call vinsert%set(knot_bad, knot_bad, knot_bad, Xc_bad)
        call vinsert%insert_knots(1, [0.5_rk], [3])

        call vremove%set(knot_bad, knot_bad, knot_bad, Xc_bad)
        call vremove%remove_knots(1, [0.5_rk], [1])

        call velev%set(knot_bad, knot_bad, knot_bad, Xc_bad)
        call velev%elevate_degree(1, -1)
        call vinsert%err%print()
        call vremove%err%print()
        call velev%err%print()

        call ut%test(ti)%check(&
            name     = "volume rejects negative degree elevation",&
            res      = velev%err%ok,&
            expected = .false.,&
            msg      = "volume negative degree elevation",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0134


    subroutine forcad_nurbs_volume_0135(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)

        ivals = vempty%get_multiplicity(1)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume multiplicity getter",&
            res      = size(ivals),&
            expected = 0,&
            msg      = "volume unset mult",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0135


    subroutine forcad_nurbs_volume_0136(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume continuity getter",&
            res      = size(ivals),&
            expected = 0,&
            msg      = "volume unset cont",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0136


    subroutine forcad_nurbs_volume_0137(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid-direction multiplicity getter",&
            res      = size(ivals),&
            expected = 0,&
            msg      = "volume bad mult dir",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0137


    subroutine forcad_nurbs_volume_0138(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid-direction control-count getter",&
            res      = vempty%get_nc(4),&
            expected = 0,&
            msg      = "volume bad nc dir",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0138


    subroutine forcad_nurbs_volume_0139(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid-direction degree getter",&
            res      = vempty%get_degree(4),&
            expected = 0,&
            msg      = "volume bad deg dir",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0139


    subroutine forcad_nurbs_volume_0140(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume control-net getter",&
            res      = size(rmat),&
            expected = 0,&
            msg      = "volume unset Xc",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0140


    subroutine forcad_nurbs_volume_0141(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume control-point getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "volume unset Xci",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0141


    subroutine forcad_nurbs_volume_0142(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume control-coordinate getter",&
            res      = vempty%get_Xc(1, 1),&
            expected = 0.0_rk,&
            tol      = 1.0e-12_rk,&
            msg      = "volume unset Xcid",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0142


    subroutine forcad_nurbs_volume_0143(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume geometry getter",&
            res      = size(rmat),&
            expected = 0,&
            msg      = "volume unset Xg",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0143


    subroutine forcad_nurbs_volume_0144(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume geometry-point getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "volume unset Xgi",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0144


    subroutine forcad_nurbs_volume_0145(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume geometry-coordinate getter",&
            res      = vempty%get_Xg(1, 1),&
            expected = 0.0_rk,&
            tol      = 1.0e-12_rk,&
            msg      = "volume unset Xgid",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0145


    subroutine forcad_nurbs_volume_0146(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume weight-vector getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "volume unset Wc",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0146


    subroutine forcad_nurbs_volume_0147(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume weight getter",&
            res      = vempty%get_Wc(1),&
            expected = 0.0_rk,&
            tol      = 1.0e-12_rk,&
            msg      = "volume unset Wci",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0147


    subroutine forcad_nurbs_volume_0148(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume parameter-vector getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "volume unset Xt",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0148


    subroutine forcad_nurbs_volume_0149(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid-direction parameter getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "volume bad Xt dir",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0149


    subroutine forcad_nurbs_volume_0150(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        rvals = vempty%get_knot(1)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume knot-vector getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "volume unset knot",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0150


    subroutine forcad_nurbs_volume_0151(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        rvals = vempty%get_knot(1)
        rvals = vempty%get_knot(4)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid-direction knot-vector getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "volume bad knot dir",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0151


    subroutine forcad_nurbs_volume_0152(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        rvals = vempty%get_knot(1)
        rvals = vempty%get_knot(4)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume knot getter",&
            res      = vempty%get_knot(1, 1),&
            expected = 0.0_rk,&
            tol      = 1.0e-12_rk,&
            msg      = "volume unset knoti",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0152


    subroutine forcad_nurbs_volume_0153(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        rvals = vempty%get_knot(1)
        rvals = vempty%get_knot(4)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid-direction knot getter",&
            res      = vempty%get_knot(4, 1),&
            expected = 0.0_rk,&
            tol      = 1.0e-12_rk,&
            msg      = "volume bad knoti dir",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0153


    subroutine forcad_nurbs_volume_0154(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        rvals = vempty%get_knot(1)
        rvals = vempty%get_knot(4)
        rvals = vempty%cmp_Xg([0.0_rk, 0.0_rk, 0.0_rk])
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume geometry computation",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "volume unset cmp_Xg",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0154


    subroutine forcad_nurbs_volume_0155(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        rvals = vempty%get_knot(1)
        rvals = vempty%get_knot(4)
        rvals = vempty%cmp_Xg([0.0_rk, 0.0_rk, 0.0_rk])
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume geometry-grid shape",&
            res      = vempty%get_ng(),&
            expected = [0, 0, 0],&
            msg      = "volume unset ng",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0155


    subroutine forcad_nurbs_volume_0156(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        rvals = vempty%get_knot(1)
        rvals = vempty%get_knot(4)
        rvals = vempty%cmp_Xg([0.0_rk, 0.0_rk, 0.0_rk])
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume control-net shape",&
            res      = vempty%get_nc(),&
            expected = [0, 0, 0],&
            msg      = "volume unset nc",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0156


    subroutine forcad_nurbs_volume_0157(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        rvals = vempty%get_knot(1)
        rvals = vempty%get_knot(4)
        rvals = vempty%cmp_Xg([0.0_rk, 0.0_rk, 0.0_rk])
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume degree vector",&
            res      = vempty%get_degree(),&
            expected = [0, 0, 0],&
            msg      = "volume unset degree",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0157


    subroutine forcad_nurbs_volume_0158(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        rvals = vempty%get_knot(1)
        rvals = vempty%get_knot(4)
        rvals = vempty%cmp_Xg([0.0_rk, 0.0_rk, 0.0_rk])
        ivals = vempty%cmp_elemFace(1, 7)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid-face element connectivity",&
            res      = size(ivals),&
            expected = 0,&
            msg      = "volume bad elem face",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0158


    subroutine forcad_nurbs_volume_0159(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        rvals = vempty%get_knot(1)
        rvals = vempty%get_knot(4)
        rvals = vempty%cmp_Xg([0.0_rk, 0.0_rk, 0.0_rk])
        ivals = vempty%cmp_elemFace(1, 7)
        ivals = vempty%cmp_elemFace_Xc_vis(1, 7)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid-face control connectivity",&
            res      = size(ivals),&
            expected = 0,&
            msg      = "volume bad Xc face",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0159


    subroutine forcad_nurbs_volume_0160(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        rvals = vempty%get_knot(1)
        rvals = vempty%get_knot(4)
        rvals = vempty%cmp_Xg([0.0_rk, 0.0_rk, 0.0_rk])
        ivals = vempty%cmp_elemFace(1, 7)
        ivals = vempty%cmp_elemFace_Xc_vis(1, 7)
        ivals = vempty%cmp_elemFace_Xg_vis(1, 7)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid-face geometry connectivity",&
            res      = size(ivals),&
            expected = 0,&
            msg      = "volume bad Xg face",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0160


    subroutine forcad_nurbs_volume_0161(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = vempty%get_multiplicity(1)
        ivals = vempty%get_continuity(2)
        ivals = vempty%get_multiplicity(4)

        rmat = vempty%get_Xc()
        rvals = vempty%get_Xc(1)
        rmat = vempty%get_Xg()
        rvals = vempty%get_Xg(1)
        rvals = vempty%get_Wc()
        rvals = vempty%get_Xt(1)
        rvals = vempty%get_Xt(4)
        rvals = vempty%get_knot(1)
        rvals = vempty%get_knot(4)
        rvals = vempty%cmp_Xg([0.0_rk, 0.0_rk, 0.0_rk])
        ivals = vempty%cmp_elemFace(1, 7)
        ivals = vempty%cmp_elemFace_Xc_vis(1, 7)
        ivals = vempty%cmp_elemFace_Xg_vis(1, 7)
        call vempty%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid-face degree",&
            res      = vempty%cmp_degreeFace(7),&
            expected = [0, 0, 0],&
            msg      = "volume bad face degree",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0161


    subroutine forcad_nurbs_volume_0162(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_derivative_order_221(3) = [2,2,1]
        type(nurbs_volume) :: rational_order
        type(nurbs_volume) :: bspline_order
        real(rk) :: Xc_order(8,3), Wc_order(8), expected_order(8), magnitude
        real(rk), allocatable :: Dgc_order(:)

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk,4.0_rk,8.0_rk,12.0_rk,24.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 192.0_rk/((1.25_rk**3)*(1.8_rk**3)*(1.6_rk**2))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude,magnitude,-magnitude,-magnitude,magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order)
        call rational_order%err%print()
        call bspline_order%err%print()

        call ut%test(ti)%check(&
            name     = "rational volume mixed derivative",&
            res      = Dgc_order,&
            expected = expected_order,&
            tol      = 1024.0_rk*epsilon(1.0_rk)*magnitude,&
            msg      = "volume rational mixed deriv",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0162


    subroutine forcad_nurbs_volume_0163(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_derivative_order_221(3) = [2,2,1]
        type(nurbs_volume) :: rational_order
        type(nurbs_volume) :: bspline_order
        real(rk) :: Xc_order(8,3), Wc_order(8), expected_order(8), magnitude
        real(rk), allocatable :: Dgc_order(:)
        real(rk), allocatable :: Dgc_order_vector(:,:)

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk,4.0_rk,8.0_rk,12.0_rk,24.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 192.0_rk/((1.25_rk**3)*(1.8_rk**3)*(1.6_rk**2))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude,magnitude,-magnitude,-magnitude,magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk],&
            Xt3   = [0.2_rk])
        call rational_order%err%print()
        call bspline_order%err%print()

        call ut%test(ti)%check(&
            name     = "rational volume vector derivative",&
            res      = Dgc_order_vector(1,:),&
            expected = expected_order,&
            tol      = 1024.0_rk*epsilon(1.0_rk)*magnitude,&
            msg      = "volume rational vector deriv",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0163


    subroutine forcad_nurbs_volume_0164(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_derivative_order_221(3) = [2,2,1]
        type(nurbs_volume) :: rational_order
        type(nurbs_volume) :: bspline_order
        real(rk) :: Xc_order(8,3), Wc_order(8), expected_order(8), magnitude
        real(rk), allocatable :: Dgc_order(:)
        real(rk), allocatable :: Dgc_order_vector(:,:)
        real(rk), allocatable :: Dgc_active(:)
        integer :: first_active(3)

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk,4.0_rk,8.0_rk,12.0_rk,24.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 192.0_rk/((1.25_rk**3)*(1.8_rk**3)*(1.6_rk**2))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude,magnitude,-magnitude,-magnitude,magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk],&
            Xt3   = [0.2_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk,0.2_rk],&
            order        = volume_derivative_order_221,&
            first_active = first_active,&
            Dgc          = Dgc_active)
        call rational_order%err%print()
        call bspline_order%err%print()

        call ut%test(ti)%check(&
            name     = "active rational volume mixed derivative",&
            res      = Dgc_active,&
            expected = expected_order,&
            tol      = 1024.0_rk*epsilon(1.0_rk)*magnitude,&
            msg      = "volume active rational deriv",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0164


    subroutine forcad_nurbs_volume_0165(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_derivative_order_221(3) = [2,2,1]
        type(nurbs_volume) :: rational_order
        type(nurbs_volume) :: bspline_order
        real(rk) :: Xc_order(8,3), Wc_order(8), expected_order(8), magnitude
        real(rk), allocatable :: Dgc_order(:)
        real(rk), allocatable :: Dgc_order_vector(:,:)
        real(rk), allocatable :: Dgc_active(:)
        integer :: first_active(3)

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk,4.0_rk,8.0_rk,12.0_rk,24.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 192.0_rk/((1.25_rk**3)*(1.8_rk**3)*(1.6_rk**2))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude,magnitude,-magnitude,-magnitude,magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk],&
            Xt3   = [0.2_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk,0.2_rk],&
            order        = volume_derivative_order_221,&
            first_active = first_active,&
            Dgc          = Dgc_active)
        call rational_order%err%print()
        call bspline_order%err%print()

        call ut%test(ti)%check(&
            name     = "active volume first basis indices",&
            res      = first_active,&
            expected = [1,1,1],&
            msg      = "volume active first indices",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0165


    subroutine forcad_nurbs_volume_0166(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_derivative_order_221(3) = [2,2,1]
        type(nurbs_volume) :: rational_order
        type(nurbs_volume) :: bspline_order
        real(rk) :: Xc_order(8,3), Wc_order(8), expected_order(8), magnitude
        real(rk), allocatable :: Dgc_order(:)
        real(rk), allocatable :: Dgc_order_vector(:,:)
        real(rk), allocatable :: Dgc_active(:)
        real(rk), allocatable :: Dgc_active_vector(:,:)
        integer, allocatable :: first_vector(:,:)
        integer :: first_active(3)

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk,4.0_rk,8.0_rk,12.0_rk,24.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 192.0_rk/((1.25_rk**3)*(1.8_rk**3)*(1.6_rk**2))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude,magnitude,-magnitude,-magnitude,magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk],&
            Xt3   = [0.2_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk,0.2_rk],&
            order        = volume_derivative_order_221,&
            first_active = first_active,&
            Dgc          = Dgc_active)

        call rational_order%derivative_order_active(&
            Xt1          = [0.25_rk,0.5_rk],&
            Xt2          = [0.4_rk],&
            Xt3          = [0.2_rk],&
            order        = volume_derivative_order_221,&
            first_active = first_vector,&
            Dgc          = Dgc_active_vector)
        call rational_order%err%print()
        call bspline_order%err%print()

        call ut%test(ti)%check(&
            name     = "active volume vector derivative ordering",&
            res      = Dgc_active_vector(:,1),&
            expected = expected_order,&
            tol      = 1024.0_rk*epsilon(1.0_rk)*magnitude,&
            msg      = "volume active vector ordering",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0166


    subroutine forcad_nurbs_volume_0167(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_derivative_order_221(3) = [2,2,1]
        type(nurbs_volume) :: rational_order
        type(nurbs_volume) :: bspline_order
        real(rk) :: Xc_order(8,3), Wc_order(8), expected_order(8), magnitude
        real(rk), allocatable :: Dgc_order(:), Dgc_order_vector(:,:), Dgc_bspline(:)
        real(rk), allocatable :: Dgc_active(:)
        real(rk), allocatable :: Dgc_active_vector(:,:)
        integer, allocatable :: first_vector(:,:)
        integer :: first_active(3)

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk,4.0_rk,8.0_rk,12.0_rk,24.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 192.0_rk/((1.25_rk**3)*(1.8_rk**3)*(1.6_rk**2))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude,magnitude,-magnitude,-magnitude,magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk],&
            Xt3   = [0.2_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk,0.2_rk],&
            order        = volume_derivative_order_221,&
            first_active = first_active,&
            Dgc          = Dgc_active)

        call rational_order%derivative_order_active(&
            Xt1          = [0.25_rk,0.5_rk],&
            Xt2          = [0.4_rk],&
            Xt3          = [0.2_rk],&
            order        = volume_derivative_order_221,&
            first_active = first_vector,&
            Dgc          = Dgc_active_vector)

        call bspline_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_bspline)
        call rational_order%err%print()
        call bspline_order%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline volume derivative above degree",&
            res      = Dgc_bspline,&
            expected = [0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk],&
            tol      = 32.0_rk*epsilon(1.0_rk),&
            msg      = "volume bspline high deriv",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0167


    subroutine forcad_nurbs_volume_0168(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_derivative_order_101(3) = [1,0,1]
        integer, parameter :: volume_derivative_order_221(3) = [2,2,1]
        type(nurbs_volume) :: rational_order
        type(nurbs_volume) :: bspline_order
        type(nurbs_volume) :: active_volume
        real(rk) :: Xc_order(8,3), Wc_order(8), expected_order(8), magnitude
        real(rk) :: Xc_active(12,3), Wc_active(12)
        real(rk), allocatable :: Dgc_order(:), Dgc_order_vector(:,:), Dgc_bspline(:)
        real(rk), allocatable :: Dgc_active(:)
        real(rk), allocatable :: Dgc_active_vector(:,:)
        real(rk), allocatable :: Dgc_dense_active(:)
        real(rk), allocatable :: Dgc_reconstructed(:)
        integer, allocatable :: first_vector(:,:)
        integer :: first_active(3), l1, l2, l3, i1, i2, i3

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk,4.0_rk,8.0_rk,12.0_rk,24.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 192.0_rk/((1.25_rk**3)*(1.8_rk**3)*(1.6_rk**2))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude,magnitude,-magnitude,-magnitude,magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk],&
            Xt3   = [0.2_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk,0.2_rk],&
            order        = volume_derivative_order_221,&
            first_active = first_active,&
            Dgc          = Dgc_active)

        call rational_order%derivative_order_active(&
            Xt1          = [0.25_rk,0.5_rk],&
            Xt2          = [0.4_rk],&
            Xt3          = [0.2_rk],&
            order        = volume_derivative_order_221,&
            first_active = first_vector,&
            Dgc          = Dgc_active_vector)

        call bspline_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_bspline)

        Xc_active = 0.0_rk
        Wc_active = [1.0_rk,1.1_rk,0.9_rk,1.3_rk,0.8_rk,1.2_rk,0.7_rk,1.4_rk,1.05_rk,0.95_rk,1.25_rk,0.85_rk]
        call active_volume%set(&
            knot1 = [0.0_rk,0.0_rk,0.35_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_active,&
            Wc    = Wc_active)
        call active_volume%derivative_order(&
            Xt    = [0.8_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_101,&
            Dgc   = Dgc_dense_active)
        call active_volume%derivative_order_active(&
            Xt           = [0.8_rk,0.4_rk,0.2_rk],&
            order        = volume_derivative_order_101,&
            first_active = first_active,&
            Dgc          = Dgc_active)
        allocate(Dgc_reconstructed(size(Dgc_dense_active)), source=0.0_rk)
        do l3 = 0, 1
            do l2 = 0, 1
                do l1 = 0, 1
                    i1 = first_active(1) + l1
                    i2 = first_active(2) + l2
                    i3 = first_active(3) + l3
                    Dgc_reconstructed(i1+3*((i2-1)+2*(i3-1))) = Dgc_active(l1+1+2*(l2+2*l3))
                end do
            end do
        end do
        call rational_order%err%print()
        call bspline_order%err%print()
        call active_volume%err%print()

        call ut%test(ti)%check(&
            name     = "active volume derivative reconstruction",&
            res      = Dgc_reconstructed,&
            expected = Dgc_dense_active,&
            tol      = 512.0_rk*epsilon(1.0_rk)*max(1.0_rk,maxval(abs(Dgc_dense_active))),&
            msg      = "volume active dense reconstruction",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0168


    subroutine forcad_nurbs_volume_0169(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_derivative_order_101(3) = [1,0,1]
        integer, parameter :: volume_derivative_order_221(3) = [2,2,1]
        type(nurbs_volume) :: rational_order
        type(nurbs_volume) :: bspline_order
        type(nurbs_volume) :: active_volume
        real(rk) :: Xc_order(8,3), Wc_order(8), expected_order(8), magnitude
        real(rk) :: Xc_active(12,3), Wc_active(12)
        real(rk), allocatable :: Dgc_order(:), Dgc_order_vector(:,:), Dgc_bspline(:)
        real(rk), allocatable :: Dgc_active(:)
        real(rk), allocatable :: Dgc_active_vector(:,:)
        real(rk), allocatable :: Dgc_dense_active(:)
        real(rk), allocatable :: Dgc_reconstructed(:)
        integer, allocatable :: first_vector(:,:)
        integer :: first_active(3), l1, l2, l3, i1, i2, i3

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk,4.0_rk,8.0_rk,12.0_rk,24.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 192.0_rk/((1.25_rk**3)*(1.8_rk**3)*(1.6_rk**2))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude,magnitude,-magnitude,-magnitude,magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk],&
            Xt3   = [0.2_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk,0.2_rk],&
            order        = volume_derivative_order_221,&
            first_active = first_active,&
            Dgc          = Dgc_active)

        call rational_order%derivative_order_active(&
            Xt1          = [0.25_rk,0.5_rk],&
            Xt2          = [0.4_rk],&
            Xt3          = [0.2_rk],&
            order        = volume_derivative_order_221,&
            first_active = first_vector,&
            Dgc          = Dgc_active_vector)

        call bspline_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_bspline)

        Xc_active = 0.0_rk
        Wc_active = [1.0_rk,1.1_rk,0.9_rk,1.3_rk,0.8_rk,1.2_rk,0.7_rk,1.4_rk,1.05_rk,0.95_rk,1.25_rk,0.85_rk]
        call active_volume%set(&
            knot1 = [0.0_rk,0.0_rk,0.35_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_active,&
            Wc    = Wc_active)
        call active_volume%derivative_order(&
            Xt    = [0.8_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_101,&
            Dgc   = Dgc_dense_active)
        call active_volume%derivative_order_active(&
            Xt           = [0.8_rk,0.4_rk,0.2_rk],&
            order        = volume_derivative_order_101,&
            first_active = first_active,&
            Dgc          = Dgc_active)
        allocate(Dgc_reconstructed(size(Dgc_dense_active)), source=0.0_rk)
        do l3 = 0, 1
            do l2 = 0, 1
                do l1 = 0, 1
                    i1 = first_active(1) + l1
                    i2 = first_active(2) + l2
                    i3 = first_active(3) + l3
                    Dgc_reconstructed(i1+3*((i2-1)+2*(i3-1))) = Dgc_active(l1+1+2*(l2+2*l3))
                end do
            end do
        end do
        call rational_order%err%print()
        call bspline_order%err%print()
        call active_volume%err%print()

        call ut%test(ti)%check(&
            name     = "active volume derivative storage reduction",&
            res      = size(Dgc_active) < size(Dgc_dense_active),&
            expected = .true.,&
            msg      = "volume active lower memory",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0169


    subroutine forcad_nurbs_volume_0170(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_derivative_order_000(3) = [0,0,0]
        integer, parameter :: volume_derivative_order_101(3) = [1,0,1]
        integer, parameter :: volume_derivative_order_221(3) = [2,2,1]
        type(nurbs_volume) :: rational_order, bspline_order, active_volume, empty_active
        real(rk) :: Xc_order(8,3), Wc_order(8), expected_order(8), magnitude
        real(rk) :: Xc_active(12,3), Wc_active(12)
        real(rk), allocatable :: Dgc_order(:), Dgc_order_vector(:,:), Dgc_bspline(:)
        real(rk), allocatable :: Dgc_active(:), Dgc_active_vector(:,:), Dgc_dense_active(:), Dgc_reconstructed(:), Dgc_empty(:)
        integer, allocatable :: first_vector(:,:)
        integer :: first_active(3), l1, l2, l3, i1, i2, i3

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk,4.0_rk,8.0_rk,12.0_rk,24.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 192.0_rk/((1.25_rk**3)*(1.8_rk**3)*(1.6_rk**2))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude,magnitude,-magnitude,-magnitude,magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk],&
            Xt3   = [0.2_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk,0.2_rk],&
            order        = volume_derivative_order_221,&
            first_active = first_active,&
            Dgc          = Dgc_active)

        call rational_order%derivative_order_active(&
            Xt1          = [0.25_rk,0.5_rk],&
            Xt2          = [0.4_rk],&
            Xt3          = [0.2_rk],&
            order        = volume_derivative_order_221,&
            first_active = first_vector,&
            Dgc          = Dgc_active_vector)

        call bspline_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_221,&
            Dgc   = Dgc_bspline)

        Xc_active = 0.0_rk
        Wc_active = [1.0_rk,1.1_rk,0.9_rk,1.3_rk,0.8_rk,1.2_rk,0.7_rk,1.4_rk,1.05_rk,0.95_rk,1.25_rk,0.85_rk]
        call active_volume%set(&
            knot1 = [0.0_rk,0.0_rk,0.35_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot3 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_active,&
            Wc    = Wc_active)
        call active_volume%derivative_order(&
            Xt    = [0.8_rk,0.4_rk,0.2_rk],&
            order = volume_derivative_order_101,&
            Dgc   = Dgc_dense_active)
        call active_volume%derivative_order_active(&
            Xt           = [0.8_rk,0.4_rk,0.2_rk],&
            order        = volume_derivative_order_101,&
            first_active = first_active,&
            Dgc          = Dgc_active)
        allocate(Dgc_reconstructed(size(Dgc_dense_active)), source=0.0_rk)
        do l3 = 0, 1
            do l2 = 0, 1
                do l1 = 0, 1
                    i1 = first_active(1) + l1
                    i2 = first_active(2) + l2
                    i3 = first_active(3) + l3
                    Dgc_reconstructed(i1+3*((i2-1)+2*(i3-1))) = Dgc_active(l1+1+2*(l2+2*l3))
                end do
            end do
        end do

        call empty_active%derivative_order_active(&
            Xt           = [0.5_rk,0.5_rk,0.5_rk],&
            order        = volume_derivative_order_000,&
            first_active = first_active,&
            Dgc          = Dgc_empty)
        call rational_order%err%print()
        call bspline_order%err%print()
        call active_volume%err%print()
        call empty_active%err%print()

        call ut%test(ti)%check(&
            name     = "unset volume active derivative",&
            res      = size(Dgc_empty),&
            expected = 0,&
            msg      = "volume unset active derivative",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0170


    subroutine forcad_nurbs_volume_0171(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian diagnostic",&
            res      = volume_hessian%err%ok,&
            expected = .true.,&
            msg      = "Ansatz Hessian diagnostic is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0171


    subroutine forcad_nurbs_volume_0172(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian xx scaling",&
            res      = d2T_dX2(1,1,1),&
            expected = 1.0_rk/32.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "The physical xx derivative must include x-coordinate scaling",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0172


    subroutine forcad_nurbs_volume_0173(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian xy scaling",&
            res      = d2T_dX2(1,1,2),&
            expected = 1.0_rk/24.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "The physical xy derivative must include both coordinate scales",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0173


    subroutine forcad_nurbs_volume_0174(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian xz scaling",&
            res      = d2T_dX2(1,1,3),&
            expected = 1.0_rk/32.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "The physical xz derivative must include both coordinate scales",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0174


    subroutine forcad_nurbs_volume_0175(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian yy scaling",&
            res      = d2T_dX2(1,2,2),&
            expected = 1.0_rk/72.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "The physical yy derivative must include y-coordinate scaling",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0175


    subroutine forcad_nurbs_volume_0176(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian yz scaling",&
            res      = d2T_dX2(1,2,3),&
            expected = 1.0_rk/48.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "The physical yz derivative must include both coordinate scales",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0176


    subroutine forcad_nurbs_volume_0177(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian zz scaling",&
            res      = d2T_dX2(1,3,3),&
            expected = 1.0_rk/128.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "The physical zz derivative must include z-coordinate scaling",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0177


    subroutine forcad_nurbs_volume_0178(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian partition",&
            res      = maxval(abs(sum(d2T_dX2,dim=1))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Physical volume second derivatives must sum to zero",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0178


    subroutine forcad_nurbs_volume_0179(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian symmetry",&
            res      = maxval(abs(d2T_dX2(:,1,2)-d2T_dX2(:,2,1))) + maxval(abs(d2T_dX2(:,1,3)-d2T_dX2(:,3,1))) + &
                maxval(abs(d2T_dX2(:,2,3)-d2T_dX2(:,3,2))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Mixed physical volume derivatives must be symmetric",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0179


    subroutine forcad_nurbs_volume_0180(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz positive measure",&
            res      = dV > 0.0_rk,&
            expected = .true.,&
            msg      = "A regular volume element must have a positive quadrature measure",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0180


    subroutine forcad_nurbs_volume_0181(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)

        call volume_hessian%finalize()
        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            Wc     = hessian_Wc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "rational ansatz partition",&
            res      = sum(T),&
            expected = 1.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Rational volume shape functions must form a partition of unity",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0181


    subroutine forcad_nurbs_volume_0182(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)

        call volume_hessian%finalize()
        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            Wc     = hessian_Wc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "rational ansatz gradient partition",&
            res      = maxval(abs(sum(dT_dX,dim=1))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Rational physical volume gradients must sum to zero",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0182


    subroutine forcad_nurbs_volume_0183(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)

        call volume_hessian%finalize()
        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            Wc     = hessian_Wc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "rational ansatz Hessian partition",&
            res      = maxval(abs(sum(d2T_dX2,dim=1))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Rational physical volume second derivatives must sum to zero",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0183


    subroutine forcad_nurbs_volume_0184(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)

        call volume_hessian%finalize()
        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            Wc     = hessian_Wc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "rational ansatz Hessian symmetry",&
            res      = maxval(abs(d2T_dX2(:,1,2)-d2T_dX2(:,2,1))) + maxval(abs(d2T_dX2(:,1,3)-d2T_dX2(:,3,1))) + &
                maxval(abs(d2T_dX2(:,2,3)-d2T_dX2(:,3,2))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Rational mixed physical volume derivatives must be symmetric",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0184


    subroutine forcad_nurbs_volume_0185(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6)
        real(rk) :: hessian_Xc(27,3)
        real(rk) :: hessian_Wc(27)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)

        call volume_hessian%finalize()
        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            Wc     = hessian_Wc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "rational ansatz chain rule",&
            res      = maxval(abs([ sum(d2T_dX2(:,1,1)*hessian_Xc(:,1)), sum(d2T_dX2(:,1,2)*hessian_Xc(:,1)),&
                sum(d2T_dX2(:,1,3)*hessian_Xc(:,1)), sum(d2T_dX2(:,2,2)*hessian_Xc(:,1)), sum(d2T_dX2(:,2,3)*hessian_Xc(:,1)),&
                sum(d2T_dX2(:,3,3)*hessian_Xc(:,1)), sum(d2T_dX2(:,1,1)*hessian_Xc(:,2)), sum(d2T_dX2(:,1,2)*hessian_Xc(:,2)),&
                sum(d2T_dX2(:,1,3)*hessian_Xc(:,2)), sum(d2T_dX2(:,2,2)*hessian_Xc(:,2)), sum(d2T_dX2(:,2,3)*hessian_Xc(:,2)),&
                sum(d2T_dX2(:,3,3)*hessian_Xc(:,2)), sum(d2T_dX2(:,1,1)*hessian_Xc(:,3)), sum(d2T_dX2(:,1,2)*hessian_Xc(:,3)),&
                sum(d2T_dX2(:,1,3)*hessian_Xc(:,3)), sum(d2T_dX2(:,2,2)*hessian_Xc(:,3)), sum(d2T_dX2(:,2,3)*hessian_Xc(:,3)),&
                sum(d2T_dX2(:,3,3)*hessian_Xc(:,3))])),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Rational ansatz chain rule is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0185


    subroutine forcad_nurbs_volume_0186(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6), hessian_Xc(27,3), hessian_Wc(27), dV, u, v, w
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)

        call volume_hessian%finalize()
        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            Wc     = hessian_Wc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)

        do kk = 1, 3
            w = 0.5_rk*real(kk-1,rk)
            do jj = 1, 3
                v = 0.5_rk*real(jj-1,rk)
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    u = 0.5_rk*real(ii-1,rk)
                    hessian_Xc(point,1) = u
                    if (ii == 3) hessian_Xc(point,1) = hessian_Xc(point,1) + 0.25_rk
                    hessian_Xc(point,2) = v + 0.2_rk*u*v
                    hessian_Xc(point,3) = w + 0.15_rk*v*w
                end do
            end do
        end do
        call volume_hessian%finalize()
        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "curved ansatz chain rule",&
            res      = maxval(abs([ sum(d2T_dX2(:,1,1)*hessian_Xc(:,1)), sum(d2T_dX2(:,1,2)*hessian_Xc(:,1)),&
                sum(d2T_dX2(:,1,3)*hessian_Xc(:,1)), sum(d2T_dX2(:,2,2)*hessian_Xc(:,1)), sum(d2T_dX2(:,2,3)*hessian_Xc(:,1)),&
                sum(d2T_dX2(:,3,3)*hessian_Xc(:,1)), sum(d2T_dX2(:,1,1)*hessian_Xc(:,2)), sum(d2T_dX2(:,1,2)*hessian_Xc(:,2)),&
                sum(d2T_dX2(:,1,3)*hessian_Xc(:,2)), sum(d2T_dX2(:,2,2)*hessian_Xc(:,2)), sum(d2T_dX2(:,2,3)*hessian_Xc(:,2)),&
                sum(d2T_dX2(:,3,3)*hessian_Xc(:,2)), sum(d2T_dX2(:,1,1)*hessian_Xc(:,3)), sum(d2T_dX2(:,1,2)*hessian_Xc(:,3)),&
                sum(d2T_dX2(:,1,3)*hessian_Xc(:,3)), sum(d2T_dX2(:,2,2)*hessian_Xc(:,3)), sum(d2T_dX2(:,2,3)*hessian_Xc(:,3)),&
                sum(d2T_dX2(:,3,3)*hessian_Xc(:,3))])),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "The physical Hessian must reproduce a curved coordinate map",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0186


    subroutine forcad_nurbs_volume_0187(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_101(3) = [1,0,1]
        integer, parameter :: volume_quadrature_order_111(3) = [1,1,1]
        type(nurbs_volume) :: volume_hessian
        real(rk) :: hessian_knot(6), hessian_Xc(27,3), hessian_Wc(27), dV, u, v, w
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, kk, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do kk = 1, 3
            do jj = 1, 3
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    hessian_Xc(point,:) = [&
                        real(ii-1,rk),&
                        1.5_rk*real(jj-1,rk),&
                        2.0_rk*real(kk-1,rk)]
                    hessian_Wc(point) = 1.0_rk + 0.01_rk*real(ii-1,rk) + &
                        0.02_rk*real(jj-1,rk) + 0.03_rk*real(kk-1,rk)
                end do
            end do
        end do

        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)

        call volume_hessian%finalize()
        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            Wc     = hessian_Wc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)

        do kk = 1, 3
            w = 0.5_rk*real(kk-1,rk)
            do jj = 1, 3
                v = 0.5_rk*real(jj-1,rk)
                do ii = 1, 3
                    point = (kk-1)*9 + (jj-1)*3 + ii
                    u = 0.5_rk*real(ii-1,rk)
                    hessian_Xc(point,1) = u
                    if (ii == 3) hessian_Xc(point,1) = hessian_Xc(point,1) + 0.25_rk
                    hessian_Xc(point,2) = v + 0.2_rk*u*v
                    hessian_Xc(point,3) = w + 0.15_rk*v*w
                end do
            end do
        end do
        call volume_hessian%finalize()
        call volume_hessian%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            knot3  = hessian_knot,&
            Xc     = hessian_Xc,&
            degree = [2,2,2])
        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_111)

        call volume_hessian%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dV          = dV,&
            ngauss      = volume_quadrature_order_101)
        call volume_hessian%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz invalid quadrature",&
            res      = .not. volume_hessian%err%ok .and. size(T) == 0 .and. size(dT_dX) == 0 .and. size(d2T_dX2) == 0 .and. &
                abs(dV) <= tiny(1.0_rk),&
            expected = .true.,&
            msg      = "Invalid quadrature must return empty outputs and a diagnostic",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0187


    subroutine forcad_nurbs_volume_0188(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: family_volume
        real(rk) :: family_knot1(7), family_knot2(7), family_knot3(7), family_Xc(64,3)
        real(rk), allocatable :: Xt1(:)
        real(rk), allocatable :: Xt2(:)
        real(rk), allocatable :: Xt3(:)
        real(rk), allocatable :: Xt_grid(:,:)
        real(rk), allocatable :: Xg_reference(:,:)
        integer :: ii
        integer :: jj
        integer :: kk
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        family_knot2 = [0.0_rk, 0.4_rk, 1.1_rk, 2.0_rk, 3.2_rk, 4.0_rk, 4.8_rk]
        family_knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.4_rk, 2.2_rk, 3.0_rk]
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    point = (kk-1)*16 + (jj-1)*4 + ii
                    family_Xc(point,:) = [real(ii-1,rk),real(jj-1,rk),real(kk-1,rk)]
                end do
            end do
        end do

        call family_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call family_volume%create(res1=5, res2=5, res3=5)
        Xt1 = family_volume%get_Xt(1)
        Xt2 = family_volume%get_Xt(2)
        Xt3 = family_volume%get_Xt(3)
        call ndgrid(Xt1, Xt2, Xt3, Xt_grid)
        Xg_reference = family_volume%get_Xg()
        call family_volume%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families finite geometry",&
            res      = all(ieee_is_finite(Xg_reference)),&
            expected = .true.,&
            msg      = "Mixed knot families finite geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0188


    subroutine forcad_nurbs_volume_0189(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: family_volume
        real(rk) :: family_knot1(7), family_knot2(7), family_knot3(7), family_Xc(64,3)
        real(rk), allocatable :: Xt1(:)
        real(rk), allocatable :: Xt2(:)
        real(rk), allocatable :: Xt3(:)
        real(rk), allocatable :: Xt_grid(:,:)
        real(rk), allocatable :: Xg_reference(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: kk
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        family_knot2 = [0.0_rk, 0.4_rk, 1.1_rk, 2.0_rk, 3.2_rk, 4.0_rk, 4.8_rk]
        family_knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.4_rk, 2.2_rk, 3.0_rk]
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    point = (kk-1)*16 + (jj-1)*4 + ii
                    family_Xc(point,:) = [real(ii-1,rk),real(jj-1,rk),real(kk-1,rk)]
                end do
            end do
        end do

        call family_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call family_volume%create(res1=5, res2=5, res3=5)
        Xt1 = family_volume%get_Xt(1)
        Xt2 = family_volume%get_Xt(2)
        Xt3 = family_volume%get_Xt(3)
        call ndgrid(Xt1, Xt2, Xt3, Xt_grid)
        Xg_reference = family_volume%get_Xg()
        elem = family_volume%get_elem_Xg_vis()
        call family_volume%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families create connectivity",&
            res      = size(elem,1),&
            expected = 64,&
            msg      = "Mixed knot families create connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0189


    subroutine forcad_nurbs_volume_0190(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: family_volume
        real(rk) :: family_knot1(7), family_knot2(7), family_knot3(7), family_Xc(64,3)
        real(rk), allocatable :: Xt1(:)
        real(rk), allocatable :: Xt2(:)
        real(rk), allocatable :: Xt3(:)
        real(rk), allocatable :: Xt_grid(:,:)
        real(rk), allocatable :: Xg_reference(:,:)
        real(rk), allocatable :: Tgrid(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: kk
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        family_knot2 = [0.0_rk, 0.4_rk, 1.1_rk, 2.0_rk, 3.2_rk, 4.0_rk, 4.8_rk]
        family_knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.4_rk, 2.2_rk, 3.0_rk]
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    point = (kk-1)*16 + (jj-1)*4 + ii
                    family_Xc(point,:) = [real(ii-1,rk),real(jj-1,rk),real(kk-1,rk)]
                end do
            end do
        end do

        call family_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call family_volume%create(res1=5, res2=5, res3=5)
        Xt1 = family_volume%get_Xt(1)
        Xt2 = family_volume%get_Xt(2)
        Xt3 = family_volume%get_Xt(3)
        call ndgrid(Xt1, Xt2, Xt3, Xt_grid)
        Xg_reference = family_volume%get_Xg()
        elem = family_volume%get_elem_Xg_vis()
        call family_volume%basis(res1=3, res2=3, res3=3, Tgc=Tgrid)
        elem = family_volume%get_elem_Xg_vis()
        call family_volume%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families basis connectivity",&
            res      = size(elem,1),&
            expected = 64,&
            msg      = "Mixed knot families basis connectivity is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0190


    subroutine forcad_nurbs_volume_0191(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: INVARIANT_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: family_volume
        real(rk) :: family_knot1(7), family_knot2(7), family_knot3(7), family_Xc(64,3)
        real(rk) :: xi(3)
        real(rk), allocatable :: Xt1(:)
        real(rk), allocatable :: Xt2(:)
        real(rk), allocatable :: Xt3(:)
        real(rk), allocatable :: Xt_grid(:,:)
        real(rk), allocatable :: Xg_reference(:,:)
        real(rk), allocatable :: active1(:)
        real(rk), allocatable :: active2(:)
        real(rk), allocatable :: active3(:)
        real(rk), allocatable :: T(:)
        real(rk), allocatable :: Tgrid(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: kk
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        family_knot2 = [0.0_rk, 0.4_rk, 1.1_rk, 2.0_rk, 3.2_rk, 4.0_rk, 4.8_rk]
        family_knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.4_rk, 2.2_rk, 3.0_rk]
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    point = (kk-1)*16 + (jj-1)*4 + ii
                    family_Xc(point,:) = [real(ii-1,rk),real(jj-1,rk),real(kk-1,rk)]
                end do
            end do
        end do

        call family_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call family_volume%create(res1=5, res2=5, res3=5)
        Xt1 = family_volume%get_Xt(1)
        Xt2 = family_volume%get_Xt(2)
        Xt3 = family_volume%get_Xt(3)
        call ndgrid(Xt1, Xt2, Xt3, Xt_grid)
        Xg_reference = family_volume%get_Xg()
        elem = family_volume%get_elem_Xg_vis()
        call family_volume%basis(res1=3, res2=3, res3=3, Tgc=Tgrid)
        elem = family_volume%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,4,2)
        active2 = active_knots(family_knot2,4,2)
        active3 = active_knots(family_knot3,4,2)
        xi = [&
            0.5_rk*(active1(1)+active1(2)),&
            0.5_rk*(active2(1)+active2(2)),&
            0.5_rk*(active3(1)+active3(2))]
        call family_volume%basis(xi, T)
        call family_volume%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families partition",&
            res      = sum(T),&
            expected = 1.0_rk,&
            tol      = INVARIANT_TOL,&
            msg      = "Mixed knot families partition is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0191


    subroutine forcad_nurbs_volume_0192(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: family_volume
        real(rk) :: family_knot1(7), family_knot2(7), family_knot3(7), family_Xc(64,3)
        real(rk) :: xi(3)
        real(rk), allocatable :: Xt1(:)
        real(rk), allocatable :: Xt2(:)
        real(rk), allocatable :: Xt3(:)
        real(rk), allocatable :: Xt_grid(:,:)
        real(rk), allocatable :: Xg_reference(:,:)
        real(rk), allocatable :: active1(:)
        real(rk), allocatable :: active2(:)
        real(rk), allocatable :: active3(:)
        real(rk), allocatable :: T(:)
        real(rk), allocatable :: Tgrid(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: kk
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        family_knot2 = [0.0_rk, 0.4_rk, 1.1_rk, 2.0_rk, 3.2_rk, 4.0_rk, 4.8_rk]
        family_knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.4_rk, 2.2_rk, 3.0_rk]
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    point = (kk-1)*16 + (jj-1)*4 + ii
                    family_Xc(point,:) = [real(ii-1,rk),real(jj-1,rk),real(kk-1,rk)]
                end do
            end do
        end do

        call family_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call family_volume%create(res1=5, res2=5, res3=5)
        Xt1 = family_volume%get_Xt(1)
        Xt2 = family_volume%get_Xt(2)
        Xt3 = family_volume%get_Xt(3)
        call ndgrid(Xt1, Xt2, Xt3, Xt_grid)
        Xg_reference = family_volume%get_Xg()
        elem = family_volume%get_elem_Xg_vis()
        call family_volume%basis(res1=3, res2=3, res3=3, Tgc=Tgrid)
        elem = family_volume%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,4,2)
        active2 = active_knots(family_knot2,4,2)
        active3 = active_knots(family_knot3,4,2)
        xi = [&
            0.5_rk*(active1(1)+active1(2)),&
            0.5_rk*(active2(1)+active2(2)),&
            0.5_rk*(active3(1)+active3(2))]
        call family_volume%basis(xi, T)
        elem = family_volume%cmp_elem()
        call family_volume%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families element count",&
            res      = size(elem,1),&
            expected = active_span_count(family_knot1,4,2)*active_span_count(family_knot2,4,2)* &
                active_span_count(family_knot3,4,2),&
            msg      = "Mixed knot families element count is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0192


    subroutine forcad_nurbs_volume_0193(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_333(3) = [3,3,3]
        type(nurbs_volume) :: family_volume
        real(rk) :: family_knot1(7), family_knot2(7), family_knot3(7), family_Xc(64,3)
        real(rk) :: xi(3)
        real(rk) :: dV
        real(rk), allocatable :: Xt1(:)
        real(rk), allocatable :: Xt2(:)
        real(rk), allocatable :: Xt3(:)
        real(rk), allocatable :: Xt_grid(:,:)
        real(rk), allocatable :: Xg_reference(:,:)
        real(rk), allocatable :: active1(:), active2(:), active3(:), T(:), Tgrid(:,:), dT_dX(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: kk
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        family_knot2 = [0.0_rk, 0.4_rk, 1.1_rk, 2.0_rk, 3.2_rk, 4.0_rk, 4.8_rk]
        family_knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.4_rk, 2.2_rk, 3.0_rk]
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    point = (kk-1)*16 + (jj-1)*4 + ii
                    family_Xc(point,:) = [real(ii-1,rk),real(jj-1,rk),real(kk-1,rk)]
                end do
            end do
        end do

        call family_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call family_volume%create(res1=5, res2=5, res3=5)
        Xt1 = family_volume%get_Xt(1)
        Xt2 = family_volume%get_Xt(2)
        Xt3 = family_volume%get_Xt(3)
        call ndgrid(Xt1, Xt2, Xt3, Xt_grid)
        Xg_reference = family_volume%get_Xg()
        elem = family_volume%get_elem_Xg_vis()
        call family_volume%basis(res1=3, res2=3, res3=3, Tgc=Tgrid)
        elem = family_volume%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,4,2)
        active2 = active_knots(family_knot2,4,2)
        active3 = active_knots(family_knot3,4,2)
        xi = [&
            0.5_rk*(active1(1)+active1(2)),&
            0.5_rk*(active2(1)+active2(2)),&
            0.5_rk*(active3(1)+active3(2))]
        call family_volume%basis(xi, T)
        elem = family_volume%cmp_elem()
        call family_volume%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = T,&
            dTgc_dXg = dT_dX,&
            dV       = dV,&
            ngauss   = volume_quadrature_order_333)
        call family_volume%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families ansatz",&
            res      = dV > 0.0_rk,&
            expected = .true.,&
            msg      = "A regular mixed-family volume element must have positive measure",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0193


    subroutine forcad_nurbs_volume_0194(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_333(3) = [3,3,3]
        real(rk), parameter :: REFINEMENT_TOL = 3.0e-10_rk
        type(nurbs_volume) :: family_volume, refined_volume
        real(rk) :: family_knot1(7), family_knot2(7), family_knot3(7), family_Xc(64,3)
        real(rk) :: xi(3)
        real(rk) :: dV
        real(rk), allocatable :: Xt1(:), Xt2(:), Xt3(:), Xt_grid(:,:), Xg_reference(:,:), Xg_family(:,:)
        real(rk), allocatable :: active1(:), active2(:), active3(:), T(:), Tgrid(:,:), dT_dX(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: kk
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        family_knot2 = [0.0_rk, 0.4_rk, 1.1_rk, 2.0_rk, 3.2_rk, 4.0_rk, 4.8_rk]
        family_knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.4_rk, 2.2_rk, 3.0_rk]
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    point = (kk-1)*16 + (jj-1)*4 + ii
                    family_Xc(point,:) = [real(ii-1,rk),real(jj-1,rk),real(kk-1,rk)]
                end do
            end do
        end do

        call family_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call family_volume%create(res1=5, res2=5, res3=5)
        Xt1 = family_volume%get_Xt(1)
        Xt2 = family_volume%get_Xt(2)
        Xt3 = family_volume%get_Xt(3)
        call ndgrid(Xt1, Xt2, Xt3, Xt_grid)
        Xg_reference = family_volume%get_Xg()
        elem = family_volume%get_elem_Xg_vis()
        call family_volume%basis(res1=3, res2=3, res3=3, Tgc=Tgrid)
        elem = family_volume%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,4,2)
        active2 = active_knots(family_knot2,4,2)
        active3 = active_knots(family_knot3,4,2)
        xi = [&
            0.5_rk*(active1(1)+active1(2)),&
            0.5_rk*(active2(1)+active2(2)),&
            0.5_rk*(active3(1)+active3(2))]
        call family_volume%basis(xi, T)
        elem = family_volume%cmp_elem()
        call family_volume%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = T,&
            dTgc_dXg = dT_dX,&
            dV       = dV,&
            ngauss   = volume_quadrature_order_333)

        call refined_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call refined_volume%insert_knots(1, [xi(1)], [1])
        call refined_volume%insert_knots(2, [xi(2)], [1])
        call refined_volume%insert_knots(3, [xi(3)], [1])
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()
        call family_volume%err%print()
        call refined_volume%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families insertion",&
            res      = Xg_family,&
            expected = Xg_reference,&
            tol      = REFINEMENT_TOL,&
            msg      = "Mixed knot families insertion is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0194


    subroutine forcad_nurbs_volume_0195(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_333(3) = [3,3,3]
        real(rk), parameter :: REFINEMENT_TOL = 3.0e-10_rk
        type(nurbs_volume) :: family_volume, refined_volume
        real(rk) :: family_knot1(7), family_knot2(7), family_knot3(7), family_Xc(64,3)
        real(rk) :: xi(3)
        real(rk) :: dV
        real(rk), allocatable :: Xt1(:), Xt2(:), Xt3(:), Xt_grid(:,:), Xg_reference(:,:), Xg_family(:,:)
        real(rk), allocatable :: active1(:), active2(:), active3(:), T(:), Tgrid(:,:), dT_dX(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: kk
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        family_knot2 = [0.0_rk, 0.4_rk, 1.1_rk, 2.0_rk, 3.2_rk, 4.0_rk, 4.8_rk]
        family_knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.4_rk, 2.2_rk, 3.0_rk]
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    point = (kk-1)*16 + (jj-1)*4 + ii
                    family_Xc(point,:) = [real(ii-1,rk),real(jj-1,rk),real(kk-1,rk)]
                end do
            end do
        end do

        call family_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call family_volume%create(res1=5, res2=5, res3=5)
        Xt1 = family_volume%get_Xt(1)
        Xt2 = family_volume%get_Xt(2)
        Xt3 = family_volume%get_Xt(3)
        call ndgrid(Xt1, Xt2, Xt3, Xt_grid)
        Xg_reference = family_volume%get_Xg()
        elem = family_volume%get_elem_Xg_vis()
        call family_volume%basis(res1=3, res2=3, res3=3, Tgc=Tgrid)
        elem = family_volume%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,4,2)
        active2 = active_knots(family_knot2,4,2)
        active3 = active_knots(family_knot3,4,2)
        xi = [&
            0.5_rk*(active1(1)+active1(2)),&
            0.5_rk*(active2(1)+active2(2)),&
            0.5_rk*(active3(1)+active3(2))]
        call family_volume%basis(xi, T)
        elem = family_volume%cmp_elem()
        call family_volume%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = T,&
            dTgc_dXg = dT_dX,&
            dV       = dV,&
            ngauss   = volume_quadrature_order_333)

        call refined_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call refined_volume%insert_knots(1, [xi(1)], [1])
        call refined_volume%insert_knots(2, [xi(2)], [1])
        call refined_volume%insert_knots(3, [xi(3)], [1])
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()
        call refined_volume%remove_knots(1, [xi(1)], [1])
        call refined_volume%remove_knots(2, [xi(2)], [1])
        call refined_volume%remove_knots(3, [xi(3)], [1])
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()
        call family_volume%err%print()
        call refined_volume%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families removal",&
            res      = Xg_family,&
            expected = Xg_reference,&
            tol      = REFINEMENT_TOL,&
            msg      = "Mixed knot families removal is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0195


    subroutine forcad_nurbs_volume_0196(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_333(3) = [3,3,3]
        real(rk), parameter :: REFINEMENT_TOL = 3.0e-10_rk
        type(nurbs_volume) :: family_volume, refined_volume
        real(rk) :: family_knot1(7), family_knot2(7), family_knot3(7), family_Xc(64,3)
        real(rk) :: xi(3)
        real(rk) :: dV
        real(rk), allocatable :: Xt1(:), Xt2(:), Xt3(:), Xt_grid(:,:), Xg_reference(:,:), Xg_family(:,:)
        real(rk), allocatable :: active1(:), active2(:), active3(:), T(:), Tgrid(:,:), dT_dX(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: kk
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        family_knot2 = [0.0_rk, 0.4_rk, 1.1_rk, 2.0_rk, 3.2_rk, 4.0_rk, 4.8_rk]
        family_knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.4_rk, 2.2_rk, 3.0_rk]
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    point = (kk-1)*16 + (jj-1)*4 + ii
                    family_Xc(point,:) = [real(ii-1,rk),real(jj-1,rk),real(kk-1,rk)]
                end do
            end do
        end do

        call family_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call family_volume%create(res1=5, res2=5, res3=5)
        Xt1 = family_volume%get_Xt(1)
        Xt2 = family_volume%get_Xt(2)
        Xt3 = family_volume%get_Xt(3)
        call ndgrid(Xt1, Xt2, Xt3, Xt_grid)
        Xg_reference = family_volume%get_Xg()
        elem = family_volume%get_elem_Xg_vis()
        call family_volume%basis(res1=3, res2=3, res3=3, Tgc=Tgrid)
        elem = family_volume%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,4,2)
        active2 = active_knots(family_knot2,4,2)
        active3 = active_knots(family_knot3,4,2)
        xi = [&
            0.5_rk*(active1(1)+active1(2)),&
            0.5_rk*(active2(1)+active2(2)),&
            0.5_rk*(active3(1)+active3(2))]
        call family_volume%basis(xi, T)
        elem = family_volume%cmp_elem()
        call family_volume%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = T,&
            dTgc_dXg = dT_dX,&
            dV       = dV,&
            ngauss   = volume_quadrature_order_333)

        call refined_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call refined_volume%insert_knots(1, [xi(1)], [1])
        call refined_volume%insert_knots(2, [xi(2)], [1])
        call refined_volume%insert_knots(3, [xi(3)], [1])
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()
        call refined_volume%remove_knots(1, [xi(1)], [1])
        call refined_volume%remove_knots(2, [xi(2)], [1])
        call refined_volume%remove_knots(3, [xi(3)], [1])
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()

        call refined_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call refined_volume%elevate_degree(1,1)
        call refined_volume%elevate_degree(2,1)
        call refined_volume%elevate_degree(3,1)
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()
        call family_volume%err%print()
        call refined_volume%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families elevation",&
            res      = Xg_family,&
            expected = Xg_reference,&
            tol      = REFINEMENT_TOL,&
            msg      = "Mixed knot families elevation is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0196


    subroutine forcad_nurbs_volume_0197(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_333(3) = [3,3,3]
        real(rk), parameter :: FIT_TOL = 5.0e-6_rk
        type(nurbs_volume) :: family_volume, refined_volume
        real(rk) :: family_knot1(7), family_knot2(7), family_knot3(7), family_Xc(64,3)
        real(rk) :: xi(3), dV, fit_error
        real(rk), allocatable :: Xt1(:), Xt2(:), Xt3(:), Xt_grid(:,:), Xg_reference(:,:), Xg_family(:,:)
        real(rk), allocatable :: active1(:), active2(:), active3(:), T(:), Tgrid(:,:), dT_dX(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii, jj, kk, point, ndata(3)

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        family_knot2 = [0.0_rk, 0.4_rk, 1.1_rk, 2.0_rk, 3.2_rk, 4.0_rk, 4.8_rk]
        family_knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.4_rk, 2.2_rk, 3.0_rk]
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    point = (kk-1)*16 + (jj-1)*4 + ii
                    family_Xc(point,:) = [real(ii-1,rk),real(jj-1,rk),real(kk-1,rk)]
                end do
            end do
        end do

        call family_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call family_volume%create(res1=5, res2=5, res3=5)
        Xt1 = family_volume%get_Xt(1)
        Xt2 = family_volume%get_Xt(2)
        Xt3 = family_volume%get_Xt(3)
        call ndgrid(Xt1, Xt2, Xt3, Xt_grid)
        Xg_reference = family_volume%get_Xg()
        elem = family_volume%get_elem_Xg_vis()
        call family_volume%basis(res1=3, res2=3, res3=3, Tgc=Tgrid)
        elem = family_volume%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,4,2)
        active2 = active_knots(family_knot2,4,2)
        active3 = active_knots(family_knot3,4,2)
        xi = [&
            0.5_rk*(active1(1)+active1(2)),&
            0.5_rk*(active2(1)+active2(2)),&
            0.5_rk*(active3(1)+active3(2))]
        call family_volume%basis(xi, T)
        elem = family_volume%cmp_elem()
        call family_volume%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = T,&
            dTgc_dXg = dT_dX,&
            dV       = dV,&
            ngauss   = volume_quadrature_order_333)

        call refined_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call refined_volume%insert_knots(1, [xi(1)], [1])
        call refined_volume%insert_knots(2, [xi(2)], [1])
        call refined_volume%insert_knots(3, [xi(3)], [1])
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()
        call refined_volume%remove_knots(1, [xi(1)], [1])
        call refined_volume%remove_knots(2, [xi(2)], [1])
        call refined_volume%remove_knots(3, [xi(3)], [1])
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()

        call refined_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call refined_volume%elevate_degree(1,1)
        call refined_volume%elevate_degree(2,1)
        call refined_volume%elevate_degree(3,1)
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()

        ndata = [size(Xt1),size(Xt2),size(Xt3)]
        call refined_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call refined_volume%lsq_fit_bspline(Xt_grid, Xg_reference, ndata)
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()
        fit_error = norm2(Xg_family-Xg_reference)/max(1.0_rk,norm2(Xg_reference))
        call family_volume%err%print()
        call refined_volume%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families B-spline fit",&
            res      = fit_error,&
            expected = 0.0_rk,&
            tol      = FIT_TOL,&
            msg      = "Mixed knot families B-spline fit is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0197


    subroutine forcad_nurbs_volume_0198(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: volume_quadrature_order_333(3) = [3,3,3]
        real(rk), parameter :: FIT_TOL = 5.0e-6_rk
        type(nurbs_volume) :: family_volume, refined_volume
        real(rk) :: family_knot1(7), family_knot2(7), family_knot3(7), family_Xc(64,3)
        real(rk) :: xi(3), dV, fit_error
        real(rk), allocatable :: Xt1(:), Xt2(:), Xt3(:), Xt_grid(:,:), Xg_reference(:,:), Xg_family(:,:)
        real(rk), allocatable :: active1(:), active2(:), active3(:), T(:), Tgrid(:,:), dT_dX(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii, jj, kk, point, ndata(3)

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        family_knot2 = [0.0_rk, 0.4_rk, 1.1_rk, 2.0_rk, 3.2_rk, 4.0_rk, 4.8_rk]
        family_knot3 = [0.0_rk, 0.0_rk, 0.0_rk, 0.6_rk, 1.4_rk, 2.2_rk, 3.0_rk]
        do kk = 1, 4
            do jj = 1, 4
                do ii = 1, 4
                    point = (kk-1)*16 + (jj-1)*4 + ii
                    family_Xc(point,:) = [real(ii-1,rk),real(jj-1,rk),real(kk-1,rk)]
                end do
            end do
        end do

        call family_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call family_volume%create(res1=5, res2=5, res3=5)
        Xt1 = family_volume%get_Xt(1)
        Xt2 = family_volume%get_Xt(2)
        Xt3 = family_volume%get_Xt(3)
        call ndgrid(Xt1, Xt2, Xt3, Xt_grid)
        Xg_reference = family_volume%get_Xg()
        elem = family_volume%get_elem_Xg_vis()
        call family_volume%basis(res1=3, res2=3, res3=3, Tgc=Tgrid)
        elem = family_volume%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,4,2)
        active2 = active_knots(family_knot2,4,2)
        active3 = active_knots(family_knot3,4,2)
        xi = [&
            0.5_rk*(active1(1)+active1(2)),&
            0.5_rk*(active2(1)+active2(2)),&
            0.5_rk*(active3(1)+active3(2))]
        call family_volume%basis(xi, T)
        elem = family_volume%cmp_elem()
        call family_volume%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = T,&
            dTgc_dXg = dT_dX,&
            dV       = dV,&
            ngauss   = volume_quadrature_order_333)

        call refined_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call refined_volume%insert_knots(1, [xi(1)], [1])
        call refined_volume%insert_knots(2, [xi(2)], [1])
        call refined_volume%insert_knots(3, [xi(3)], [1])
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()
        call refined_volume%remove_knots(1, [xi(1)], [1])
        call refined_volume%remove_knots(2, [xi(2)], [1])
        call refined_volume%remove_knots(3, [xi(3)], [1])
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()

        call refined_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call refined_volume%elevate_degree(1,1)
        call refined_volume%elevate_degree(2,1)
        call refined_volume%elevate_degree(3,1)
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()

        ndata = [size(Xt1),size(Xt2),size(Xt3)]
        call refined_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call refined_volume%lsq_fit_bspline(Xt_grid, Xg_reference, ndata)
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()
        fit_error = norm2(Xg_family-Xg_reference)/max(1.0_rk,norm2(Xg_reference))

        call refined_volume%set(family_knot1, family_knot2, family_knot3, family_Xc, degree=[2,2,2])
        call refined_volume%lsq_fit_nurbs(&
            Xt    = Xt_grid,&
            Xdata = Xg_reference,&
            ndata = ndata,&
            maxit = 30,&
            tol   = sqrt(epsilon(1.0_rk)))
        call refined_volume%create(Xt=Xt_grid)
        Xg_family = refined_volume%get_Xg()
        fit_error = norm2(Xg_family-Xg_reference)/max(1.0_rk,norm2(Xg_reference))
        call family_volume%err%print()
        call refined_volume%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families NURBS fit",&
            res      = fit_error,&
            expected = 0.0_rk,&
            tol      = FIT_TOL,&
            msg      = "Mixed knot families NURBS fit is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0198


    subroutine forcad_nurbs_volume_0199(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: full_break_volume
        real(rk) :: knot1_full(9), knot2_full(4), knot3_full(4), Xc_full(24,3), Wc_full(24), Xt_eval(20,3)
        real(rk) :: xline(6), zline(6), wline(6), u_eval(5), v_eval(2), w_eval(2)
        real(rk), allocatable :: Xg_before(:,:), Xg_after(:,:)
        integer :: ii, jj, kk, point

        knot1_full = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 0.5_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot2_full = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot3_full = knot2_full
        xline = [0.0_rk, 0.2_rk, 0.5_rk, 0.5_rk, 0.8_rk, 1.0_rk]
        zline = [0.0_rk, 0.8_rk, 0.0_rk, 2.0_rk, 2.8_rk, 2.0_rk]
        wline = [1.0_rk, 0.7_rk, 1.2_rk, 0.9_rk, 1.3_rk, 1.0_rk]
        u_eval = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk]
        v_eval = [0.2_rk, 0.8_rk]
        w_eval = [0.3_rk, 0.7_rk]
        do kk = 1, 2
            do jj = 1, 2
                do ii = 1, 6
                    point = (kk-1)*12 + (jj-1)*6 + ii
                    Xc_full(point,:) = [&
                        xline(ii),&
                        real(jj-1,rk)+0.1_rk*real(kk-1,rk),&
                        zline(ii)+real(kk-1,rk)]
                    Wc_full(point) = wline(ii)*(1.0_rk+0.04_rk*real(jj+kk-2,rk))
                end do
            end do
        end do
        do kk = 1, 2
            do jj = 1, 2
                do ii = 1, 5
                    point = (kk-1)*10 + (jj-1)*5 + ii
                    Xt_eval(point,:) = [u_eval(ii),v_eval(jj),w_eval(kk)]
                end do
            end do
        end do

        call full_break_volume%set(knot1_full, knot2_full, knot3_full, Xc_full, Wc_full, degree=[2,1,1])
        call full_break_volume%create(Xt=Xt_eval)
        Xg_before = full_break_volume%get_Xg()
        call full_break_volume%elevate_degree(1,1)
        call full_break_volume%create(Xt=Xt_eval)
        Xg_after = full_break_volume%get_Xg()
        call full_break_volume%err%print()

        call ut%test(ti)%check(&
            name     = "full-break elevation geometry",&
            res      = Xg_after,&
            expected = Xg_before,&
            tol      = 3.0e-11_rk,&
            msg      = "Full-break elevation geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0199


    subroutine forcad_nurbs_volume_0200(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: full_break_volume
        real(rk) :: knot1_full(9), knot2_full(4), knot3_full(4), Xc_full(24,3), Wc_full(24), Xt_eval(20,3)
        real(rk) :: xline(6), zline(6), wline(6), u_eval(5), v_eval(2), w_eval(2)
        real(rk), allocatable :: Xg_before(:,:), Xg_after(:,:)
        integer :: ii, jj, kk, point

        knot1_full = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 0.5_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot2_full = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot3_full = knot2_full
        xline = [0.0_rk, 0.2_rk, 0.5_rk, 0.5_rk, 0.8_rk, 1.0_rk]
        zline = [0.0_rk, 0.8_rk, 0.0_rk, 2.0_rk, 2.8_rk, 2.0_rk]
        wline = [1.0_rk, 0.7_rk, 1.2_rk, 0.9_rk, 1.3_rk, 1.0_rk]
        u_eval = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk]
        v_eval = [0.2_rk, 0.8_rk]
        w_eval = [0.3_rk, 0.7_rk]
        do kk = 1, 2
            do jj = 1, 2
                do ii = 1, 6
                    point = (kk-1)*12 + (jj-1)*6 + ii
                    Xc_full(point,:) = [&
                        xline(ii),&
                        real(jj-1,rk)+0.1_rk*real(kk-1,rk),&
                        zline(ii)+real(kk-1,rk)]
                    Wc_full(point) = wline(ii)*(1.0_rk+0.04_rk*real(jj+kk-2,rk))
                end do
            end do
        end do
        do kk = 1, 2
            do jj = 1, 2
                do ii = 1, 5
                    point = (kk-1)*10 + (jj-1)*5 + ii
                    Xt_eval(point,:) = [u_eval(ii),v_eval(jj),w_eval(kk)]
                end do
            end do
        end do

        call full_break_volume%set(knot1_full, knot2_full, knot3_full, Xc_full, Wc_full, degree=[2,1,1])
        call full_break_volume%create(Xt=Xt_eval)
        Xg_before = full_break_volume%get_Xg()
        call full_break_volume%elevate_degree(1,1)
        call full_break_volume%create(Xt=Xt_eval)
        Xg_after = full_break_volume%get_Xg()
        call full_break_volume%err%print()

        call ut%test(ti)%check(&
            name     = "full-break elevation degree",&
            res      = full_break_volume%get_degree(),&
            expected = [3,1,1],&
            msg      = "Full-break elevation degree is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0200


    subroutine forcad_nurbs_volume_0201(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NGAUSS(3) = [2,2,2]
        real(rk), parameter :: H = 1.0e-10_rk
        type(nurbs_volume) :: tiny_volume
        real(rk) :: linear_knot(4)
        real(rk) :: Xc_tiny(8,3)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:)
        integer :: i1, i2, i3, point

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        do i3 = 1, 2
            do i2 = 1, 2
                do i1 = 1, 2
                    point = (i3-1)*4 + (i2-1)*2 + i1
                    Xc_tiny(point,:) = [&
                        real(i1-1,rk)*H,&
                        real(i2-1,rk)*H,&
                        real(i3-1,rk)*H]
                end do
            end do
        end do
        call tiny_volume%set(linear_knot, linear_knot, linear_knot, Xc_tiny, degree=[1,1,1])
        call tiny_volume%ansatz(1, 1, T, dT_dX, dV, ngauss=NGAUSS)
        call tiny_volume%err%print()

        call ut%test(ti)%check(&
            name     = "tiny volume ansatz diagnostic",&
            res      = tiny_volume%err%ok,&
            expected = .true.,&
            msg      = "A regular volume remains valid when its physical scale is small",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0201


    subroutine forcad_nurbs_volume_0202(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NGAUSS(3) = [2,2,2]
        real(rk), parameter :: H = 1.0e-10_rk
        type(nurbs_volume) :: tiny_volume
        real(rk) :: linear_knot(4)
        real(rk) :: Xc_tiny(8,3)
        real(rk) :: dV
        real(rk), allocatable :: T(:), dT_dX(:,:)
        integer :: i1, i2, i3, point

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        do i3 = 1, 2
            do i2 = 1, 2
                do i1 = 1, 2
                    point = (i3-1)*4 + (i2-1)*2 + i1
                    Xc_tiny(point,:) = [&
                        real(i1-1,rk)*H,&
                        real(i2-1,rk)*H,&
                        real(i3-1,rk)*H]
                end do
            end do
        end do
        call tiny_volume%set(linear_knot, linear_knot, linear_knot, Xc_tiny, degree=[1,1,1])
        call tiny_volume%ansatz(1, 1, T, dT_dX, dV, ngauss=NGAUSS)
        call tiny_volume%err%print()

        call ut%test(ti)%check(&
            name     = "tiny volume ansatz measure",&
            res      = dV > 0.0_rk,&
            expected = .true.,&
            msg      = "A small regular volume must retain positive measure",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0202


    subroutine forcad_nurbs_volume_0203(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NGAUSS(3) = [2,2,2]
        real(rk), parameter :: H = 1.0e-10_rk
        type(nurbs_volume) :: tiny_volume
        real(rk) :: linear_knot(4), Xc_tiny(8,3), volume_measure, dV
        real(rk), allocatable :: T(:), dT_dX(:,:)
        integer :: i1, i2, i3, point

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        do i3 = 1, 2
            do i2 = 1, 2
                do i1 = 1, 2
                    point = (i3-1)*4 + (i2-1)*2 + i1
                    Xc_tiny(point,:) = [&
                        real(i1-1,rk)*H,&
                        real(i2-1,rk)*H,&
                        real(i3-1,rk)*H]
                end do
            end do
        end do
        call tiny_volume%set(linear_knot, linear_knot, linear_knot, Xc_tiny, degree=[1,1,1])
        call tiny_volume%ansatz(1, 1, T, dT_dX, dV, ngauss=NGAUSS)
        call tiny_volume%cmp_volume(volume_measure, ngauss=NGAUSS)
        call tiny_volume%err%print()

        call ut%test(ti)%check(&
            name     = "tiny volume measure",&
            res      = volume_measure/(H*H*H),&
            expected = 1.0_rk,&
            tol      = 3.0e-12_rk,&
            msg      = "Volume integration must be relative-scale robust",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0203


    subroutine forcad_nurbs_volume_0204(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: periodic_volume
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(24,3), periodic_Wc(24), ring_Wc(6)
        real(rk) :: radius
        integer :: i1, i2, i3, n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    periodic_Xc(n,:) = [radius*ring(i1,1), radius*ring(i1,2), real(i3-1,rk)]
                    periodic_Wc(n) = ring_Wc(i1)
                end do
            end do
        end do

        call periodic_volume%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            knot3    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1,1],&
            wrap_parameters = [.true.,.false.,.false.])
        call periodic_volume%err%print()

        call ut%test(ti)%check(&
            name     = "periodic volume parameter wrapping",&
            res      = periodic_volume%get_parameter_wrapping(1),&
            expected = .true.,&
            msg      = "Periodic volume parameter wrapping is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0204


    subroutine forcad_nurbs_volume_0205(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: periodic_volume
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(24,3), periodic_Wc(24), ring_Wc(6)
        real(rk) :: radius
        integer :: i1, i2, i3, n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    periodic_Xc(n,:) = [radius*ring(i1,1), radius*ring(i1,2), real(i3-1,rk)]
                    periodic_Wc(n) = ring_Wc(i1)
                end do
            end do
        end do

        call periodic_volume%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            knot3    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1,1],&
            wrap_parameters = [.true.,.false.,.false.])
        call periodic_volume%err%print()

        call ut%test(ti)%check(&
            name     = "bounded volume second direction",&
            res      = periodic_volume%get_parameter_wrapping(2),&
            expected = .false.,&
            msg      = "Enabling one volume direction must not wrap the second direction",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0205


    subroutine forcad_nurbs_volume_0206(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: periodic_volume
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(24,3), periodic_Wc(24), ring_Wc(6)
        real(rk) :: radius
        integer :: i1, i2, i3, n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    periodic_Xc(n,:) = [radius*ring(i1,1), radius*ring(i1,2), real(i3-1,rk)]
                    periodic_Wc(n) = ring_Wc(i1)
                end do
            end do
        end do

        call periodic_volume%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            knot3    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1,1],&
            wrap_parameters = [.true.,.false.,.false.])
        call periodic_volume%err%print()

        call ut%test(ti)%check(&
            name     = "bounded volume third direction",&
            res      = periodic_volume%get_parameter_wrapping(3),&
            expected = .false.,&
            msg      = "Enabling one volume direction must not wrap the third direction",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0206


    subroutine forcad_nurbs_volume_0207(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: PERIODIC_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: periodic_volume
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(24,3), periodic_Wc(24), ring_Wc(6)
        real(rk) :: radius
        integer :: i1, i2, i3, n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    periodic_Xc(n,:) = [radius*ring(i1,1), radius*ring(i1,2), real(i3-1,rk)]
                    periodic_Wc(n) = ring_Wc(i1)
                end do
            end do
        end do

        call periodic_volume%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            knot3    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1,1],&
            wrap_parameters = [.true.,.false.,.false.])
        call periodic_volume%err%print()

        call ut%test(ti)%check(&
            name     = "wrapped periodic NURBS volume geometry",&
            res      = periodic_volume%cmp_Xg([10.75_rk,0.35_rk,0.65_rk]),&
            expected = periodic_volume%cmp_Xg([2.75_rk,0.35_rk,0.65_rk]),&
            tol      = PERIODIC_TOL,&
            msg      = "Wrapped periodic NURBS volume geometry is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0207


    subroutine forcad_nurbs_volume_0208(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: PERIODIC_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: periodic_volume
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(24,3), periodic_Wc(24), ring_Wc(6)
        real(rk), allocatable :: T_start(:), T_wrapped(:), dT_start(:,:), dT_wrapped(:,:)
        real(rk) :: radius
        integer :: i1, i2, i3, n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    periodic_Xc(n,:) = [radius*ring(i1,1), radius*ring(i1,2), real(i3-1,rk)]
                    periodic_Wc(n) = ring_Wc(i1)
                end do
            end do
        end do

        call periodic_volume%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            knot3    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1,1],&
            wrap_parameters = [.true.,.false.,.false.])
        call periodic_volume%derivative([2.75_rk,0.35_rk,0.65_rk], dT_start, T_start)
        call periodic_volume%derivative([-1.25_rk,0.35_rk,0.65_rk], dT_wrapped, T_wrapped)
        call periodic_volume%err%print()

        call ut%test(ti)%check(&
            name     = "wrapped periodic NURBS volume basis",&
            res      = T_wrapped,&
            expected = T_start,&
            tol      = PERIODIC_TOL,&
            msg      = "Wrapped periodic NURBS volume basis is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0208


    subroutine forcad_nurbs_volume_0209(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: PERIODIC_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: periodic_volume
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(24,3), periodic_Wc(24), ring_Wc(6)
        real(rk), allocatable :: T_start(:), T_wrapped(:), dT_start(:,:), dT_wrapped(:,:)
        real(rk) :: radius
        integer :: i1, i2, i3, n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    periodic_Xc(n,:) = [radius*ring(i1,1), radius*ring(i1,2), real(i3-1,rk)]
                    periodic_Wc(n) = ring_Wc(i1)
                end do
            end do
        end do

        call periodic_volume%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            knot3    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1,1],&
            wrap_parameters = [.true.,.false.,.false.])
        call periodic_volume%derivative([2.75_rk,0.35_rk,0.65_rk], dT_start, T_start)
        call periodic_volume%derivative([-1.25_rk,0.35_rk,0.65_rk], dT_wrapped, T_wrapped)
        call periodic_volume%err%print()

        call ut%test(ti)%check(&
            name     = "wrapped periodic NURBS volume derivatives",&
            res      = dT_wrapped,&
            expected = dT_start,&
            tol      = PERIODIC_TOL,&
            msg      = "Wrapped periodic NURBS volume derivatives are incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0209


    subroutine forcad_nurbs_volume_0210(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: periodic_volume
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(24,3), periodic_Wc(24), ring_Wc(6)
        real(rk), allocatable :: T_start(:), T_wrapped(:), dT_start(:,:), dT_wrapped(:,:)
        real(rk) :: radius
        integer :: i1, i2, i3, n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    periodic_Xc(n,:) = [radius*ring(i1,1), radius*ring(i1,2), real(i3-1,rk)]
                    periodic_Wc(n) = ring_Wc(i1)
                end do
            end do
        end do

        call periodic_volume%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            knot3    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1,1],&
            wrap_parameters = [.true.,.false.,.false.])
        call periodic_volume%derivative([2.75_rk,0.35_rk,0.65_rk], dT_start, T_start)
        call periodic_volume%derivative([-1.25_rk,0.35_rk,0.65_rk], dT_wrapped, T_wrapped)
        call periodic_volume%insert_knots(1, [3.5_rk], [1])
        call periodic_volume%err%print()

        call ut%test(ti)%check(&
            name     = "parameter wrapping after volume refinement",&
            res      = periodic_volume%get_parameter_wrapping(1),&
            expected = .true.,&
            msg      = "Parameter wrapping after volume refinement is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0210


    subroutine forcad_nurbs_volume_0211(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: PERIODIC_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: periodic_volume
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(24,3), periodic_Wc(24), ring_Wc(6)
        real(rk), allocatable :: T_start(:), T_wrapped(:), dT_start(:,:), dT_wrapped(:,:)
        real(rk) :: radius
        integer :: i1, i2, i3, n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    periodic_Xc(n,:) = [radius*ring(i1,1), radius*ring(i1,2), real(i3-1,rk)]
                    periodic_Wc(n) = ring_Wc(i1)
                end do
            end do
        end do

        call periodic_volume%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            knot3    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1,1],&
            wrap_parameters = [.true.,.false.,.false.])
        call periodic_volume%derivative([2.75_rk,0.35_rk,0.65_rk], dT_start, T_start)
        call periodic_volume%derivative([-1.25_rk,0.35_rk,0.65_rk], dT_wrapped, T_wrapped)
        call periodic_volume%insert_knots(1, [3.5_rk], [1])
        call periodic_volume%err%print()

        call ut%test(ti)%check(&
            name     = "periodic volume refinement wrapping",&
            res      = periodic_volume%cmp_Xg([10.75_rk,0.35_rk,0.65_rk]),&
            expected = periodic_volume%cmp_Xg([2.75_rk,0.35_rk,0.65_rk]),&
            tol      = PERIODIC_TOL,&
            msg      = "A refined periodic volume must retain repeated evaluation",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0211


    subroutine forcad_nurbs_volume_0212(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: export_volume, copied_volume
        real(rk) :: linear_knot(4), export_Xc(8,3)

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        export_Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        export_Xc(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        export_Xc(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        export_Xc(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        export_Xc(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        export_Xc(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        export_Xc(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        export_Xc(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]

        call export_volume%set(linear_knot, linear_knot, linear_knot, export_Xc, degree=[1,1,1])
        copied_volume = export_volume
        call export_volume%modify_Xc(5.0_rk, 1, 1)
        call export_volume%err%print()
        call copied_volume%err%print()

        call ut%test(ti)%check(&
            name     = "volume intrinsic assignment deep copy",&
            res      = copied_volume%get_Xc(1,1),&
            expected = 0.0_rk,&
            tol      = 0.0_rk,&
            msg      = "Volume intrinsic assignment deep copy is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0212


    subroutine forcad_nurbs_volume_0213(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: vtkfile = "vtk/forcad_test_volume_Xth_in_Xg.vtk"
        character(len=*), parameter :: vtk_header = "# vtk DataFile Version 2.0"
        type(nurbs_volume) :: export_volume, copied_volume
        real(rk) :: linear_knot(4), export_Xc(8,3)
        character(len=len(vtk_header)) :: actual_header
        integer :: file_unit, io_status
        logical :: file_exists

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        export_Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        export_Xc(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        export_Xc(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        export_Xc(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        export_Xc(5,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        export_Xc(6,:) = [1.0_rk, 0.0_rk, 1.0_rk]
        export_Xc(7,:) = [0.0_rk, 1.0_rk, 1.0_rk]
        export_Xc(8,:) = [1.0_rk, 1.0_rk, 1.0_rk]

        call export_volume%set(linear_knot, linear_knot, linear_knot, export_Xc, degree=[1,1,1])
        copied_volume = export_volume
        call export_volume%modify_Xc(5.0_rk, 1, 1)

        open(newunit=file_unit, file=vtkfile, status="replace")
        close(file_unit, status="delete")
        call copied_volume%export_Xth_in_Xg(vtkfile, res=3)
        actual_header = ""
        io_status = 1
        inquire(file=vtkfile, exist=file_exists)
        if (file_exists) then
            open(&
                newunit = file_unit,&
                file    = vtkfile,&
                access  = "stream",&
                form    = "unformatted",&
                action  = "read",&
                status  = "old",&
                iostat  = io_status)
            if (io_status == 0) then
                read(file_unit, iostat=io_status) actual_header
                close(file_unit)
            end if
        end if
        call export_volume%err%print()
        call copied_volume%err%print()

        call ut%test(ti)%check(&
            name     = "volume export Xth in Xg",&
            res      = file_exists .and. io_status == 0 .and. actual_header == vtk_header,&
            expected = .true.,&
            msg      = "Volume export Xth in Xg is incorrect.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0213


    subroutine forcad_nurbs_volume_0214(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot1(8) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: knot2(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: knot3(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(20,3) = reshape([&
            0.0_rk,0.0_rk,1.0_rk,0.0_rk,0.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,0.0_rk,0.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,0.0_rk,0.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,0.0_rk,0.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,&
            1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,&
            1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,&
            1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,&
            1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [20,3])
        type(nurbs_volume) :: volume
        real(rk), allocatable :: knot_after(:), Xc_after(:,:)
        logical :: preserved

        call volume%set(knot1, knot2, knot3, Xc)
        call volume%remove_knots(1, [0.5_rk], [1])
        call volume%err%print()

        knot_after = volume%get_knot(1)
        Xc_after = volume%get_Xc()
        preserved = volume%err%ok .and. size(knot_after) == size(knot1) .and. &
            size(Xc_after,1) == size(Xc,1) .and. size(Xc_after,2) == size(Xc,2)
        if (preserved) preserved = maxval(abs(knot_after-knot1)) <= epsilon(1.0_rk) .and. &
            maxval(abs(Xc_after-Xc)) <= epsilon(1.0_rk)

        call ut%test(ti)%check(&
            name     = "volume rejects geometry-changing removal",&
            res      = preserved,&
            expected = .true.,&
            msg      = "Volume removal changed nonremovable geometry.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0214


    subroutine forcad_nurbs_volume_0215(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume

        call volume%set_hexahedron([1.0_rk,1.0_rk,1.0_rk], [2,2,2])
        call volume%set(&
            Xth_dir1    = [0.0_rk,0.5_rk,1.0_rk],&
            Xth_dir2    = [0.0_rk,0.5_rk,1.0_rk],&
            Xth_dir3    = [0.0_rk,0.5_rk,1.0_rk],&
            degree      = [2,2,2],&
            continuity1 = [-1,1,-1],&
            continuity2 = [-1,1,-1],&
            continuity3 = [-1,1,-1])
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "set2 clears old volume controls",&
            res      = size(volume%get_Xc(),1),&
            expected = 0,&
            msg      = "Omitting Xc must remove the previous volume controls.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0215


    subroutine forcad_nurbs_volume_0216(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume
        real(rk), parameter :: periodic_knot(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        real(rk), parameter :: linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: ring_Wc(6) = [&
            1.0_rk,1.4_rk,0.8_rk,1.2_rk,1.0_rk,1.4_rk]
        real(rk) :: ring(6,2), Xc(24,3), Wc(24), radius
        integer :: i1, i2, i3, n

        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    Xc(n,:) = [&
                        radius*ring(i1,1),radius*ring(i1,2),real(i3-1,rk)]
                    Wc(n) = ring_Wc(i1)
                end do
            end do
        end do

        call volume%set(&
            knot1  = periodic_knot,&
            knot2  = linear_knot,&
            knot3  = linear_knot,&
            Xc     = Xc,&
            Wc     = Wc,&
            degree = [2,1,1])
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "volume directional topology",&
            res      = volume%get_parameter_topology(1) == "periodic" .and. &
                volume%get_parameter_topology(2) == "bounded" .and. &
                volume%get_parameter_topology(3) == "bounded",&
            expected = .true.,&
            msg      = "Volume parameter topologies were classified incorrectly.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0216


    subroutine forcad_nurbs_volume_0217(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume
        real(rk), parameter :: periodic_knot(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        real(rk), parameter :: linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        integer, parameter :: expected_map(24) = [&
            1,2,3,4,1,2,5,6,7,8,5,6,9,10,11,12,9,10,13,14,15,16,13,14]
        real(rk) :: ring(6,2), Xc(24,3), radius
        integer, allocatable :: map(:)
        integer :: i1, i2, i3, n

        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    Xc(n,:) = [&
                        radius*ring(i1,1),radius*ring(i1,2),real(i3-1,rk)]
                end do
            end do
        end do

        call volume%set(&
            knot1  = periodic_knot,&
            knot2  = linear_knot,&
            knot3  = linear_knot,&
            Xc     = Xc,&
            degree = [2,1,1])
        map = volume%cmp_dof_map()
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "volume periodic DOF map",&
            res      = map,&
            expected = expected_map,&
            msg      = "Repeated periodic volume layers were not identified.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0217


    subroutine forcad_nurbs_volume_0218(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume
        real(rk), parameter :: periodic_knot(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        real(rk), parameter :: linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: ring_Wc(6) = [&
            1.0_rk,1.4_rk,0.8_rk,1.2_rk,1.0_rk,1.4_rk]
        real(rk) :: ring(6,2), Xc(24,3), Wc(24), radius
        integer :: i1, i2, i3, n

        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    Xc(n,:) = [&
                        radius*ring(i1,1),radius*ring(i1,2),real(i3-1,rk)]
                    Wc(n) = ring_Wc(i1)
                end do
            end do
        end do

        call volume%set(&
            knot1  = periodic_knot,&
            knot2  = linear_knot,&
            knot3  = linear_knot,&
            Xc     = Xc,&
            Wc     = Wc,&
            degree = [2,1,1])
        call volume%elevate_degree(1, 1)
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "volume elevated periodic topology",&
            res      = volume%get_parameter_topology(1) == "periodic",&
            expected = .true.,&
            msg      = "Volume degree elevation lost periodic topology.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0218


    subroutine forcad_nurbs_volume_0219(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 65536.0_rk*epsilon(1.0_rk)
        type(nurbs_volume) :: volume
        real(rk), parameter :: periodic_knot(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        real(rk), parameter :: linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: ring_Wc(6) = [&
            1.0_rk,1.4_rk,0.8_rk,1.2_rk,1.0_rk,1.4_rk]
        real(rk), parameter :: Xt(3,5) = reshape([&
            2.0_rk,0.0_rk,0.0_rk,2.2_rk,0.3_rk,0.2_rk,&
            3.1_rk,0.5_rk,0.5_rk,4.4_rk,0.7_rk,0.8_rk,&
            6.0_rk,1.0_rk,1.0_rk], [3,5])
        real(rk) :: ring(6,2), Xc(24,3), Wc(24), radius
        real(rk) :: actual(5,3), expected(5,3)
        integer :: i, i1, i2, i3, n

        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    Xc(n,:) = [&
                        radius*ring(i1,1),radius*ring(i1,2),real(i3-1,rk)]
                    Wc(n) = ring_Wc(i1)
                end do
            end do
        end do

        call volume%set(&
            knot1  = periodic_knot,&
            knot2  = linear_knot,&
            knot3  = linear_knot,&
            Xc     = Xc,&
            Wc     = Wc,&
            degree = [2,1,1])
        do i = 1, size(Xt,2)
            expected(i,:) = volume%cmp_Xg(Xt(:,i))
        end do
        call volume%elevate_degree(1, 1)
        do i = 1, size(Xt,2)
            actual(i,:) = volume%cmp_Xg(Xt(:,i))
        end do
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "volume periodic elevation geometry",&
            res      = actual,&
            expected = expected,&
            tol      = tol,&
            msg      = "Periodic degree elevation changed volume geometry.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0219


    subroutine forcad_nurbs_volume_0220(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume
        real(rk), parameter :: periodic_knot(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        real(rk), parameter :: linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: knots_to_remove(1) = [4.0_rk]
        integer, parameter :: removals(1) = [1]
        real(rk) :: ring(6,2), Xc(24,3), radius
        integer :: i1, i2, i3, n

        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        do i3 = 1, 2
            do i2 = 1, 2
                radius = 1.0_rk + 0.25_rk*real(i2-1,rk)
                do i1 = 1, 6
                    n = i1 + (i2-1)*6 + (i3-1)*12
                    Xc(n,:) = [&
                        radius*ring(i1,1),radius*ring(i1,2),real(i3-1,rk)]
                end do
            end do
        end do

        call volume%set(&
            knot1  = periodic_knot,&
            knot2  = linear_knot,&
            knot3  = linear_knot,&
            Xc     = Xc,&
            degree = [2,1,1])
        call volume%remove_knots(1, knots_to_remove, removals)
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "periodic volume removal diagnostic",&
            res      = volume%err%ok,&
            expected = .false.,&
            msg      = "Periodic removal must not destroy volume topology.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0220


    subroutine forcad_nurbs_volume_0221(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(8,3) = reshape([&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [8,3])
        integer, parameter :: degree(3) = [1,1,1]
        integer, parameter :: ngauss(3) = [1,1,1]
        type(nurbs_volume) :: volume
        real(rk), allocatable :: Tgc(:), dTgc_dXg(:,:), d2Tgc_dXg2(:,:,:)
        real(rk) :: dV

        call volume%set(&
            knot1  = knot,&
            knot2  = knot,&
            knot3  = knot,&
            Xc     = Xc,&
            degree = degree)
        call volume%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = Tgc,&
            dTgc_dXg    = dTgc_dXg,&
            dV          = dV,&
            ngauss      = ngauss,&
            d2Tgc_dXg2  = d2Tgc_dXg2,&
            strict      = .true.)
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "strict volume rejects inversion",&
            res      = .not. volume%err%ok .and. abs(dV) <= tiny(1.0_rk) .and. &
                maxval(abs(dTgc_dXg)) <= tiny(1.0_rk) .and. &
                maxval(abs(d2Tgc_dXg2)) <= tiny(1.0_rk),&
            expected = .true.,&
            msg      = "Strict ansatz accepted an inverted volume map.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0221


    subroutine forcad_nurbs_volume_0222(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(8,3) = reshape([&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [8,3])
        integer, parameter :: degree(3) = [1,1,1]
        integer, parameter :: ngauss(3) = [1,1,1]
        type(nurbs_volume) :: volume
        real(rk), allocatable :: Tgc(:), dTgc_dXg(:,:)
        real(rk) :: dV

        call volume%set(&
            knot1  = knot,&
            knot2  = knot,&
            knot3  = knot,&
            Xc     = Xc,&
            degree = degree)
        call volume%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = Tgc,&
            dTgc_dXg = dTgc_dXg,&
            dV       = dV,&
            ngauss   = ngauss)
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "volume default measures inversion",&
            res      = volume%err%ok .and. dV > 0.0_rk .and. &
                all(ieee_is_finite(dTgc_dXg)),&
            expected = .true.,&
            msg      = "Default ansatz lost absolute volume measure.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0222


    subroutine forcad_nurbs_volume_0223(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(8,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [8,3])
        integer, parameter :: degree(3) = [1,1,1]
        integer, parameter :: ngauss(3) = [1,1,1]
        type(nurbs_volume) :: volume
        real(rk), allocatable :: Tgc(:), dTgc_dXg(:,:)
        real(rk) :: dV

        call volume%set(&
            knot1  = knot,&
            knot2  = knot,&
            knot3  = knot,&
            Xc     = Xc,&
            degree = degree)
        call volume%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = Tgc,&
            dTgc_dXg = dTgc_dXg,&
            dV       = dV,&
            ngauss   = ngauss,&
            strict   = .true.)
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "strict volume accepts orientation",&
            res      = volume%err%ok .and. dV > 0.0_rk .and. &
                all(ieee_is_finite(dTgc_dXg)),&
            expected = .true.,&
            msg      = "Strict ansatz rejected a positive volume map.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0223


    subroutine forcad_nurbs_volume_0224(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 128.0_rk*epsilon(1.0_rk)
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(8,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [8,3])
        real(rk), parameter :: X(4,3) = reshape([&
            0.0_rk,1.0_rk,0.5_rk,0.25_rk,&
            0.0_rk,0.25_rk,1.0_rk,0.5_rk,&
            0.0_rk,0.5_rk,0.25_rk,1.0_rk], [4,3])
        type(nurbs_volume) :: volume
        real(rk), allocatable :: Xt1(:), Xt2(:), Xt3(:), Xg(:,:)
        integer, allocatable :: elem(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(3)
        integer :: id

        call volume%set(&
            knot1  = knot,&
            knot2  = knot,&
            knot3  = knot,&
            Xc     = Xc,&
            degree = [1,1,1])
        call volume%create(res1=3, res2=3, res3=3)
        call volume%put_to_nurbs(X)
        call volume%nearest_point(&
            point_Xg   = [1.0_rk,0.25_rk,0.5_rk],&
            nearest_Xg = nearest_Xg,&
            nearest_Xt = nearest_Xt,&
            id         = id)
        Xt1 = volume%get_Xt(1)
        Xt2 = volume%get_Xt(2)
        Xt3 = volume%get_Xt(3)
        Xg = volume%get_Xg()
        elem = volume%get_elem_Xg_vis()
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "put_to_nurbs volume cache replacement",&
            res      = volume%err%ok .and. all(volume%get_ng() == [4,1,1]) .and. &
                size(Xt1) == 0 .and. size(Xt2) == 0 .and. size(Xt3) == 0 .and. &
                all(shape(Xg) == [4,3]) .and. all(shape(elem) == [0,8]) .and. &
                id == 2 .and. maxval(abs(nearest_Xt-[1.0_rk,0.25_rk,0.5_rk])) <= tol .and. &
                maxval(abs(nearest_Xg-[1.0_rk,0.25_rk,0.5_rk])) <= tol,&
            expected = .true.,&
            msg      = "Volume parameter tuples and geometry cache are inconsistent.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0224


    subroutine forcad_nurbs_volume_0225(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 128.0_rk*epsilon(1.0_rk)
        real(rk), parameter :: knot1(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: knot2(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: u(9) = [&
            0.0_rk,0.125_rk,0.25_rk,0.375_rk,0.5_rk,&
            0.625_rk,0.75_rk,0.875_rk,1.0_rk]
        real(rk), parameter :: y(9) = [&
            -0.52613569629176171_rk,-2.7932097690057072_rk,&
            -0.12457766409497265_rk,-1.0921215668918383_rk,&
            -0.10431650515569141_rk,-1.1755376952645962_rk,&
            -2.2127487058295481_rk,0.02135077665918364_rk,&
            -0.26585387032100716_rk]
        type(nurbs_volume) :: volume
        real(rk) :: Xt(36,3), Xdata(36,3), Xc(16,3)
        real(rk), allocatable :: Wc(:)
        integer :: i, j, k, point
        logical :: fit_failed

        do k = 1, 2
            do j = 1, 2
                do i = 1, 9
                    point = (k-1)*18 + (j-1)*9 + i
                    Xt(point,:) = [u(i),real(j-1,rk),real(k-1,rk)]
                    Xdata(point,:) = [y(i),real(j-1,rk),real(k-1,rk)]
                end do
                do i = 1, 4
                    point = (k-1)*8 + (j-1)*4 + i
                    Xc(point,:) = [real(i-1,rk),real(j-1,rk),real(k-1,rk)]
                end do
            end do
        end do

        call volume%set(&
            knot1  = knot1,&
            knot2  = knot2,&
            knot3  = knot2,&
            Xc     = Xc,&
            degree = [2,1,1])
        call volume%lsq_fit_nurbs(&
            Xt       = Xt,&
            Xdata    = Xdata,&
            ndata    = [9,2,2],&
            maxit    = 2,&
            tol      = 0.0_rk,&
            mu0      = 0.0_rk,&
            reg_logw = 0.0_rk)
        fit_failed = .not. volume%err%ok
        call volume%err%print()
        call volume%err%reset()
        Wc = volume%get_Wc()

        call ut%test(ti)%check(&
            name     = "volume NURBS fit maxit rollback",&
            res      = fit_failed .and. size(Wc) == 16 .and. &
                maxval(abs(Wc-1.0_rk)) <= tol,&
            expected = .true.,&
            msg      = "Volume fit retained a pending trial after maxit.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0225


    subroutine forcad_nurbs_volume_0226(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 1.0e-12_rk
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(8,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [8,3])
        type(nurbs_volume) :: volume
        real(rk) :: nearest_Xt(3), nearest_Xg(3)

        call volume%set(&
            knot1  = knot,&
            knot2  = knot,&
            knot3  = knot,&
            Xc     = Xc,&
            degree = [1,1,1])
        call volume%nearest_point2(&
            point_Xg   = [-1.0_rk,0.37_rk,0.63_rk],&
            tol        = tol,&
            maxit      = 5,&
            nearest_Xt = nearest_Xt,&
            nearest_Xg = nearest_Xg)
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "volume projected boundary descent",&
            res      = volume%err%ok .and. &
                maxval(abs(nearest_Xt-[0.0_rk,0.37_rk,0.63_rk])) <= tol .and. &
                maxval(abs(nearest_Xg-[0.0_rk,0.37_rk,0.63_rk])) <= tol,&
            expected = .true.,&
            msg      = "An active volume bound blocked free-direction descent.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0226


    subroutine forcad_nurbs_volume_0227(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 16384.0_rk*epsilon(1.0_rk)
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(8,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [8,3])
        real(rk), parameter :: Wc(8) = [&
            1.0_rk,2.0_rk,3.0_rk,4.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk]
        type(nurbs_volume) :: reference, small_scale, large_scale
        real(rk), allocatable :: Xg_ref(:), Xg_small(:), Xg_large(:)
        real(rk), allocatable :: T_ref(:), T_small(:), T_large(:)
        real(rk), allocatable :: dT_ref(:,:), dT_small(:,:), dT_large(:,:)
        real(rk), allocatable :: d2T_ref(:,:), d2T_small(:,:), d2T_large(:,:)

        call reference%set(knot, knot, knot, Xc, Wc, degree=[1,1,1])
        call small_scale%set(knot, knot, knot, Xc, 1.0e-200_rk*Wc, degree=[1,1,1])
        call large_scale%set(knot, knot, knot, Xc, 1.0e200_rk*Wc, degree=[1,1,1])
        Xg_ref = reference%cmp_Xg([0.23_rk,0.47_rk,0.71_rk])
        Xg_small = small_scale%cmp_Xg([0.23_rk,0.47_rk,0.71_rk])
        Xg_large = large_scale%cmp_Xg([0.23_rk,0.47_rk,0.71_rk])
        call reference%derivative2([0.23_rk,0.47_rk,0.71_rk], d2T_ref, dT_ref, T_ref)
        call small_scale%derivative2([0.23_rk,0.47_rk,0.71_rk], d2T_small, dT_small, T_small)
        call large_scale%derivative2([0.23_rk,0.47_rk,0.71_rk], d2T_large, dT_large, T_large)
        call reference%err%print()
        call small_scale%err%print()
        call large_scale%err%print()

        call ut%test(ti)%check(&
            name     = "volume projective-scale invariance",&
            res      = reference%err%ok .and. small_scale%err%ok .and. &
                large_scale%err%ok .and. all(ieee_is_finite(Xg_small)) .and. &
                all(ieee_is_finite(Xg_large)) .and. &
                maxval(abs(Xg_small-Xg_ref)) <= tol .and. &
                maxval(abs(Xg_large-Xg_ref)) <= tol .and. &
                maxval(abs(T_small-T_ref)) <= tol .and. &
                maxval(abs(T_large-T_ref)) <= tol .and. &
                maxval(abs(dT_small-dT_ref)) <= tol .and. &
                maxval(abs(dT_large-dT_ref)) <= tol .and. &
                maxval(abs(d2T_small-d2T_ref)) <= tol .and. &
                maxval(abs(d2T_large-d2T_ref)) <= tol,&
            expected = .true.,&
            msg      = "Volume values changed under a common weight scale.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0227


    subroutine forcad_nurbs_volume_0228(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 16384.0_rk*epsilon(1.0_rk)
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(8,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [8,3])
        real(rk), parameter :: Wc(8) = [&
            1.0_rk,2.0_rk,3.0_rk,4.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk]
        integer, parameter :: ngauss(3) = [5,5,5]
        type(nurbs_volume) :: reference, small_scale, large_scale
        real(rk) :: volume_ref, volume_small, volume_large

        call reference%set(knot, knot, knot, Xc, Wc, degree=[1,1,1])
        call small_scale%set(knot, knot, knot, Xc, 1.0e-200_rk*Wc, degree=[1,1,1])
        call large_scale%set(knot, knot, knot, Xc, 1.0e200_rk*Wc, degree=[1,1,1])
        call reference%cmp_volume(volume_ref, ngauss=ngauss)
        call small_scale%cmp_volume(volume_small, ngauss=ngauss)
        call large_scale%cmp_volume(volume_large, ngauss=ngauss)
        call reference%err%print()
        call small_scale%err%print()
        call large_scale%err%print()

        call ut%test(ti)%check(&
            name     = "volume measure weight-scale invariance",&
            res      = reference%err%ok .and. small_scale%err%ok .and. &
                large_scale%err%ok .and. ieee_is_finite(volume_ref) .and. &
                ieee_is_finite(volume_small) .and. ieee_is_finite(volume_large) .and. &
                volume_ref > 0.0_rk .and. abs(volume_small-volume_ref) <= tol .and. &
                abs(volume_large-volume_ref) <= tol,&
            expected = .true.,&
            msg      = "Volume measure changed under a common weight scale.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0228


    subroutine forcad_nurbs_volume_0229(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 5.0e-11_rk
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xt(5,3) = reshape([&
            0.0_rk,0.21_rk,0.5_rk,0.78_rk,1.0_rk,&
            0.0_rk,0.67_rk,0.34_rk,0.91_rk,1.0_rk,&
            0.0_rk,0.43_rk,0.82_rk,0.16_rk,1.0_rk], [5,3])
        real(rk), parameter :: Xc(8,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [8,3])
        real(rk), parameter :: Wc(8) = 1.0e-200_rk*[&
            1.0_rk,2.0_rk,3.0_rk,4.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk]
        type(nurbs_volume) :: volume
        real(rk), allocatable :: Xg_ref(:,:), Xg_insert(:,:), Xg_remove(:,:), Xg_elevate(:,:)
        real(rk), allocatable :: Wc_after(:)

        call volume%set(knot, knot, knot, Xc, Wc, degree=[1,1,1])
        call volume%create(Xt=Xt)
        Xg_ref = volume%get_Xg()
        call volume%insert_knots(1, [0.5_rk], [1])
        call volume%create(Xt=Xt)
        Xg_insert = volume%get_Xg()
        call volume%remove_knots(1, [0.5_rk], [1])
        call volume%create(Xt=Xt)
        Xg_remove = volume%get_Xg()
        call volume%elevate_degree(3, 1)
        call volume%create(Xt=Xt)
        Xg_elevate = volume%get_Xg()
        Wc_after = volume%get_Wc()
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "volume refinement preserves tiny weights",&
            res      = volume%err%ok .and. all(ieee_is_finite(Xg_elevate)) .and. &
                all(ieee_is_finite(Wc_after)) .and. all(Wc_after > 0.0_rk) .and. &
                maxval(Wc_after) < 1.0e-190_rk .and. &
                maxval(abs(Xg_insert-Xg_ref)) <= tol .and. &
                maxval(abs(Xg_remove-Xg_ref)) <= tol .and. &
                maxval(abs(Xg_elevate-Xg_ref)) <= tol,&
            expected = .true.,&
            msg      = "Volume refinement changed geometry or weight gauge.",&
            group    = "forcad_nurbs_volume")
        ti = ti + 1

    end subroutine forcad_nurbs_volume_0229


    subroutine run_nurbs_volume_tests(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call forcad_nurbs_volume_0001(ut, ti)
        call forcad_nurbs_volume_0002(ut, ti)
        call forcad_nurbs_volume_0003(ut, ti)
        call forcad_nurbs_volume_0004(ut, ti)
        call forcad_nurbs_volume_0005(ut, ti)
        call forcad_nurbs_volume_0006(ut, ti)
        call forcad_nurbs_volume_0007(ut, ti)
        call forcad_nurbs_volume_0008(ut, ti)
        call forcad_nurbs_volume_0009(ut, ti)
        call forcad_nurbs_volume_0010(ut, ti)
        call forcad_nurbs_volume_0011(ut, ti)
        call forcad_nurbs_volume_0012(ut, ti)
        call forcad_nurbs_volume_0013(ut, ti)
        call forcad_nurbs_volume_0014(ut, ti)
        call forcad_nurbs_volume_0015(ut, ti)
        call forcad_nurbs_volume_0016(ut, ti)
        call forcad_nurbs_volume_0017(ut, ti)
        call forcad_nurbs_volume_0018(ut, ti)
        call forcad_nurbs_volume_0019(ut, ti)
        call forcad_nurbs_volume_0020(ut, ti)
        call forcad_nurbs_volume_0021(ut, ti)
        call forcad_nurbs_volume_0022(ut, ti)
        call forcad_nurbs_volume_0023(ut, ti)
        call forcad_nurbs_volume_0024(ut, ti)
        call forcad_nurbs_volume_0025(ut, ti)
        call forcad_nurbs_volume_0026(ut, ti)
        call forcad_nurbs_volume_0027(ut, ti)
        call forcad_nurbs_volume_0028(ut, ti)
        call forcad_nurbs_volume_0029(ut, ti)
        call forcad_nurbs_volume_0030(ut, ti)
        call forcad_nurbs_volume_0031(ut, ti)
        call forcad_nurbs_volume_0032(ut, ti)
        call forcad_nurbs_volume_0033(ut, ti)
        call forcad_nurbs_volume_0034(ut, ti)
        call forcad_nurbs_volume_0035(ut, ti)
        call forcad_nurbs_volume_0036(ut, ti)
        call forcad_nurbs_volume_0037(ut, ti)
        call forcad_nurbs_volume_0038(ut, ti)
        call forcad_nurbs_volume_0039(ut, ti)
        call forcad_nurbs_volume_0040(ut, ti)
        call forcad_nurbs_volume_0041(ut, ti)
        call forcad_nurbs_volume_0042(ut, ti)
        call forcad_nurbs_volume_0043(ut, ti)
        call forcad_nurbs_volume_0044(ut, ti)
        call forcad_nurbs_volume_0045(ut, ti)
        call forcad_nurbs_volume_0046(ut, ti)
        call forcad_nurbs_volume_0047(ut, ti)
        call forcad_nurbs_volume_0048(ut, ti)
        call forcad_nurbs_volume_0049(ut, ti)
        call forcad_nurbs_volume_0050(ut, ti)
        call forcad_nurbs_volume_0051(ut, ti)
        call forcad_nurbs_volume_0052(ut, ti)
        call forcad_nurbs_volume_0053(ut, ti)
        call forcad_nurbs_volume_0054(ut, ti)
        call forcad_nurbs_volume_0055(ut, ti)
        call forcad_nurbs_volume_0056(ut, ti)
        call forcad_nurbs_volume_0057(ut, ti)
        call forcad_nurbs_volume_0058(ut, ti)
        call forcad_nurbs_volume_0059(ut, ti)
        call forcad_nurbs_volume_0060(ut, ti)
        call forcad_nurbs_volume_0061(ut, ti)
        call forcad_nurbs_volume_0062(ut, ti)
        call forcad_nurbs_volume_0063(ut, ti)
        call forcad_nurbs_volume_0064(ut, ti)
        call forcad_nurbs_volume_0065(ut, ti)
        call forcad_nurbs_volume_0066(ut, ti)
        call forcad_nurbs_volume_0067(ut, ti)
        call forcad_nurbs_volume_0068(ut, ti)
        call forcad_nurbs_volume_0069(ut, ti)
        call forcad_nurbs_volume_0070(ut, ti)
        call forcad_nurbs_volume_0071(ut, ti)
        call forcad_nurbs_volume_0072(ut, ti)
        call forcad_nurbs_volume_0073(ut, ti)
        call forcad_nurbs_volume_0074(ut, ti)
        call forcad_nurbs_volume_0075(ut, ti)
        call forcad_nurbs_volume_0076(ut, ti)
        call forcad_nurbs_volume_0077(ut, ti)
        call forcad_nurbs_volume_0078(ut, ti)
        call forcad_nurbs_volume_0079(ut, ti)
        call forcad_nurbs_volume_0080(ut, ti)
        call forcad_nurbs_volume_0081(ut, ti)
        call forcad_nurbs_volume_0082(ut, ti)
        call forcad_nurbs_volume_0083(ut, ti)
        call forcad_nurbs_volume_0084(ut, ti)
        call forcad_nurbs_volume_0085(ut, ti)
        call forcad_nurbs_volume_0086(ut, ti)
        call forcad_nurbs_volume_0087(ut, ti)
        call forcad_nurbs_volume_0088(ut, ti)
        call forcad_nurbs_volume_0089(ut, ti)
        call forcad_nurbs_volume_0090(ut, ti)
        call forcad_nurbs_volume_0091(ut, ti)
        call forcad_nurbs_volume_0092(ut, ti)
        call forcad_nurbs_volume_0093(ut, ti)
        call forcad_nurbs_volume_0094(ut, ti)
        call forcad_nurbs_volume_0095(ut, ti)
        call forcad_nurbs_volume_0096(ut, ti)
        call forcad_nurbs_volume_0097(ut, ti)
        call forcad_nurbs_volume_0098(ut, ti)
        call forcad_nurbs_volume_0099(ut, ti)
        call forcad_nurbs_volume_0100(ut, ti)
        call forcad_nurbs_volume_0101(ut, ti)
        call forcad_nurbs_volume_0102(ut, ti)
        call forcad_nurbs_volume_0103(ut, ti)
        call forcad_nurbs_volume_0104(ut, ti)
        call forcad_nurbs_volume_0105(ut, ti)
        call forcad_nurbs_volume_0106(ut, ti)
        call forcad_nurbs_volume_0107(ut, ti)
        call forcad_nurbs_volume_0108(ut, ti)
        call forcad_nurbs_volume_0109(ut, ti)
        call forcad_nurbs_volume_0110(ut, ti)
        call forcad_nurbs_volume_0111(ut, ti)
        call forcad_nurbs_volume_0112(ut, ti)
        call forcad_nurbs_volume_0113(ut, ti)
        call forcad_nurbs_volume_0114(ut, ti)
        call forcad_nurbs_volume_0115(ut, ti)
        call forcad_nurbs_volume_0116(ut, ti)
        call forcad_nurbs_volume_0117(ut, ti)
        call forcad_nurbs_volume_0118(ut, ti)
        call forcad_nurbs_volume_0119(ut, ti)
        call forcad_nurbs_volume_0120(ut, ti)
        call forcad_nurbs_volume_0121(ut, ti)
        call forcad_nurbs_volume_0122(ut, ti)
        call forcad_nurbs_volume_0123(ut, ti)
        call forcad_nurbs_volume_0124(ut, ti)
        call forcad_nurbs_volume_0125(ut, ti)
        call forcad_nurbs_volume_0126(ut, ti)
        call forcad_nurbs_volume_0127(ut, ti)
        call forcad_nurbs_volume_0128(ut, ti)
        call forcad_nurbs_volume_0129(ut, ti)
        call forcad_nurbs_volume_0130(ut, ti)
        call forcad_nurbs_volume_0131(ut, ti)
        call forcad_nurbs_volume_0132(ut, ti)
        call forcad_nurbs_volume_0133(ut, ti)
        call forcad_nurbs_volume_0134(ut, ti)
        call forcad_nurbs_volume_0135(ut, ti)
        call forcad_nurbs_volume_0136(ut, ti)
        call forcad_nurbs_volume_0137(ut, ti)
        call forcad_nurbs_volume_0138(ut, ti)
        call forcad_nurbs_volume_0139(ut, ti)
        call forcad_nurbs_volume_0140(ut, ti)
        call forcad_nurbs_volume_0141(ut, ti)
        call forcad_nurbs_volume_0142(ut, ti)
        call forcad_nurbs_volume_0143(ut, ti)
        call forcad_nurbs_volume_0144(ut, ti)
        call forcad_nurbs_volume_0145(ut, ti)
        call forcad_nurbs_volume_0146(ut, ti)
        call forcad_nurbs_volume_0147(ut, ti)
        call forcad_nurbs_volume_0148(ut, ti)
        call forcad_nurbs_volume_0149(ut, ti)
        call forcad_nurbs_volume_0150(ut, ti)
        call forcad_nurbs_volume_0151(ut, ti)
        call forcad_nurbs_volume_0152(ut, ti)
        call forcad_nurbs_volume_0153(ut, ti)
        call forcad_nurbs_volume_0154(ut, ti)
        call forcad_nurbs_volume_0155(ut, ti)
        call forcad_nurbs_volume_0156(ut, ti)
        call forcad_nurbs_volume_0157(ut, ti)
        call forcad_nurbs_volume_0158(ut, ti)
        call forcad_nurbs_volume_0159(ut, ti)
        call forcad_nurbs_volume_0160(ut, ti)
        call forcad_nurbs_volume_0161(ut, ti)
        call forcad_nurbs_volume_0162(ut, ti)
        call forcad_nurbs_volume_0163(ut, ti)
        call forcad_nurbs_volume_0164(ut, ti)
        call forcad_nurbs_volume_0165(ut, ti)
        call forcad_nurbs_volume_0166(ut, ti)
        call forcad_nurbs_volume_0167(ut, ti)
        call forcad_nurbs_volume_0168(ut, ti)
        call forcad_nurbs_volume_0169(ut, ti)
        call forcad_nurbs_volume_0170(ut, ti)
        call forcad_nurbs_volume_0171(ut, ti)
        call forcad_nurbs_volume_0172(ut, ti)
        call forcad_nurbs_volume_0173(ut, ti)
        call forcad_nurbs_volume_0174(ut, ti)
        call forcad_nurbs_volume_0175(ut, ti)
        call forcad_nurbs_volume_0176(ut, ti)
        call forcad_nurbs_volume_0177(ut, ti)
        call forcad_nurbs_volume_0178(ut, ti)
        call forcad_nurbs_volume_0179(ut, ti)
        call forcad_nurbs_volume_0180(ut, ti)
        call forcad_nurbs_volume_0181(ut, ti)
        call forcad_nurbs_volume_0182(ut, ti)
        call forcad_nurbs_volume_0183(ut, ti)
        call forcad_nurbs_volume_0184(ut, ti)
        call forcad_nurbs_volume_0185(ut, ti)
        call forcad_nurbs_volume_0186(ut, ti)
        call forcad_nurbs_volume_0187(ut, ti)
        call forcad_nurbs_volume_0188(ut, ti)
        call forcad_nurbs_volume_0189(ut, ti)
        call forcad_nurbs_volume_0190(ut, ti)
        call forcad_nurbs_volume_0191(ut, ti)
        call forcad_nurbs_volume_0192(ut, ti)
        call forcad_nurbs_volume_0193(ut, ti)
        call forcad_nurbs_volume_0194(ut, ti)
        call forcad_nurbs_volume_0195(ut, ti)
        call forcad_nurbs_volume_0196(ut, ti)
        call forcad_nurbs_volume_0197(ut, ti)
        call forcad_nurbs_volume_0198(ut, ti)
        call forcad_nurbs_volume_0199(ut, ti)
        call forcad_nurbs_volume_0200(ut, ti)
        call forcad_nurbs_volume_0201(ut, ti)
        call forcad_nurbs_volume_0202(ut, ti)
        call forcad_nurbs_volume_0203(ut, ti)
        call forcad_nurbs_volume_0204(ut, ti)
        call forcad_nurbs_volume_0205(ut, ti)
        call forcad_nurbs_volume_0206(ut, ti)
        call forcad_nurbs_volume_0207(ut, ti)
        call forcad_nurbs_volume_0208(ut, ti)
        call forcad_nurbs_volume_0209(ut, ti)
        call forcad_nurbs_volume_0210(ut, ti)
        call forcad_nurbs_volume_0211(ut, ti)
        call forcad_nurbs_volume_0212(ut, ti)
        call forcad_nurbs_volume_0213(ut, ti)
        call forcad_nurbs_volume_0214(ut, ti)
        call forcad_nurbs_volume_0215(ut, ti)
        call forcad_nurbs_volume_0216(ut, ti)
        call forcad_nurbs_volume_0217(ut, ti)
        call forcad_nurbs_volume_0218(ut, ti)
        call forcad_nurbs_volume_0219(ut, ti)
        call forcad_nurbs_volume_0220(ut, ti)
        call forcad_nurbs_volume_0221(ut, ti)
        call forcad_nurbs_volume_0222(ut, ti)
        call forcad_nurbs_volume_0223(ut, ti)
        call forcad_nurbs_volume_0224(ut, ti)
        call forcad_nurbs_volume_0225(ut, ti)
        call forcad_nurbs_volume_0226(ut, ti)
        call forcad_nurbs_volume_0227(ut, ti)
        call forcad_nurbs_volume_0228(ut, ti)
        call forcad_nurbs_volume_0229(ut, ti)

    end subroutine run_nurbs_volume_tests

end module forcad_test_nurbs_volume
