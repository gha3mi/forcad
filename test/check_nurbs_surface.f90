module forcad_test_nurbs_surface

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use forcad_kinds, only: rk
    use forcad_nurbs_surface, only: nurbs_surface
    use forcad_utils, only: ndgrid, rotation, linspace, kron_eye, active_knots, active_span_count
    use forunittest, only: unit_tests

    implicit none

    private
    public :: run_nurbs_surface_tests

contains

    subroutine forcad_nurbs_surface_0001(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS surface area",&
            res      = area,&
            expected = 25.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "The NURBS surface area result does not match the expected value.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0001


    subroutine forcad_nurbs_surface_0002(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline surface area",&
            res      = areab,&
            expected = 25.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Surface area is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0002


    subroutine forcad_nurbs_surface_0003(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS sampled nearest-point geometry",&
            res      = nearest_Xg,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Sampled nearest-point geometry is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0003


    subroutine forcad_nurbs_surface_0004(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS sampled nearest-point parameters",&
            res      = nearest_Xt,&
            expected = [0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Sampled nearest-point parameters are incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0004


    subroutine forcad_nurbs_surface_0005(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS sampled nearest-point index",&
            res      = id,&
            expected = 1,&
            msg      = "Sampled nearest-point index is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0005


    subroutine forcad_nurbs_surface_0006(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline sampled nearest-point geometry",&
            res      = nearest_Xg,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Sampled nearest-point geometry is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0006


    subroutine forcad_nurbs_surface_0007(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline sampled nearest-point parameters",&
            res      = nearest_Xt,&
            expected = [0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Sampled nearest-point parameters are incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0007


    subroutine forcad_nurbs_surface_0008(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline sampled nearest-point index",&
            res      = id,&
            expected = 1,&
            msg      = "Sampled nearest-point index is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0008


    subroutine forcad_nurbs_surface_0009(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS iterative nearest-point geometry",&
            res      = nearest_Xg,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Iterative nearest-point geometry is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0009


    subroutine forcad_nurbs_surface_0010(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS iterative nearest-point parameters",&
            res      = nearest_Xt,&
            expected = [0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Iterative nearest-point parameters are incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0010


    subroutine forcad_nurbs_surface_0011(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline iterative nearest-point geometry",&
            res      = nearest_Xg,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Iterative nearest-point geometry is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0011


    subroutine forcad_nurbs_surface_0012(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:)
        real(rk), allocatable :: Wc(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline iterative nearest-point parameters",&
            res      = nearest_Xt,&
            expected = [0.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Iterative nearest-point parameters are incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0012


    subroutine forcad_nurbs_surface_0013(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS explicit-knot setter geometry",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Explicit-knot setter geometry is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0013


    subroutine forcad_nurbs_surface_0014(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline explicit-knot setter geometry",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Explicit-knot setter geometry is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0014


    subroutine forcad_nurbs_surface_0015(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS degree-control setter geometry",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Degree-control setter geometry is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0015


    subroutine forcad_nurbs_surface_0016(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline degree-control setter geometry",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Degree-control setter geometry is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0016


    subroutine forcad_nurbs_surface_0017(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS explicit-grid geometry recreation",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Explicit-grid geometry recreation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0017


    subroutine forcad_nurbs_surface_0018(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline explicit-grid geometry recreation",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Explicit-grid geometry recreation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0018


    subroutine forcad_nurbs_surface_0019(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS control-net getter",&
            res      = nurbs%get_Xc(),&
            expected = Xc,&
            tol      = 1e-5_rk,&
            msg      = "Control-net getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0019


    subroutine forcad_nurbs_surface_0020(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline control-net getter",&
            res      = bsp%get_Xc(),&
            expected = Xc,&
            tol      = 1e-5_rk,&
            msg      = "Control-net getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0020


    subroutine forcad_nurbs_surface_0021(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS control-point getter",&
            res      = nurbs%get_Xc(1),&
            expected = Xc(1,:),&
            tol      = 1e-5_rk,&
            msg      = "Control-point getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0021


    subroutine forcad_nurbs_surface_0022(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline control-point getter",&
            res      = bsp%get_Xc(1),&
            expected = Xc(1,:),&
            tol      = 1e-5_rk,&
            msg      = "Control-point getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0022


    subroutine forcad_nurbs_surface_0023(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS control-coordinate getter",&
            res      = nurbs%get_Xc(1,1),&
            expected = Xc(1,1),&
            tol      = 1e-5_rk,&
            msg      = "Control-coordinate getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0023


    subroutine forcad_nurbs_surface_0024(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline control-coordinate getter",&
            res      = bsp%get_Xc(1,1),&
            expected = Xc(1,1),&
            tol      = 1e-5_rk,&
            msg      = "Control-coordinate getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0024


    subroutine forcad_nurbs_surface_0025(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS geometry getter",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Geometry getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0025


    subroutine forcad_nurbs_surface_0026(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline geometry getter",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Geometry getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0026


    subroutine forcad_nurbs_surface_0027(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS geometry-point getter",&
            res      = nurbs%get_Xg(1),&
            expected = Xg(1,:),&
            tol      = 1e-5_rk,&
            msg      = "Geometry-point getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0027


    subroutine forcad_nurbs_surface_0028(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline geometry-point getter",&
            res      = bsp%get_Xg(1),&
            expected = Xgb(1,:),&
            tol      = 1e-5_rk,&
            msg      = "Geometry-point getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0028


    subroutine forcad_nurbs_surface_0029(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS geometry-coordinate getter",&
            res      = nurbs%get_Xg(1,1),&
            expected = Xg(1,1),&
            tol      = 1e-5_rk,&
            msg      = "Geometry-coordinate getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0029


    subroutine forcad_nurbs_surface_0030(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline geometry-coordinate getter",&
            res      = bsp%get_Xg(1,1),&
            expected = Xgb(1,1),&
            tol      = 1e-5_rk,&
            msg      = "Geometry-coordinate getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0030


    subroutine forcad_nurbs_surface_0031(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS weight-vector getter",&
            res      = nurbs%get_Wc(),&
            expected = Wc,&
            tol      = 1e-5_rk,&
            msg      = "Weight-vector getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0031


    subroutine forcad_nurbs_surface_0032(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS weight getter",&
            res      = nurbs%get_Wc(1),&
            expected = Wc(1),&
            tol      = 1e-5_rk,&
            msg      = "Weight getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0032


    subroutine forcad_nurbs_surface_0033(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS directional knot-vector getter",&
            res      = nurbs%get_knot(1),&
            expected = knot1,&
            tol      = 1e-5_rk,&
            msg      = "Directional knot-vector getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0033


    subroutine forcad_nurbs_surface_0034(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline directional knot-vector getter",&
            res      = bsp%get_knot(1),&
            expected = knot1,&
            tol      = 1e-5_rk,&
            msg      = "Directional knot-vector getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0034


    subroutine forcad_nurbs_surface_0035(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS knot getter",&
            res      = nurbs%get_knot(1,1),&
            expected = knot1(1),&
            tol      = 1e-5_rk,&
            msg      = "The NURBS knot getter result does not match the expected value.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0035


    subroutine forcad_nurbs_surface_0036(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline knot getter",&
            res      = bsp%get_knot(1,1),&
            expected = knot1(1),&
            tol      = 1e-5_rk,&
            msg      = "Knot getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0036


    subroutine forcad_nurbs_surface_0037(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS directional degree getter",&
            res      = nurbs%get_degree(1),&
            expected = 1,&
            msg      = "Directional degree getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0037


    subroutine forcad_nurbs_surface_0038(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline directional degree getter",&
            res      = bsp%get_degree(1),&
            expected = 1,&
            msg      = "Directional degree getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0038


    subroutine forcad_nurbs_surface_0039(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS knot multiplicity getter",&
            res      = nurbs%get_multiplicity(1),&
            expected = [2,2],&
            msg      = "Knot multiplicity getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0039


    subroutine forcad_nurbs_surface_0040(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline knot multiplicity getter",&
            res      = bsp%get_multiplicity(1),&
            expected = [2,2],&
            msg      = "Knot multiplicity getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0040


    subroutine forcad_nurbs_surface_0041(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS knot continuity getter",&
            res      = nurbs%get_continuity(1),&
            expected = [-1,-1],&
            msg      = "Knot continuity getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0041


    subroutine forcad_nurbs_surface_0042(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline knot continuity getter",&
            res      = bsp%get_continuity(1),&
            expected = [-1,-1],&
            msg      = "Knot continuity getter is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0042


    subroutine forcad_nurbs_surface_0043(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS control-net shape",&
            res      = nurbs%get_nc(),&
            expected = [2, 2],&
            msg      = "Control-net shape is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0043


    subroutine forcad_nurbs_surface_0044(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline control-net shape",&
            res      = bsp%get_nc(),&
            expected = [2, 2],&
            msg      = "Control-net shape is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0044


    subroutine forcad_nurbs_surface_0045(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS recomputed control-net shape",&
            res      = nurbs%get_nc(),&
            expected = [2, 2],&
            msg      = "Recomputed control-net shape is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0045


    subroutine forcad_nurbs_surface_0046(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline recomputed control-net shape",&
            res      = bsp%get_nc(),&
            expected = [2, 2],&
            msg      = "Recomputed control-net shape is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0046


    subroutine forcad_nurbs_surface_0047(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS custom control-net visualization connectivity",&
            res      = nurbs%get_elem_Xc_vis(),&
            expected = elemConn,&
            msg      = "Custom control-net visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0047


    subroutine forcad_nurbs_surface_0048(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0048


    subroutine forcad_nurbs_surface_0049(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline custom control-net visualization connectivity",&
            res      = bsp%get_elem_Xc_vis(),&
            expected = elemConn,&
            msg      = "Custom control-net visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0049


    subroutine forcad_nurbs_surface_0050(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0050


    subroutine forcad_nurbs_surface_0051(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS custom geometry visualization connectivity",&
            res      = nurbs%get_elem_Xg_vis(),&
            expected = elemConn,&
            msg      = "Custom geometry visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0051


    subroutine forcad_nurbs_surface_0052(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0052


    subroutine forcad_nurbs_surface_0053(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
        call bsp%set_elem_Xg_vis(elemConn)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline custom geometry visualization connectivity",&
            res      = bsp%get_elem_Xg_vis(),&
            expected = elemConn,&
            msg      = "Custom geometry visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0053


    subroutine forcad_nurbs_surface_0054(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0054


    subroutine forcad_nurbs_surface_0055(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0055


    subroutine forcad_nurbs_surface_0056(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0056


    subroutine forcad_nurbs_surface_0057(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS modified-control geometry recreation",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Modified-control geometry recreation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0057


    subroutine forcad_nurbs_surface_0058(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline modified-control geometry recreation",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Modified-control geometry recreation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0058


    subroutine forcad_nurbs_surface_0059(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:)
        real(rk), allocatable :: Tgc1(:)
        real(rk), allocatable :: Tgc1b(:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS boundary active basis values",&
            res      = Tgc1e,&
            expected = Tgc1(1:3),&
            tol      = 1e-5_rk,&
            msg      = "Boundary active basis values are incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0059


    subroutine forcad_nurbs_surface_0060(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:)
        real(rk), allocatable :: Tgc1(:)
        real(rk), allocatable :: Tgc1b(:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline boundary active basis values",&
            res      = Tgc1be,&
            expected = Tgc1b(1:3),&
            tol      = 1e-5_rk,&
            msg      = "Boundary active basis values are incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0060


    subroutine forcad_nurbs_surface_0061(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:)
        real(rk), allocatable :: Tgc1(:)
        real(rk), allocatable :: Tgc1b(:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS boundary element basis selection",&
            res      = maxval(abs(Tgc1e - [Tgc1(1), Tgc1(4)])),&
            expected = 0.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Boundary element basis selection is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0061


    subroutine forcad_nurbs_surface_0062(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:)
        real(rk), allocatable :: Tgc1(:)
        real(rk), allocatable :: Tgc1b(:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline boundary element basis selection",&
            res      = maxval(abs(Tgc1be - [Tgc1b(1), Tgc1b(4)])),&
            expected = 0.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Boundary element basis selection is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0062


    subroutine forcad_nurbs_surface_0063(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0063


    subroutine forcad_nurbs_surface_0064(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0064


    subroutine forcad_nurbs_surface_0065(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0065


    subroutine forcad_nurbs_surface_0066(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0066


    subroutine forcad_nurbs_surface_0067(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0067


    subroutine forcad_nurbs_surface_0068(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0068


    subroutine forcad_nurbs_surface_0069(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0069


    subroutine forcad_nurbs_surface_0070(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0070


    subroutine forcad_nurbs_surface_0071(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0071


    subroutine forcad_nurbs_surface_0072(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0072


    subroutine forcad_nurbs_surface_0073(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0073


    subroutine forcad_nurbs_surface_0074(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0074


    subroutine forcad_nurbs_surface_0075(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0075


    subroutine forcad_nurbs_surface_0076(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0076


    subroutine forcad_nurbs_surface_0077(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0077


    subroutine forcad_nurbs_surface_0078(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0078


    subroutine forcad_nurbs_surface_0079(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS knot insertion geometry preservation",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Knot insertion geometry preservation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0079


    subroutine forcad_nurbs_surface_0080(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline knot insertion geometry preservation",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Knot insertion geometry preservation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0080


    subroutine forcad_nurbs_surface_0081(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%elevate_degree(1, 2)
        call nurbs%elevate_degree(2, 2)

        call bsp%elevate_degree(1, 2)
        call bsp%elevate_degree(2, 2)

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS degree elevation geometry preservation",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Degree elevation geometry preservation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0081


    subroutine forcad_nurbs_surface_0082(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%elevate_degree(1, 2)
        call nurbs%elevate_degree(2, 2)

        call bsp%elevate_degree(1, 2)
        call bsp%elevate_degree(2, 2)

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline degree elevation geometry preservation",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Degree elevation geometry preservation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0082


    subroutine forcad_nurbs_surface_0083(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%elevate_degree(1, 2)
        call nurbs%elevate_degree(2, 2)

        call bsp%elevate_degree(1, 2)
        call bsp%elevate_degree(2, 2)

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call bsp%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call bsp%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS knot removal geometry preservation",&
            res      = nurbs%get_Xg(),&
            expected = Xg,&
            tol      = 1e-5_rk,&
            msg      = "Knot removal geometry preservation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0083


    subroutine forcad_nurbs_surface_0084(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%elevate_degree(1, 2)
        call nurbs%elevate_degree(2, 2)

        call bsp%elevate_degree(1, 2)
        call bsp%elevate_degree(2, 2)

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call bsp%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call bsp%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline knot removal geometry preservation",&
            res      = bsp%get_Xg(),&
            expected = Xgb,&
            tol      = 1e-5_rk,&
            msg      = "Knot removal geometry preservation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0084


    subroutine forcad_nurbs_surface_0085(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%elevate_degree(1, 2)
        call nurbs%elevate_degree(2, 2)

        call bsp%elevate_degree(1, 2)
        call bsp%elevate_degree(2, 2)

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call bsp%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call bsp%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%set_tetragon([2.0_rk, 2.0_rk], [2,2])
        call bsp%set_tetragon([2.0_rk, 2.0_rk], [2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk])
        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)

        call nurbs%set_tetragon([2.0_rk, 3.0_rk], [3,2])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "tetragon surface control-net shape",&
            res      = nurbs%get_nc(),&
            expected = [3,2],&
            msg      = "Tetragon surface control-net shape is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0085


    subroutine forcad_nurbs_surface_0086(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%elevate_degree(1, 2)
        call nurbs%elevate_degree(2, 2)

        call bsp%elevate_degree(1, 2)
        call bsp%elevate_degree(2, 2)

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call bsp%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call bsp%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%set_tetragon([2.0_rk, 2.0_rk], [2,2])
        call bsp%set_tetragon([2.0_rk, 2.0_rk], [2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk])
        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)

        call nurbs%set_tetragon([2.0_rk, 3.0_rk], [3,2])
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "tetragon surface final control point",&
            res      = nurbs%get_Xc(6),&
            expected = [2.0_rk, 3.0_rk, 0.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Tetragon surface final control point is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0086


    subroutine forcad_nurbs_surface_0087(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%elevate_degree(1, 2)
        call nurbs%elevate_degree(2, 2)

        call bsp%elevate_degree(1, 2)
        call bsp%elevate_degree(2, 2)

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call bsp%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call bsp%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%set_tetragon([2.0_rk, 2.0_rk], [2,2])
        call bsp%set_tetragon([2.0_rk, 2.0_rk], [2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk])
        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)

        call nurbs%set_tetragon([2.0_rk, 3.0_rk], [3,2])

        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "ring surface control-net shape",&
            res      = nurbs%get_nc(),&
            expected = [7,2],&
            msg      = "Ring surface control-net shape is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0087


    subroutine forcad_nurbs_surface_0088(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%elevate_degree(1, 2)
        call nurbs%elevate_degree(2, 2)

        call bsp%elevate_degree(1, 2)
        call bsp%elevate_degree(2, 2)

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call bsp%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call bsp%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%set_tetragon([2.0_rk, 2.0_rk], [2,2])
        call bsp%set_tetragon([2.0_rk, 2.0_rk], [2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk])
        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)

        call nurbs%set_tetragon([2.0_rk, 3.0_rk], [3,2])

        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "ring surface rational state",&
            res      = nurbs%is_rational(),&
            expected = .true.,&
            msg      = "Ring surface rational state is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0088


    subroutine forcad_nurbs_surface_0089(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: i
        type(nurbs_surface) :: nurbs, bsp
        real(rk) :: knot1(4), knot2(4), area, areab
        real(rk), allocatable :: Xc(:,:), Wc(:), Xg(:,:), Xgb(:,:)
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Tgc(:,:), dTgc(:,:,:), Tgcb(:,:), dTgcb(:,:,:), d2Tgc(:,:,:), d2Tgcb(:,:,:)
        real(rk), allocatable :: Tgc1(:), dTgc1(:,:), Tgc1b(:), dTgc1b(:,:), d2Tgc1(:,:), d2Tgc1b(:,:)
        real(rk), allocatable :: Tgc1e(:), Tgc1be(:)
        real(rk) :: nearest_Xg(3), nearest_Xt(2)
        integer :: id

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

        call nurbs%create(20, 20)
        call bsp%create(20, 20)

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_area(area)
        call bsp%cmp_area(areab)

        call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)

        call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1e-10_rk, 10, nearest_Xt, nearest_Xg)

        Xg = nurbs%get_Xg()
        Xgb = bsp%get_Xg()

        call nurbs%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc, Wc)
        call bsp%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1,1], [-1, -1], [-1, -1], Xc)

        call nurbs%set([2,2], Xc, Wc)
        call bsp%set([2,2], Xc)

        call nurbs%create(Xt1 = nurbs%get_Xt(1), Xt2 = nurbs%get_Xt(2))
        call bsp%create(Xt1 = bsp%get_Xt(1), Xt2 = bsp%get_Xt(2))

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%cmp_nc()
        call bsp%cmp_nc()

        elemConn = nurbs%cmp_elem_Xc_vis([1,1])
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xc_vis()
        call nurbs%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis([1,1])
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xc_vis()
        call bsp%set_elem_Xc_vis(elemConn)
        deallocate(elemConn)

        elemConn = nurbs%cmp_elem_Xg_vis([1,1])
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = nurbs%cmp_elem_Xg_vis()
        call nurbs%set_elem_Xg_vis(elemConn)
        deallocate(elemConn)
        elemConn = bsp%cmp_elem_Xg_vis([1,1])
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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%basis(res1=20, res2=20, Tgc=Tgc)
        call bsp%basis(res1=20, res2=20, Tgc=Tgc)

        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1e, elem=[1,2,3])
        call bsp%basis(Xt=[0.0_rk, 0.0_rk], Tgc=Tgc1be, elem=[1,2,3])
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1)
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1b)
        call nurbs%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1e, elem=[1,4])
        call bsp%basis(Xt=[0.25_rk, 0.5_rk], Tgc=Tgc1be, elem=[1,4])

        call nurbs%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)
        call bsp%basis(Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Tgc=Tgc)

        call nurbs%derivative(res1=20, res2=20, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(res1=20, res2=20, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
        call bsp%derivative(Xt=[0.0_rk,0.0_rk], dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])

        call nurbs%derivative2(res1=20, res2=20, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(res1=20, res2=20, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
        call bsp%derivative2(&
            Xt1=[(real(i-1, rk) / real(20-1, rk), i=1, 20)], Xt2=[(real(i-1, rk) / real(20-1, rk), i=1, 20)],&
            d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

        call nurbs%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
        call bsp%derivative2(Xt=[0.0_rk,0.0_rk], d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)

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

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%elevate_degree(1, 2)
        call nurbs%elevate_degree(2, 2)

        call bsp%elevate_degree(1, 2)
        call bsp%elevate_degree(2, 2)

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call nurbs%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call bsp%remove_knots(1, [0.25_rk, 0.75_rk], [2,1])
        call bsp%remove_knots(2, [0.25_rk, 0.75_rk], [2,1])

        call nurbs%create()
        call bsp%create()

        call nurbs%export_Xc("vtk/test_nurbs_surface_Xc.vtk")
        call bsp%export_Xc("vtk/test_bsp_surface_Xc.vtk")

        call nurbs%export_Xg("vtk/test_nurbs_surface_Xg.vtk")
        call bsp%export_Xg("vtk/test_bsp_surface_Xg.vtk")

        call nurbs%export_Xth("vtk/test_nurbs_surface_Xth.vtk")
        call bsp%export_Xth("vtk/test_bsp_surface_Xth.vtk")

        call nurbs%export_iges("iges/test_nurbs_surface.iges")
        call bsp%export_iges("iges/test_bsp_surface.iges")

        call nurbs%set_tetragon([2.0_rk, 2.0_rk], [2,2])
        call bsp%set_tetragon([2.0_rk, 2.0_rk], [2,2], [1.0_rk,1.0_rk,0.9_rk,0.9_rk])
        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)

        call nurbs%set_tetragon([2.0_rk, 3.0_rk], [3,2])

        call nurbs%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)

        call nurbs%set_C([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk)
        call nurbs%err%print()
        call bsp%err%print()

        call ut%test(ti)%check(&
            name     = "C-shaped surface control-net shape",&
            res      = nurbs%get_nc(),&
            expected = [5,2],&
            msg      = "C-shaped surface control-net shape is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0089


    subroutine forcad_nurbs_surface_0090(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: bsp_fit
        integer :: j, n(2), ndata
        real(rk), parameter :: pi = acos(-1.0_rk)
        real(rk), allocatable :: Xt1(:), Xt2(:), Xt(:,:), Xdata(:,:), Xg_eval(:,:)
        real(rk) :: err1, err2, err3, rms

        n = [6,6]
        allocate(Xt1(n(1)), Xt2(n(2)))
        do concurrent (j = 1: n(1))
            Xt1(j) = real(j-1, rk)/real(n(1)-1, rk)
        end do
        do concurrent (j = 1: n(2))
            Xt2(j) = real(j-1, rk)/real(n(2)-1, rk)
        end do
        call ndgrid(Xt1, Xt2, Xt)

        ndata = n(1)*n(2)
        allocate(Xdata(ndata, 3))
        do j = 1, ndata
            Xdata(j,1) = Xt(j,1)
            Xdata(j,2) = Xt(j,2)
            Xdata(j,3) = 0.1_rk * sin(2.0_rk*pi*Xt(j,1)) * cos(2.0_rk*pi*Xt(j,2))
        end do

        call bsp_fit%set(&
            degree      = [2, 2],&
            Xth_dir1    = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk],&
            Xth_dir2    = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk],&
            continuity1 = [ -1   ,   1    ,   1   ,   1    ,  -1   ],&
            continuity2 = [ -1   ,   1    ,   1   ,   1    ,  -1   ])

        call bsp_fit%lsq_fit_bspline(Xt, Xdata, n)
        call bsp_fit%create(n(1), n(2))
        Xg_eval = bsp_fit%get_Xg()

        err1 = norm2(Xg_eval(:,1) - Xdata(:,1)) / norm2(Xdata(:,1))
        err2 = norm2(Xg_eval(:,2) - Xdata(:,2)) / norm2(Xdata(:,2))
        err3 = norm2(Xg_eval(:,3) - Xdata(:,3)) / norm2(Xdata(:,3))
        rms  = sqrt((err1**2 + err2**2 + err3**2)/3.0_rk)
        call bsp_fit%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline surface least-squares fit",&
            res      = rms,&
            expected = 0.0_rk,&
            tol      = 1e-6_rk,&
            msg      = "Surface least-squares fit is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0090


    subroutine forcad_nurbs_surface_0091(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nonopen
        real(rk), allocatable :: Xc_nonopen(:,:)
        real(rk), allocatable :: Xt_nonopen(:)
        real(rk) :: knot_nonopen(7)
        integer :: ii, jj, idx

        allocate(Xc_nonopen(16, 3), source=0.0_rk)
        do jj = 1, 4
            do ii = 1, 4
                idx = (jj - 1)*4 + ii
                Xc_nonopen(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), 0.0_rk]
            end do
        end do
        knot_nonopen = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call nonopen%set(knot_nonopen, knot_nonopen, Xc_nonopen, degree=[2, 2])
        call nonopen%create(5, 5)
        Xt_nonopen = nonopen%get_Xt(1)
        call nonopen%err%print()

        call ut%test(ti)%check(&
            name     = "unclamped surface degree",&
            res      = nonopen%get_degree(),&
            expected = [2, 2],&
            msg      = "Unclamped surface degree is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0091


    subroutine forcad_nurbs_surface_0092(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nonopen
        real(rk), allocatable :: Xc_nonopen(:,:)
        real(rk), allocatable :: Xt_nonopen(:)
        real(rk) :: knot_nonopen(7)
        integer :: ii, jj, idx

        allocate(Xc_nonopen(16, 3), source=0.0_rk)
        do jj = 1, 4
            do ii = 1, 4
                idx = (jj - 1)*4 + ii
                Xc_nonopen(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), 0.0_rk]
            end do
        end do
        knot_nonopen = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call nonopen%set(knot_nonopen, knot_nonopen, Xc_nonopen, degree=[2, 2])
        call nonopen%create(5, 5)
        Xt_nonopen = nonopen%get_Xt(1)
        call nonopen%err%print()

        call ut%test(ti)%check(&
            name     = "unclamped surface control-net shape",&
            res      = nonopen%get_nc(),&
            expected = [4, 4],&
            msg      = "Unclamped surface control-net shape is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0092


    subroutine forcad_nurbs_surface_0093(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nonopen
        real(rk), allocatable :: Xc_nonopen(:,:)
        real(rk), allocatable :: Xt_nonopen(:)
        real(rk) :: knot_nonopen(7)
        integer :: ii, jj, idx

        allocate(Xc_nonopen(16, 3), source=0.0_rk)
        do jj = 1, 4
            do ii = 1, 4
                idx = (jj - 1)*4 + ii
                Xc_nonopen(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), 0.0_rk]
            end do
        end do
        knot_nonopen = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call nonopen%set(knot_nonopen, knot_nonopen, Xc_nonopen, degree=[2, 2])
        call nonopen%create(5, 5)
        Xt_nonopen = nonopen%get_Xt(1)
        call nonopen%err%print()

        call ut%test(ti)%check(&
            name     = "unclamped surface active domain",&
            res      = [Xt_nonopen(1), Xt_nonopen(size(Xt_nonopen))],&
            expected = [2.0_rk, 4.0_rk],&
            tol      = 1e-5_rk,&
            msg      = "Unclamped surface active domain is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0093


    subroutine forcad_nurbs_surface_0094(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nonopen
        real(rk), allocatable :: Xc_nonopen(:,:)
        real(rk), allocatable :: Xt_nonopen(:)
        real(rk), allocatable :: Tgc_nonopen(:)
        real(rk) :: knot_nonopen(7), pt(2)
        integer :: ii, jj, idx

        allocate(Xc_nonopen(16, 3), source=0.0_rk)
        do jj = 1, 4
            do ii = 1, 4
                idx = (jj - 1)*4 + ii
                Xc_nonopen(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), 0.0_rk]
            end do
        end do
        knot_nonopen = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call nonopen%set(knot_nonopen, knot_nonopen, Xc_nonopen, degree=[2, 2])
        call nonopen%create(5, 5)
        Xt_nonopen = nonopen%get_Xt(1)

        pt = [2.5_rk, 3.5_rk]
        call nonopen%basis(pt, Tgc_nonopen)
        call nonopen%err%print()

        call ut%test(ti)%check(&
            name     = "unclamped surface partition of unity",&
            res      = sum(Tgc_nonopen),&
            expected = 1.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Unclamped surface partition of unity is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0094


    subroutine forcad_nurbs_surface_0095(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nonopen
        real(rk), allocatable :: Xc_nonopen(:,:)
        real(rk), allocatable :: Xt_nonopen(:)
        real(rk), allocatable :: Tgc_nonopen(:)
        real(rk) :: knot_nonopen(7), pt(2)
        integer :: ii, jj, idx

        allocate(Xc_nonopen(16, 3), source=0.0_rk)
        do jj = 1, 4
            do ii = 1, 4
                idx = (jj - 1)*4 + ii
                Xc_nonopen(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), 0.0_rk]
            end do
        end do
        knot_nonopen = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call nonopen%set(knot_nonopen, knot_nonopen, Xc_nonopen, degree=[2, 2])
        call nonopen%create(5, 5)
        Xt_nonopen = nonopen%get_Xt(1)

        pt = [2.5_rk, 3.5_rk]
        call nonopen%basis(pt, Tgc_nonopen)

        pt = [1.5_rk, 2.5_rk]
        call nonopen%basis(pt, Tgc_nonopen)
        call nonopen%err%print()

        call ut%test(ti)%check(&
            name     = "inactive unclamped surface span basis",&
            res      = sum(Tgc_nonopen),&
            expected = 0.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Inactive unclamped surface span basis is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0095


    subroutine forcad_nurbs_surface_0096(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: nonopen
        real(rk), allocatable :: Xc_nonopen(:,:), Xt_nonopen(:), Tgc_nonopen(:), dTgc_nonopen(:,:)
        real(rk) :: knot_nonopen(7), pt(2)
        integer :: ii, jj, idx

        allocate(Xc_nonopen(16, 3), source=0.0_rk)
        do jj = 1, 4
            do ii = 1, 4
                idx = (jj - 1)*4 + ii
                Xc_nonopen(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), 0.0_rk]
            end do
        end do
        knot_nonopen = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call nonopen%set(knot_nonopen, knot_nonopen, Xc_nonopen, degree=[2, 2])
        call nonopen%create(5, 5)
        Xt_nonopen = nonopen%get_Xt(1)

        pt = [2.5_rk, 3.5_rk]
        call nonopen%basis(pt, Tgc_nonopen)

        pt = [1.5_rk, 2.5_rk]
        call nonopen%basis(pt, Tgc_nonopen)

        pt = [2.5_rk, 3.5_rk]
        call nonopen%derivative(pt, dTgc_nonopen, Tgc_nonopen, elem=[0, 17])
        call nonopen%err%print()

        call ut%test(ti)%check(&
            name     = "out-of-domain unclamped surface derivatives",&
            res      = max(maxval(abs(Tgc_nonopen)), maxval(abs(dTgc_nonopen))),&
            expected = 0.0_rk,&
            tol      = 1e-5_rk,&
            msg      = "Out-of-domain unclamped surface derivatives are incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0096


    subroutine forcad_nurbs_surface_0097(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: rational
        real(rk), allocatable :: Xc_out(:,:), Wc_out(:), Xg_out(:), Tgc_out(:)
        real(rk) :: kout(4)

        allocate(Xc_out(4, 3), Wc_out(4))
        Xc_out(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xc_out(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xc_out(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xc_out(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Wc_out = [1.0_rk, 0.9_rk, 1.0_rk, 1.0_rk]
        kout = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

        call rational%set(kout, kout, Xc_out, Wc_out)
        Xg_out = rational%cmp_Xg([-0.25_rk, 0.0_rk])
        call rational%basis([-0.25_rk, 0.0_rk], Tgc_out)
        call rational%err%print()

        call ut%test(ti)%check(&
            name     = "out-of-domain surface geometry",&
            res      = all(abs(Xg_out) <= 1e-5_rk),&
            expected = .true.,&
            msg      = "Out-of-domain surface geometry is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0097


    subroutine forcad_nurbs_surface_0098(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: rational
        real(rk), allocatable :: Xc_out(:,:), Wc_out(:), Xg_out(:), Tgc_out(:)
        real(rk) :: kout(4)

        allocate(Xc_out(4, 3), Wc_out(4))
        Xc_out(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xc_out(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xc_out(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xc_out(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Wc_out = [1.0_rk, 0.9_rk, 1.0_rk, 1.0_rk]
        kout = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

        call rational%set(kout, kout, Xc_out, Wc_out)
        Xg_out = rational%cmp_Xg([-0.25_rk, 0.0_rk])
        call rational%basis([-0.25_rk, 0.0_rk], Tgc_out)
        call rational%err%print()

        call ut%test(ti)%check(&
            name     = "out-of-domain surface basis",&
            res      = all(abs(Tgc_out) <= 1e-5_rk),&
            expected = .true.,&
            msg      = "Out-of-domain surface basis is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0098


    subroutine forcad_nurbs_surface_0099(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: repeated
        real(rk), allocatable :: Xc_repeated(:,:)
        real(rk), allocatable :: Xt_repeated(:)
        real(rk) :: knot_repeated(9)
        integer :: ii, jj, idx

        allocate(Xc_repeated(36, 3), source=0.0_rk)
        do jj = 1, 6
            do ii = 1, 6
                idx = ii + (jj - 1)*6
                Xc_repeated(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), 0.1_rk*real((ii - 1)*(jj - 1), rk)]
            end do
        end do
        knot_repeated = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 3.0_rk, 3.0_rk, 4.0_rk, 6.0_rk, 9.0_rk]

        call repeated%set(knot_repeated, knot_repeated, Xc_repeated, degree=[2, 2])
        call repeated%create(5, 5)
        Xt_repeated = repeated%get_Xt(1)
        call repeated%err%print()

        call ut%test(ti)%check(&
            name     = "repeated-knot surface active domain",&
            res      = [Xt_repeated(1), Xt_repeated(size(Xt_repeated))],&
            expected = [2.0_rk, 4.0_rk],&
            tol      = 1.0e-10_rk,&
            msg      = "Repeated-knot surface active domain is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0099


    subroutine forcad_nurbs_surface_0100(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: repeated
        real(rk), allocatable :: Xc_repeated(:,:), Xt_repeated(:), Tgc_repeated(:)
        real(rk) :: knot_repeated(9), pt(2)
        integer :: ii, jj, idx

        allocate(Xc_repeated(36, 3), source=0.0_rk)
        do jj = 1, 6
            do ii = 1, 6
                idx = ii + (jj - 1)*6
                Xc_repeated(idx,:) = [real(ii - 1, rk), real(jj - 1, rk), 0.1_rk*real((ii - 1)*(jj - 1), rk)]
            end do
        end do
        knot_repeated = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 3.0_rk, 3.0_rk, 4.0_rk, 6.0_rk, 9.0_rk]

        call repeated%set(knot_repeated, knot_repeated, Xc_repeated, degree=[2, 2])
        call repeated%create(5, 5)
        Xt_repeated = repeated%get_Xt(1)

        pt = [2.5_rk, 3.25_rk]
        call repeated%basis(pt, Tgc_repeated)
        call repeated%err%print()

        call ut%test(ti)%check(&
            name     = "repeated-knot surface partition of unity",&
            res      = sum(Tgc_repeated),&
            expected = 1.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "Repeated-knot surface partition of unity is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0100


    subroutine forcad_nurbs_surface_0101(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface_ref, surface_scaled
        real(rk) :: knot_scale(4), pt(2)
        real(rk), allocatable :: Xc_scale(:,:), W_scale(:)
        real(rk), allocatable :: Xg_ref(:), Xg_scaled(:), T_ref(:), T_scaled(:), dT_ref(:,:), dT_scaled(:,:)

        allocate(Xc_scale(4, 3), W_scale(4))
        Xc_scale(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xc_scale(2,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xc_scale(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xc_scale(4,:) = [2.0_rk, 1.0_rk, 0.5_rk]
        W_scale = [1.0_rk, 0.7_rk, 1.3_rk, 0.9_rk]
        knot_scale = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        pt = [0.37_rk, 0.61_rk]

        call surface_ref%set(knot_scale, knot_scale, Xc_scale, W_scale)
        call surface_scaled%set(knot_scale, knot_scale, Xc_scale, 1.0e-20_rk*W_scale)
        Xg_ref = surface_ref%cmp_Xg(pt)
        Xg_scaled = surface_scaled%cmp_Xg(pt)
        call surface_ref%derivative(pt, dT_ref, T_ref)
        call surface_scaled%derivative(pt, dT_scaled, T_scaled)
        call surface_ref%err%print()
        call surface_scaled%err%print()

        call ut%test(ti)%check(&
            name     = "affine knot scaling surface invariance",&
            res      = maxval(abs(Xg_ref - Xg_scaled)) <= 1.0e-10_rk .and. maxval(abs(T_ref - T_scaled)) <= 1.0e-10_rk .and. &
                maxval(abs(dT_ref - dT_scaled)) <= 1.0e-10_rk,&
            expected = .true.,&
            msg      = "Affine knot scaling surface invariance is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0101


    subroutine forcad_nurbs_surface_0102(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface_fd
        real(rk), allocatable :: Wc_fd(:), T0(:), Tm(:), Tp(:), dT0(:,:), dTm(:,:), dTp(:,:)
        real(rk), allocatable :: d2T0(:,:), d2Tm(:,:), d2Tp(:,:), cfd(:,:), cfd2(:,:)
        real(rk) :: h, pt(2), pm(2), pp(2)
        integer :: dir

        Wc_fd = [1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, &
            1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk]
        h = 1.0e-6_rk
        pt = [0.5_rk, 0.3_rk]

        call surface_fd%set_tetragon(L=[5.0_rk, 8.0_rk], nc=[4,4], Wc=Wc_fd)
        call surface_fd%derivative2(Xt=pt, d2Tgc=d2T0, dTgc=dT0, Tgc=T0)
        allocate(cfd(size(T0),2), cfd2(size(d2T0,1),size(d2T0,2)))
        do dir = 1, 2
            pm = pt
            pp = pt
            pm(dir) = pm(dir) - h
            pp(dir) = pp(dir) + h
            call surface_fd%derivative2(Xt=pm, d2Tgc=d2Tm, dTgc=dTm, Tgc=Tm)
            call surface_fd%derivative2(Xt=pp, d2Tgc=d2Tp, dTgc=dTp, Tgc=Tp)
            cfd(:,dir) = (Tp - Tm)/(2.0_rk*h)
            cfd2(:,dir) = reshape((dTp - dTm)/(2.0_rk*h), [size(d2T0,1)])
        end do
        call surface_fd%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS surface first-derivative finite difference",&
            res      = norm2(cfd - dT0),&
            expected = 0.0_rk,&
            tol      = 1.0e-7_rk,&
            msg      = "Surface first-derivative finite difference is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0102


    subroutine forcad_nurbs_surface_0103(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface_fd
        real(rk), allocatable :: Wc_fd(:), T0(:), Tm(:), Tp(:), dT0(:,:), dTm(:,:), dTp(:,:)
        real(rk), allocatable :: d2T0(:,:), d2Tm(:,:), d2Tp(:,:), cfd(:,:), cfd2(:,:)
        real(rk) :: h, pt(2), pm(2), pp(2)
        integer :: dir

        Wc_fd = [1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, &
            1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk]
        h = 1.0e-6_rk
        pt = [0.5_rk, 0.3_rk]

        call surface_fd%set_tetragon(L=[5.0_rk, 8.0_rk], nc=[4,4], Wc=Wc_fd)
        call surface_fd%derivative2(Xt=pt, d2Tgc=d2T0, dTgc=dT0, Tgc=T0)
        allocate(cfd(size(T0),2), cfd2(size(d2T0,1),size(d2T0,2)))
        do dir = 1, 2
            pm = pt
            pp = pt
            pm(dir) = pm(dir) - h
            pp(dir) = pp(dir) + h
            call surface_fd%derivative2(Xt=pm, d2Tgc=d2Tm, dTgc=dTm, Tgc=Tm)
            call surface_fd%derivative2(Xt=pp, d2Tgc=d2Tp, dTgc=dTp, Tgc=Tp)
            cfd(:,dir) = (Tp - Tm)/(2.0_rk*h)
            cfd2(:,dir) = reshape((dTp - dTm)/(2.0_rk*h), [size(d2T0,1)])
        end do
        call surface_fd%err%print()

        call ut%test(ti)%check(&
            name     = "NURBS surface second-derivative finite difference",&
            res      = norm2(cfd2 - d2T0),&
            expected = 0.0_rk,&
            tol      = 1.0e-6_rk,&
            msg      = "Surface second-derivative finite difference is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0103


    subroutine forcad_nurbs_surface_0104(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface_fd
        real(rk), allocatable :: Wc_fd(:), T0(:), Tm(:), Tp(:), dT0(:,:), dTm(:,:), dTp(:,:)
        real(rk), allocatable :: d2T0(:,:), d2Tm(:,:), d2Tp(:,:), cfd(:,:), cfd2(:,:)
        real(rk) :: h, pt(2), pm(2), pp(2)
        integer :: dir

        Wc_fd = [1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, &
            1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk]
        h = 1.0e-6_rk
        pt = [0.5_rk, 0.3_rk]

        call surface_fd%set_tetragon(L=[5.0_rk, 8.0_rk], nc=[4,4], Wc=Wc_fd)
        call surface_fd%derivative2(Xt=pt, d2Tgc=d2T0, dTgc=dT0, Tgc=T0)
        allocate(cfd(size(T0),2), cfd2(size(d2T0,1),size(d2T0,2)))
        do dir = 1, 2
            pm = pt
            pp = pt
            pm(dir) = pm(dir) - h
            pp(dir) = pp(dir) + h
            call surface_fd%derivative2(Xt=pm, d2Tgc=d2Tm, dTgc=dTm, Tgc=Tm)
            call surface_fd%derivative2(Xt=pp, d2Tgc=d2Tp, dTgc=dTp, Tgc=Tp)
            cfd(:,dir) = (Tp - Tm)/(2.0_rk*h)
            cfd2(:,dir) = reshape((dTp - dTm)/(2.0_rk*h), [size(d2T0,1)])
        end do

        call surface_fd%finalize()
        call surface_fd%set_tetragon(L=[5.0_rk, 8.0_rk], nc=[4,4])
        call surface_fd%derivative2(Xt=pt, d2Tgc=d2T0, dTgc=dT0, Tgc=T0)
        do dir = 1, 2
            pm = pt
            pp = pt
            pm(dir) = pm(dir) - h
            pp(dir) = pp(dir) + h
            call surface_fd%derivative2(Xt=pm, d2Tgc=d2Tm, dTgc=dTm, Tgc=Tm)
            call surface_fd%derivative2(Xt=pp, d2Tgc=d2Tp, dTgc=dTp, Tgc=Tp)
            cfd(:,dir) = (Tp - Tm)/(2.0_rk*h)
            cfd2(:,dir) = reshape((dTp - dTm)/(2.0_rk*h), [size(d2T0,1)])
        end do
        call surface_fd%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline surface first-derivative finite difference",&
            res      = norm2(cfd - dT0),&
            expected = 0.0_rk,&
            tol      = 1.0e-7_rk,&
            msg      = "Surface first-derivative finite difference is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0104


    subroutine forcad_nurbs_surface_0105(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface_fd
        real(rk), allocatable :: Wc_fd(:), T0(:), Tm(:), Tp(:), dT0(:,:), dTm(:,:), dTp(:,:)
        real(rk), allocatable :: d2T0(:,:), d2Tm(:,:), d2Tp(:,:), cfd(:,:), cfd2(:,:)
        real(rk) :: h, pt(2), pm(2), pp(2)
        integer :: dir

        Wc_fd = [1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, &
            1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk]
        h = 1.0e-6_rk
        pt = [0.5_rk, 0.3_rk]

        call surface_fd%set_tetragon(L=[5.0_rk, 8.0_rk], nc=[4,4], Wc=Wc_fd)
        call surface_fd%derivative2(Xt=pt, d2Tgc=d2T0, dTgc=dT0, Tgc=T0)
        allocate(cfd(size(T0),2), cfd2(size(d2T0,1),size(d2T0,2)))
        do dir = 1, 2
            pm = pt
            pp = pt
            pm(dir) = pm(dir) - h
            pp(dir) = pp(dir) + h
            call surface_fd%derivative2(Xt=pm, d2Tgc=d2Tm, dTgc=dTm, Tgc=Tm)
            call surface_fd%derivative2(Xt=pp, d2Tgc=d2Tp, dTgc=dTp, Tgc=Tp)
            cfd(:,dir) = (Tp - Tm)/(2.0_rk*h)
            cfd2(:,dir) = reshape((dTp - dTm)/(2.0_rk*h), [size(d2T0,1)])
        end do

        call surface_fd%finalize()
        call surface_fd%set_tetragon(L=[5.0_rk, 8.0_rk], nc=[4,4])
        call surface_fd%derivative2(Xt=pt, d2Tgc=d2T0, dTgc=dT0, Tgc=T0)
        do dir = 1, 2
            pm = pt
            pp = pt
            pm(dir) = pm(dir) - h
            pp(dir) = pp(dir) + h
            call surface_fd%derivative2(Xt=pm, d2Tgc=d2Tm, dTgc=dTm, Tgc=Tm)
            call surface_fd%derivative2(Xt=pp, d2Tgc=d2Tp, dTgc=dTp, Tgc=Tp)
            cfd(:,dir) = (Tp - Tm)/(2.0_rk*h)
            cfd2(:,dir) = reshape((dTp - dTm)/(2.0_rk*h), [size(d2T0,1)])
        end do
        call surface_fd%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline surface second-derivative finite difference",&
            res      = norm2(cfd2 - d2T0),&
            expected = 0.0_rk,&
            tol      = 1.0e-6_rk,&
            msg      = "Surface second-derivative finite difference is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0105


    subroutine forcad_nurbs_surface_0106(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface_ref
        real(rk), allocatable :: Xc_ref(:,:)
        real(rk), allocatable :: Xg_ref(:,:)
        real(rk) :: knot_ref(7), Xt_ref(4)
        integer :: ii, jj, idx

        allocate(Xc_ref(16, 3), source=0.0_rk)
        do jj = 1, 4
            do ii = 1, 4
                idx = ii + (jj - 1)*4
                Xc_ref(idx, 1) = real(ii - 1, rk)
                Xc_ref(idx, 2) = real(jj - 1, rk)
                Xc_ref(idx, 3) = 0.1_rk*real((ii - 1)*(jj - 1), rk)
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 2.4_rk, 3.1_rk, 4.0_rk]

        call surface_ref%set(knot_ref, knot_ref, Xc_ref, degree=[2, 2])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)
        Xg_ref = surface_ref%get_Xg()

        call surface_ref%insert_knots(1, [2.5_rk], [1])
        call surface_ref%insert_knots(2, [3.5_rk], [1])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)
        call surface_ref%err%print()

        call ut%test(ti)%check(&
            name     = "surface knot insertion geometry preservation",&
            res      = maxval(abs(surface_ref%get_Xg() - Xg_ref)),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "Surface knot insertion geometry preservation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0106


    subroutine forcad_nurbs_surface_0107(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface_ref
        real(rk), allocatable :: Xc_ref(:,:)
        real(rk), allocatable :: Xg_ref(:,:)
        real(rk) :: knot_ref(7), Xt_ref(4)
        integer :: ii, jj, idx

        allocate(Xc_ref(16, 3), source=0.0_rk)
        do jj = 1, 4
            do ii = 1, 4
                idx = ii + (jj - 1)*4
                Xc_ref(idx, 1) = real(ii - 1, rk)
                Xc_ref(idx, 2) = real(jj - 1, rk)
                Xc_ref(idx, 3) = 0.1_rk*real((ii - 1)*(jj - 1), rk)
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 2.4_rk, 3.1_rk, 4.0_rk]

        call surface_ref%set(knot_ref, knot_ref, Xc_ref, degree=[2, 2])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)
        Xg_ref = surface_ref%get_Xg()

        call surface_ref%insert_knots(1, [2.5_rk], [1])
        call surface_ref%insert_knots(2, [3.5_rk], [1])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)

        call surface_ref%remove_knots(1, [2.5_rk], [1])
        call surface_ref%remove_knots(2, [3.5_rk], [1])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)
        call surface_ref%err%print()

        call ut%test(ti)%check(&
            name     = "surface insertion-removal geometry preservation",&
            res      = maxval(abs(surface_ref%get_Xg() - Xg_ref)),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "Surface insertion-removal geometry preservation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0107


    subroutine forcad_nurbs_surface_0108(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface_ref
        real(rk), allocatable :: Xc_ref(:,:), Xg_ref(:,:), Xt1_after(:), Xt2_after(:)
        real(rk) :: knot_ref(7), Xt_ref(4)
        integer :: ii, jj, idx

        allocate(Xc_ref(16, 3), source=0.0_rk)
        do jj = 1, 4
            do ii = 1, 4
                idx = ii + (jj - 1)*4
                Xc_ref(idx, 1) = real(ii - 1, rk)
                Xc_ref(idx, 2) = real(jj - 1, rk)
                Xc_ref(idx, 3) = 0.1_rk*real((ii - 1)*(jj - 1), rk)
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 2.4_rk, 3.1_rk, 4.0_rk]

        call surface_ref%set(knot_ref, knot_ref, Xc_ref, degree=[2, 2])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)
        Xg_ref = surface_ref%get_Xg()

        call surface_ref%insert_knots(1, [2.5_rk], [1])
        call surface_ref%insert_knots(2, [3.5_rk], [1])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)

        call surface_ref%remove_knots(1, [2.5_rk], [1])
        call surface_ref%remove_knots(2, [3.5_rk], [1])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)

        call surface_ref%elevate_degree(1, 1)
        call surface_ref%elevate_degree(2, 1)
        call surface_ref%create(res1=4, res2=4)
        Xt1_after = surface_ref%get_Xt(1)
        Xt2_after = surface_ref%get_Xt(2)
        call surface_ref%err%print()

        call ut%test(ti)%check(&
            name     = "elevated surface direction-one domain",&
            res      = [Xt1_after(1), Xt1_after(size(Xt1_after))],&
            expected = [2.0_rk, 4.0_rk],&
            tol      = 1.0e-10_rk,&
            msg      = "Elevated surface direction-one domain is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0108


    subroutine forcad_nurbs_surface_0109(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface_ref
        real(rk), allocatable :: Xc_ref(:,:), Xg_ref(:,:), Xt1_after(:), Xt2_after(:)
        real(rk) :: knot_ref(7), Xt_ref(4)
        integer :: ii, jj, idx

        allocate(Xc_ref(16, 3), source=0.0_rk)
        do jj = 1, 4
            do ii = 1, 4
                idx = ii + (jj - 1)*4
                Xc_ref(idx, 1) = real(ii - 1, rk)
                Xc_ref(idx, 2) = real(jj - 1, rk)
                Xc_ref(idx, 3) = 0.1_rk*real((ii - 1)*(jj - 1), rk)
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 2.4_rk, 3.1_rk, 4.0_rk]

        call surface_ref%set(knot_ref, knot_ref, Xc_ref, degree=[2, 2])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)
        Xg_ref = surface_ref%get_Xg()

        call surface_ref%insert_knots(1, [2.5_rk], [1])
        call surface_ref%insert_knots(2, [3.5_rk], [1])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)

        call surface_ref%remove_knots(1, [2.5_rk], [1])
        call surface_ref%remove_knots(2, [3.5_rk], [1])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)

        call surface_ref%elevate_degree(1, 1)
        call surface_ref%elevate_degree(2, 1)
        call surface_ref%create(res1=4, res2=4)
        Xt1_after = surface_ref%get_Xt(1)
        Xt2_after = surface_ref%get_Xt(2)
        call surface_ref%err%print()

        call ut%test(ti)%check(&
            name     = "elevated surface direction-two domain",&
            res      = [Xt2_after(1), Xt2_after(size(Xt2_after))],&
            expected = [2.0_rk, 4.0_rk],&
            tol      = 1.0e-10_rk,&
            msg      = "Elevated surface direction-two domain is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0109


    subroutine forcad_nurbs_surface_0110(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface_ref
        real(rk), allocatable :: Xc_ref(:,:), Xg_ref(:,:), Xt1_after(:), Xt2_after(:)
        real(rk) :: knot_ref(7), Xt_ref(4)
        integer :: ii, jj, idx

        allocate(Xc_ref(16, 3), source=0.0_rk)
        do jj = 1, 4
            do ii = 1, 4
                idx = ii + (jj - 1)*4
                Xc_ref(idx, 1) = real(ii - 1, rk)
                Xc_ref(idx, 2) = real(jj - 1, rk)
                Xc_ref(idx, 3) = 0.1_rk*real((ii - 1)*(jj - 1), rk)
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 2.4_rk, 3.1_rk, 4.0_rk]

        call surface_ref%set(knot_ref, knot_ref, Xc_ref, degree=[2, 2])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)
        Xg_ref = surface_ref%get_Xg()

        call surface_ref%insert_knots(1, [2.5_rk], [1])
        call surface_ref%insert_knots(2, [3.5_rk], [1])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)

        call surface_ref%remove_knots(1, [2.5_rk], [1])
        call surface_ref%remove_knots(2, [3.5_rk], [1])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)

        call surface_ref%elevate_degree(1, 1)
        call surface_ref%elevate_degree(2, 1)
        call surface_ref%create(res1=4, res2=4)
        Xt1_after = surface_ref%get_Xt(1)
        Xt2_after = surface_ref%get_Xt(2)
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)
        call surface_ref%err%print()

        call ut%test(ti)%check(&
            name     = "surface degree elevation geometry preservation",&
            res      = maxval(abs(surface_ref%get_Xg() - Xg_ref)),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "Surface degree elevation geometry preservation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0110


    subroutine forcad_nurbs_surface_0111(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface_ref
        real(rk), allocatable :: Xc_ref(:,:), Xg_ref(:,:), Xt1_after(:), Xt2_after(:)
        real(rk) :: knot_ref(7), Xt_ref(4)
        integer :: ii, jj, idx

        allocate(Xc_ref(16, 3), source=0.0_rk)
        do jj = 1, 4
            do ii = 1, 4
                idx = ii + (jj - 1)*4
                Xc_ref(idx, 1) = real(ii - 1, rk)
                Xc_ref(idx, 2) = real(jj - 1, rk)
                Xc_ref(idx, 3) = 0.1_rk*real((ii - 1)*(jj - 1), rk)
            end do
        end do
        knot_ref = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        Xt_ref = [2.0_rk, 2.4_rk, 3.1_rk, 4.0_rk]

        call surface_ref%set(knot_ref, knot_ref, Xc_ref, degree=[2, 2])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)
        Xg_ref = surface_ref%get_Xg()

        call surface_ref%insert_knots(1, [2.5_rk], [1])
        call surface_ref%insert_knots(2, [3.5_rk], [1])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)

        call surface_ref%remove_knots(1, [2.5_rk], [1])
        call surface_ref%remove_knots(2, [3.5_rk], [1])
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)

        call surface_ref%elevate_degree(1, 1)
        call surface_ref%elevate_degree(2, 1)
        call surface_ref%create(res1=4, res2=4)
        Xt1_after = surface_ref%get_Xt(1)
        Xt2_after = surface_ref%get_Xt(2)
        call surface_ref%create(Xt1=Xt_ref, Xt2=Xt_ref)
        call surface_ref%err%print()

        call ut%test(ti)%check(&
            name     = "elevated surface degree",&
            res      = surface_ref%get_degree(),&
            expected = [3, 3],&
            msg      = "Elevated surface degree is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0111


    subroutine forcad_nurbs_surface_0112(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sh0
        type(nurbs_surface) :: shr
        type(nurbs_surface) :: shb
        real(rk), allocatable :: Xc_s(:,:)
        real(rk), allocatable :: Wc_s(:)
        real(rk), allocatable :: Xc0(:,:)
        real(rk), allocatable :: knot1_s(:)
        real(rk), allocatable :: knot2_s(:)
        real(rk), allocatable :: S1(:,:)
        real(rk), allocatable :: Bdirect(:,:)
        integer :: ncoord
        integer :: nc0
        integer :: ndof_old

        allocate(Xc_s(4,3), Wc_s(4), knot1_s(4), knot2_s(4))
        Xc_s(1,:) = [ 2.5_rk, -2.5_rk, 0.0_rk]
        Xc_s(2,:) = [-2.5_rk, -2.5_rk, 0.0_rk]
        Xc_s(3,:) = [ 2.5_rk,  2.5_rk, 0.0_rk]
        Xc_s(4,:) = [-2.5_rk,  2.5_rk, 0.0_rk]
        Wc_s = [1.0_rk, 0.5_rk, 1.0_rk, 1.0_rk]
        knot1_s = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot2_s = knot1_s

        call shr%set(knot1_s, knot2_s, Xc_s, Wc_s)
        sh0 = shr
        Xc0 = sh0%get_Xc()
        ncoord = size(Xc0,2)
        nc0 = size(Xc0,1)
        ndof_old = nc0*ncoord

        call shr%elevate_degree(1, 3, Bs=S1)
        call shb%set(sh0%get_knot(1), sh0%get_knot(2), Xc0, sh0%get_Wc())
        call shb%elevate_degree(1, 3, B=Bdirect)
        call sh0%err%print()
        call shr%err%print()
        call shb%err%print()

        call ut%test(ti)%check(&
            name     = "surface direction-one degree-elevation map",&
            res      = maxval(abs(Bdirect - kron_eye(S1, ncoord))),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "surface elevate dir1 B-only map",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0112


    subroutine forcad_nurbs_surface_0113(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sh0
        type(nurbs_surface) :: shr
        type(nurbs_surface) :: shb
        real(rk), allocatable :: Xc_s(:,:)
        real(rk), allocatable :: Wc_s(:)
        real(rk), allocatable :: Xc0(:,:)
        real(rk), allocatable :: knot1_s(:)
        real(rk), allocatable :: knot2_s(:)
        real(rk), allocatable :: S1(:,:)
        real(rk), allocatable :: S2(:,:)
        real(rk), allocatable :: Bdirect(:,:)
        integer :: ncoord
        integer :: nc0
        integer :: ndof_old

        allocate(Xc_s(4,3), Wc_s(4), knot1_s(4), knot2_s(4))
        Xc_s(1,:) = [ 2.5_rk, -2.5_rk, 0.0_rk]
        Xc_s(2,:) = [-2.5_rk, -2.5_rk, 0.0_rk]
        Xc_s(3,:) = [ 2.5_rk,  2.5_rk, 0.0_rk]
        Xc_s(4,:) = [-2.5_rk,  2.5_rk, 0.0_rk]
        Wc_s = [1.0_rk, 0.5_rk, 1.0_rk, 1.0_rk]
        knot1_s = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot2_s = knot1_s

        call shr%set(knot1_s, knot2_s, Xc_s, Wc_s)
        sh0 = shr
        Xc0 = sh0%get_Xc()
        ncoord = size(Xc0,2)
        nc0 = size(Xc0,1)
        ndof_old = nc0*ncoord

        call shr%elevate_degree(1, 3, Bs=S1)
        call shb%set(sh0%get_knot(1), sh0%get_knot(2), Xc0, sh0%get_Wc())
        call shb%elevate_degree(1, 3, B=Bdirect)
        call shr%elevate_degree(2, 2, Bs=S2)
        call shb%elevate_degree(2, 2, B=Bdirect)
        call sh0%err%print()
        call shr%err%print()
        call shb%err%print()

        call ut%test(ti)%check(&
            name     = "surface direction-two degree-elevation map",&
            res      = maxval(abs(Bdirect - kron_eye(S2, ncoord))),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "surface elevate dir2 B-only map",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0113


    subroutine forcad_nurbs_surface_0114(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sh0
        type(nurbs_surface) :: shr
        type(nurbs_surface) :: shb
        real(rk), allocatable :: Xc_s(:,:)
        real(rk), allocatable :: Wc_s(:)
        real(rk), allocatable :: Xc0(:,:)
        real(rk), allocatable :: knot1_s(:)
        real(rk), allocatable :: knot2_s(:)
        real(rk), allocatable :: S1(:,:)
        real(rk), allocatable :: S2(:,:)
        real(rk), allocatable :: S3(:,:)
        real(rk), allocatable :: Bdirect(:,:)
        real(rk), allocatable :: u1(:), u2(:)
        integer, allocatable :: r1(:), r2(:)
        integer :: ncoord
        integer :: nc0
        integer :: ndof_old

        allocate(Xc_s(4,3), Wc_s(4), knot1_s(4), knot2_s(4))
        Xc_s(1,:) = [ 2.5_rk, -2.5_rk, 0.0_rk]
        Xc_s(2,:) = [-2.5_rk, -2.5_rk, 0.0_rk]
        Xc_s(3,:) = [ 2.5_rk,  2.5_rk, 0.0_rk]
        Xc_s(4,:) = [-2.5_rk,  2.5_rk, 0.0_rk]
        Wc_s = [1.0_rk, 0.5_rk, 1.0_rk, 1.0_rk]
        knot1_s = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot2_s = knot1_s

        call shr%set(knot1_s, knot2_s, Xc_s, Wc_s)
        sh0 = shr
        Xc0 = sh0%get_Xc()
        ncoord = size(Xc0,2)
        nc0 = size(Xc0,1)
        ndof_old = nc0*ncoord

        call shr%elevate_degree(1, 3, Bs=S1)
        call shb%set(sh0%get_knot(1), sh0%get_knot(2), Xc0, sh0%get_Wc())
        call shb%elevate_degree(1, 3, B=Bdirect)
        call shr%elevate_degree(2, 2, Bs=S2)
        call shb%elevate_degree(2, 2, B=Bdirect)
        u1 = linspace(0.0_rk, 1.0_rk, 7)
        u1 = u1(2:6)
        u2 = linspace(0.0_rk, 1.0_rk, 6)
        u2 = u2(2:5)
        allocate(r1(size(u1)), source=2)
        allocate(r2(size(u2)), source=1)
        call shr%insert_knots(1, u1, r1, Bs=S3)
        call shb%insert_knots(1, u1, r1, B=Bdirect)
        call sh0%err%print()
        call shr%err%print()
        call shb%err%print()

        call ut%test(ti)%check(&
            name     = "surface direction-one knot-insertion map",&
            res      = maxval(abs(Bdirect - kron_eye(S3, ncoord))),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "surface insert dir1 B-only map",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0114


    subroutine forcad_nurbs_surface_0115(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sh0
        type(nurbs_surface) :: shr
        type(nurbs_surface) :: shb
        real(rk), allocatable :: Xc_s(:,:)
        real(rk), allocatable :: Wc_s(:)
        real(rk), allocatable :: Xc0(:,:)
        real(rk), allocatable :: knot1_s(:)
        real(rk), allocatable :: knot2_s(:)
        real(rk), allocatable :: S1(:,:)
        real(rk), allocatable :: S2(:,:)
        real(rk), allocatable :: S3(:,:)
        real(rk), allocatable :: S4(:,:)
        real(rk), allocatable :: Bdirect(:,:)
        real(rk), allocatable :: u1(:), u2(:)
        integer, allocatable :: r1(:), r2(:)
        integer :: ncoord
        integer :: nc0
        integer :: ndof_old

        allocate(Xc_s(4,3), Wc_s(4), knot1_s(4), knot2_s(4))
        Xc_s(1,:) = [ 2.5_rk, -2.5_rk, 0.0_rk]
        Xc_s(2,:) = [-2.5_rk, -2.5_rk, 0.0_rk]
        Xc_s(3,:) = [ 2.5_rk,  2.5_rk, 0.0_rk]
        Xc_s(4,:) = [-2.5_rk,  2.5_rk, 0.0_rk]
        Wc_s = [1.0_rk, 0.5_rk, 1.0_rk, 1.0_rk]
        knot1_s = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot2_s = knot1_s

        call shr%set(knot1_s, knot2_s, Xc_s, Wc_s)
        sh0 = shr
        Xc0 = sh0%get_Xc()
        ncoord = size(Xc0,2)
        nc0 = size(Xc0,1)
        ndof_old = nc0*ncoord

        call shr%elevate_degree(1, 3, Bs=S1)
        call shb%set(sh0%get_knot(1), sh0%get_knot(2), Xc0, sh0%get_Wc())
        call shb%elevate_degree(1, 3, B=Bdirect)
        call shr%elevate_degree(2, 2, Bs=S2)
        call shb%elevate_degree(2, 2, B=Bdirect)
        u1 = linspace(0.0_rk, 1.0_rk, 7)
        u1 = u1(2:6)
        u2 = linspace(0.0_rk, 1.0_rk, 6)
        u2 = u2(2:5)
        allocate(r1(size(u1)), source=2)
        allocate(r2(size(u2)), source=1)
        call shr%insert_knots(1, u1, r1, Bs=S3)
        call shb%insert_knots(1, u1, r1, B=Bdirect)
        call shr%insert_knots(2, u2, r2, Bs=S4)
        call shb%insert_knots(2, u2, r2, B=Bdirect)
        call sh0%err%print()
        call shr%err%print()
        call shb%err%print()

        call ut%test(ti)%check(&
            name     = "surface direction-two knot-insertion map",&
            res      = maxval(abs(Bdirect - kron_eye(S4, ncoord))),&
            expected = 0.0_rk,&
            tol      = 1.0e-10_rk,&
            msg      = "surface insert dir2 B-only map",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0115


    subroutine forcad_nurbs_surface_0116(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sh0, shr, shfd, shb
        real(rk), allocatable :: Xc_s(:,:), Wc_s(:), Xc0(:,:), Xp(:,:), Xm(:,:)
        real(rk), allocatable :: Xcp_vec(:), Xcm_vec(:), knot1_s(:), knot2_s(:)
        real(rk), allocatable :: S1(:,:), S2(:,:), S3(:,:), S4(:,:), Bs(:,:), B(:,:), Bdirect(:,:), Bfd(:,:)
        real(rk), allocatable :: u1(:), u2(:)
        integer, allocatable :: r1(:), r2(:)
        real(rk) :: rel_err
        integer :: ncoord, nc0, ndof_old, ndof_new, idx, ic, dc

        allocate(Xc_s(4,3), Wc_s(4), knot1_s(4), knot2_s(4))
        Xc_s(1,:) = [ 2.5_rk, -2.5_rk, 0.0_rk]
        Xc_s(2,:) = [-2.5_rk, -2.5_rk, 0.0_rk]
        Xc_s(3,:) = [ 2.5_rk,  2.5_rk, 0.0_rk]
        Xc_s(4,:) = [-2.5_rk,  2.5_rk, 0.0_rk]
        Wc_s = [1.0_rk, 0.5_rk, 1.0_rk, 1.0_rk]
        knot1_s = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot2_s = knot1_s

        call shr%set(knot1_s, knot2_s, Xc_s, Wc_s)
        sh0 = shr
        Xc0 = sh0%get_Xc()
        ncoord = size(Xc0,2)
        nc0 = size(Xc0,1)
        ndof_old = nc0*ncoord

        call shr%elevate_degree(1, 3, Bs=S1)
        call shb%set(sh0%get_knot(1), sh0%get_knot(2), Xc0, sh0%get_Wc())
        call shb%elevate_degree(1, 3, B=Bdirect)
        call shr%elevate_degree(2, 2, Bs=S2)
        call shb%elevate_degree(2, 2, B=Bdirect)
        u1 = linspace(0.0_rk, 1.0_rk, 7)
        u1 = u1(2:6)
        u2 = linspace(0.0_rk, 1.0_rk, 6)
        u2 = u2(2:5)
        allocate(r1(size(u1)), source=2)
        allocate(r2(size(u2)), source=1)
        call shr%insert_knots(1, u1, r1, Bs=S3)
        call shb%insert_knots(1, u1, r1, B=Bdirect)
        call shr%insert_knots(2, u2, r2, Bs=S4)
        call shb%insert_knots(2, u2, r2, B=Bdirect)
        Bs = matmul(S4, matmul(S3, matmul(S2, S1)))
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

            call shfd%set(sh0%get_knot(1), sh0%get_knot(2), Xp, sh0%get_Wc())
            call shfd%elevate_degree(1, 3)
            call shfd%elevate_degree(2, 2)
            call shfd%insert_knots(1, u1, r1)
            call shfd%insert_knots(2, u2, r2)
            Xcp_vec = reshape(transpose(shfd%get_Xc()), [ndof_new])

            call shfd%set(sh0%get_knot(1), sh0%get_knot(2), Xm, sh0%get_Wc())
            call shfd%elevate_degree(1, 3)
            call shfd%elevate_degree(2, 2)
            call shfd%insert_knots(1, u1, r1)
            call shfd%insert_knots(2, u2, r2)
            Xcm_vec = reshape(transpose(shfd%get_Xc()), [ndof_new])

            Bfd(:,idx) = (Xcp_vec - Xcm_vec) * 5.0e4_rk
        end do
        rel_err = norm2(Bfd - B)/norm2(Bfd)
        call sh0%err%print()
        call shr%err%print()
        call shfd%err%print()
        call shb%err%print()

        call ut%test(ti)%check(&
            name     = "surface refinement-map finite difference",&
            res      = rel_err,&
            expected = 0.0_rk,&
            tol      = 1.0e-7_rk,&
            msg      = "Surface refinement-map finite difference is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0116


    subroutine forcad_nurbs_surface_0117(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: plate
        real(rk), allocatable :: Xp(:,:)

        call plate%set_tetragon([1.0_rk, 1.0_rk], [2,2])
        Xp = plate%get_Xc()
        call plate%put_to_nurbs(Xp, reshape([1, 2, 3, 4], [1,4]))
        call plate%err%print()

        call ut%test(ti)%check(&
            name     = "surface control-to-geometry projection",&
            res      = maxval(abs(plate%get_Xg() - Xp)),&
            expected = 0.0_rk,&
            tol      = 1e-10_rk,&
            msg      = "Surface control-to-geometry projection is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0117


    subroutine forcad_nurbs_surface_0118(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: plate
        real(rk), allocatable :: Xp(:,:)

        call plate%set_tetragon([1.0_rk, 1.0_rk], [2,2])
        Xp = plate%get_Xc()
        call plate%put_to_nurbs(Xp, reshape([1, 2, 3, 4], [1,4]))
        call plate%err%print()

        call ut%test(ti)%check(&
            name     = "surface projected visualization connectivity",&
            res      = shape(plate%get_elem_Xg_vis()),&
            expected = [1,4],&
            msg      = "Surface projected visualization connectivity is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0118


    subroutine forcad_nurbs_surface_0119(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sgen
        real(rk) :: Xgen(6,3)
        character(len=*), parameter :: fgen = "vtk/forcad_surface_general.iges"
        character(len=512) :: line
        integer :: unit, ios
        logical :: found

        Xgen(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xgen(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xgen(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xgen(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xgen(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xgen(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        call sgen%set([1, 1], [3, 2], Xgen)
        call sgen%export_iges(fgen)
        found = .false.
        open(newunit=unit, file=fgen, status="old", action="read", iostat=ios)
        if (ios == 0) then
            do
                read(unit, "(A)", iostat=ios) line
                if (ios /= 0) exit
                if (index(line, "128,2,1,1,1,0,0,1") > 0) then
                    found = .true.
                    exit
                end if
            end do
            close(unit)
        end if
        call sgen%err%print()

        call ut%test(ti)%check(&
            name     = "surface IGES control-net dimensions",&
            res      = found,&
            expected = .true.,&
            msg      = "Surface IGES control-net dimensions are incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0119


    subroutine forcad_nurbs_surface_0120(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sproj
        real(rk) :: Xproj(4,3), Xt_out(2), Xg_out(3)

        Xproj(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xproj(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xproj(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xproj(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        call sproj%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xproj)
        call sproj%nearest_point2([0.2_rk, 0.3_rk, 0.4_rk], 1.0e-12_rk, 0, Xt_out, Xg_out)
        call sproj%err%print()

        call ut%test(ti)%check(&
            name     = "zero-iteration surface projection outputs",&
            res      = all(Xt_out >= 0.0_rk .and. Xt_out <= 1.0_rk) .and. all(ieee_is_finite(Xg_out)),&
            expected = .true.,&
            msg      = "Zero-iteration surface projection outputs are incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0120


    subroutine forcad_nurbs_surface_0121(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sbad
        real(rk) :: Xbad(4,3)
        real(rk) :: Xnew(6,3)
        real(rk) :: Wbad(4)
        real(rk) :: Wbad6(6)
        real(rk) :: Wgood(4)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad6 = 1.0_rk
        Wbad6(2) = 0.0_rk
        Wgood = 1.0_rk
        call sbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)
        call sbad%err%print()

        call ut%test(ti)%check(&
            name     = "surface rejects nonpositive weights",&
            res      = sbad%err%ok,&
            expected = .false.,&
            msg      = "A surface must reject nonpositive weights.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0121


    subroutine forcad_nurbs_surface_0122(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sbad
        real(rk) :: Xbad(4,3)
        real(rk) :: Xnew(6,3)
        real(rk) :: Wbad(4)
        real(rk) :: Wbad6(6)
        real(rk) :: Wgood(4)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad6 = 1.0_rk
        Wbad6(2) = 0.0_rk
        Wgood = 1.0_rk
        call sbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)
        call sbad%err%print()

        call ut%test(ti)%check(&
            name     = "invalid-weight surface rational state",&
            res      = sbad%is_rational(),&
            expected = .false.,&
            msg      = "Invalid-weight surface rational state is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0122


    subroutine forcad_nurbs_surface_0123(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sbad
        type(nurbs_surface) :: smod
        real(rk) :: Xbad(4,3)
        real(rk) :: Xnew(6,3)
        real(rk) :: Wbad(4)
        real(rk) :: Wbad6(6)
        real(rk) :: Wgood(4)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad6 = 1.0_rk
        Wbad6(2) = 0.0_rk
        Wgood = 1.0_rk
        call sbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call smod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call smod%modify_Wc(0.5_rk, 1)
        call sbad%err%print()
        call smod%err%print()

        call ut%test(ti)%check(&
            name     = "modified surface rational state",&
            res      = smod%is_rational(),&
            expected = .true.,&
            msg      = "surface rational cache update",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0123


    subroutine forcad_nurbs_surface_0124(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sbad
        type(nurbs_surface) :: smod
        real(rk) :: Xbad(4,3)
        real(rk) :: Xnew(6,3)
        real(rk) :: Wbad(4)
        real(rk) :: Wbad6(6)
        real(rk) :: Wgood(4)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad6 = 1.0_rk
        Wbad6(2) = 0.0_rk
        Wgood = 1.0_rk
        call sbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call smod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call smod%modify_Wc(0.5_rk, 1)
        call smod%modify_Wc(1.0_rk, 1)
        call sbad%err%print()
        call smod%err%print()

        call ut%test(ti)%check(&
            name     = "restored uniform surface rational state",&
            res      = smod%is_rational(),&
            expected = .false.,&
            msg      = "surface uniform cache restore",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0124


    subroutine forcad_nurbs_surface_0125(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sbad
        type(nurbs_surface) :: smod
        real(rk) :: Xbad(4,3)
        real(rk) :: Xnew(6,3)
        real(rk) :: Wbad(4)
        real(rk) :: Wbad6(6)
        real(rk) :: Wgood(4)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad6 = 1.0_rk
        Wbad6(2) = 0.0_rk
        Wgood = 1.0_rk
        call sbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call smod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call smod%modify_Wc(0.5_rk, 1)
        call smod%modify_Wc(1.0_rk, 1)
        call smod%modify_Wc(-1.0_rk, 1)
        call sbad%err%print()
        call smod%err%print()

        call ut%test(ti)%check(&
            name     = "surface rejects negative weight modification",&
            res      = smod%err%ok,&
            expected = .false.,&
            msg      = "Surface rejects negative weight modification is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0125


    subroutine forcad_nurbs_surface_0126(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sbad
        type(nurbs_surface) :: smod
        type(nurbs_surface) :: sstate
        real(rk) :: Xbad(4,3), Xnew(6,3), Wbad(4), Wbad6(6), Wgood(4), Wstate(4)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad6 = 1.0_rk
        Wbad6(2) = 0.0_rk
        Wgood = 1.0_rk
        call sbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call smod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call smod%modify_Wc(0.5_rk, 1)
        call smod%modify_Wc(1.0_rk, 1)
        call smod%modify_Wc(-1.0_rk, 1)

        Wstate = 1.0_rk
        Wstate(2) = 0.5_rk
        call sstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wstate)
        call sstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)
        call sbad%err%print()
        call smod%err%print()
        call sstate%err%print()

        call ut%test(ti)%check(&
            name     = "unweighted surface reset clears rational state",&
            res      = sstate%is_rational(),&
            expected = .false.,&
            msg      = "Unweighted surface reset clears rational state is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0126


    subroutine forcad_nurbs_surface_0127(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sbad
        type(nurbs_surface) :: smod
        type(nurbs_surface) :: sstate
        type(nurbs_surface) :: sdegree1
        real(rk) :: Xbad(4,3), Xnew(6,3), Wbad(4), Wbad6(6), Wgood(4), Wstate(4)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad6 = 1.0_rk
        Wbad6(2) = 0.0_rk
        Wgood = 1.0_rk
        call sbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call smod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call smod%modify_Wc(0.5_rk, 1)
        call smod%modify_Wc(1.0_rk, 1)
        call smod%modify_Wc(-1.0_rk, 1)

        Wstate = 1.0_rk
        Wstate(2) = 0.5_rk
        call sstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wstate)
        call sstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)

        call sdegree1%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, degree=[1])
        call sbad%err%print()
        call smod%err%print()
        call sstate%err%print()
        call sdegree1%err%print()

        call ut%test(ti)%check(&
            name     = "surface rejects incomplete degree vector",&
            res      = sdegree1%err%ok,&
            expected = .false.,&
            msg      = "Surface rejects incomplete degree vector is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0127


    subroutine forcad_nurbs_surface_0128(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sbad
        type(nurbs_surface) :: smod
        type(nurbs_surface) :: sstate
        type(nurbs_surface) :: sdegree1
        type(nurbs_surface) :: sdegree2
        real(rk) :: Xbad(4,3), Xnew(6,3), Wbad(4), Wbad6(6), Wgood(4), Wstate(4)

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad6 = 1.0_rk
        Wbad6(2) = 0.0_rk
        Wgood = 1.0_rk
        call sbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call smod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call smod%modify_Wc(0.5_rk, 1)
        call smod%modify_Wc(1.0_rk, 1)
        call smod%modify_Wc(-1.0_rk, 1)

        Wstate = 1.0_rk
        Wstate(2) = 0.5_rk
        call sstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wstate)
        call sstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)

        call sdegree1%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, degree=[1])

        call sdegree2%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1], [-1, -1], [-1, -1])
        call sbad%err%print()
        call smod%err%print()
        call sstate%err%print()
        call sdegree1%err%print()
        call sdegree2%err%print()

        call ut%test(ti)%check(&
            name     = "surface rejects inconsistent parameter metadata",&
            res      = sdegree2%err%ok,&
            expected = .false.,&
            msg      = "Surface rejects inconsistent parameter metadata is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0128


    subroutine forcad_nurbs_surface_0129(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sbad, smod, sstate, sdegree1, sdegree2, satomic
        real(rk) :: Xbad(4,3), Xnew(6,3), Wbad(4), Wbad6(6), Wgood(4), Wstate(4)
        real(rk), allocatable :: Xkeep(:,:)
        logical :: rejected
        integer :: preserved

        Xbad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xbad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xbad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xbad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xnew(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xnew(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
        Xnew(4,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xnew(5,:) = [1.0_rk, 1.0_rk, 0.0_rk]
        Xnew(6,:) = [2.0_rk, 1.0_rk, 0.0_rk]
        Wbad = 1.0_rk
        Wbad(2) = 0.0_rk
        Wbad6 = 1.0_rk
        Wbad6(2) = 0.0_rk
        Wgood = 1.0_rk
        call sbad%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wbad)

        call smod%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wgood)
        call smod%modify_Wc(0.5_rk, 1)
        call smod%modify_Wc(1.0_rk, 1)
        call smod%modify_Wc(-1.0_rk, 1)

        Wstate = 1.0_rk
        Wstate(2) = 0.5_rk
        call sstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wstate)
        call sstate%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad)

        call sdegree1%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, degree=[1])

        call sdegree2%set([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1], [-1, -1], [-1, -1])

        call satomic%set([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xbad, Wstate)
        call satomic%set([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xnew, Wbad6, degree=[2, 1])
        rejected = .not. satomic%err%ok
        call satomic%err%print()
        call satomic%err%reset()
        preserved = 0
        if (rejected .and. all(satomic%get_nc() == [2, 2]) .and. all(satomic%get_degree() == [1, 1])) then
            Xkeep = satomic%get_Xc()
            if (size(Xkeep, 1) == size(Xbad, 1) .and. size(Xkeep, 2) == size(Xbad, 2)) then
                if (maxval(abs(Xkeep - Xbad)) <= 1.0e-12_rk .and. satomic%is_rational()) preserved = 1
            end if
        end if
        call sbad%err%print()
        call smod%err%print()
        call sstate%err%print()
        call sdegree1%err%print()
        call sdegree2%err%print()
        call satomic%err%print()

        call ut%test(ti)%check(&
            name     = "failed surface reset preserves object state",&
            res      = preserved,&
            expected = 1,&
            msg      = "Failed surface reset preserves object state is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0129


    subroutine forcad_nurbs_surface_0130(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sinsert
        real(rk) :: knot_bad(4), Xc_bad(4,3)

        knot_bad = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        Xc_bad = 0.0_rk
        Xc_bad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xc_bad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xc_bad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xc_bad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]

        call sinsert%set(knot_bad, knot_bad, Xc_bad)
        call sinsert%insert_knots(1, [0.5_rk], [3])
        call sinsert%err%print()

        call ut%test(ti)%check(&
            name     = "surface rejects excessive knot insertion",&
            res      = sinsert%err%ok,&
            expected = .false.,&
            msg      = "surface bad insert mult",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0130


    subroutine forcad_nurbs_surface_0131(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sinsert
        type(nurbs_surface) :: sremove
        real(rk) :: knot_bad(4), Xc_bad(4,3)

        knot_bad = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        Xc_bad = 0.0_rk
        Xc_bad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xc_bad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xc_bad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xc_bad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]

        call sinsert%set(knot_bad, knot_bad, Xc_bad)
        call sinsert%insert_knots(1, [0.5_rk], [3])

        call sremove%set(knot_bad, knot_bad, Xc_bad)
        call sremove%remove_knots(1, [0.5_rk], [1])
        call sinsert%err%print()
        call sremove%err%print()

        call ut%test(ti)%check(&
            name     = "surface rejects absent knot removal",&
            res      = sremove%err%ok,&
            expected = .false.,&
            msg      = "surface absent remove knot",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0131


    subroutine forcad_nurbs_surface_0132(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sinsert, sremove, selev
        real(rk) :: knot_bad(4), Xc_bad(4,3)

        knot_bad = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        Xc_bad = 0.0_rk
        Xc_bad(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        Xc_bad(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        Xc_bad(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        Xc_bad(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]

        call sinsert%set(knot_bad, knot_bad, Xc_bad)
        call sinsert%insert_knots(1, [0.5_rk], [3])

        call sremove%set(knot_bad, knot_bad, Xc_bad)
        call sremove%remove_knots(1, [0.5_rk], [1])

        call selev%set(knot_bad, knot_bad, Xc_bad)
        call selev%elevate_degree(1, -1)
        call sinsert%err%print()
        call sremove%err%print()
        call selev%err%print()

        call ut%test(ti)%check(&
            name     = "surface rejects negative degree elevation",&
            res      = selev%err%ok,&
            expected = .false.,&
            msg      = "surface negative degree elevation",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0132


    subroutine forcad_nurbs_surface_0133(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)

        ivals = sempty%get_multiplicity(1)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface multiplicity getter",&
            res      = size(ivals),&
            expected = 0,&
            msg      = "surface unset mult",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0133


    subroutine forcad_nurbs_surface_0134(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface continuity getter",&
            res      = size(ivals),&
            expected = 0,&
            msg      = "surface unset cont",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0134


    subroutine forcad_nurbs_surface_0135(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "surface invalid-direction multiplicity getter",&
            res      = size(ivals),&
            expected = 0,&
            msg      = "surface bad mult dir",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0135


    subroutine forcad_nurbs_surface_0136(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "surface invalid-direction control-count getter",&
            res      = sempty%get_nc(3),&
            expected = 0,&
            msg      = "surface bad nc dir",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0136


    subroutine forcad_nurbs_surface_0137(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "surface invalid-direction degree getter",&
            res      = sempty%get_degree(3),&
            expected = 0,&
            msg      = "surface bad deg dir",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0137


    subroutine forcad_nurbs_surface_0138(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface control-net getter",&
            res      = size(rmat),&
            expected = 0,&
            msg      = "surface unset Xc",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0138


    subroutine forcad_nurbs_surface_0139(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface control-point getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "surface unset Xci",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0139


    subroutine forcad_nurbs_surface_0140(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface control-coordinate getter",&
            res      = sempty%get_Xc(1, 1),&
            expected = 0.0_rk,&
            tol      = 1.0e-12_rk,&
            msg      = "surface unset Xcid",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0140


    subroutine forcad_nurbs_surface_0141(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface geometry getter",&
            res      = size(rmat),&
            expected = 0,&
            msg      = "surface unset Xg",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0141


    subroutine forcad_nurbs_surface_0142(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface geometry-point getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "surface unset Xgi",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0142


    subroutine forcad_nurbs_surface_0143(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface geometry-coordinate getter",&
            res      = sempty%get_Xg(1, 1),&
            expected = 0.0_rk,&
            tol      = 1.0e-12_rk,&
            msg      = "surface unset Xgid",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0143


    subroutine forcad_nurbs_surface_0144(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        rvals = sempty%get_Wc()
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface weight-vector getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "surface unset Wc",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0144


    subroutine forcad_nurbs_surface_0145(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        rvals = sempty%get_Wc()
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface weight getter",&
            res      = sempty%get_Wc(1),&
            expected = 0.0_rk,&
            tol      = 1.0e-12_rk,&
            msg      = "surface unset Wci",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0145


    subroutine forcad_nurbs_surface_0146(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        rvals = sempty%get_Wc()
        rvals = sempty%get_Xt(1)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface parameter-vector getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "surface unset Xt",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0146


    subroutine forcad_nurbs_surface_0147(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        rvals = sempty%get_Wc()
        rvals = sempty%get_Xt(1)
        rvals = sempty%get_Xt(3)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "surface invalid-direction parameter getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "surface bad Xt dir",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0147


    subroutine forcad_nurbs_surface_0148(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        rvals = sempty%get_Wc()
        rvals = sempty%get_Xt(1)
        rvals = sempty%get_Xt(3)
        rvals = sempty%get_knot(1)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface knot-vector getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "surface unset knot",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0148


    subroutine forcad_nurbs_surface_0149(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        rvals = sempty%get_Wc()
        rvals = sempty%get_Xt(1)
        rvals = sempty%get_Xt(3)
        rvals = sempty%get_knot(1)
        rvals = sempty%get_knot(3)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "surface invalid-direction knot-vector getter",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "surface bad knot dir",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0149


    subroutine forcad_nurbs_surface_0150(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        rvals = sempty%get_Wc()
        rvals = sempty%get_Xt(1)
        rvals = sempty%get_Xt(3)
        rvals = sempty%get_knot(1)
        rvals = sempty%get_knot(3)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface knot getter",&
            res      = sempty%get_knot(1, 1),&
            expected = 0.0_rk,&
            tol      = 1.0e-12_rk,&
            msg      = "surface unset knoti",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0150


    subroutine forcad_nurbs_surface_0151(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        rvals = sempty%get_Wc()
        rvals = sempty%get_Xt(1)
        rvals = sempty%get_Xt(3)
        rvals = sempty%get_knot(1)
        rvals = sempty%get_knot(3)
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "surface invalid-direction knot getter",&
            res      = sempty%get_knot(3, 1),&
            expected = 0.0_rk,&
            tol      = 1.0e-12_rk,&
            msg      = "surface bad knoti dir",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0151


    subroutine forcad_nurbs_surface_0152(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        rvals = sempty%get_Wc()
        rvals = sempty%get_Xt(1)
        rvals = sempty%get_Xt(3)
        rvals = sempty%get_knot(1)
        rvals = sempty%get_knot(3)
        rvals = sempty%cmp_Xg([0.0_rk, 0.0_rk])
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface geometry computation",&
            res      = size(rvals),&
            expected = 0,&
            msg      = "surface unset cmp_Xg",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0152


    subroutine forcad_nurbs_surface_0153(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        rvals = sempty%get_Wc()
        rvals = sempty%get_Xt(1)
        rvals = sempty%get_Xt(3)
        rvals = sempty%get_knot(1)
        rvals = sempty%get_knot(3)
        rvals = sempty%cmp_Xg([0.0_rk, 0.0_rk])
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface geometry-grid shape",&
            res      = sempty%get_ng(),&
            expected = [0, 0],&
            msg      = "surface unset ng",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0153


    subroutine forcad_nurbs_surface_0154(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        rvals = sempty%get_Wc()
        rvals = sempty%get_Xt(1)
        rvals = sempty%get_Xt(3)
        rvals = sempty%get_knot(1)
        rvals = sempty%get_knot(3)
        rvals = sempty%cmp_Xg([0.0_rk, 0.0_rk])
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface control-net shape",&
            res      = sempty%get_nc(),&
            expected = [0, 0],&
            msg      = "surface unset nc",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0154


    subroutine forcad_nurbs_surface_0155(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sempty
        integer, allocatable :: ivals(:)
        real(rk), allocatable :: rvals(:), rmat(:,:)

        ivals = sempty%get_multiplicity(1)
        ivals = sempty%get_continuity(2)
        ivals = sempty%get_multiplicity(3)

        rmat = sempty%get_Xc()
        rvals = sempty%get_Xc(1)
        rmat = sempty%get_Xg()
        rvals = sempty%get_Xg(1)
        rvals = sempty%get_Wc()
        rvals = sempty%get_Xt(1)
        rvals = sempty%get_Xt(3)
        rvals = sempty%get_knot(1)
        rvals = sempty%get_knot(3)
        rvals = sempty%cmp_Xg([0.0_rk, 0.0_rk])
        call sempty%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface degree vector",&
            res      = sempty%get_degree(),&
            expected = [0, 0],&
            msg      = "surface unset degree",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0155


    subroutine forcad_nurbs_surface_0156(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_derivative_order_32(2) = [3,2]
        type(nurbs_surface) :: rational_order
        type(nurbs_surface) :: bspline_order
        real(rk) :: Xc_order(4,3), Wc_order(4), expected_order(4), magnitude
        real(rk), allocatable :: Dgc_order(:)

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 144.0_rk/((1.25_rk**4)*(1.8_rk**3))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order)
        call rational_order%err%print()
        call bspline_order%err%print()

        call ut%test(ti)%check(&
            name     = "rational surface mixed derivative",&
            res      = Dgc_order,&
            expected = expected_order,&
            tol      = 512.0_rk*epsilon(1.0_rk)*magnitude,&
            msg      = "surface rational mixed deriv",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0156


    subroutine forcad_nurbs_surface_0157(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_derivative_order_32(2) = [3,2]
        type(nurbs_surface) :: rational_order
        type(nurbs_surface) :: bspline_order
        real(rk) :: Xc_order(4,3), Wc_order(4), expected_order(4), magnitude
        real(rk), allocatable :: Dgc_order(:)
        real(rk), allocatable :: Dgc_order_vector(:,:)

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 144.0_rk/((1.25_rk**4)*(1.8_rk**3))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk])
        call rational_order%err%print()
        call bspline_order%err%print()

        call ut%test(ti)%check(&
            name     = "rational surface vector derivative",&
            res      = Dgc_order_vector(1,:),&
            expected = expected_order,&
            tol      = 512.0_rk*epsilon(1.0_rk)*magnitude,&
            msg      = "surface rational vector deriv",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0157


    subroutine forcad_nurbs_surface_0158(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_derivative_order_32(2) = [3,2]
        type(nurbs_surface) :: rational_order
        type(nurbs_surface) :: bspline_order
        real(rk) :: Xc_order(4,3), Wc_order(4), expected_order(4), magnitude
        real(rk), allocatable :: Dgc_order(:)
        real(rk), allocatable :: Dgc_order_vector(:,:)
        real(rk), allocatable :: Dgc_active(:)
        integer :: first_active(2)

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 144.0_rk/((1.25_rk**4)*(1.8_rk**3))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk],&
            order        = surface_derivative_order_32,&
            first_active = first_active,&
            Dgc          = Dgc_active)
        call rational_order%err%print()
        call bspline_order%err%print()

        call ut%test(ti)%check(&
            name     = "active rational surface mixed derivative",&
            res      = Dgc_active,&
            expected = expected_order,&
            tol      = 512.0_rk*epsilon(1.0_rk)*magnitude,&
            msg      = "surface active rational deriv",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0158


    subroutine forcad_nurbs_surface_0159(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_derivative_order_32(2) = [3,2]
        type(nurbs_surface) :: rational_order
        type(nurbs_surface) :: bspline_order
        real(rk) :: Xc_order(4,3), Wc_order(4), expected_order(4), magnitude
        real(rk), allocatable :: Dgc_order(:)
        real(rk), allocatable :: Dgc_order_vector(:,:)
        real(rk), allocatable :: Dgc_active(:)
        integer :: first_active(2)

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 144.0_rk/((1.25_rk**4)*(1.8_rk**3))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk],&
            order        = surface_derivative_order_32,&
            first_active = first_active,&
            Dgc          = Dgc_active)
        call rational_order%err%print()
        call bspline_order%err%print()

        call ut%test(ti)%check(&
            name     = "active surface first basis indices",&
            res      = first_active,&
            expected = [1,1],&
            msg      = "surface active first indices",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0159


    subroutine forcad_nurbs_surface_0160(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_derivative_order_32(2) = [3,2]
        type(nurbs_surface) :: rational_order
        type(nurbs_surface) :: bspline_order
        real(rk) :: Xc_order(4,3), Wc_order(4), expected_order(4), magnitude
        real(rk), allocatable :: Dgc_order(:)
        real(rk), allocatable :: Dgc_order_vector(:,:)
        real(rk), allocatable :: Dgc_active(:)
        real(rk), allocatable :: Dgc_active_vector(:,:)
        integer, allocatable :: first_vector(:,:)
        integer :: first_active(2)

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 144.0_rk/((1.25_rk**4)*(1.8_rk**3))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk],&
            order        = surface_derivative_order_32,&
            first_active = first_active,&
            Dgc          = Dgc_active)

        call rational_order%derivative_order_active(&
            Xt1          = [0.25_rk,0.5_rk],&
            Xt2          = [0.4_rk],&
            order        = surface_derivative_order_32,&
            first_active = first_vector,&
            Dgc          = Dgc_active_vector)
        call rational_order%err%print()
        call bspline_order%err%print()

        call ut%test(ti)%check(&
            name     = "active surface vector derivative ordering",&
            res      = Dgc_active_vector(:,1),&
            expected = expected_order,&
            tol      = 512.0_rk*epsilon(1.0_rk)*magnitude,&
            msg      = "surface active vector ordering",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0160


    subroutine forcad_nurbs_surface_0161(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_derivative_order_32(2) = [3,2]
        type(nurbs_surface) :: rational_order
        type(nurbs_surface) :: bspline_order
        real(rk) :: Xc_order(4,3), Wc_order(4), expected_order(4), magnitude
        real(rk), allocatable :: Dgc_order(:), Dgc_order_vector(:,:), Dgc_bspline(:)
        real(rk), allocatable :: Dgc_active(:)
        real(rk), allocatable :: Dgc_active_vector(:,:)
        integer, allocatable :: first_vector(:,:)
        integer :: first_active(2)

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 144.0_rk/((1.25_rk**4)*(1.8_rk**3))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk],&
            order        = surface_derivative_order_32,&
            first_active = first_active,&
            Dgc          = Dgc_active)

        call rational_order%derivative_order_active(&
            Xt1          = [0.25_rk,0.5_rk],&
            Xt2          = [0.4_rk],&
            order        = surface_derivative_order_32,&
            first_active = first_vector,&
            Dgc          = Dgc_active_vector)

        call bspline_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_bspline)
        call rational_order%err%print()
        call bspline_order%err%print()

        call ut%test(ti)%check(&
            name     = "B-spline surface derivative above degree",&
            res      = Dgc_bspline,&
            expected = [0.0_rk,0.0_rk,0.0_rk,0.0_rk],&
            tol      = 32.0_rk*epsilon(1.0_rk),&
            msg      = "surface bspline high deriv",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0161


    subroutine forcad_nurbs_surface_0162(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_derivative_order_11(2) = [1,1]
        integer, parameter :: surface_derivative_order_32(2) = [3,2]
        type(nurbs_surface) :: rational_order
        type(nurbs_surface) :: bspline_order
        type(nurbs_surface) :: active_surface
        real(rk) :: Xc_order(4,3), Wc_order(4), expected_order(4), magnitude
        real(rk) :: Xc_active(9,3), Wc_active(9)
        real(rk), allocatable :: Dgc_order(:), Dgc_order_vector(:,:), Dgc_bspline(:)
        real(rk), allocatable :: Dgc_active(:)
        real(rk), allocatable :: Dgc_active_vector(:,:)
        real(rk), allocatable :: Dgc_dense_active(:)
        real(rk), allocatable :: Dgc_reconstructed(:)
        integer, allocatable :: first_vector(:,:)
        integer :: first_active(2), l1, l2, i1, i2

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 144.0_rk/((1.25_rk**4)*(1.8_rk**3))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk],&
            order        = surface_derivative_order_32,&
            first_active = first_active,&
            Dgc          = Dgc_active)

        call rational_order%derivative_order_active(&
            Xt1          = [0.25_rk,0.5_rk],&
            Xt2          = [0.4_rk],&
            order        = surface_derivative_order_32,&
            first_active = first_vector,&
            Dgc          = Dgc_active_vector)

        call bspline_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_bspline)

        Xc_active = 0.0_rk
        Wc_active = [1.0_rk,1.1_rk,0.9_rk,1.3_rk,0.8_rk,1.2_rk,0.7_rk,1.4_rk,1.05_rk]
        call active_surface%set(&
            knot1 = [0.0_rk,0.0_rk,0.35_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,0.7_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_active,&
            Wc    = Wc_active)
        call active_surface%derivative_order(&
            Xt    = [0.8_rk,0.2_rk],&
            order = surface_derivative_order_11,&
            Dgc   = Dgc_dense_active)
        call active_surface%derivative_order_active(&
            Xt           = [0.8_rk,0.2_rk],&
            order        = surface_derivative_order_11,&
            first_active = first_active,&
            Dgc          = Dgc_active)
        allocate(Dgc_reconstructed(size(Dgc_dense_active)), source=0.0_rk)
        do l2 = 0, 1
            do l1 = 0, 1
                i1 = first_active(1) + l1
                i2 = first_active(2) + l2
                Dgc_reconstructed(i1+(i2-1)*3) = Dgc_active(l1+1+2*l2)
            end do
        end do
        call rational_order%err%print()
        call bspline_order%err%print()
        call active_surface%err%print()

        call ut%test(ti)%check(&
            name     = "active surface derivative reconstruction",&
            res      = Dgc_reconstructed,&
            expected = Dgc_dense_active,&
            tol      = 256.0_rk*epsilon(1.0_rk)*max(1.0_rk,maxval(abs(Dgc_dense_active))),&
            msg      = "surface active dense reconstruction",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0162


    subroutine forcad_nurbs_surface_0163(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_derivative_order_11(2) = [1,1]
        integer, parameter :: surface_derivative_order_32(2) = [3,2]
        type(nurbs_surface) :: rational_order
        type(nurbs_surface) :: bspline_order
        type(nurbs_surface) :: active_surface
        real(rk) :: Xc_order(4,3), Wc_order(4), expected_order(4), magnitude
        real(rk) :: Xc_active(9,3), Wc_active(9)
        real(rk), allocatable :: Dgc_order(:), Dgc_order_vector(:,:), Dgc_bspline(:)
        real(rk), allocatable :: Dgc_active(:)
        real(rk), allocatable :: Dgc_active_vector(:,:)
        real(rk), allocatable :: Dgc_dense_active(:)
        real(rk), allocatable :: Dgc_reconstructed(:)
        integer, allocatable :: first_vector(:,:)
        integer :: first_active(2), l1, l2, i1, i2

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 144.0_rk/((1.25_rk**4)*(1.8_rk**3))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk],&
            order        = surface_derivative_order_32,&
            first_active = first_active,&
            Dgc          = Dgc_active)

        call rational_order%derivative_order_active(&
            Xt1          = [0.25_rk,0.5_rk],&
            Xt2          = [0.4_rk],&
            order        = surface_derivative_order_32,&
            first_active = first_vector,&
            Dgc          = Dgc_active_vector)

        call bspline_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_bspline)

        Xc_active = 0.0_rk
        Wc_active = [1.0_rk,1.1_rk,0.9_rk,1.3_rk,0.8_rk,1.2_rk,0.7_rk,1.4_rk,1.05_rk]
        call active_surface%set(&
            knot1 = [0.0_rk,0.0_rk,0.35_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,0.7_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_active,&
            Wc    = Wc_active)
        call active_surface%derivative_order(&
            Xt    = [0.8_rk,0.2_rk],&
            order = surface_derivative_order_11,&
            Dgc   = Dgc_dense_active)
        call active_surface%derivative_order_active(&
            Xt           = [0.8_rk,0.2_rk],&
            order        = surface_derivative_order_11,&
            first_active = first_active,&
            Dgc          = Dgc_active)
        allocate(Dgc_reconstructed(size(Dgc_dense_active)), source=0.0_rk)
        do l2 = 0, 1
            do l1 = 0, 1
                i1 = first_active(1) + l1
                i2 = first_active(2) + l2
                Dgc_reconstructed(i1+(i2-1)*3) = Dgc_active(l1+1+2*l2)
            end do
        end do
        call rational_order%err%print()
        call bspline_order%err%print()
        call active_surface%err%print()

        call ut%test(ti)%check(&
            name     = "active surface derivative storage reduction",&
            res      = size(Dgc_active) < size(Dgc_dense_active),&
            expected = .true.,&
            msg      = "surface active lower memory",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0163


    subroutine forcad_nurbs_surface_0164(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_derivative_order_00(2) = [0,0]
        integer, parameter :: surface_derivative_order_11(2) = [1,1]
        integer, parameter :: surface_derivative_order_32(2) = [3,2]
        type(nurbs_surface) :: rational_order, bspline_order, active_surface, empty_active
        real(rk) :: Xc_order(4,3), Wc_order(4), expected_order(4), magnitude
        real(rk) :: Xc_active(9,3), Wc_active(9)
        real(rk), allocatable :: Dgc_order(:), Dgc_order_vector(:,:), Dgc_bspline(:)
        real(rk), allocatable :: Dgc_active(:), Dgc_active_vector(:,:), Dgc_dense_active(:), Dgc_reconstructed(:), Dgc_empty(:)
        integer, allocatable :: first_vector(:,:)
        integer :: first_active(2), l1, l2, i1, i2

        Xc_order = 0.0_rk
        Wc_order = [1.0_rk,2.0_rk,3.0_rk,6.0_rk]
        call rational_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order,&
            Wc    = Wc_order)
        call bspline_order%set(&
            knot1 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_order)
        magnitude = 144.0_rk/((1.25_rk**4)*(1.8_rk**3))
        expected_order = [-magnitude,magnitude,magnitude,-magnitude]

        call rational_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order)

        call rational_order%derivative_order(&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_order_vector,&
            Xt1   = [0.25_rk,0.5_rk],&
            Xt2   = [0.4_rk])

        call rational_order%derivative_order_active(&
            Xt           = [0.25_rk,0.4_rk],&
            order        = surface_derivative_order_32,&
            first_active = first_active,&
            Dgc          = Dgc_active)

        call rational_order%derivative_order_active(&
            Xt1          = [0.25_rk,0.5_rk],&
            Xt2          = [0.4_rk],&
            order        = surface_derivative_order_32,&
            first_active = first_vector,&
            Dgc          = Dgc_active_vector)

        call bspline_order%derivative_order(&
            Xt    = [0.25_rk,0.4_rk],&
            order = surface_derivative_order_32,&
            Dgc   = Dgc_bspline)

        Xc_active = 0.0_rk
        Wc_active = [1.0_rk,1.1_rk,0.9_rk,1.3_rk,0.8_rk,1.2_rk,0.7_rk,1.4_rk,1.05_rk]
        call active_surface%set(&
            knot1 = [0.0_rk,0.0_rk,0.35_rk,1.0_rk,1.0_rk],&
            knot2 = [0.0_rk,0.0_rk,0.7_rk,1.0_rk,1.0_rk],&
            Xc    = Xc_active,&
            Wc    = Wc_active)
        call active_surface%derivative_order(&
            Xt    = [0.8_rk,0.2_rk],&
            order = surface_derivative_order_11,&
            Dgc   = Dgc_dense_active)
        call active_surface%derivative_order_active(&
            Xt           = [0.8_rk,0.2_rk],&
            order        = surface_derivative_order_11,&
            first_active = first_active,&
            Dgc          = Dgc_active)
        allocate(Dgc_reconstructed(size(Dgc_dense_active)), source=0.0_rk)
        do l2 = 0, 1
            do l1 = 0, 1
                i1 = first_active(1) + l1
                i2 = first_active(2) + l2
                Dgc_reconstructed(i1+(i2-1)*3) = Dgc_active(l1+1+2*l2)
            end do
        end do

        call empty_active%derivative_order_active(&
            Xt           = [0.5_rk,0.5_rk],&
            order        = surface_derivative_order_00,&
            first_active = first_active,&
            Dgc          = Dgc_empty)
        call rational_order%err%print()
        call bspline_order%err%print()
        call active_surface%err%print()
        call empty_active%err%print()

        call ut%test(ti)%check(&
            name     = "unset surface active derivative",&
            res      = size(Dgc_empty),&
            expected = 0,&
            msg      = "surface unset active derivative",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0164


    subroutine forcad_nurbs_surface_0165(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian diagnostic",&
            res      = surface%err%ok,&
            expected = .true.,&
            msg      = "Ansatz Hessian diagnostic is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0165


    subroutine forcad_nurbs_surface_0166(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian xx scaling",&
            res      = d2T_dX2(1,1,1),&
            expected = 0.125_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "The physical xx derivative must include the x-coordinate scaling",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0166


    subroutine forcad_nurbs_surface_0167(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian xy scaling",&
            res      = d2T_dX2(1,1,2),&
            expected = 1.0_rk/6.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Ansatz Hessian xy scaling is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0167


    subroutine forcad_nurbs_surface_0168(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian yy scaling",&
            res      = d2T_dX2(1,2,2),&
            expected = 1.0_rk/18.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "The physical yy derivative must include the y-coordinate scaling",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0168


    subroutine forcad_nurbs_surface_0169(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian symmetry",&
            res      = maxval(abs(d2T_dX2(:,1,2)-d2T_dX2(:,2,1))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Mixed physical surface derivatives must be symmetric",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0169


    subroutine forcad_nurbs_surface_0170(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz Hessian partition",&
            res      = maxval(abs(sum(d2T_dX2,dim=1))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Physical surface second derivatives must sum to zero",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0170


    subroutine forcad_nurbs_surface_0171(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz positive measure",&
            res      = dA > 0.0_rk,&
            expected = .true.,&
            msg      = "Ansatz positive measure is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0171


    subroutine forcad_nurbs_surface_0172(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            Wc     = hessian_Wc,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "rational ansatz partition",&
            res      = sum(T),&
            expected = 1.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Rational surface shape functions must form a partition of unity",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0172


    subroutine forcad_nurbs_surface_0173(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            Wc     = hessian_Wc,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "rational ansatz gradient partition",&
            res      = maxval(abs(sum(dT_dX,dim=1))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Rational physical surface gradients must sum to zero",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0173


    subroutine forcad_nurbs_surface_0174(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            Wc     = hessian_Wc,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "rational ansatz Hessian partition",&
            res      = maxval(abs(sum(d2T_dX2,dim=1))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Rational physical surface second derivatives must sum to zero",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0174


    subroutine forcad_nurbs_surface_0175(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            Wc     = hessian_Wc,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "rational ansatz Hessian symmetry",&
            res      = maxval(abs(d2T_dX2(:,1,2)-d2T_dX2(:,2,1))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Rational mixed physical derivatives must be symmetric",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0175


    subroutine forcad_nurbs_surface_0176(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            Wc     = hessian_Wc,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "rational ansatz x chain rule",&
            res      = sum(d2T_dX2(:,1,1)*Xc2(:,1)),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "The rational Hessian must reproduce the affine x coordinate",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0176


    subroutine forcad_nurbs_surface_0177(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6)
        real(rk) :: Xc2(9,2)
        real(rk) :: hessian_Wc(9)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            Wc     = hessian_Wc,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "rational ansatz y chain rule",&
            res      = sum(d2T_dX2(:,2,2)*Xc2(:,2)),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "The rational Hessian must reproduce the affine y coordinate",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0177


    subroutine forcad_nurbs_surface_0178(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6), Xc2(9,2), Xc3(9,3), hessian_Wc(9), dA, u, v
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            Wc     = hessian_Wc,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        do jj = 1, 3
            v = 0.5_rk*real(jj-1,rk)
            do ii = 1, 3
                point = (jj-1)*3 + ii
                u = 0.5_rk*real(ii-1,rk)
                Xc3(point,1) = u
                if (ii == 3) Xc3(point,1) = Xc3(point,1) + 0.25_rk
                Xc3(point,2) = v + 0.2_rk*u*v
                Xc3(point,3) = 0.0_rk
            end do
        end do
        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc3,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "embedded ansatz Hessian symmetry",&
            res      = maxval(abs(d2T_dX2(:,1,2)-d2T_dX2(:,2,1))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Embedded-surface mixed physical derivatives must be symmetric",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0178


    subroutine forcad_nurbs_surface_0179(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6), Xc2(9,2), Xc3(9,3), hessian_Wc(9), dA, u, v
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            Wc     = hessian_Wc,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        do jj = 1, 3
            v = 0.5_rk*real(jj-1,rk)
            do ii = 1, 3
                point = (jj-1)*3 + ii
                u = 0.5_rk*real(ii-1,rk)
                Xc3(point,1) = u
                if (ii == 3) Xc3(point,1) = Xc3(point,1) + 0.25_rk
                Xc3(point,2) = v + 0.2_rk*u*v
                Xc3(point,3) = 0.0_rk
            end do
        end do
        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc3,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "embedded ansatz normal component",&
            res      = maxval(abs(d2T_dX2(:,:,3))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Embedded ansatz normal component is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0179


    subroutine forcad_nurbs_surface_0180(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6), Xc2(9,2), Xc3(9,3), hessian_Wc(9), dA, u, v
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            Wc     = hessian_Wc,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        do jj = 1, 3
            v = 0.5_rk*real(jj-1,rk)
            do ii = 1, 3
                point = (jj-1)*3 + ii
                u = 0.5_rk*real(ii-1,rk)
                Xc3(point,1) = u
                if (ii == 3) Xc3(point,1) = Xc3(point,1) + 0.25_rk
                Xc3(point,2) = v + 0.2_rk*u*v
                Xc3(point,3) = 0.0_rk
            end do
        end do
        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc3,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "embedded ansatz normal symmetry",&
            res      = maxval(abs(d2T_dX2(:,3,:))),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Normal Hessian rows and columns must both vanish",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0180


    subroutine forcad_nurbs_surface_0181(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        real(rk), parameter :: HESSIAN_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6), Xc2(9,2), Xc3(9,3), hessian_Wc(9), dA, u, v
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            Wc     = hessian_Wc,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        do jj = 1, 3
            v = 0.5_rk*real(jj-1,rk)
            do ii = 1, 3
                point = (jj-1)*3 + ii
                u = 0.5_rk*real(ii-1,rk)
                Xc3(point,1) = u
                if (ii == 3) Xc3(point,1) = Xc3(point,1) + 0.25_rk
                Xc3(point,2) = v + 0.2_rk*u*v
                Xc3(point,3) = 0.0_rk
            end do
        end do
        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc3,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "embedded ansatz chain rule",&
            res      = maxval(abs([ sum(d2T_dX2(:,1,1)*Xc3(:,1)), sum(d2T_dX2(:,1,2)*Xc3(:,1)), sum(d2T_dX2(:,2,2)*Xc3(:,1)),&
                sum(d2T_dX2(:,1,1)*Xc3(:,2)), sum(d2T_dX2(:,1,2)*Xc3(:,2)), sum(d2T_dX2(:,2,2)*Xc3(:,2))])),&
            expected = 0.0_rk,&
            tol      = HESSIAN_TOL,&
            msg      = "Embedded ansatz chain rule is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0181


    subroutine forcad_nurbs_surface_0182(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_01(2) = [0,1]
        integer, parameter :: surface_quadrature_order_11(2) = [1,1]
        type(nurbs_surface) :: surface
        real(rk) :: hessian_knot(6), Xc2(9,2), Xc3(9,3), hessian_Wc(9), dA, u, v
        real(rk), allocatable :: T(:), dT_dX(:,:), d2T_dX2(:,:,:)
        integer :: ii, jj, point

        hessian_knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        do jj = 1, 3
            do ii = 1, 3
                point = (jj-1)*3 + ii
                Xc2(point,:) = [real(ii-1,rk), 1.5_rk*real(jj-1,rk)]
                hessian_Wc(point) = 1.0_rk + 0.02_rk*real(ii-1,rk) + 0.03_rk*real(jj-1,rk)
            end do
        end do

        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc2,&
            Wc     = hessian_Wc,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        do jj = 1, 3
            v = 0.5_rk*real(jj-1,rk)
            do ii = 1, 3
                point = (jj-1)*3 + ii
                u = 0.5_rk*real(ii-1,rk)
                Xc3(point,1) = u
                if (ii == 3) Xc3(point,1) = Xc3(point,1) + 0.25_rk
                Xc3(point,2) = v + 0.2_rk*u*v
                Xc3(point,3) = 0.0_rk
            end do
        end do
        call surface%finalize()
        call surface%set(&
            knot1  = hessian_knot,&
            knot2  = hessian_knot,&
            Xc     = Xc3,&
            degree = [2,2])
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_11)

        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = T,&
            dTgc_dXg    = dT_dX,&
            d2Tgc_dXg2  = d2T_dX2,&
            dA          = dA,&
            ngauss      = surface_quadrature_order_01)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "ansatz invalid quadrature",&
            res      = .not. surface%err%ok .and. size(T) == 0 .and. size(dT_dX) == 0 .and. size(d2T_dX2) == 0 .and. abs(dA) <= &
                tiny(1.0_rk),&
            expected = .true.,&
            msg      = "Invalid quadrature must return empty outputs and a diagnostic",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0182


    subroutine forcad_nurbs_surface_0183(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: family_surface
        real(rk) :: family_knot1(9)
        real(rk) :: family_knot2(9)
        real(rk) :: family_Xc(36,3)
        real(rk), allocatable :: Xt1(:)
        real(rk), allocatable :: Xt2(:)
        real(rk), allocatable :: Xt_grid(:,:)
        real(rk), allocatable :: Xg_reference(:,:)
        integer :: ii
        integer :: jj
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk]
        family_knot2 = [-1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk]
        do jj = 1, 6
            do ii = 1, 6
                point = (jj-1)*6 + ii
                family_Xc(point,:) = [&
                    real(ii-1,rk),&
                    real(jj-1,rk),&
                    0.02_rk*real((ii-1)*(jj-1),rk)]
            end do
        end do

        call family_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call family_surface%create(res1=7, res2=7)
        Xt1 = family_surface%get_Xt(1)
        Xt2 = family_surface%get_Xt(2)
        call ndgrid(Xt1, Xt2, Xt_grid)
        Xg_reference = family_surface%get_Xg()
        call family_surface%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families finite geometry",&
            res      = all(ieee_is_finite(Xg_reference)),&
            expected = .true.,&
            msg      = "Mixed knot families finite geometry is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0183


    subroutine forcad_nurbs_surface_0184(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: family_surface
        real(rk) :: family_knot1(9)
        real(rk) :: family_knot2(9)
        real(rk) :: family_Xc(36,3)
        real(rk), allocatable :: Xt1(:)
        real(rk), allocatable :: Xt2(:)
        real(rk), allocatable :: Xt_grid(:,:)
        real(rk), allocatable :: Xg_reference(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk]
        family_knot2 = [-1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk]
        do jj = 1, 6
            do ii = 1, 6
                point = (jj-1)*6 + ii
                family_Xc(point,:) = [&
                    real(ii-1,rk),&
                    real(jj-1,rk),&
                    0.02_rk*real((ii-1)*(jj-1),rk)]
            end do
        end do

        call family_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call family_surface%create(res1=7, res2=7)
        Xt1 = family_surface%get_Xt(1)
        Xt2 = family_surface%get_Xt(2)
        call ndgrid(Xt1, Xt2, Xt_grid)
        Xg_reference = family_surface%get_Xg()
        elem = family_surface%get_elem_Xg_vis()
        call family_surface%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families create connectivity",&
            res      = size(elem,1),&
            expected = 36,&
            msg      = "A 7 by 7 visualization grid must produce 36 quadrilateral cells",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0184


    subroutine forcad_nurbs_surface_0185(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: family_surface
        real(rk) :: family_knot1(9)
        real(rk) :: family_knot2(9)
        real(rk) :: family_Xc(36,3)
        real(rk), allocatable :: Xt1(:)
        real(rk), allocatable :: Xt2(:)
        real(rk), allocatable :: Xt_grid(:,:)
        real(rk), allocatable :: Xg_reference(:,:)
        real(rk), allocatable :: Tgrid(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk]
        family_knot2 = [-1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk]
        do jj = 1, 6
            do ii = 1, 6
                point = (jj-1)*6 + ii
                family_Xc(point,:) = [&
                    real(ii-1,rk),&
                    real(jj-1,rk),&
                    0.02_rk*real((ii-1)*(jj-1),rk)]
            end do
        end do

        call family_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call family_surface%create(res1=7, res2=7)
        Xt1 = family_surface%get_Xt(1)
        Xt2 = family_surface%get_Xt(2)
        call ndgrid(Xt1, Xt2, Xt_grid)
        Xg_reference = family_surface%get_Xg()
        elem = family_surface%get_elem_Xg_vis()
        call family_surface%basis(res1=3, res2=3, Tgc=Tgrid)
        elem = family_surface%get_elem_Xg_vis()
        call family_surface%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families basis connectivity",&
            res      = size(elem,1),&
            expected = 36,&
            msg      = "Mixed knot families basis connectivity is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0185


    subroutine forcad_nurbs_surface_0186(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: INVARIANT_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: family_surface
        real(rk) :: family_knot1(9)
        real(rk) :: family_knot2(9)
        real(rk) :: family_Xc(36,3)
        real(rk) :: xi(2)
        real(rk), allocatable :: Xt1(:)
        real(rk), allocatable :: Xt2(:)
        real(rk), allocatable :: Xt_grid(:,:)
        real(rk), allocatable :: Xg_reference(:,:)
        real(rk), allocatable :: active1(:)
        real(rk), allocatable :: active2(:)
        real(rk), allocatable :: T(:)
        real(rk), allocatable :: Tgrid(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk]
        family_knot2 = [-1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk]
        do jj = 1, 6
            do ii = 1, 6
                point = (jj-1)*6 + ii
                family_Xc(point,:) = [&
                    real(ii-1,rk),&
                    real(jj-1,rk),&
                    0.02_rk*real((ii-1)*(jj-1),rk)]
            end do
        end do

        call family_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call family_surface%create(res1=7, res2=7)
        Xt1 = family_surface%get_Xt(1)
        Xt2 = family_surface%get_Xt(2)
        call ndgrid(Xt1, Xt2, Xt_grid)
        Xg_reference = family_surface%get_Xg()
        elem = family_surface%get_elem_Xg_vis()
        call family_surface%basis(res1=3, res2=3, Tgc=Tgrid)
        elem = family_surface%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,6,2)
        active2 = active_knots(family_knot2,6,2)
        xi = [0.5_rk*(active1(1)+active1(2)), 0.5_rk*(active2(1)+active2(2))]
        call family_surface%basis(xi, T)
        call family_surface%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families partition",&
            res      = sum(T),&
            expected = 1.0_rk,&
            tol      = INVARIANT_TOL,&
            msg      = "Mixed knot families partition is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0186


    subroutine forcad_nurbs_surface_0187(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: family_surface
        real(rk) :: family_knot1(9)
        real(rk) :: family_knot2(9)
        real(rk) :: family_Xc(36,3)
        real(rk) :: xi(2)
        real(rk), allocatable :: Xt1(:)
        real(rk), allocatable :: Xt2(:)
        real(rk), allocatable :: Xt_grid(:,:)
        real(rk), allocatable :: Xg_reference(:,:)
        real(rk), allocatable :: active1(:)
        real(rk), allocatable :: active2(:)
        real(rk), allocatable :: T(:)
        real(rk), allocatable :: Tgrid(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk]
        family_knot2 = [-1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk]
        do jj = 1, 6
            do ii = 1, 6
                point = (jj-1)*6 + ii
                family_Xc(point,:) = [&
                    real(ii-1,rk),&
                    real(jj-1,rk),&
                    0.02_rk*real((ii-1)*(jj-1),rk)]
            end do
        end do

        call family_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call family_surface%create(res1=7, res2=7)
        Xt1 = family_surface%get_Xt(1)
        Xt2 = family_surface%get_Xt(2)
        call ndgrid(Xt1, Xt2, Xt_grid)
        Xg_reference = family_surface%get_Xg()
        elem = family_surface%get_elem_Xg_vis()
        call family_surface%basis(res1=3, res2=3, Tgc=Tgrid)
        elem = family_surface%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,6,2)
        active2 = active_knots(family_knot2,6,2)
        xi = [0.5_rk*(active1(1)+active1(2)), 0.5_rk*(active2(1)+active2(2))]
        call family_surface%basis(xi, T)
        elem = family_surface%cmp_elem()
        call family_surface%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families element count",&
            res      = size(elem,1),&
            expected = active_span_count(family_knot1,6,2)*active_span_count(family_knot2,6,2),&
            msg      = "Mixed knot families element count is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0187


    subroutine forcad_nurbs_surface_0188(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_33(2) = [3,3]
        type(nurbs_surface) :: family_surface
        real(rk) :: family_knot1(9)
        real(rk) :: family_knot2(9)
        real(rk) :: family_Xc(36,3)
        real(rk) :: xi(2)
        real(rk) :: dA
        real(rk), allocatable :: Xt1(:)
        real(rk), allocatable :: Xt2(:)
        real(rk), allocatable :: Xt_grid(:,:)
        real(rk), allocatable :: Xg_reference(:,:)
        real(rk), allocatable :: active1(:), active2(:), T(:), Tgrid(:,:), dT_dX(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk]
        family_knot2 = [-1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk]
        do jj = 1, 6
            do ii = 1, 6
                point = (jj-1)*6 + ii
                family_Xc(point,:) = [&
                    real(ii-1,rk),&
                    real(jj-1,rk),&
                    0.02_rk*real((ii-1)*(jj-1),rk)]
            end do
        end do

        call family_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call family_surface%create(res1=7, res2=7)
        Xt1 = family_surface%get_Xt(1)
        Xt2 = family_surface%get_Xt(2)
        call ndgrid(Xt1, Xt2, Xt_grid)
        Xg_reference = family_surface%get_Xg()
        elem = family_surface%get_elem_Xg_vis()
        call family_surface%basis(res1=3, res2=3, Tgc=Tgrid)
        elem = family_surface%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,6,2)
        active2 = active_knots(family_knot2,6,2)
        xi = [0.5_rk*(active1(1)+active1(2)), 0.5_rk*(active2(1)+active2(2))]
        call family_surface%basis(xi, T)
        elem = family_surface%cmp_elem()
        call family_surface%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = T,&
            dTgc_dXg = dT_dX,&
            dA       = dA,&
            ngauss   = surface_quadrature_order_33)
        call family_surface%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families ansatz",&
            res      = dA > 0.0_rk,&
            expected = .true.,&
            msg      = "Mixed knot families ansatz is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0188


    subroutine forcad_nurbs_surface_0189(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_33(2) = [3,3]
        real(rk), parameter :: REFINEMENT_TOL = 2.0e-10_rk
        type(nurbs_surface) :: family_surface, refined_surface
        real(rk) :: family_knot1(9)
        real(rk) :: family_knot2(9)
        real(rk) :: family_Xc(36,3)
        real(rk) :: xi(2)
        real(rk) :: dA
        real(rk), allocatable :: Xt1(:), Xt2(:), Xt_grid(:,:), Xg_reference(:,:), Xg_family(:,:)
        real(rk), allocatable :: active1(:), active2(:), T(:), Tgrid(:,:), dT_dX(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk]
        family_knot2 = [-1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk]
        do jj = 1, 6
            do ii = 1, 6
                point = (jj-1)*6 + ii
                family_Xc(point,:) = [&
                    real(ii-1,rk),&
                    real(jj-1,rk),&
                    0.02_rk*real((ii-1)*(jj-1),rk)]
            end do
        end do

        call family_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call family_surface%create(res1=7, res2=7)
        Xt1 = family_surface%get_Xt(1)
        Xt2 = family_surface%get_Xt(2)
        call ndgrid(Xt1, Xt2, Xt_grid)
        Xg_reference = family_surface%get_Xg()
        elem = family_surface%get_elem_Xg_vis()
        call family_surface%basis(res1=3, res2=3, Tgc=Tgrid)
        elem = family_surface%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,6,2)
        active2 = active_knots(family_knot2,6,2)
        xi = [0.5_rk*(active1(1)+active1(2)), 0.5_rk*(active2(1)+active2(2))]
        call family_surface%basis(xi, T)
        elem = family_surface%cmp_elem()
        call family_surface%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = T,&
            dTgc_dXg = dT_dX,&
            dA       = dA,&
            ngauss   = surface_quadrature_order_33)

        call refined_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call refined_surface%insert_knots(1, [xi(1)], [1])
        call refined_surface%insert_knots(2, [xi(2)], [1])
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()
        call family_surface%err%print()
        call refined_surface%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families insertion",&
            res      = Xg_family,&
            expected = Xg_reference,&
            tol      = REFINEMENT_TOL,&
            msg      = "Mixed knot families insertion is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0189


    subroutine forcad_nurbs_surface_0190(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_33(2) = [3,3]
        real(rk), parameter :: REFINEMENT_TOL = 2.0e-10_rk
        type(nurbs_surface) :: family_surface, refined_surface
        real(rk) :: family_knot1(9)
        real(rk) :: family_knot2(9)
        real(rk) :: family_Xc(36,3)
        real(rk) :: xi(2)
        real(rk) :: dA
        real(rk), allocatable :: Xt1(:), Xt2(:), Xt_grid(:,:), Xg_reference(:,:), Xg_family(:,:)
        real(rk), allocatable :: active1(:), active2(:), T(:), Tgrid(:,:), dT_dX(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk]
        family_knot2 = [-1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk]
        do jj = 1, 6
            do ii = 1, 6
                point = (jj-1)*6 + ii
                family_Xc(point,:) = [&
                    real(ii-1,rk),&
                    real(jj-1,rk),&
                    0.02_rk*real((ii-1)*(jj-1),rk)]
            end do
        end do

        call family_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call family_surface%create(res1=7, res2=7)
        Xt1 = family_surface%get_Xt(1)
        Xt2 = family_surface%get_Xt(2)
        call ndgrid(Xt1, Xt2, Xt_grid)
        Xg_reference = family_surface%get_Xg()
        elem = family_surface%get_elem_Xg_vis()
        call family_surface%basis(res1=3, res2=3, Tgc=Tgrid)
        elem = family_surface%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,6,2)
        active2 = active_knots(family_knot2,6,2)
        xi = [0.5_rk*(active1(1)+active1(2)), 0.5_rk*(active2(1)+active2(2))]
        call family_surface%basis(xi, T)
        elem = family_surface%cmp_elem()
        call family_surface%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = T,&
            dTgc_dXg = dT_dX,&
            dA       = dA,&
            ngauss   = surface_quadrature_order_33)

        call refined_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call refined_surface%insert_knots(1, [xi(1)], [1])
        call refined_surface%insert_knots(2, [xi(2)], [1])
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()
        call refined_surface%remove_knots(1, [xi(1)], [1])
        call refined_surface%remove_knots(2, [xi(2)], [1])
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()
        call family_surface%err%print()
        call refined_surface%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families removal",&
            res      = Xg_family,&
            expected = Xg_reference,&
            tol      = REFINEMENT_TOL,&
            msg      = "Mixed knot families removal is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0190


    subroutine forcad_nurbs_surface_0191(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_33(2) = [3,3]
        real(rk), parameter :: REFINEMENT_TOL = 2.0e-10_rk
        type(nurbs_surface) :: family_surface, refined_surface
        real(rk) :: family_knot1(9)
        real(rk) :: family_knot2(9)
        real(rk) :: family_Xc(36,3)
        real(rk) :: xi(2)
        real(rk) :: dA
        real(rk), allocatable :: Xt1(:), Xt2(:), Xt_grid(:,:), Xg_reference(:,:), Xg_family(:,:)
        real(rk), allocatable :: active1(:), active2(:), T(:), Tgrid(:,:), dT_dX(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii
        integer :: jj
        integer :: point

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk]
        family_knot2 = [-1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk]
        do jj = 1, 6
            do ii = 1, 6
                point = (jj-1)*6 + ii
                family_Xc(point,:) = [&
                    real(ii-1,rk),&
                    real(jj-1,rk),&
                    0.02_rk*real((ii-1)*(jj-1),rk)]
            end do
        end do

        call family_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call family_surface%create(res1=7, res2=7)
        Xt1 = family_surface%get_Xt(1)
        Xt2 = family_surface%get_Xt(2)
        call ndgrid(Xt1, Xt2, Xt_grid)
        Xg_reference = family_surface%get_Xg()
        elem = family_surface%get_elem_Xg_vis()
        call family_surface%basis(res1=3, res2=3, Tgc=Tgrid)
        elem = family_surface%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,6,2)
        active2 = active_knots(family_knot2,6,2)
        xi = [0.5_rk*(active1(1)+active1(2)), 0.5_rk*(active2(1)+active2(2))]
        call family_surface%basis(xi, T)
        elem = family_surface%cmp_elem()
        call family_surface%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = T,&
            dTgc_dXg = dT_dX,&
            dA       = dA,&
            ngauss   = surface_quadrature_order_33)

        call refined_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call refined_surface%insert_knots(1, [xi(1)], [1])
        call refined_surface%insert_knots(2, [xi(2)], [1])
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()
        call refined_surface%remove_knots(1, [xi(1)], [1])
        call refined_surface%remove_knots(2, [xi(2)], [1])
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()

        call refined_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call refined_surface%elevate_degree(1,1)
        call refined_surface%elevate_degree(2,1)
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()
        call family_surface%err%print()
        call refined_surface%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families elevation",&
            res      = Xg_family,&
            expected = Xg_reference,&
            tol      = REFINEMENT_TOL,&
            msg      = "Mixed knot families elevation is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0191


    subroutine forcad_nurbs_surface_0192(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_33(2) = [3,3]
        real(rk), parameter :: FIT_TOL = 3.0e-6_rk
        type(nurbs_surface) :: family_surface, refined_surface
        real(rk) :: family_knot1(9), family_knot2(9), family_Xc(36,3), xi(2), dA, fit_error
        real(rk), allocatable :: Xt1(:), Xt2(:), Xt_grid(:,:), Xg_reference(:,:), Xg_family(:,:)
        real(rk), allocatable :: active1(:), active2(:), T(:), Tgrid(:,:), dT_dX(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii, jj, point, ndata(2)

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk]
        family_knot2 = [-1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk]
        do jj = 1, 6
            do ii = 1, 6
                point = (jj-1)*6 + ii
                family_Xc(point,:) = [&
                    real(ii-1,rk),&
                    real(jj-1,rk),&
                    0.02_rk*real((ii-1)*(jj-1),rk)]
            end do
        end do

        call family_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call family_surface%create(res1=7, res2=7)
        Xt1 = family_surface%get_Xt(1)
        Xt2 = family_surface%get_Xt(2)
        call ndgrid(Xt1, Xt2, Xt_grid)
        Xg_reference = family_surface%get_Xg()
        elem = family_surface%get_elem_Xg_vis()
        call family_surface%basis(res1=3, res2=3, Tgc=Tgrid)
        elem = family_surface%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,6,2)
        active2 = active_knots(family_knot2,6,2)
        xi = [0.5_rk*(active1(1)+active1(2)), 0.5_rk*(active2(1)+active2(2))]
        call family_surface%basis(xi, T)
        elem = family_surface%cmp_elem()
        call family_surface%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = T,&
            dTgc_dXg = dT_dX,&
            dA       = dA,&
            ngauss   = surface_quadrature_order_33)

        call refined_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call refined_surface%insert_knots(1, [xi(1)], [1])
        call refined_surface%insert_knots(2, [xi(2)], [1])
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()
        call refined_surface%remove_knots(1, [xi(1)], [1])
        call refined_surface%remove_knots(2, [xi(2)], [1])
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()

        call refined_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call refined_surface%elevate_degree(1,1)
        call refined_surface%elevate_degree(2,1)
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()

        ndata = [size(Xt1),size(Xt2)]
        call refined_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call refined_surface%lsq_fit_bspline(Xt_grid, Xg_reference, ndata)
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()
        fit_error = norm2(Xg_family-Xg_reference)/max(1.0_rk,norm2(Xg_reference))
        call family_surface%err%print()
        call refined_surface%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families B-spline fit",&
            res      = fit_error,&
            expected = 0.0_rk,&
            tol      = FIT_TOL,&
            msg      = "Mixed knot families B-spline fit is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0192


    subroutine forcad_nurbs_surface_0193(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: surface_quadrature_order_33(2) = [3,3]
        real(rk), parameter :: FIT_TOL = 3.0e-6_rk
        type(nurbs_surface) :: family_surface, refined_surface
        real(rk) :: family_knot1(9), family_knot2(9), family_Xc(36,3), xi(2), dA, fit_error
        real(rk), allocatable :: Xt1(:), Xt2(:), Xt_grid(:,:), Xg_reference(:,:), Xg_family(:,:)
        real(rk), allocatable :: active1(:), active2(:), T(:), Tgrid(:,:), dT_dX(:,:)
        integer, allocatable :: elem(:,:)
        integer :: ii, jj, point, ndata(2)

        family_knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk]
        family_knot2 = [-1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk]
        do jj = 1, 6
            do ii = 1, 6
                point = (jj-1)*6 + ii
                family_Xc(point,:) = [&
                    real(ii-1,rk),&
                    real(jj-1,rk),&
                    0.02_rk*real((ii-1)*(jj-1),rk)]
            end do
        end do

        call family_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call family_surface%create(res1=7, res2=7)
        Xt1 = family_surface%get_Xt(1)
        Xt2 = family_surface%get_Xt(2)
        call ndgrid(Xt1, Xt2, Xt_grid)
        Xg_reference = family_surface%get_Xg()
        elem = family_surface%get_elem_Xg_vis()
        call family_surface%basis(res1=3, res2=3, Tgc=Tgrid)
        elem = family_surface%get_elem_Xg_vis()

        active1 = active_knots(family_knot1,6,2)
        active2 = active_knots(family_knot2,6,2)
        xi = [0.5_rk*(active1(1)+active1(2)), 0.5_rk*(active2(1)+active2(2))]
        call family_surface%basis(xi, T)
        elem = family_surface%cmp_elem()
        call family_surface%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = T,&
            dTgc_dXg = dT_dX,&
            dA       = dA,&
            ngauss   = surface_quadrature_order_33)

        call refined_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call refined_surface%insert_knots(1, [xi(1)], [1])
        call refined_surface%insert_knots(2, [xi(2)], [1])
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()
        call refined_surface%remove_knots(1, [xi(1)], [1])
        call refined_surface%remove_knots(2, [xi(2)], [1])
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()

        call refined_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call refined_surface%elevate_degree(1,1)
        call refined_surface%elevate_degree(2,1)
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()

        ndata = [size(Xt1),size(Xt2)]
        call refined_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call refined_surface%lsq_fit_bspline(Xt_grid, Xg_reference, ndata)
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()
        fit_error = norm2(Xg_family-Xg_reference)/max(1.0_rk,norm2(Xg_reference))

        call refined_surface%set(family_knot1, family_knot2, family_Xc, degree=[2,2])
        call refined_surface%lsq_fit_nurbs(&
            Xt    = Xt_grid,&
            Xdata = Xg_reference,&
            ndata = ndata,&
            maxit = 30,&
            tol   = sqrt(epsilon(1.0_rk)))
        call refined_surface%create(Xt=Xt_grid)
        Xg_family = refined_surface%get_Xg()
        fit_error = norm2(Xg_family-Xg_reference)/max(1.0_rk,norm2(Xg_reference))
        call family_surface%err%print()
        call refined_surface%err%print()

        call ut%test(ti)%check(&
            name     = "mixed knot families NURBS fit",&
            res      = fit_error,&
            expected = 0.0_rk,&
            tol      = FIT_TOL,&
            msg      = "Mixed knot families NURBS fit is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0193


    subroutine forcad_nurbs_surface_0194(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: full_break_surface
        real(rk) :: knot1_full(9), knot2_full(4), Xc_full(12,3), Wc_full(12), Xt_eval(10,2)
        real(rk) :: xline(6), zline(6), wline(6), u_eval(5), v_eval(2)
        real(rk), allocatable :: Xg_before(:,:), Xg_after(:,:)
        integer :: ii, jj, point

        knot1_full = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 0.5_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot2_full = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xline = [0.0_rk, 0.2_rk, 0.5_rk, 0.5_rk, 0.8_rk, 1.0_rk]
        zline = [0.0_rk, 0.8_rk, 0.0_rk, 2.0_rk, 2.8_rk, 2.0_rk]
        wline = [1.0_rk, 0.7_rk, 1.2_rk, 0.9_rk, 1.3_rk, 1.0_rk]
        u_eval = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk]
        v_eval = [0.2_rk, 0.8_rk]
        do jj = 1, 2
            do ii = 1, 6
                point = (jj-1)*6 + ii
                Xc_full(point,:) = [xline(ii),real(jj-1,rk),zline(ii)+0.1_rk*real(jj-1,rk)]
                Wc_full(point) = wline(ii)*(1.0_rk+0.05_rk*real(jj-1,rk))
            end do
        end do
        do jj = 1, 2
            do ii = 1, 5
                point = (jj-1)*5 + ii
                Xt_eval(point,:) = [u_eval(ii),v_eval(jj)]
            end do
        end do

        call full_break_surface%set(knot1_full, knot2_full, Xc_full, Wc_full, degree=[2,1])
        call full_break_surface%create(Xt=Xt_eval)
        Xg_before = full_break_surface%get_Xg()
        call full_break_surface%elevate_degree(1,1)
        call full_break_surface%create(Xt=Xt_eval)
        Xg_after = full_break_surface%get_Xg()
        call full_break_surface%err%print()

        call ut%test(ti)%check(&
            name     = "full-break elevation geometry",&
            res      = Xg_after,&
            expected = Xg_before,&
            tol      = 2.0e-11_rk,&
            msg      = "Full-break elevation geometry is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0194


    subroutine forcad_nurbs_surface_0195(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: full_break_surface
        real(rk) :: knot1_full(9), knot2_full(4), Xc_full(12,3), Wc_full(12), Xt_eval(10,2)
        real(rk) :: xline(6), zline(6), wline(6), u_eval(5), v_eval(2)
        real(rk), allocatable :: Xg_before(:,:), Xg_after(:,:)
        integer :: ii, jj, point

        knot1_full = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 0.5_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot2_full = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xline = [0.0_rk, 0.2_rk, 0.5_rk, 0.5_rk, 0.8_rk, 1.0_rk]
        zline = [0.0_rk, 0.8_rk, 0.0_rk, 2.0_rk, 2.8_rk, 2.0_rk]
        wline = [1.0_rk, 0.7_rk, 1.2_rk, 0.9_rk, 1.3_rk, 1.0_rk]
        u_eval = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk]
        v_eval = [0.2_rk, 0.8_rk]
        do jj = 1, 2
            do ii = 1, 6
                point = (jj-1)*6 + ii
                Xc_full(point,:) = [xline(ii),real(jj-1,rk),zline(ii)+0.1_rk*real(jj-1,rk)]
                Wc_full(point) = wline(ii)*(1.0_rk+0.05_rk*real(jj-1,rk))
            end do
        end do
        do jj = 1, 2
            do ii = 1, 5
                point = (jj-1)*5 + ii
                Xt_eval(point,:) = [u_eval(ii),v_eval(jj)]
            end do
        end do

        call full_break_surface%set(knot1_full, knot2_full, Xc_full, Wc_full, degree=[2,1])
        call full_break_surface%create(Xt=Xt_eval)
        Xg_before = full_break_surface%get_Xg()
        call full_break_surface%elevate_degree(1,1)
        call full_break_surface%create(Xt=Xt_eval)
        Xg_after = full_break_surface%get_Xg()
        call full_break_surface%err%print()

        call ut%test(ti)%check(&
            name     = "full-break elevation degree",&
            res      = full_break_surface%get_degree(),&
            expected = [3,1],&
            msg      = "Full-break elevation degree is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0195


    subroutine forcad_nurbs_surface_0196(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NGAUSS(2) = [2,2]
        real(rk), parameter :: H = 1.0e-10_rk
        type(nurbs_surface) :: tiny_surface_2d
        real(rk) :: linear_knot(4)
        real(rk) :: Xc2_tiny(4,2)
        real(rk) :: Xc3_tiny(4,3)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:)
        integer :: i1, i2, point

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        do i2 = 1, 2
            do i1 = 1, 2
                point = (i2-1)*2 + i1
                Xc2_tiny(point,:) = [real(i1-1,rk)*H,real(i2-1,rk)*H]
                Xc3_tiny(point,:) = [real(i1-1,rk)*H,real(i2-1,rk)*H,0.0_rk]
            end do
        end do

        call tiny_surface_2d%set(linear_knot, linear_knot, Xc2_tiny, degree=[1,1])
        call tiny_surface_2d%ansatz(1, 1, T, dT_dX, dA, ngauss=NGAUSS)
        call tiny_surface_2d%err%print()

        call ut%test(ti)%check(&
            name     = "tiny planar ansatz diagnostic",&
            res      = tiny_surface_2d%err%ok,&
            expected = .true.,&
            msg      = "A regular surface remains valid when its physical scale is small",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0196


    subroutine forcad_nurbs_surface_0197(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NGAUSS(2) = [2,2]
        real(rk), parameter :: H = 1.0e-10_rk
        type(nurbs_surface) :: tiny_surface_2d
        real(rk) :: linear_knot(4)
        real(rk) :: Xc2_tiny(4,2)
        real(rk) :: Xc3_tiny(4,3)
        real(rk) :: dA
        real(rk), allocatable :: T(:), dT_dX(:,:)
        integer :: i1, i2, point

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        do i2 = 1, 2
            do i1 = 1, 2
                point = (i2-1)*2 + i1
                Xc2_tiny(point,:) = [real(i1-1,rk)*H,real(i2-1,rk)*H]
                Xc3_tiny(point,:) = [real(i1-1,rk)*H,real(i2-1,rk)*H,0.0_rk]
            end do
        end do

        call tiny_surface_2d%set(linear_knot, linear_knot, Xc2_tiny, degree=[1,1])
        call tiny_surface_2d%ansatz(1, 1, T, dT_dX, dA, ngauss=NGAUSS)
        call tiny_surface_2d%err%print()

        call ut%test(ti)%check(&
            name     = "tiny planar ansatz measure",&
            res      = dA > 0.0_rk,&
            expected = .true.,&
            msg      = "A small regular planar surface must retain positive measure",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0197


    subroutine forcad_nurbs_surface_0198(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NGAUSS(2) = [2,2]
        real(rk), parameter :: H = 1.0e-10_rk
        type(nurbs_surface) :: tiny_surface_2d
        real(rk) :: linear_knot(4), Xc2_tiny(4,2), Xc3_tiny(4,3), area, dA
        real(rk), allocatable :: T(:), dT_dX(:,:)
        integer :: i1, i2, point

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        do i2 = 1, 2
            do i1 = 1, 2
                point = (i2-1)*2 + i1
                Xc2_tiny(point,:) = [real(i1-1,rk)*H,real(i2-1,rk)*H]
                Xc3_tiny(point,:) = [real(i1-1,rk)*H,real(i2-1,rk)*H,0.0_rk]
            end do
        end do

        call tiny_surface_2d%set(linear_knot, linear_knot, Xc2_tiny, degree=[1,1])
        call tiny_surface_2d%ansatz(1, 1, T, dT_dX, dA, ngauss=NGAUSS)
        call tiny_surface_2d%cmp_area(area, ngauss=NGAUSS)
        call tiny_surface_2d%err%print()

        call ut%test(ti)%check(&
            name     = "tiny planar area",&
            res      = area/(H*H),&
            expected = 1.0_rk,&
            tol      = 2.0e-12_rk,&
            msg      = "Area integration must be relative-scale robust",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0198


    subroutine forcad_nurbs_surface_0199(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NGAUSS(2) = [2,2]
        real(rk), parameter :: H = 1.0e-10_rk
        type(nurbs_surface) :: tiny_surface_2d, tiny_surface_3d
        real(rk) :: linear_knot(4), Xc2_tiny(4,2), Xc3_tiny(4,3), area, dA
        real(rk), allocatable :: T(:), dT_dX(:,:)
        integer :: i1, i2, point

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        do i2 = 1, 2
            do i1 = 1, 2
                point = (i2-1)*2 + i1
                Xc2_tiny(point,:) = [real(i1-1,rk)*H,real(i2-1,rk)*H]
                Xc3_tiny(point,:) = [real(i1-1,rk)*H,real(i2-1,rk)*H,0.0_rk]
            end do
        end do

        call tiny_surface_2d%set(linear_knot, linear_knot, Xc2_tiny, degree=[1,1])
        call tiny_surface_2d%ansatz(1, 1, T, dT_dX, dA, ngauss=NGAUSS)
        call tiny_surface_2d%cmp_area(area, ngauss=NGAUSS)

        call tiny_surface_3d%set(linear_knot, linear_knot, Xc3_tiny, degree=[1,1])
        call tiny_surface_3d%ansatz(1, 1, T, dT_dX, dA, ngauss=NGAUSS)
        call tiny_surface_2d%err%print()
        call tiny_surface_3d%err%print()

        call ut%test(ti)%check(&
            name     = "tiny embedded ansatz diagnostic",&
            res      = tiny_surface_3d%err%ok,&
            expected = .true.,&
            msg      = "A small regular embedded surface must remain valid",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0199


    subroutine forcad_nurbs_surface_0200(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NGAUSS(2) = [2,2]
        real(rk), parameter :: H = 1.0e-10_rk
        type(nurbs_surface) :: tiny_surface_2d, tiny_surface_3d
        real(rk) :: linear_knot(4), Xc2_tiny(4,2), Xc3_tiny(4,3), area, dA
        real(rk), allocatable :: T(:), dT_dX(:,:)
        integer :: i1, i2, point

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        do i2 = 1, 2
            do i1 = 1, 2
                point = (i2-1)*2 + i1
                Xc2_tiny(point,:) = [real(i1-1,rk)*H,real(i2-1,rk)*H]
                Xc3_tiny(point,:) = [real(i1-1,rk)*H,real(i2-1,rk)*H,0.0_rk]
            end do
        end do

        call tiny_surface_2d%set(linear_knot, linear_knot, Xc2_tiny, degree=[1,1])
        call tiny_surface_2d%ansatz(1, 1, T, dT_dX, dA, ngauss=NGAUSS)
        call tiny_surface_2d%cmp_area(area, ngauss=NGAUSS)

        call tiny_surface_3d%set(linear_knot, linear_knot, Xc3_tiny, degree=[1,1])
        call tiny_surface_3d%ansatz(1, 1, T, dT_dX, dA, ngauss=NGAUSS)
        call tiny_surface_2d%err%print()
        call tiny_surface_3d%err%print()

        call ut%test(ti)%check(&
            name     = "tiny embedded ansatz measure",&
            res      = dA > 0.0_rk,&
            expected = .true.,&
            msg      = "A small regular embedded surface must retain positive measure",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0200


    subroutine forcad_nurbs_surface_0201(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NGAUSS(2) = [2,2]
        real(rk), parameter :: H = 1.0e-10_rk
        type(nurbs_surface) :: tiny_surface_2d, tiny_surface_3d
        real(rk) :: linear_knot(4), Xc2_tiny(4,2), Xc3_tiny(4,3), area, dA
        real(rk), allocatable :: T(:), dT_dX(:,:)
        integer :: i1, i2, point

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        do i2 = 1, 2
            do i1 = 1, 2
                point = (i2-1)*2 + i1
                Xc2_tiny(point,:) = [real(i1-1,rk)*H,real(i2-1,rk)*H]
                Xc3_tiny(point,:) = [real(i1-1,rk)*H,real(i2-1,rk)*H,0.0_rk]
            end do
        end do

        call tiny_surface_2d%set(linear_knot, linear_knot, Xc2_tiny, degree=[1,1])
        call tiny_surface_2d%ansatz(1, 1, T, dT_dX, dA, ngauss=NGAUSS)
        call tiny_surface_2d%cmp_area(area, ngauss=NGAUSS)

        call tiny_surface_3d%set(linear_knot, linear_knot, Xc3_tiny, degree=[1,1])
        call tiny_surface_3d%ansatz(1, 1, T, dT_dX, dA, ngauss=NGAUSS)
        call tiny_surface_3d%cmp_area(area, ngauss=NGAUSS)
        call tiny_surface_2d%err%print()
        call tiny_surface_3d%err%print()

        call ut%test(ti)%check(&
            name     = "tiny embedded area",&
            res      = area/(H*H),&
            expected = 1.0_rk,&
            tol      = 2.0e-12_rk,&
            msg      = "Embedded area integration must be relative-scale robust",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0201


    subroutine forcad_nurbs_surface_0202(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: periodic_surface
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(12,3), periodic_Wc(12), ring_Wc(6)
        integer :: i1
        integer :: i2
        integer :: n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                periodic_Xc(n,:) = [ring(i1,1), ring(i1,2), real(i2-1,rk)]
                periodic_Wc(n) = ring_Wc(i1)
            end do
        end do

        call periodic_surface%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1],&
            wrap_parameters = [.true.,.false.])
        call periodic_surface%err%print()

        call ut%test(ti)%check(&
            name     = "periodic surface parameter wrapping",&
            res      = periodic_surface%get_parameter_wrapping(1),&
            expected = .true.,&
            msg      = "Periodic surface parameter wrapping is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0202


    subroutine forcad_nurbs_surface_0203(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: periodic_surface
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(12,3), periodic_Wc(12), ring_Wc(6)
        integer :: i1
        integer :: i2
        integer :: n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                periodic_Xc(n,:) = [ring(i1,1), ring(i1,2), real(i2-1,rk)]
                periodic_Wc(n) = ring_Wc(i1)
            end do
        end do

        call periodic_surface%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1],&
            wrap_parameters = [.true.,.false.])
        call periodic_surface%err%print()

        call ut%test(ti)%check(&
            name     = "bounded surface second direction",&
            res      = periodic_surface%get_parameter_wrapping(2),&
            expected = .false.,&
            msg      = "Enabling one surface direction must not wrap another direction",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0203


    subroutine forcad_nurbs_surface_0204(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: PERIODIC_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: periodic_surface
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(12,3), periodic_Wc(12), ring_Wc(6)
        integer :: i1
        integer :: i2
        integer :: n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                periodic_Xc(n,:) = [ring(i1,1), ring(i1,2), real(i2-1,rk)]
                periodic_Wc(n) = ring_Wc(i1)
            end do
        end do

        call periodic_surface%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1],&
            wrap_parameters = [.true.,.false.])
        call periodic_surface%err%print()

        call ut%test(ti)%check(&
            name     = "wrapped periodic NURBS surface geometry",&
            res      = periodic_surface%cmp_Xg([10.75_rk,0.35_rk]),&
            expected = periodic_surface%cmp_Xg([2.75_rk,0.35_rk]),&
            tol      = PERIODIC_TOL,&
            msg      = "Wrapped periodic NURBS surface geometry is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0204


    subroutine forcad_nurbs_surface_0205(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: PERIODIC_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: periodic_surface
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(12,3), periodic_Wc(12), ring_Wc(6)
        real(rk), allocatable :: T_start(:), T_wrapped(:), dT_start(:,:), dT_wrapped(:,:)
        integer :: i1
        integer :: i2
        integer :: n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                periodic_Xc(n,:) = [ring(i1,1), ring(i1,2), real(i2-1,rk)]
                periodic_Wc(n) = ring_Wc(i1)
            end do
        end do

        call periodic_surface%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1],&
            wrap_parameters = [.true.,.false.])
        call periodic_surface%derivative([2.75_rk,0.35_rk], dT_start, T_start)
        call periodic_surface%derivative([-1.25_rk,0.35_rk], dT_wrapped, T_wrapped)
        call periodic_surface%err%print()

        call ut%test(ti)%check(&
            name     = "wrapped periodic NURBS surface basis",&
            res      = T_wrapped,&
            expected = T_start,&
            tol      = PERIODIC_TOL,&
            msg      = "Wrapped periodic NURBS surface basis is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0205


    subroutine forcad_nurbs_surface_0206(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: PERIODIC_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: periodic_surface
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(12,3), periodic_Wc(12), ring_Wc(6)
        real(rk), allocatable :: T_start(:), T_wrapped(:), dT_start(:,:), dT_wrapped(:,:)
        integer :: i1
        integer :: i2
        integer :: n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                periodic_Xc(n,:) = [ring(i1,1), ring(i1,2), real(i2-1,rk)]
                periodic_Wc(n) = ring_Wc(i1)
            end do
        end do

        call periodic_surface%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1],&
            wrap_parameters = [.true.,.false.])
        call periodic_surface%derivative([2.75_rk,0.35_rk], dT_start, T_start)
        call periodic_surface%derivative([-1.25_rk,0.35_rk], dT_wrapped, T_wrapped)
        call periodic_surface%err%print()

        call ut%test(ti)%check(&
            name     = "wrapped periodic NURBS surface derivatives",&
            res      = dT_wrapped,&
            expected = dT_start,&
            tol      = PERIODIC_TOL,&
            msg      = "Wrapped periodic NURBS surface derivatives are incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0206


    subroutine forcad_nurbs_surface_0207(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: periodic_surface
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(12,3), periodic_Wc(12), ring_Wc(6)
        real(rk), allocatable :: T_start(:), T_wrapped(:), dT_start(:,:), dT_wrapped(:,:)
        integer :: i1
        integer :: i2
        integer :: n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                periodic_Xc(n,:) = [ring(i1,1), ring(i1,2), real(i2-1,rk)]
                periodic_Wc(n) = ring_Wc(i1)
            end do
        end do

        call periodic_surface%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1],&
            wrap_parameters = [.true.,.false.])
        call periodic_surface%derivative([2.75_rk,0.35_rk], dT_start, T_start)
        call periodic_surface%derivative([-1.25_rk,0.35_rk], dT_wrapped, T_wrapped)
        call periodic_surface%insert_knots(1, [3.5_rk], [1])
        call periodic_surface%err%print()

        call ut%test(ti)%check(&
            name     = "parameter wrapping after surface refinement",&
            res      = periodic_surface%get_parameter_wrapping(1),&
            expected = .true.,&
            msg      = "Parameter wrapping after surface refinement is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0207


    subroutine forcad_nurbs_surface_0208(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: PERIODIC_TOL = 4096.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: periodic_surface
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(12,3), periodic_Wc(12), ring_Wc(6)
        real(rk), allocatable :: T_start(:), T_wrapped(:), dT_start(:,:), dT_wrapped(:,:)
        integer :: i1
        integer :: i2
        integer :: n

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                periodic_Xc(n,:) = [ring(i1,1), ring(i1,2), real(i2-1,rk)]
                periodic_Wc(n) = ring_Wc(i1)
            end do
        end do

        call periodic_surface%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1],&
            wrap_parameters = [.true.,.false.])
        call periodic_surface%derivative([2.75_rk,0.35_rk], dT_start, T_start)
        call periodic_surface%derivative([-1.25_rk,0.35_rk], dT_wrapped, T_wrapped)
        call periodic_surface%insert_knots(1, [3.5_rk], [1])
        call periodic_surface%err%print()

        call ut%test(ti)%check(&
            name     = "periodic surface refinement wrapping",&
            res      = periodic_surface%cmp_Xg([10.75_rk,0.35_rk]),&
            expected = periodic_surface%cmp_Xg([2.75_rk,0.35_rk]),&
            tol      = PERIODIC_TOL,&
            msg      = "A refined periodic surface must retain repeated evaluation",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0208


    subroutine forcad_nurbs_surface_0209(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: periodic_surface
        real(rk) :: periodic_knot(9), linear_knot(4), ring(6,2), periodic_Xc(12,3), periodic_Wc(12), ring_Wc(6)
        real(rk), allocatable :: T_start(:), T_wrapped(:), dT_start(:,:), dT_wrapped(:,:)
        character(len=512) :: iges_line
        integer :: i1, i2, n, iges_unit, ios
        logical :: found

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        ring_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                periodic_Xc(n,:) = [ring(i1,1), ring(i1,2), real(i2-1,rk)]
                periodic_Wc(n) = ring_Wc(i1)
            end do
        end do

        call periodic_surface%set(&
            knot1    = periodic_knot,&
            knot2    = linear_knot,&
            Xc       = periodic_Xc,&
            Wc       = periodic_Wc,&
            degree   = [2,1],&
            wrap_parameters = [.true.,.false.])
        call periodic_surface%derivative([2.75_rk,0.35_rk], dT_start, T_start)
        call periodic_surface%derivative([-1.25_rk,0.35_rk], dT_wrapped, T_wrapped)
        call periodic_surface%export_iges("vtk/forcad_periodic_surface.iges")
        found = .false.
        open(newunit=iges_unit, file="vtk/forcad_periodic_surface.iges", status="old", action="read", iostat=ios)
        if (ios == 0) then
            do
                read(iges_unit, "(A)", iostat=ios) iges_line
                if (ios /= 0) exit
                if (index(iges_line, "128,5,1,2,1,1,0,0,1,0") > 0) then
                    found = .true.
                    exit
                end if
            end do
            close(iges_unit)
        end if
        call periodic_surface%err%print()

        call ut%test(ti)%check(&
            name     = "periodic surface IGES properties",&
            res      = found,&
            expected = .true.,&
            msg      = "Periodic surface IGES properties are incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0209


    subroutine forcad_nurbs_surface_0210(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: export_surface, copied_surface
        real(rk) :: linear_knot(4), export_Xc(4,3)

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        export_Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        export_Xc(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        export_Xc(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        export_Xc(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]

        call export_surface%set(linear_knot, linear_knot, export_Xc, degree=[1,1])
        copied_surface = export_surface
        call export_surface%modify_Xc(5.0_rk, 1, 1)
        call export_surface%err%print()
        call copied_surface%err%print()

        call ut%test(ti)%check(&
            name     = "surface intrinsic assignment deep copy",&
            res      = copied_surface%get_Xc(1,1),&
            expected = 0.0_rk,&
            tol      = 0.0_rk,&
            msg      = "Surface intrinsic assignment deep copy is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0210


    subroutine forcad_nurbs_surface_0211(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: vtkfile = "vtk/forcad_test_surface_Xth_in_Xg.vtk"
        character(len=*), parameter :: vtk_header = "# vtk DataFile Version 2.0"
        type(nurbs_surface) :: export_surface, copied_surface
        real(rk) :: linear_knot(4), export_Xc(4,3)
        character(len=len(vtk_header)) :: actual_header
        integer :: file_unit, io_status
        logical :: file_exists

        linear_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        export_Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
        export_Xc(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        export_Xc(3,:) = [0.0_rk, 1.0_rk, 0.0_rk]
        export_Xc(4,:) = [1.0_rk, 1.0_rk, 0.0_rk]

        call export_surface%set(linear_knot, linear_knot, export_Xc, degree=[1,1])
        copied_surface = export_surface
        call export_surface%modify_Xc(5.0_rk, 1, 1)

        open(newunit=file_unit, file=vtkfile, status="replace")
        close(file_unit, status="delete")
        call copied_surface%export_Xth_in_Xg(vtkfile, res=3)
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
        call export_surface%err%print()
        call copied_surface%err%print()

        call ut%test(ti)%check(&
            name     = "surface export Xth in Xg",&
            res      = file_exists .and. io_status == 0 .and. actual_header == vtk_header,&
            expected = .true.,&
            msg      = "Surface export Xth in Xg is incorrect.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0211


    subroutine forcad_nurbs_surface_0212(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot1(8) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: knot2(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(10,3) = reshape([&
            0.0_rk,0.0_rk,1.0_rk,0.0_rk,0.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,0.0_rk,0.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,&
            1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [10,3])
        type(nurbs_surface) :: surface
        real(rk), allocatable :: knot_after(:), Xc_after(:,:)
        logical :: preserved

        call surface%set(knot1, knot2, Xc)
        call surface%remove_knots(1, [0.5_rk], [1])
        call surface%err%print()

        knot_after = surface%get_knot(1)
        Xc_after = surface%get_Xc()
        preserved = surface%err%ok .and. size(knot_after) == size(knot1) .and. &
            size(Xc_after,1) == size(Xc,1) .and. size(Xc_after,2) == size(Xc,2)
        if (preserved) preserved = maxval(abs(knot_after-knot1)) <= epsilon(1.0_rk) .and. &
            maxval(abs(Xc_after-Xc)) <= epsilon(1.0_rk)

        call ut%test(ti)%check(&
            name     = "surface rejects geometry-changing removal",&
            res      = preserved,&
            expected = .true.,&
            msg      = "Surface removal changed nonremovable geometry.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0212


    subroutine forcad_nurbs_surface_0213(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface

        call surface%set_tetragon([1.0_rk,1.0_rk], [2,2])
        call surface%set(&
            Xth_dir1    = [0.0_rk,0.5_rk,1.0_rk],&
            Xth_dir2    = [0.0_rk,0.5_rk,1.0_rk],&
            degree      = [2,2],&
            continuity1 = [-1,1,-1],&
            continuity2 = [-1,1,-1])
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "set2 clears old surface controls",&
            res      = size(surface%get_Xc(),1),&
            expected = 0,&
            msg      = "Omitting Xc must remove the previous surface controls.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0213


    subroutine forcad_nurbs_surface_0214(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: filename = "vtk/forcad_wrapped_open_surface.iges"
        type(nurbs_surface) :: surface
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(4,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk], [4,3])
        character(len=512) :: iges_line
        integer :: iges_unit, ios
        logical :: found

        call surface%set(&
            knot1           = knot,&
            knot2           = knot,&
            Xc              = Xc,&
            degree          = [1,1],&
            wrap_parameters = [.true.,.true.])
        call surface%export_iges(filename)
        found = .false.
        open(newunit=iges_unit, file=filename, status="old", action="read", iostat=ios)
        if (ios == 0) then
            do
                read(iges_unit, "(A)", iostat=ios) iges_line
                if (ios /= 0) exit
                if (index(iges_line, "128,1,1,1,1,0,0,1,0,0") > 0) then
                    found = .true.
                    exit
                end if
            end do
            close(iges_unit)
        end if
        open(newunit=iges_unit, file=filename, status="old", iostat=ios)
        if (ios == 0) close(iges_unit, status="delete")
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "wrapped open surface IGES properties",&
            res      = found,&
            expected = .true.,&
            msg      = "Wrapping must not mark open surface directions closed.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0214


    subroutine forcad_nurbs_surface_0215(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: filename = "vtk/forcad_closed_ring.iges"
        type(nurbs_surface) :: surface
        character(len=512) :: iges_line
        integer :: iges_unit, ios
        logical :: found

        call surface%set_ring([0.0_rk,0.0_rk,0.0_rk], 1.0_rk, 2.0_rk)
        call surface%export_iges(filename)
        found = .false.
        open(newunit=iges_unit, file=filename, status="old", action="read", iostat=ios)
        if (ios == 0) then
            do
                read(iges_unit, "(A)", iostat=ios) iges_line
                if (ios /= 0) exit
                if (index(iges_line, "128,6,1,2,1,1,0,0,0,0") > 0) then
                    found = .true.
                    exit
                end if
            end do
            close(iges_unit)
        end if
        open(newunit=iges_unit, file=filename, status="old", iostat=ios)
        if (ios == 0) close(iges_unit, status="delete")
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "clamped ring IGES properties",&
            res      = found,&
            expected = .true.,&
            msg      = "A clamped ring must be closed only in its first direction.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0215


    subroutine forcad_nurbs_surface_0216(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface
        real(rk), parameter :: periodic_knot(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        real(rk), parameter :: linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: ring_Wc(6) = [1.0_rk,1.4_rk,0.8_rk,1.2_rk,1.0_rk,1.4_rk]
        real(rk) :: ring(6,2), Xc(12,3), Wc(12)
        integer :: i1, i2, n

        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                Xc(n,:) = [ring(i1,1),ring(i1,2),real(i2-1,rk)]
                Wc(n) = ring_Wc(i1)
            end do
        end do

        call surface%set(&
            knot1  = periodic_knot,&
            knot2  = linear_knot,&
            Xc     = Xc,&
            Wc     = Wc,&
            degree = [2,1])
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "surface directional topology",&
            res      = surface%get_parameter_topology(1) == "periodic" .and. &
                surface%get_parameter_topology(2) == "bounded",&
            expected = .true.,&
            msg      = "Surface parameter topologies were classified incorrectly.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0216


    subroutine forcad_nurbs_surface_0217(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface
        real(rk), parameter :: periodic_knot(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        real(rk), parameter :: linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        integer, parameter :: expected_map(12) = [1,2,3,4,1,2,5,6,7,8,5,6]
        real(rk) :: ring(6,2), Xc(12,3)
        integer, allocatable :: map(:)
        integer :: i1, i2, n

        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                Xc(n,:) = [ring(i1,1),ring(i1,2),real(i2-1,rk)]
            end do
        end do

        call surface%set(&
            knot1  = periodic_knot,&
            knot2  = linear_knot,&
            Xc     = Xc,&
            degree = [2,1])
        map = surface%cmp_dof_map()
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "surface periodic DOF map",&
            res      = map,&
            expected = expected_map,&
            msg      = "Repeated periodic surface layers were not identified.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0217


    subroutine forcad_nurbs_surface_0218(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface
        real(rk), parameter :: periodic_knot(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        real(rk), parameter :: linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: ring_Wc(6) = [1.0_rk,1.4_rk,0.8_rk,1.2_rk,1.0_rk,1.4_rk]
        real(rk) :: ring(6,2), Xc(12,3), Wc(12)
        integer :: i1, i2, n

        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                Xc(n,:) = [ring(i1,1),ring(i1,2),real(i2-1,rk)]
                Wc(n) = ring_Wc(i1)
            end do
        end do

        call surface%set(&
            knot1  = periodic_knot,&
            knot2  = linear_knot,&
            Xc     = Xc,&
            Wc     = Wc,&
            degree = [2,1])
        call surface%elevate_degree(1, 1)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "surface elevated periodic topology",&
            res      = surface%get_parameter_topology(1) == "periodic",&
            expected = .true.,&
            msg      = "Surface degree elevation lost periodic topology.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0218


    subroutine forcad_nurbs_surface_0219(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 32768.0_rk*epsilon(1.0_rk)
        type(nurbs_surface) :: surface
        real(rk), parameter :: periodic_knot(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        real(rk), parameter :: linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: ring_Wc(6) = [1.0_rk,1.4_rk,0.8_rk,1.2_rk,1.0_rk,1.4_rk]
        real(rk), parameter :: Xt(2,5) = reshape([&
            2.0_rk,0.0_rk,2.2_rk,0.3_rk,3.1_rk,0.5_rk,&
            4.4_rk,0.7_rk,6.0_rk,1.0_rk], [2,5])
        real(rk) :: ring(6,2), Xc(12,3), Wc(12)
        real(rk) :: actual(5,3), expected(5,3)
        integer :: i, i1, i2, n

        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                Xc(n,:) = [ring(i1,1),ring(i1,2),real(i2-1,rk)]
                Wc(n) = ring_Wc(i1)
            end do
        end do
        call surface%set(&
            knot1  = periodic_knot,&
            knot2  = linear_knot,&
            Xc     = Xc,&
            Wc     = Wc,&
            degree = [2,1])
        do i = 1, size(Xt,2)
            expected(i,:) = surface%cmp_Xg(Xt(:,i))
        end do
        call surface%elevate_degree(1, 1)
        do i = 1, size(Xt,2)
            actual(i,:) = surface%cmp_Xg(Xt(:,i))
        end do
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "surface periodic elevation geometry",&
            res      = actual,&
            expected = expected,&
            tol      = tol,&
            msg      = "Periodic degree elevation changed surface geometry.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0219


    subroutine forcad_nurbs_surface_0220(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface
        real(rk), parameter :: periodic_knot(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        real(rk), parameter :: linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk) :: ring(6,2), Xc(12,3)
        integer :: i1, i2, n

        ring(1,:) = [ 1.0_rk, 0.0_rk]
        ring(2,:) = [ 0.0_rk, 1.0_rk]
        ring(3,:) = [-1.0_rk, 0.0_rk]
        ring(4,:) = [ 0.0_rk,-1.0_rk]
        ring(5,:) = ring(1,:)
        ring(6,:) = ring(2,:)
        do i2 = 1, 2
            do i1 = 1, 6
                n = i1 + (i2-1)*6
                Xc(n,:) = [ring(i1,1),ring(i1,2),real(i2-1,rk)]
            end do
        end do

        call surface%set(&
            knot1  = periodic_knot,&
            knot2  = linear_knot,&
            Xc     = Xc,&
            degree = [2,1])
        call surface%remove_knots(1, [4.0_rk], [1])
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "periodic surface removal diagnostic",&
            res      = surface%err%ok,&
            expected = .false.,&
            msg      = "Periodic removal must not destroy surface topology.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0220


    subroutine forcad_nurbs_surface_0221(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(4,2) = reshape([&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk], [4,2])
        integer, parameter :: degree(2) = [1,1]
        integer, parameter :: ngauss(2) = [1,1]
        type(nurbs_surface) :: surface
        real(rk), allocatable :: Tgc(:), dTgc_dXg(:,:), d2Tgc_dXg2(:,:,:)
        real(rk) :: dA

        call surface%set(&
            knot1  = knot,&
            knot2  = knot,&
            Xc     = Xc,&
            degree = degree)
        call surface%ansatz(&
            ie          = 1,&
            ig          = 1,&
            Tgc         = Tgc,&
            dTgc_dXg    = dTgc_dXg,&
            dA          = dA,&
            ngauss      = ngauss,&
            d2Tgc_dXg2  = d2Tgc_dXg2,&
            strict      = .true.)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "strict surface rejects inversion",&
            res      = .not. surface%err%ok .and. abs(dA) <= tiny(1.0_rk) .and. &
                maxval(abs(dTgc_dXg)) <= tiny(1.0_rk) .and. &
                maxval(abs(d2Tgc_dXg2)) <= tiny(1.0_rk),&
            expected = .true.,&
            msg      = "Strict ansatz accepted an inverted surface map.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0221


    subroutine forcad_nurbs_surface_0222(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(4,2) = reshape([&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk], [4,2])
        integer, parameter :: degree(2) = [1,1]
        integer, parameter :: ngauss(2) = [1,1]
        type(nurbs_surface) :: surface
        real(rk), allocatable :: Tgc(:), dTgc_dXg(:,:)
        real(rk) :: dA

        call surface%set(&
            knot1  = knot,&
            knot2  = knot,&
            Xc     = Xc,&
            degree = degree)
        call surface%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = Tgc,&
            dTgc_dXg = dTgc_dXg,&
            dA       = dA,&
            ngauss   = ngauss)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "surface default measures inversion",&
            res      = surface%err%ok .and. dA > 0.0_rk .and. &
                all(ieee_is_finite(dTgc_dXg)),&
            expected = .true.,&
            msg      = "Default ansatz lost absolute surface measure.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0222


    subroutine forcad_nurbs_surface_0223(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk], [4,2])
        integer, parameter :: degree(2) = [1,1]
        integer, parameter :: ngauss(2) = [1,1]
        type(nurbs_surface) :: surface
        real(rk), allocatable :: Tgc(:), dTgc_dXg(:,:)
        real(rk) :: dA

        call surface%set(&
            knot1  = knot,&
            knot2  = knot,&
            Xc     = Xc,&
            degree = degree)
        call surface%ansatz(&
            ie       = 1,&
            ig       = 1,&
            Tgc      = Tgc,&
            dTgc_dXg = dTgc_dXg,&
            dA       = dA,&
            ngauss   = ngauss,&
            strict   = .true.)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "strict surface accepts orientation",&
            res      = surface%err%ok .and. dA > 0.0_rk .and. &
                all(ieee_is_finite(dTgc_dXg)),&
            expected = .true.,&
            msg      = "Strict ansatz rejected a positive surface map.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0223


    subroutine forcad_nurbs_surface_0224(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 128.0_rk*epsilon(1.0_rk)
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk], [4,2])
        real(rk), parameter :: X(3,2) = reshape([&
            0.0_rk,1.0_rk,0.5_rk,&
            0.0_rk,0.25_rk,1.0_rk], [3,2])
        type(nurbs_surface) :: surface
        real(rk), allocatable :: Xt1(:), Xt2(:), Xg(:,:)
        integer, allocatable :: elem(:,:)
        real(rk) :: nearest_Xg(2), nearest_Xt(2)
        integer :: id

        call surface%set(&
            knot1  = knot,&
            knot2  = knot,&
            Xc     = Xc,&
            degree = [1,1])
        call surface%create(res1=4, res2=3)
        call surface%put_to_nurbs(X)
        call surface%nearest_point(&
            point_Xg   = [1.0_rk,0.25_rk],&
            nearest_Xg = nearest_Xg,&
            nearest_Xt = nearest_Xt,&
            id         = id)
        Xt1 = surface%get_Xt(1)
        Xt2 = surface%get_Xt(2)
        Xg = surface%get_Xg()
        elem = surface%get_elem_Xg_vis()
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "put_to_nurbs surface cache replacement",&
            res      = surface%err%ok .and. all(surface%get_ng() == [3,1]) .and. &
                size(Xt1) == 0 .and. size(Xt2) == 0 .and. &
                all(shape(Xg) == [3,2]) .and. all(shape(elem) == [0,4]) .and. &
                id == 2 .and. maxval(abs(nearest_Xt-[1.0_rk,0.25_rk])) <= tol .and. &
                maxval(abs(nearest_Xg-[1.0_rk,0.25_rk])) <= tol,&
            expected = .true.,&
            msg      = "Surface parameter tuples and geometry cache are inconsistent.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0224


    subroutine forcad_nurbs_surface_0225(ut, ti)
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
        type(nurbs_surface) :: surface
        real(rk) :: Xt(18,2), Xdata(18,2), Xc(8,2)
        real(rk), allocatable :: Wc(:)
        integer :: i, j, point
        logical :: fit_failed

        do j = 1, 2
            do i = 1, 9
                point = (j-1)*9 + i
                Xt(point,:) = [u(i),real(j-1,rk)]
                Xdata(point,:) = [y(i),real(j-1,rk)]
            end do
            do i = 1, 4
                point = (j-1)*4 + i
                Xc(point,:) = [real(i-1,rk),real(j-1,rk)]
            end do
        end do

        call surface%set(&
            knot1  = knot1,&
            knot2  = knot2,&
            Xc     = Xc,&
            degree = [2,1])
        call surface%lsq_fit_nurbs(&
            Xt       = Xt,&
            Xdata    = Xdata,&
            ndata    = [9,2],&
            maxit    = 2,&
            tol      = 0.0_rk,&
            mu0      = 0.0_rk,&
            reg_logw = 0.0_rk)
        fit_failed = .not. surface%err%ok
        call surface%err%print()
        call surface%err%reset()
        Wc = surface%get_Wc()

        call ut%test(ti)%check(&
            name     = "surface NURBS fit maxit rollback",&
            res      = fit_failed .and. size(Wc) == 8 .and. &
                maxval(abs(Wc-1.0_rk)) <= tol,&
            expected = .true.,&
            msg      = "Surface fit retained a pending trial after maxit.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0225


    subroutine forcad_nurbs_surface_0226(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 1.0e-12_rk
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk], [4,2])
        type(nurbs_surface) :: surface
        real(rk) :: nearest_Xt(2), nearest_Xg(2)

        call surface%set(&
            knot1  = knot,&
            knot2  = knot,&
            Xc     = Xc,&
            degree = [1,1])
        call surface%nearest_point2(&
            point_Xg   = [-1.0_rk,0.37_rk],&
            tol        = tol,&
            maxit      = 5,&
            nearest_Xt = nearest_Xt,&
            nearest_Xg = nearest_Xg)
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "surface projected boundary descent",&
            res      = surface%err%ok .and. &
                maxval(abs(nearest_Xt-[0.0_rk,0.37_rk])) <= tol .and. &
                maxval(abs(nearest_Xg-[0.0_rk,0.37_rk])) <= tol,&
            expected = .true.,&
            msg      = "An active surface bound blocked free-direction descent.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0226


    subroutine forcad_nurbs_surface_0227(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 8192.0_rk*epsilon(1.0_rk)
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(4,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.2_rk,-0.1_rk,0.3_rk], [4,3])
        real(rk), parameter :: Wc(4) = [1.0_rk,2.0_rk,3.0_rk,4.0_rk]
        type(nurbs_surface) :: reference, small_scale, large_scale
        real(rk), allocatable :: Xg_ref(:), Xg_small(:), Xg_large(:)
        real(rk), allocatable :: T_ref(:), T_small(:), T_large(:)
        real(rk), allocatable :: dT_ref(:,:), dT_small(:,:), dT_large(:,:)
        real(rk), allocatable :: d2T_ref(:,:), d2T_small(:,:), d2T_large(:,:)

        call reference%set(knot, knot, Xc, Wc, degree=[1,1])
        call small_scale%set(knot, knot, Xc, 1.0e-200_rk*Wc, degree=[1,1])
        call large_scale%set(knot, knot, Xc, 1.0e200_rk*Wc, degree=[1,1])
        Xg_ref = reference%cmp_Xg([0.37_rk,0.61_rk])
        Xg_small = small_scale%cmp_Xg([0.37_rk,0.61_rk])
        Xg_large = large_scale%cmp_Xg([0.37_rk,0.61_rk])
        call reference%derivative2([0.37_rk,0.61_rk], d2T_ref, dT_ref, T_ref)
        call small_scale%derivative2([0.37_rk,0.61_rk], d2T_small, dT_small, T_small)
        call large_scale%derivative2([0.37_rk,0.61_rk], d2T_large, dT_large, T_large)
        call reference%err%print()
        call small_scale%err%print()
        call large_scale%err%print()

        call ut%test(ti)%check(&
            name     = "surface projective-scale invariance",&
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
            msg      = "Surface values changed under a common weight scale.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0227


    subroutine forcad_nurbs_surface_0228(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 8192.0_rk*epsilon(1.0_rk)
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xc(4,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.2_rk,-0.1_rk,0.3_rk], [4,3])
        real(rk), parameter :: Wc(4) = [1.0_rk,2.0_rk,3.0_rk,4.0_rk]
        integer, parameter :: ngauss(2) = [8,8]
        type(nurbs_surface) :: reference, small_scale, large_scale
        real(rk) :: area_ref, area_small, area_large

        call reference%set(knot, knot, Xc, Wc, degree=[1,1])
        call small_scale%set(knot, knot, Xc, 1.0e-200_rk*Wc, degree=[1,1])
        call large_scale%set(knot, knot, Xc, 1.0e200_rk*Wc, degree=[1,1])
        call reference%cmp_area(area_ref, ngauss=ngauss)
        call small_scale%cmp_area(area_small, ngauss=ngauss)
        call large_scale%cmp_area(area_large, ngauss=ngauss)
        call reference%err%print()
        call small_scale%err%print()
        call large_scale%err%print()

        call ut%test(ti)%check(&
            name     = "surface area weight-scale invariance",&
            res      = reference%err%ok .and. small_scale%err%ok .and. &
                large_scale%err%ok .and. ieee_is_finite(area_ref) .and. &
                ieee_is_finite(area_small) .and. ieee_is_finite(area_large) .and. &
                area_ref > 0.0_rk .and. abs(area_small-area_ref) <= tol .and. &
                abs(area_large-area_ref) <= tol,&
            expected = .true.,&
            msg      = "Surface area changed under a common weight scale.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0228


    subroutine forcad_nurbs_surface_0229(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 3.0e-11_rk
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xt(5,2) = reshape([&
            0.0_rk,0.21_rk,0.5_rk,0.78_rk,1.0_rk,&
            0.0_rk,0.67_rk,0.34_rk,0.91_rk,1.0_rk], [5,2])
        real(rk), parameter :: Xc(4,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.2_rk,-0.1_rk,0.3_rk], [4,3])
        real(rk), parameter :: Wc(4) = &
            1.0e-200_rk*[1.0_rk,2.0_rk,3.0_rk,4.0_rk]
        type(nurbs_surface) :: surface
        real(rk), allocatable :: Xg_ref(:,:), Xg_insert(:,:), Xg_remove(:,:), Xg_elevate(:,:)
        real(rk), allocatable :: Wc_after(:)

        call surface%set(knot, knot, Xc, Wc, degree=[1,1])
        call surface%create(Xt=Xt)
        Xg_ref = surface%get_Xg()
        call surface%insert_knots(1, [0.5_rk], [1])
        call surface%create(Xt=Xt)
        Xg_insert = surface%get_Xg()
        call surface%remove_knots(1, [0.5_rk], [1])
        call surface%create(Xt=Xt)
        Xg_remove = surface%get_Xg()
        call surface%elevate_degree(2, 1)
        call surface%create(Xt=Xt)
        Xg_elevate = surface%get_Xg()
        Wc_after = surface%get_Wc()
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "surface refinement preserves tiny weights",&
            res      = surface%err%ok .and. all(ieee_is_finite(Xg_elevate)) .and. &
                all(ieee_is_finite(Wc_after)) .and. all(Wc_after > 0.0_rk) .and. &
                maxval(Wc_after) < 1.0e-190_rk .and. &
                maxval(abs(Xg_insert-Xg_ref)) <= tol .and. &
                maxval(abs(Xg_remove-Xg_ref)) <= tol .and. &
                maxval(abs(Xg_elevate-Xg_ref)) <= tol,&
            expected = .true.,&
            msg      = "Surface refinement changed geometry or weight gauge.",&
            group    = "forcad_nurbs_surface")
        ti = ti + 1

    end subroutine forcad_nurbs_surface_0229


    subroutine run_nurbs_surface_tests(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call forcad_nurbs_surface_0001(ut, ti)
        call forcad_nurbs_surface_0002(ut, ti)
        call forcad_nurbs_surface_0003(ut, ti)
        call forcad_nurbs_surface_0004(ut, ti)
        call forcad_nurbs_surface_0005(ut, ti)
        call forcad_nurbs_surface_0006(ut, ti)
        call forcad_nurbs_surface_0007(ut, ti)
        call forcad_nurbs_surface_0008(ut, ti)
        call forcad_nurbs_surface_0009(ut, ti)
        call forcad_nurbs_surface_0010(ut, ti)
        call forcad_nurbs_surface_0011(ut, ti)
        call forcad_nurbs_surface_0012(ut, ti)
        call forcad_nurbs_surface_0013(ut, ti)
        call forcad_nurbs_surface_0014(ut, ti)
        call forcad_nurbs_surface_0015(ut, ti)
        call forcad_nurbs_surface_0016(ut, ti)
        call forcad_nurbs_surface_0017(ut, ti)
        call forcad_nurbs_surface_0018(ut, ti)
        call forcad_nurbs_surface_0019(ut, ti)
        call forcad_nurbs_surface_0020(ut, ti)
        call forcad_nurbs_surface_0021(ut, ti)
        call forcad_nurbs_surface_0022(ut, ti)
        call forcad_nurbs_surface_0023(ut, ti)
        call forcad_nurbs_surface_0024(ut, ti)
        call forcad_nurbs_surface_0025(ut, ti)
        call forcad_nurbs_surface_0026(ut, ti)
        call forcad_nurbs_surface_0027(ut, ti)
        call forcad_nurbs_surface_0028(ut, ti)
        call forcad_nurbs_surface_0029(ut, ti)
        call forcad_nurbs_surface_0030(ut, ti)
        call forcad_nurbs_surface_0031(ut, ti)
        call forcad_nurbs_surface_0032(ut, ti)
        call forcad_nurbs_surface_0033(ut, ti)
        call forcad_nurbs_surface_0034(ut, ti)
        call forcad_nurbs_surface_0035(ut, ti)
        call forcad_nurbs_surface_0036(ut, ti)
        call forcad_nurbs_surface_0037(ut, ti)
        call forcad_nurbs_surface_0038(ut, ti)
        call forcad_nurbs_surface_0039(ut, ti)
        call forcad_nurbs_surface_0040(ut, ti)
        call forcad_nurbs_surface_0041(ut, ti)
        call forcad_nurbs_surface_0042(ut, ti)
        call forcad_nurbs_surface_0043(ut, ti)
        call forcad_nurbs_surface_0044(ut, ti)
        call forcad_nurbs_surface_0045(ut, ti)
        call forcad_nurbs_surface_0046(ut, ti)
        call forcad_nurbs_surface_0047(ut, ti)
        call forcad_nurbs_surface_0048(ut, ti)
        call forcad_nurbs_surface_0049(ut, ti)
        call forcad_nurbs_surface_0050(ut, ti)
        call forcad_nurbs_surface_0051(ut, ti)
        call forcad_nurbs_surface_0052(ut, ti)
        call forcad_nurbs_surface_0053(ut, ti)
        call forcad_nurbs_surface_0054(ut, ti)
        call forcad_nurbs_surface_0055(ut, ti)
        call forcad_nurbs_surface_0056(ut, ti)
        call forcad_nurbs_surface_0057(ut, ti)
        call forcad_nurbs_surface_0058(ut, ti)
        call forcad_nurbs_surface_0059(ut, ti)
        call forcad_nurbs_surface_0060(ut, ti)
        call forcad_nurbs_surface_0061(ut, ti)
        call forcad_nurbs_surface_0062(ut, ti)
        call forcad_nurbs_surface_0063(ut, ti)
        call forcad_nurbs_surface_0064(ut, ti)
        call forcad_nurbs_surface_0065(ut, ti)
        call forcad_nurbs_surface_0066(ut, ti)
        call forcad_nurbs_surface_0067(ut, ti)
        call forcad_nurbs_surface_0068(ut, ti)
        call forcad_nurbs_surface_0069(ut, ti)
        call forcad_nurbs_surface_0070(ut, ti)
        call forcad_nurbs_surface_0071(ut, ti)
        call forcad_nurbs_surface_0072(ut, ti)
        call forcad_nurbs_surface_0073(ut, ti)
        call forcad_nurbs_surface_0074(ut, ti)
        call forcad_nurbs_surface_0075(ut, ti)
        call forcad_nurbs_surface_0076(ut, ti)
        call forcad_nurbs_surface_0077(ut, ti)
        call forcad_nurbs_surface_0078(ut, ti)
        call forcad_nurbs_surface_0079(ut, ti)
        call forcad_nurbs_surface_0080(ut, ti)
        call forcad_nurbs_surface_0081(ut, ti)
        call forcad_nurbs_surface_0082(ut, ti)
        call forcad_nurbs_surface_0083(ut, ti)
        call forcad_nurbs_surface_0084(ut, ti)
        call forcad_nurbs_surface_0085(ut, ti)
        call forcad_nurbs_surface_0086(ut, ti)
        call forcad_nurbs_surface_0087(ut, ti)
        call forcad_nurbs_surface_0088(ut, ti)
        call forcad_nurbs_surface_0089(ut, ti)
        call forcad_nurbs_surface_0090(ut, ti)
        call forcad_nurbs_surface_0091(ut, ti)
        call forcad_nurbs_surface_0092(ut, ti)
        call forcad_nurbs_surface_0093(ut, ti)
        call forcad_nurbs_surface_0094(ut, ti)
        call forcad_nurbs_surface_0095(ut, ti)
        call forcad_nurbs_surface_0096(ut, ti)
        call forcad_nurbs_surface_0097(ut, ti)
        call forcad_nurbs_surface_0098(ut, ti)
        call forcad_nurbs_surface_0099(ut, ti)
        call forcad_nurbs_surface_0100(ut, ti)
        call forcad_nurbs_surface_0101(ut, ti)
        call forcad_nurbs_surface_0102(ut, ti)
        call forcad_nurbs_surface_0103(ut, ti)
        call forcad_nurbs_surface_0104(ut, ti)
        call forcad_nurbs_surface_0105(ut, ti)
        call forcad_nurbs_surface_0106(ut, ti)
        call forcad_nurbs_surface_0107(ut, ti)
        call forcad_nurbs_surface_0108(ut, ti)
        call forcad_nurbs_surface_0109(ut, ti)
        call forcad_nurbs_surface_0110(ut, ti)
        call forcad_nurbs_surface_0111(ut, ti)
        call forcad_nurbs_surface_0112(ut, ti)
        call forcad_nurbs_surface_0113(ut, ti)
        call forcad_nurbs_surface_0114(ut, ti)
        call forcad_nurbs_surface_0115(ut, ti)
        call forcad_nurbs_surface_0116(ut, ti)
        call forcad_nurbs_surface_0117(ut, ti)
        call forcad_nurbs_surface_0118(ut, ti)
        call forcad_nurbs_surface_0119(ut, ti)
        call forcad_nurbs_surface_0120(ut, ti)
        call forcad_nurbs_surface_0121(ut, ti)
        call forcad_nurbs_surface_0122(ut, ti)
        call forcad_nurbs_surface_0123(ut, ti)
        call forcad_nurbs_surface_0124(ut, ti)
        call forcad_nurbs_surface_0125(ut, ti)
        call forcad_nurbs_surface_0126(ut, ti)
        call forcad_nurbs_surface_0127(ut, ti)
        call forcad_nurbs_surface_0128(ut, ti)
        call forcad_nurbs_surface_0129(ut, ti)
        call forcad_nurbs_surface_0130(ut, ti)
        call forcad_nurbs_surface_0131(ut, ti)
        call forcad_nurbs_surface_0132(ut, ti)
        call forcad_nurbs_surface_0133(ut, ti)
        call forcad_nurbs_surface_0134(ut, ti)
        call forcad_nurbs_surface_0135(ut, ti)
        call forcad_nurbs_surface_0136(ut, ti)
        call forcad_nurbs_surface_0137(ut, ti)
        call forcad_nurbs_surface_0138(ut, ti)
        call forcad_nurbs_surface_0139(ut, ti)
        call forcad_nurbs_surface_0140(ut, ti)
        call forcad_nurbs_surface_0141(ut, ti)
        call forcad_nurbs_surface_0142(ut, ti)
        call forcad_nurbs_surface_0143(ut, ti)
        call forcad_nurbs_surface_0144(ut, ti)
        call forcad_nurbs_surface_0145(ut, ti)
        call forcad_nurbs_surface_0146(ut, ti)
        call forcad_nurbs_surface_0147(ut, ti)
        call forcad_nurbs_surface_0148(ut, ti)
        call forcad_nurbs_surface_0149(ut, ti)
        call forcad_nurbs_surface_0150(ut, ti)
        call forcad_nurbs_surface_0151(ut, ti)
        call forcad_nurbs_surface_0152(ut, ti)
        call forcad_nurbs_surface_0153(ut, ti)
        call forcad_nurbs_surface_0154(ut, ti)
        call forcad_nurbs_surface_0155(ut, ti)
        call forcad_nurbs_surface_0156(ut, ti)
        call forcad_nurbs_surface_0157(ut, ti)
        call forcad_nurbs_surface_0158(ut, ti)
        call forcad_nurbs_surface_0159(ut, ti)
        call forcad_nurbs_surface_0160(ut, ti)
        call forcad_nurbs_surface_0161(ut, ti)
        call forcad_nurbs_surface_0162(ut, ti)
        call forcad_nurbs_surface_0163(ut, ti)
        call forcad_nurbs_surface_0164(ut, ti)
        call forcad_nurbs_surface_0165(ut, ti)
        call forcad_nurbs_surface_0166(ut, ti)
        call forcad_nurbs_surface_0167(ut, ti)
        call forcad_nurbs_surface_0168(ut, ti)
        call forcad_nurbs_surface_0169(ut, ti)
        call forcad_nurbs_surface_0170(ut, ti)
        call forcad_nurbs_surface_0171(ut, ti)
        call forcad_nurbs_surface_0172(ut, ti)
        call forcad_nurbs_surface_0173(ut, ti)
        call forcad_nurbs_surface_0174(ut, ti)
        call forcad_nurbs_surface_0175(ut, ti)
        call forcad_nurbs_surface_0176(ut, ti)
        call forcad_nurbs_surface_0177(ut, ti)
        call forcad_nurbs_surface_0178(ut, ti)
        call forcad_nurbs_surface_0179(ut, ti)
        call forcad_nurbs_surface_0180(ut, ti)
        call forcad_nurbs_surface_0181(ut, ti)
        call forcad_nurbs_surface_0182(ut, ti)
        call forcad_nurbs_surface_0183(ut, ti)
        call forcad_nurbs_surface_0184(ut, ti)
        call forcad_nurbs_surface_0185(ut, ti)
        call forcad_nurbs_surface_0186(ut, ti)
        call forcad_nurbs_surface_0187(ut, ti)
        call forcad_nurbs_surface_0188(ut, ti)
        call forcad_nurbs_surface_0189(ut, ti)
        call forcad_nurbs_surface_0190(ut, ti)
        call forcad_nurbs_surface_0191(ut, ti)
        call forcad_nurbs_surface_0192(ut, ti)
        call forcad_nurbs_surface_0193(ut, ti)
        call forcad_nurbs_surface_0194(ut, ti)
        call forcad_nurbs_surface_0195(ut, ti)
        call forcad_nurbs_surface_0196(ut, ti)
        call forcad_nurbs_surface_0197(ut, ti)
        call forcad_nurbs_surface_0198(ut, ti)
        call forcad_nurbs_surface_0199(ut, ti)
        call forcad_nurbs_surface_0200(ut, ti)
        call forcad_nurbs_surface_0201(ut, ti)
        call forcad_nurbs_surface_0202(ut, ti)
        call forcad_nurbs_surface_0203(ut, ti)
        call forcad_nurbs_surface_0204(ut, ti)
        call forcad_nurbs_surface_0205(ut, ti)
        call forcad_nurbs_surface_0206(ut, ti)
        call forcad_nurbs_surface_0207(ut, ti)
        call forcad_nurbs_surface_0208(ut, ti)
        call forcad_nurbs_surface_0209(ut, ti)
        call forcad_nurbs_surface_0210(ut, ti)
        call forcad_nurbs_surface_0211(ut, ti)
        call forcad_nurbs_surface_0212(ut, ti)
        call forcad_nurbs_surface_0213(ut, ti)
        call forcad_nurbs_surface_0214(ut, ti)
        call forcad_nurbs_surface_0215(ut, ti)
        call forcad_nurbs_surface_0216(ut, ti)
        call forcad_nurbs_surface_0217(ut, ti)
        call forcad_nurbs_surface_0218(ut, ti)
        call forcad_nurbs_surface_0219(ut, ti)
        call forcad_nurbs_surface_0220(ut, ti)
        call forcad_nurbs_surface_0221(ut, ti)
        call forcad_nurbs_surface_0222(ut, ti)
        call forcad_nurbs_surface_0223(ut, ti)
        call forcad_nurbs_surface_0224(ut, ti)
        call forcad_nurbs_surface_0225(ut, ti)
        call forcad_nurbs_surface_0226(ut, ti)
        call forcad_nurbs_surface_0227(ut, ti)
        call forcad_nurbs_surface_0228(ut, ti)
        call forcad_nurbs_surface_0229(ut, ti)

    end subroutine run_nurbs_surface_tests

end module forcad_test_nurbs_surface
