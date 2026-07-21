module forcad_test_geometry

    use forcad_geometry, only: extrude, revolve, sweep, loft
    use forcad_kinds, only: rk
    use forcad_nurbs_curve, only: nurbs_curve
    use forcad_nurbs_surface, only: nurbs_surface
    use forcad_nurbs_volume, only: nurbs_volume
    use forunittest, only: unit_tests

    implicit none

    private
    public :: run_geometry_tests

contains

    subroutine forcad_geometry_0001(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        call curve%err%print()
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve extrusion diagnostic",&
            res      = surface%err%ok,&
            expected = .true.,&
            msg      = "Valid curve extrusion returned an error.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0001


    subroutine forcad_geometry_0002(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        call curve%err%print()
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve extrusion degree",&
            res      = surface%get_degree(),&
            expected = [2,1],&
            msg      = "Curve extrusion did not preserve the source degree.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0002


    subroutine forcad_geometry_0003(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call curve%err%print()
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve extrusion control shape",&
            res      = surface%get_nc(),&
            expected = [4,2],&
            msg      = "Curve extrusion produced the wrong control-net shape.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0003


    subroutine forcad_geometry_0004(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_tol = 1.0e3_rk*epsilon(1.0_rk)
        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call curve%err%print()
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve extrusion knot",&
            res      = surface%get_knot(1),&
            expected = geometry_knot,&
            tol      = geometry_tol,&
            msg      = "Curve extrusion changed the source knot vector.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0004


    subroutine forcad_geometry_0005(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_tol = 1.0e3_rk*epsilon(1.0_rk)
        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call curve%err%print()
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve extrusion weights",&
            res      = surface%get_Wc(),&
            expected = [geometry_Wc,geometry_Wc],&
            tol      = geometry_tol,&
            msg      = "Curve extrusion did not preserve rational weights.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0005


    subroutine forcad_geometry_0006(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_tol = 1.0e3_rk*epsilon(1.0_rk)
        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface
        real(rk), allocatable :: curve_point(:)
        real(rk), allocatable :: surface_point(:)

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])

        curve_point = curve%cmp_Xg(0.37_rk)
        surface_point = surface%cmp_Xg([0.37_rk,0.40_rk])
        call curve%err%print()
        call surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve extrusion geometry",&
            res      = surface_point,&
            expected = [curve_point,0.8_rk],&
            tol      = geometry_tol,&
            msg      = "Extruded surface does not equal curve plus the linear extrusion.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0006


    subroutine forcad_geometry_0007(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface
        type(nurbs_volume) :: volume
        real(rk), allocatable :: curve_point(:)
        real(rk), allocatable :: surface_point(:)

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])

        curve_point = curve%cmp_Xg(0.37_rk)
        surface_point = surface%cmp_Xg([0.37_rk,0.40_rk])

        volume = extrude(surface, [0.5_rk,-0.25_rk,1.0_rk])
        call curve%err%print()
        call surface%err%print()
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface extrusion degree",&
            res      = volume%get_degree(),&
            expected = [2,1,1],&
            msg      = "Surface extrusion did not preserve source degrees.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0007


    subroutine forcad_geometry_0008(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface
        type(nurbs_volume) :: volume
        real(rk), allocatable :: curve_point(:)
        real(rk), allocatable :: surface_point(:)

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])

        curve_point = curve%cmp_Xg(0.37_rk)
        surface_point = surface%cmp_Xg([0.37_rk,0.40_rk])

        volume = extrude(surface, [0.5_rk,-0.25_rk,1.0_rk])
        call curve%err%print()
        call surface%err%print()
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface extrusion rationality",&
            res      = volume%is_rational(),&
            expected = .true.,&
            msg      = "Surface extrusion lost rational weights.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0008


    subroutine forcad_geometry_0009(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_tol = 1.0e3_rk*epsilon(1.0_rk)
        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface
        type(nurbs_volume) :: volume
        real(rk), allocatable :: curve_point(:), surface_point(:), volume_point(:)

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])

        curve_point = curve%cmp_Xg(0.37_rk)
        surface_point = surface%cmp_Xg([0.37_rk,0.40_rk])

        volume = extrude(surface, [0.5_rk,-0.25_rk,1.0_rk])
        surface_point = surface%cmp_Xg([0.63_rk,0.25_rk])
        volume_point = volume%cmp_Xg([0.63_rk,0.25_rk,0.60_rk])
        call curve%err%print()
        call surface%err%print()
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface extrusion geometry",&
            res      = volume_point,&
            expected = surface_point + 0.60_rk*[0.5_rk,-0.25_rk,1.0_rk],&
            tol      = geometry_tol,&
            msg      = "Surface extrusion geometry is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0009


    subroutine forcad_geometry_0010(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface, invalid_surface
        type(nurbs_volume) :: volume
        real(rk), allocatable :: curve_point(:), surface_point(:), volume_point(:)

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])

        curve_point = curve%cmp_Xg(0.37_rk)
        surface_point = surface%cmp_Xg([0.37_rk,0.40_rk])

        volume = extrude(surface, [0.5_rk,-0.25_rk,1.0_rk])
        surface_point = surface%cmp_Xg([0.63_rk,0.25_rk])
        volume_point = volume%cmp_Xg([0.63_rk,0.25_rk,0.60_rk])

        invalid_surface = extrude(curve, [0.0_rk,0.0_rk,0.0_rk])
        call curve%err%print()
        call surface%err%print()
        call invalid_surface%err%print()
        call volume%err%print()

        call ut%test(ti)%check(&
            name     = "zero extrusion diagnostic",&
            res      = invalid_surface%err%ok,&
            expected = .false.,&
            msg      = "A zero extrusion vector must return a recoverable diagnostic.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0010


    subroutine forcad_geometry_0011(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve
        type(nurbs_surface) :: revolved_surface

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        revolved_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = 2.0_rk*geometry_pi)
        call axis_curve%err%print()
        call revolved_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve revolution diagnostic",&
            res      = revolved_surface%err%ok,&
            expected = .true.,&
            msg      = "Valid curve revolution returned an error.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0011


    subroutine forcad_geometry_0012(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve
        type(nurbs_surface) :: revolved_surface

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        revolved_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = 2.0_rk*geometry_pi)
        call axis_curve%err%print()
        call revolved_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve revolution degree",&
            res      = revolved_surface%get_degree(),&
            expected = [1,2],&
            msg      = "Curve revolution did not preserve the profile degree.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0012


    subroutine forcad_geometry_0013(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve
        type(nurbs_surface) :: revolved_surface

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        revolved_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = 2.0_rk*geometry_pi)
        call axis_curve%err%print()
        call revolved_surface%err%print()

        call ut%test(ti)%check(&
            name     = "full revolution control shape",&
            res      = revolved_surface%get_nc(),&
            expected = [2,9],&
            msg      = "Full revolution must use four exact quadratic arc segments.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0013


    subroutine forcad_geometry_0014(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve
        type(nurbs_surface) :: revolved_surface

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        revolved_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = 2.0_rk*geometry_pi)
        call axis_curve%err%print()
        call revolved_surface%err%print()

        call ut%test(ti)%check(&
            name     = "full revolution rationality",&
            res      = revolved_surface%is_rational(),&
            expected = .true.,&
            msg      = "Exact circular revolution must retain rational arc weights.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0014


    subroutine forcad_geometry_0015(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_tol = 1.0e3_rk*epsilon(1.0_rk)
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve
        type(nurbs_surface) :: revolved_surface
        real(rk), allocatable :: surface_point(:)

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        revolved_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = 2.0_rk*geometry_pi)
        surface_point = revolved_surface%cmp_Xg([0.25_rk,0.25_rk])
        call axis_curve%err%print()
        call revolved_surface%err%print()

        call ut%test(ti)%check(&
            name     = "quarter revolution geometry",&
            res      = surface_point,&
            expected = [0.0_rk,2.0_rk,-0.5_rk],&
            tol      = geometry_tol,&
            msg      = "A quarter revolution did not produce the exact circular point.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0015


    subroutine forcad_geometry_0016(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_tol = 1.0e3_rk*epsilon(1.0_rk)
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve
        type(nurbs_surface) :: revolved_surface
        real(rk), allocatable :: surface_point(:)

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        revolved_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = 2.0_rk*geometry_pi)
        surface_point = revolved_surface%cmp_Xg([0.25_rk,0.25_rk])
        surface_point = revolved_surface%cmp_Xg([0.25_rk,0.125_rk])
        call axis_curve%err%print()
        call revolved_surface%err%print()

        call ut%test(ti)%check(&
            name     = "mid-arc revolution geometry",&
            res      = surface_point,&
            expected = [sqrt(2.0_rk),sqrt(2.0_rk),-0.5_rk],&
            tol      = geometry_tol,&
            msg      = "The rational quadratic arc is not exact at 45 degrees.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0016


    subroutine forcad_geometry_0017(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_tol = 1.0e3_rk*epsilon(1.0_rk)
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve
        type(nurbs_surface) :: radial_surface
        type(nurbs_surface) :: revolved_surface
        type(nurbs_volume) :: revolved_volume
        real(rk), allocatable :: surface_point(:), volume_point(:)

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        revolved_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = 2.0_rk*geometry_pi)
        surface_point = revolved_surface%cmp_Xg([0.25_rk,0.25_rk])
        surface_point = revolved_surface%cmp_Xg([0.25_rk,0.125_rk])

        radial_surface = extrude(axis_curve, [1.0_rk,0.0_rk,0.0_rk])
        revolved_volume = revolve(&
            surface        = radial_surface,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = geometry_pi)
        volume_point = revolved_volume%cmp_Xg([0.25_rk,0.40_rk,0.50_rk])
        call axis_curve%err%print()
        call radial_surface%err%print()
        call revolved_surface%err%print()
        call revolved_volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface revolution geometry",&
            res      = volume_point,&
            expected = [0.0_rk,2.4_rk,-0.5_rk],&
            tol      = geometry_tol,&
            msg      = "Surface revolution geometry is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0017


    subroutine forcad_geometry_0018(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve
        type(nurbs_surface) :: radial_surface, revolved_surface, invalid_surface
        type(nurbs_volume) :: revolved_volume
        real(rk), allocatable :: surface_point(:), volume_point(:)

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        revolved_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = 2.0_rk*geometry_pi)
        surface_point = revolved_surface%cmp_Xg([0.25_rk,0.25_rk])
        surface_point = revolved_surface%cmp_Xg([0.25_rk,0.125_rk])

        radial_surface = extrude(axis_curve, [1.0_rk,0.0_rk,0.0_rk])
        revolved_volume = revolve(&
            surface        = radial_surface,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = geometry_pi)
        volume_point = revolved_volume%cmp_Xg([0.25_rk,0.40_rk,0.50_rk])

        invalid_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,0.0_rk],&
            angle          = geometry_pi)
        call axis_curve%err%print()
        call radial_surface%err%print()
        call revolved_surface%err%print()
        call invalid_surface%err%print()
        call revolved_volume%err%print()

        call ut%test(ti)%check(&
            name     = "zero revolution axis diagnostic",&
            res      = invalid_surface%err%ok,&
            expected = .false.,&
            msg      = "A zero revolution axis must return a recoverable diagnostic.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0018


    subroutine forcad_geometry_0019(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_tol = 1.0e3_rk*epsilon(1.0_rk)
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve
        type(nurbs_surface) :: radial_surface, revolved_surface, invalid_surface
        type(nurbs_volume) :: revolved_volume
        real(rk), allocatable :: surface_point(:), volume_point(:)

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        revolved_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = 2.0_rk*geometry_pi)
        surface_point = revolved_surface%cmp_Xg([0.25_rk,0.25_rk])
        surface_point = revolved_surface%cmp_Xg([0.25_rk,0.125_rk])

        radial_surface = extrude(axis_curve, [1.0_rk,0.0_rk,0.0_rk])
        revolved_volume = revolve(&
            surface        = radial_surface,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = geometry_pi)
        volume_point = revolved_volume%cmp_Xg([0.25_rk,0.40_rk,0.50_rk])

        invalid_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,0.0_rk],&
            angle          = geometry_pi)

        revolved_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,huge(1.0_rk)],&
            angle          = 0.5_rk*geometry_pi)
        surface_point = revolved_surface%cmp_Xg([0.25_rk,1.0_rk])
        call axis_curve%err%print()
        call radial_surface%err%print()
        call revolved_surface%err%print()
        call invalid_surface%err%print()
        call revolved_volume%err%print()

        call ut%test(ti)%check(&
            name     = "scaled revolution axis",&
            res      = surface_point,&
            expected = [0.0_rk,2.0_rk,-0.5_rk],&
            tol      = geometry_tol,&
            msg      = "Scaled revolution axis is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0019


    subroutine forcad_geometry_0020(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_spine_Xc(4,3) = reshape([&
            0.0_rk,0.4_rk,1.0_rk,1.5_rk,0.0_rk,-0.5_rk,0.5_rk,0.0_rk,&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk], [4,3])
        real(rk), parameter :: geometry_spine_Wc(4) = [1.0_rk,1.2_rk,0.8_rk,1.0_rk]
        type(nurbs_curve) :: curve, spine
        type(nurbs_surface) :: surface
        type(nurbs_surface) :: swept_surface

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call spine%set(&
            knot   = geometry_knot,&
            Xc     = geometry_spine_Xc,&
            Wc     = geometry_spine_Wc,&
            degree = 2)
        swept_surface = sweep(&
            profile = curve,&
            spine   = spine,&
            origin  = [0.0_rk,0.0_rk,0.0_rk])
        call curve%err%print()
        call spine%err%print()
        call surface%err%print()
        call swept_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve sweep diagnostic",&
            res      = swept_surface%err%ok,&
            expected = .true.,&
            msg      = "Valid curve sweep returned an error.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0020


    subroutine forcad_geometry_0021(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_spine_Xc(4,3) = reshape([&
            0.0_rk,0.4_rk,1.0_rk,1.5_rk,0.0_rk,-0.5_rk,0.5_rk,0.0_rk,&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk], [4,3])
        real(rk), parameter :: geometry_spine_Wc(4) = [1.0_rk,1.2_rk,0.8_rk,1.0_rk]
        type(nurbs_curve) :: curve, spine
        type(nurbs_surface) :: surface
        type(nurbs_surface) :: swept_surface

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call spine%set(&
            knot   = geometry_knot,&
            Xc     = geometry_spine_Xc,&
            Wc     = geometry_spine_Wc,&
            degree = 2)
        swept_surface = sweep(&
            profile = curve,&
            spine   = spine,&
            origin  = [0.0_rk,0.0_rk,0.0_rk])
        call curve%err%print()
        call spine%err%print()
        call surface%err%print()
        call swept_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve sweep degree",&
            res      = swept_surface%get_degree(),&
            expected = [2,2],&
            msg      = "Curve sweep did not preserve profile and spine degrees.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0021


    subroutine forcad_geometry_0022(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_spine_Xc(4,3) = reshape([&
            0.0_rk,0.4_rk,1.0_rk,1.5_rk,0.0_rk,-0.5_rk,0.5_rk,0.0_rk,&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk], [4,3])
        real(rk), parameter :: geometry_spine_Wc(4) = [1.0_rk,1.2_rk,0.8_rk,1.0_rk]
        type(nurbs_curve) :: curve, spine
        type(nurbs_surface) :: surface
        type(nurbs_surface) :: swept_surface

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call spine%set(&
            knot   = geometry_knot,&
            Xc     = geometry_spine_Xc,&
            Wc     = geometry_spine_Wc,&
            degree = 2)
        swept_surface = sweep(&
            profile = curve,&
            spine   = spine,&
            origin  = [0.0_rk,0.0_rk,0.0_rk])
        call curve%err%print()
        call spine%err%print()
        call surface%err%print()
        call swept_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve sweep control shape",&
            res      = swept_surface%get_nc(),&
            expected = [4,4],&
            msg      = "Curve sweep produced the wrong tensor-product control shape.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0022


    subroutine forcad_geometry_0023(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_spine_Xc(4,3) = reshape([&
            0.0_rk,0.4_rk,1.0_rk,1.5_rk,0.0_rk,-0.5_rk,0.5_rk,0.0_rk,&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk], [4,3])
        real(rk), parameter :: geometry_spine_Wc(4) = [1.0_rk,1.2_rk,0.8_rk,1.0_rk]
        type(nurbs_curve) :: curve, spine
        type(nurbs_surface) :: surface
        type(nurbs_surface) :: swept_surface

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call spine%set(&
            knot   = geometry_knot,&
            Xc     = geometry_spine_Xc,&
            Wc     = geometry_spine_Wc,&
            degree = 2)
        swept_surface = sweep(&
            profile = curve,&
            spine   = spine,&
            origin  = [0.0_rk,0.0_rk,0.0_rk])
        call curve%err%print()
        call spine%err%print()
        call surface%err%print()
        call swept_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve sweep rationality",&
            res      = swept_surface%is_rational(),&
            expected = .true.,&
            msg      = "Curve sweep lost the rational product weights.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0023


    subroutine forcad_geometry_0024(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_tol = 1.0e3_rk*epsilon(1.0_rk)
        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_spine_Xc(4,3) = reshape([&
            0.0_rk,0.4_rk,1.0_rk,1.5_rk,0.0_rk,-0.5_rk,0.5_rk,0.0_rk,&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk], [4,3])
        real(rk), parameter :: geometry_spine_Wc(4) = [1.0_rk,1.2_rk,0.8_rk,1.0_rk]
        type(nurbs_curve) :: curve, spine
        type(nurbs_surface) :: surface
        type(nurbs_surface) :: swept_surface
        real(rk), allocatable :: curve_point(:)
        real(rk), allocatable :: spine_point(:)
        real(rk), allocatable :: surface_point(:)

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call spine%set(&
            knot   = geometry_knot,&
            Xc     = geometry_spine_Xc,&
            Wc     = geometry_spine_Wc,&
            degree = 2)
        swept_surface = sweep(&
            profile = curve,&
            spine   = spine,&
            origin  = [0.0_rk,0.0_rk,0.0_rk])
        curve_point = curve%cmp_Xg(0.31_rk)
        spine_point = spine%cmp_Xg(0.67_rk)
        surface_point = swept_surface%cmp_Xg([0.31_rk,0.67_rk])
        call curve%err%print()
        call spine%err%print()
        call surface%err%print()
        call swept_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve sweep exact geometry",&
            res      = surface_point,&
            expected = [curve_point(1)+spine_point(1),curve_point(2)+spine_point(2),spine_point(3)],&
            tol      = geometry_tol,&
            msg      = "Curve sweep exact geometry is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0024


    subroutine forcad_geometry_0025(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_spine_Xc(4,3) = reshape([&
            0.0_rk,0.4_rk,1.0_rk,1.5_rk,0.0_rk,-0.5_rk,0.5_rk,0.0_rk,&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk], [4,3])
        real(rk), parameter :: geometry_spine_Wc(4) = [1.0_rk,1.2_rk,0.8_rk,1.0_rk]
        type(nurbs_curve) :: curve, spine
        type(nurbs_surface) :: surface
        type(nurbs_surface) :: swept_surface
        type(nurbs_volume) :: swept_volume
        real(rk), allocatable :: curve_point(:)
        real(rk), allocatable :: spine_point(:)
        real(rk), allocatable :: surface_point(:)

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call spine%set(&
            knot   = geometry_knot,&
            Xc     = geometry_spine_Xc,&
            Wc     = geometry_spine_Wc,&
            degree = 2)
        swept_surface = sweep(&
            profile = curve,&
            spine   = spine,&
            origin  = [0.0_rk,0.0_rk,0.0_rk])
        curve_point = curve%cmp_Xg(0.31_rk)
        spine_point = spine%cmp_Xg(0.67_rk)
        surface_point = swept_surface%cmp_Xg([0.31_rk,0.67_rk])

        swept_volume = sweep(&
            profile = surface,&
            spine   = spine,&
            origin  = [0.0_rk,0.0_rk,0.0_rk])
        call curve%err%print()
        call spine%err%print()
        call surface%err%print()
        call swept_surface%err%print()
        call swept_volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface sweep degree",&
            res      = swept_volume%get_degree(),&
            expected = [2,1,2],&
            msg      = "Surface sweep did not preserve profile and spine degrees.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0025


    subroutine forcad_geometry_0026(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_tol = 1.0e3_rk*epsilon(1.0_rk)
        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_spine_Xc(4,3) = reshape([&
            0.0_rk,0.4_rk,1.0_rk,1.5_rk,0.0_rk,-0.5_rk,0.5_rk,0.0_rk,&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk], [4,3])
        real(rk), parameter :: geometry_spine_Wc(4) = [1.0_rk,1.2_rk,0.8_rk,1.0_rk]
        type(nurbs_curve) :: curve, spine
        type(nurbs_surface) :: surface
        type(nurbs_surface) :: swept_surface
        type(nurbs_volume) :: swept_volume
        real(rk), allocatable :: curve_point(:), spine_point(:), surface_point(:), volume_point(:)

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call spine%set(&
            knot   = geometry_knot,&
            Xc     = geometry_spine_Xc,&
            Wc     = geometry_spine_Wc,&
            degree = 2)
        swept_surface = sweep(&
            profile = curve,&
            spine   = spine,&
            origin  = [0.0_rk,0.0_rk,0.0_rk])
        curve_point = curve%cmp_Xg(0.31_rk)
        spine_point = spine%cmp_Xg(0.67_rk)
        surface_point = swept_surface%cmp_Xg([0.31_rk,0.67_rk])

        swept_volume = sweep(&
            profile = surface,&
            spine   = spine,&
            origin  = [0.0_rk,0.0_rk,0.0_rk])
        surface_point = surface%cmp_Xg([0.42_rk,0.36_rk])
        spine_point = spine%cmp_Xg(0.58_rk)
        volume_point = swept_volume%cmp_Xg([0.42_rk,0.36_rk,0.58_rk])
        call curve%err%print()
        call spine%err%print()
        call surface%err%print()
        call swept_surface%err%print()
        call swept_volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface sweep exact geometry",&
            res      = volume_point,&
            expected = surface_point + spine_point,&
            tol      = geometry_tol,&
            msg      = "Surface sweep exact geometry is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0026


    subroutine forcad_geometry_0027(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_spine_Xc(4,3) = reshape([&
            0.0_rk,0.4_rk,1.0_rk,1.5_rk,0.0_rk,-0.5_rk,0.5_rk,0.0_rk,&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk], [4,3])
        real(rk), parameter :: geometry_spine_Wc(4) = [1.0_rk,1.2_rk,0.8_rk,1.0_rk]
        type(nurbs_curve) :: curve, spine
        type(nurbs_surface) :: surface, swept_surface, invalid_surface
        type(nurbs_volume) :: swept_volume
        real(rk), allocatable :: curve_point(:), spine_point(:), surface_point(:), volume_point(:)

        call curve%set(&
            knot   = geometry_knot,&
            Xc     = geometry_Xc,&
            Wc     = geometry_Wc,&
            degree = 2)
        surface = extrude(curve, [0.0_rk,0.0_rk,2.0_rk])
        call spine%set(&
            knot   = geometry_knot,&
            Xc     = geometry_spine_Xc,&
            Wc     = geometry_spine_Wc,&
            degree = 2)
        swept_surface = sweep(&
            profile = curve,&
            spine   = spine,&
            origin  = [0.0_rk,0.0_rk,0.0_rk])
        curve_point = curve%cmp_Xg(0.31_rk)
        spine_point = spine%cmp_Xg(0.67_rk)
        surface_point = swept_surface%cmp_Xg([0.31_rk,0.67_rk])

        swept_volume = sweep(&
            profile = surface,&
            spine   = spine,&
            origin  = [0.0_rk,0.0_rk,0.0_rk])
        surface_point = surface%cmp_Xg([0.42_rk,0.36_rk])
        spine_point = spine%cmp_Xg(0.58_rk)
        volume_point = swept_volume%cmp_Xg([0.42_rk,0.36_rk,0.58_rk])

        invalid_surface = sweep(&
            profile = curve,&
            spine   = spine,&
            origin  = [0.0_rk,0.0_rk])
        call curve%err%print()
        call spine%err%print()
        call surface%err%print()
        call swept_surface%err%print()
        call invalid_surface%err%print()
        call swept_volume%err%print()

        call ut%test(ti)%check(&
            name     = "sweep origin diagnostic",&
            res      = invalid_surface%err%ok,&
            expected = .false.,&
            msg      = "A mismatched sweep origin must return a recoverable diagnostic.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0027


    subroutine forcad_geometry_0028(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        real(rk), parameter :: geometry_loft_parameters(4) = [0.0_rk,0.18_rk,0.64_rk,1.0_rk]
        type(nurbs_curve) :: axis_curve, curve_sections(4)
        type(nurbs_surface) :: lofted_surface
        real(rk) :: section_Xc(4,3)
        real(rk) :: xi
        integer :: section

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        do section = 1, size(curve_sections)
            xi = real(section-1,rk)/real(size(curve_sections)-1,rk)
            section_Xc(:,1) = (1.0_rk+0.12_rk*xi)*geometry_Xc(:,1) + 0.20_rk*xi
            section_Xc(:,2) = geometry_Xc(:,2) + 0.45_rk*sin(geometry_pi*xi)*[0.0_rk,1.0_rk,1.0_rk,0.0_rk]
            section_Xc(:,3) = 3.0_rk*xi + 0.15_rk*xi*geometry_Xc(:,1)
            call curve_sections(section)%set(&
                knot   = geometry_knot,&
                Xc     = section_Xc,&
                Wc     = geometry_Wc,&
                degree = 2)
        end do
        lofted_surface = loft(&
            sections   = curve_sections,&
            parameters = geometry_loft_parameters,&
            degree     = 3)
        call axis_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call lofted_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve loft diagnostic",&
            res      = lofted_surface%err%ok,&
            expected = .true.,&
            msg      = "Valid compatible curve loft returned an error.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0028


    subroutine forcad_geometry_0029(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        real(rk), parameter :: geometry_loft_parameters(4) = [0.0_rk,0.18_rk,0.64_rk,1.0_rk]
        type(nurbs_curve) :: axis_curve, curve_sections(4)
        type(nurbs_surface) :: lofted_surface
        real(rk) :: section_Xc(4,3)
        real(rk) :: xi
        integer :: section

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        do section = 1, size(curve_sections)
            xi = real(section-1,rk)/real(size(curve_sections)-1,rk)
            section_Xc(:,1) = (1.0_rk+0.12_rk*xi)*geometry_Xc(:,1) + 0.20_rk*xi
            section_Xc(:,2) = geometry_Xc(:,2) + 0.45_rk*sin(geometry_pi*xi)*[0.0_rk,1.0_rk,1.0_rk,0.0_rk]
            section_Xc(:,3) = 3.0_rk*xi + 0.15_rk*xi*geometry_Xc(:,1)
            call curve_sections(section)%set(&
                knot   = geometry_knot,&
                Xc     = section_Xc,&
                Wc     = geometry_Wc,&
                degree = 2)
        end do
        lofted_surface = loft(&
            sections   = curve_sections,&
            parameters = geometry_loft_parameters,&
            degree     = 3)
        call axis_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call lofted_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve loft degree",&
            res      = lofted_surface%get_degree(),&
            expected = [2,3],&
            msg      = "Curve loft degree is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0029


    subroutine forcad_geometry_0030(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        real(rk), parameter :: geometry_loft_parameters(4) = [0.0_rk,0.18_rk,0.64_rk,1.0_rk]
        type(nurbs_curve) :: axis_curve, curve_sections(4)
        type(nurbs_surface) :: lofted_surface
        real(rk) :: section_Xc(4,3)
        real(rk) :: xi
        integer :: section

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        do section = 1, size(curve_sections)
            xi = real(section-1,rk)/real(size(curve_sections)-1,rk)
            section_Xc(:,1) = (1.0_rk+0.12_rk*xi)*geometry_Xc(:,1) + 0.20_rk*xi
            section_Xc(:,2) = geometry_Xc(:,2) + 0.45_rk*sin(geometry_pi*xi)*[0.0_rk,1.0_rk,1.0_rk,0.0_rk]
            section_Xc(:,3) = 3.0_rk*xi + 0.15_rk*xi*geometry_Xc(:,1)
            call curve_sections(section)%set(&
                knot   = geometry_knot,&
                Xc     = section_Xc,&
                Wc     = geometry_Wc,&
                degree = 2)
        end do
        lofted_surface = loft(&
            sections   = curve_sections,&
            parameters = geometry_loft_parameters,&
            degree     = 3)
        call axis_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call lofted_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve loft control shape",&
            res      = lofted_surface%get_nc(),&
            expected = [4,4],&
            msg      = "Curve loft produced the wrong tensor-product control shape.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0030


    subroutine forcad_geometry_0031(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        real(rk), parameter :: geometry_loft_parameters(4) = [0.0_rk,0.18_rk,0.64_rk,1.0_rk]
        type(nurbs_curve) :: axis_curve, curve_sections(4)
        type(nurbs_surface) :: lofted_surface
        real(rk) :: section_Xc(4,3)
        real(rk) :: xi
        integer :: section

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        do section = 1, size(curve_sections)
            xi = real(section-1,rk)/real(size(curve_sections)-1,rk)
            section_Xc(:,1) = (1.0_rk+0.12_rk*xi)*geometry_Xc(:,1) + 0.20_rk*xi
            section_Xc(:,2) = geometry_Xc(:,2) + 0.45_rk*sin(geometry_pi*xi)*[0.0_rk,1.0_rk,1.0_rk,0.0_rk]
            section_Xc(:,3) = 3.0_rk*xi + 0.15_rk*xi*geometry_Xc(:,1)
            call curve_sections(section)%set(&
                knot   = geometry_knot,&
                Xc     = section_Xc,&
                Wc     = geometry_Wc,&
                degree = 2)
        end do
        lofted_surface = loft(&
            sections   = curve_sections,&
            parameters = geometry_loft_parameters,&
            degree     = 3)
        call axis_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call lofted_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve loft rationality",&
            res      = lofted_surface%is_rational(),&
            expected = .true.,&
            msg      = "Curve loft lost rational section weights.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0031


    subroutine forcad_geometry_0032(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_tol = 1.0e3_rk*epsilon(1.0_rk)
        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        real(rk), parameter :: geometry_loft_parameters(4) = [0.0_rk,0.18_rk,0.64_rk,1.0_rk]
        type(nurbs_curve) :: axis_curve, curve_sections(4)
        type(nurbs_surface) :: lofted_surface
        real(rk), allocatable :: curve_point(:)
        real(rk), allocatable :: surface_point(:)
        real(rk) :: section_Xc(4,3), loft_error, xi
        integer :: section, sample

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        do section = 1, size(curve_sections)
            xi = real(section-1,rk)/real(size(curve_sections)-1,rk)
            section_Xc(:,1) = (1.0_rk+0.12_rk*xi)*geometry_Xc(:,1) + 0.20_rk*xi
            section_Xc(:,2) = geometry_Xc(:,2) + 0.45_rk*sin(geometry_pi*xi)*[0.0_rk,1.0_rk,1.0_rk,0.0_rk]
            section_Xc(:,3) = 3.0_rk*xi + 0.15_rk*xi*geometry_Xc(:,1)
            call curve_sections(section)%set(&
                knot   = geometry_knot,&
                Xc     = section_Xc,&
                Wc     = geometry_Wc,&
                degree = 2)
        end do
        lofted_surface = loft(&
            sections   = curve_sections,&
            parameters = geometry_loft_parameters,&
            degree     = 3)
        loft_error = 0.0_rk
        do section = 1, size(curve_sections)
            do sample = 1, 3
                xi = real(2*sample-1,rk)/6.0_rk
                curve_point = curve_sections(section)%cmp_Xg(xi)
                surface_point = lofted_surface%cmp_Xg([xi,geometry_loft_parameters(section)])
                loft_error = max(loft_error,maxval(abs(surface_point-curve_point)))
            end do
        end do
        call axis_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call lofted_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve loft section interpolation",&
            res      = loft_error,&
            expected = 0.0_rk,&
            tol      = geometry_tol,&
            msg      = "Curve loft does not interpolate every supplied rational section.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0032


    subroutine forcad_geometry_0033(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        real(rk), parameter :: geometry_loft_parameters(4) = [0.0_rk,0.18_rk,0.64_rk,1.0_rk]
        real(rk), parameter :: geometry_surface_loft_parameters(3) = [0.0_rk,0.40_rk,1.0_rk]
        type(nurbs_curve) :: axis_curve, curve_sections(4)
        type(nurbs_surface) :: lofted_surface
        type(nurbs_surface) :: surface_sections(3)
        type(nurbs_volume) :: lofted_volume
        real(rk), allocatable :: curve_point(:)
        real(rk), allocatable :: surface_point(:)
        real(rk) :: section_Xc(4,3), loft_error, xi
        integer :: section, sample

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        do section = 1, size(curve_sections)
            xi = real(section-1,rk)/real(size(curve_sections)-1,rk)
            section_Xc(:,1) = (1.0_rk+0.12_rk*xi)*geometry_Xc(:,1) + 0.20_rk*xi
            section_Xc(:,2) = geometry_Xc(:,2) + 0.45_rk*sin(geometry_pi*xi)*[0.0_rk,1.0_rk,1.0_rk,0.0_rk]
            section_Xc(:,3) = 3.0_rk*xi + 0.15_rk*xi*geometry_Xc(:,1)
            call curve_sections(section)%set(&
                knot   = geometry_knot,&
                Xc     = section_Xc,&
                Wc     = geometry_Wc,&
                degree = 2)
        end do
        lofted_surface = loft(&
            sections   = curve_sections,&
            parameters = geometry_loft_parameters,&
            degree     = 3)
        loft_error = 0.0_rk
        do section = 1, size(curve_sections)
            do sample = 1, 3
                xi = real(2*sample-1,rk)/6.0_rk
                curve_point = curve_sections(section)%cmp_Xg(xi)
                surface_point = lofted_surface%cmp_Xg([xi,geometry_loft_parameters(section)])
                loft_error = max(loft_error,maxval(abs(surface_point-curve_point)))
            end do
        end do

        do section = 1, size(surface_sections)
            surface_sections(section) = extrude(&
                curve  = curve_sections(section),&
                vector = [0.10_rk*real(section-2,rk),0.15_rk,0.75_rk+0.10_rk*real(section,rk)])
        end do
        lofted_volume = loft(&
            sections   = surface_sections,&
            parameters = geometry_surface_loft_parameters,&
            degree     = 2)
        call axis_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call lofted_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call lofted_volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface loft diagnostic",&
            res      = lofted_volume%err%ok,&
            expected = .true.,&
            msg      = "Valid compatible surface loft returned an error.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0033


    subroutine forcad_geometry_0034(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        real(rk), parameter :: geometry_loft_parameters(4) = [0.0_rk,0.18_rk,0.64_rk,1.0_rk]
        real(rk), parameter :: geometry_surface_loft_parameters(3) = [0.0_rk,0.40_rk,1.0_rk]
        type(nurbs_curve) :: axis_curve, curve_sections(4)
        type(nurbs_surface) :: lofted_surface
        type(nurbs_surface) :: surface_sections(3)
        type(nurbs_volume) :: lofted_volume
        real(rk), allocatable :: curve_point(:)
        real(rk), allocatable :: surface_point(:)
        real(rk) :: section_Xc(4,3), loft_error, xi
        integer :: section, sample

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        do section = 1, size(curve_sections)
            xi = real(section-1,rk)/real(size(curve_sections)-1,rk)
            section_Xc(:,1) = (1.0_rk+0.12_rk*xi)*geometry_Xc(:,1) + 0.20_rk*xi
            section_Xc(:,2) = geometry_Xc(:,2) + 0.45_rk*sin(geometry_pi*xi)*[0.0_rk,1.0_rk,1.0_rk,0.0_rk]
            section_Xc(:,3) = 3.0_rk*xi + 0.15_rk*xi*geometry_Xc(:,1)
            call curve_sections(section)%set(&
                knot   = geometry_knot,&
                Xc     = section_Xc,&
                Wc     = geometry_Wc,&
                degree = 2)
        end do
        lofted_surface = loft(&
            sections   = curve_sections,&
            parameters = geometry_loft_parameters,&
            degree     = 3)
        loft_error = 0.0_rk
        do section = 1, size(curve_sections)
            do sample = 1, 3
                xi = real(2*sample-1,rk)/6.0_rk
                curve_point = curve_sections(section)%cmp_Xg(xi)
                surface_point = lofted_surface%cmp_Xg([xi,geometry_loft_parameters(section)])
                loft_error = max(loft_error,maxval(abs(surface_point-curve_point)))
            end do
        end do

        do section = 1, size(surface_sections)
            surface_sections(section) = extrude(&
                curve  = curve_sections(section),&
                vector = [0.10_rk*real(section-2,rk),0.15_rk,0.75_rk+0.10_rk*real(section,rk)])
        end do
        lofted_volume = loft(&
            sections   = surface_sections,&
            parameters = geometry_surface_loft_parameters,&
            degree     = 2)
        call axis_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call lofted_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call lofted_volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface loft degree",&
            res      = lofted_volume%get_degree(),&
            expected = [2,1,2],&
            msg      = "Surface loft degree is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0034


    subroutine forcad_geometry_0035(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_tol = 1.0e3_rk*epsilon(1.0_rk)
        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        real(rk), parameter :: geometry_loft_parameters(4) = [0.0_rk,0.18_rk,0.64_rk,1.0_rk]
        real(rk), parameter :: geometry_surface_loft_parameters(3) = [0.0_rk,0.40_rk,1.0_rk]
        type(nurbs_curve) :: axis_curve, curve_sections(4)
        type(nurbs_surface) :: lofted_surface
        type(nurbs_surface) :: surface_sections(3)
        type(nurbs_volume) :: lofted_volume
        real(rk), allocatable :: curve_point(:), surface_point(:), volume_point(:)
        real(rk) :: section_Xc(4,3), loft_error, xi
        integer :: section, sample

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        do section = 1, size(curve_sections)
            xi = real(section-1,rk)/real(size(curve_sections)-1,rk)
            section_Xc(:,1) = (1.0_rk+0.12_rk*xi)*geometry_Xc(:,1) + 0.20_rk*xi
            section_Xc(:,2) = geometry_Xc(:,2) + 0.45_rk*sin(geometry_pi*xi)*[0.0_rk,1.0_rk,1.0_rk,0.0_rk]
            section_Xc(:,3) = 3.0_rk*xi + 0.15_rk*xi*geometry_Xc(:,1)
            call curve_sections(section)%set(&
                knot   = geometry_knot,&
                Xc     = section_Xc,&
                Wc     = geometry_Wc,&
                degree = 2)
        end do
        lofted_surface = loft(&
            sections   = curve_sections,&
            parameters = geometry_loft_parameters,&
            degree     = 3)
        loft_error = 0.0_rk
        do section = 1, size(curve_sections)
            do sample = 1, 3
                xi = real(2*sample-1,rk)/6.0_rk
                curve_point = curve_sections(section)%cmp_Xg(xi)
                surface_point = lofted_surface%cmp_Xg([xi,geometry_loft_parameters(section)])
                loft_error = max(loft_error,maxval(abs(surface_point-curve_point)))
            end do
        end do

        do section = 1, size(surface_sections)
            surface_sections(section) = extrude(&
                curve  = curve_sections(section),&
                vector = [0.10_rk*real(section-2,rk),0.15_rk,0.75_rk+0.10_rk*real(section,rk)])
        end do
        lofted_volume = loft(&
            sections   = surface_sections,&
            parameters = geometry_surface_loft_parameters,&
            degree     = 2)
        loft_error = 0.0_rk
        do section = 1, size(surface_sections)
            do sample = 1, 3
                xi = real(2*sample-1,rk)/6.0_rk
                surface_point = surface_sections(section)%cmp_Xg([xi,1.0_rk-xi])
                volume_point = lofted_volume%cmp_Xg(&
                    [xi,1.0_rk-xi,geometry_surface_loft_parameters(section)])
                loft_error = max(loft_error,maxval(abs(volume_point-surface_point)))
            end do
        end do
        call axis_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call lofted_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call lofted_volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface loft section interpolation",&
            res      = loft_error,&
            expected = 0.0_rk,&
            tol      = geometry_tol,&
            msg      = "Surface loft section interpolation is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0035


    subroutine forcad_geometry_0036(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        real(rk), parameter :: geometry_loft_parameters(4) = [0.0_rk,0.18_rk,0.64_rk,1.0_rk]
        real(rk), parameter :: geometry_surface_loft_parameters(3) = [0.0_rk,0.40_rk,1.0_rk]
        type(nurbs_curve) :: axis_curve, curve_sections(4)
        type(nurbs_surface) :: invalid_surface, lofted_surface, surface_sections(3)
        type(nurbs_volume) :: lofted_volume
        real(rk), allocatable :: curve_point(:), surface_point(:), volume_point(:)
        real(rk) :: section_Xc(4,3), loft_error, xi
        integer :: section, sample

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        do section = 1, size(curve_sections)
            xi = real(section-1,rk)/real(size(curve_sections)-1,rk)
            section_Xc(:,1) = (1.0_rk+0.12_rk*xi)*geometry_Xc(:,1) + 0.20_rk*xi
            section_Xc(:,2) = geometry_Xc(:,2) + 0.45_rk*sin(geometry_pi*xi)*[0.0_rk,1.0_rk,1.0_rk,0.0_rk]
            section_Xc(:,3) = 3.0_rk*xi + 0.15_rk*xi*geometry_Xc(:,1)
            call curve_sections(section)%set(&
                knot   = geometry_knot,&
                Xc     = section_Xc,&
                Wc     = geometry_Wc,&
                degree = 2)
        end do
        lofted_surface = loft(&
            sections   = curve_sections,&
            parameters = geometry_loft_parameters,&
            degree     = 3)
        loft_error = 0.0_rk
        do section = 1, size(curve_sections)
            do sample = 1, 3
                xi = real(2*sample-1,rk)/6.0_rk
                curve_point = curve_sections(section)%cmp_Xg(xi)
                surface_point = lofted_surface%cmp_Xg([xi,geometry_loft_parameters(section)])
                loft_error = max(loft_error,maxval(abs(surface_point-curve_point)))
            end do
        end do

        do section = 1, size(surface_sections)
            surface_sections(section) = extrude(&
                curve  = curve_sections(section),&
                vector = [0.10_rk*real(section-2,rk),0.15_rk,0.75_rk+0.10_rk*real(section,rk)])
        end do
        lofted_volume = loft(&
            sections   = surface_sections,&
            parameters = geometry_surface_loft_parameters,&
            degree     = 2)
        loft_error = 0.0_rk
        do section = 1, size(surface_sections)
            do sample = 1, 3
                xi = real(2*sample-1,rk)/6.0_rk
                surface_point = surface_sections(section)%cmp_Xg([xi,1.0_rk-xi])
                volume_point = lofted_volume%cmp_Xg(&
                    [xi,1.0_rk-xi,geometry_surface_loft_parameters(section)])
                loft_error = max(loft_error,maxval(abs(volume_point-surface_point)))
            end do
        end do

        invalid_surface = loft(&
            sections = [curve_sections(1),axis_curve],&
            degree   = 1)
        call axis_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call invalid_surface%err%print()
        call lofted_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call lofted_volume%err%print()

        call ut%test(ti)%check(&
            name     = "incompatible loft section diagnostic",&
            res      = invalid_surface%err%ok,&
            expected = .false.,&
            msg      = "Incompatible loft section diagnostic is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0036


    subroutine forcad_geometry_0037(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_knot(7) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_Xc(4,2) = reshape([&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,0.0_rk,1.0_rk,-0.5_rk,0.0_rk], [4,2])
        real(rk), parameter :: geometry_Wc(4) = [1.0_rk,0.7_rk,1.3_rk,1.0_rk]
        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        real(rk), parameter :: geometry_loft_parameters(4) = [0.0_rk,0.18_rk,0.64_rk,1.0_rk]
        real(rk), parameter :: geometry_surface_loft_parameters(3) = [0.0_rk,0.40_rk,1.0_rk]
        type(nurbs_curve) :: axis_curve, curve_sections(4)
        type(nurbs_surface) :: invalid_surface, lofted_surface, surface_sections(3)
        type(nurbs_volume) :: lofted_volume
        real(rk), allocatable :: curve_point(:), surface_point(:), volume_point(:)
        real(rk) :: section_Xc(4,3), loft_error, xi
        integer :: section, sample

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        do section = 1, size(curve_sections)
            xi = real(section-1,rk)/real(size(curve_sections)-1,rk)
            section_Xc(:,1) = (1.0_rk+0.12_rk*xi)*geometry_Xc(:,1) + 0.20_rk*xi
            section_Xc(:,2) = geometry_Xc(:,2) + 0.45_rk*sin(geometry_pi*xi)*[0.0_rk,1.0_rk,1.0_rk,0.0_rk]
            section_Xc(:,3) = 3.0_rk*xi + 0.15_rk*xi*geometry_Xc(:,1)
            call curve_sections(section)%set(&
                knot   = geometry_knot,&
                Xc     = section_Xc,&
                Wc     = geometry_Wc,&
                degree = 2)
        end do
        lofted_surface = loft(&
            sections   = curve_sections,&
            parameters = geometry_loft_parameters,&
            degree     = 3)
        loft_error = 0.0_rk
        do section = 1, size(curve_sections)
            do sample = 1, 3
                xi = real(2*sample-1,rk)/6.0_rk
                curve_point = curve_sections(section)%cmp_Xg(xi)
                surface_point = lofted_surface%cmp_Xg([xi,geometry_loft_parameters(section)])
                loft_error = max(loft_error,maxval(abs(surface_point-curve_point)))
            end do
        end do

        do section = 1, size(surface_sections)
            surface_sections(section) = extrude(&
                curve  = curve_sections(section),&
                vector = [0.10_rk*real(section-2,rk),0.15_rk,0.75_rk+0.10_rk*real(section,rk)])
        end do
        lofted_volume = loft(&
            sections   = surface_sections,&
            parameters = geometry_surface_loft_parameters,&
            degree     = 2)
        loft_error = 0.0_rk
        do section = 1, size(surface_sections)
            do sample = 1, 3
                xi = real(2*sample-1,rk)/6.0_rk
                surface_point = surface_sections(section)%cmp_Xg([xi,1.0_rk-xi])
                volume_point = lofted_volume%cmp_Xg(&
                    [xi,1.0_rk-xi,geometry_surface_loft_parameters(section)])
                loft_error = max(loft_error,maxval(abs(volume_point-surface_point)))
            end do
        end do

        invalid_surface = loft(&
            sections = [curve_sections(1),axis_curve],&
            degree   = 1)
        invalid_surface = loft(&
            sections   = curve_sections,&
            parameters = [0.0_rk,0.5_rk,0.4_rk,1.0_rk],&
            degree     = 2)
        call axis_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call invalid_surface%err%print()
        call lofted_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call lofted_volume%err%print()

        call ut%test(ti)%check(&
            name     = "loft parameter diagnostic",&
            res      = invalid_surface%err%ok,&
            expected = .false.,&
            msg      = "Loft parameter diagnostic is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0037


    subroutine forcad_geometry_0038(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_curve) :: periodic_curve
        type(nurbs_surface) :: periodic_surface
        real(rk) :: periodic_knot(9), periodic_Xc(6,3), periodic_Wc(6)

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        periodic_Xc(1,:) = [ 1.0_rk, 0.0_rk, 0.0_rk]
        periodic_Xc(2,:) = [ 0.0_rk, 1.0_rk, 0.0_rk]
        periodic_Xc(3,:) = [-1.0_rk, 0.0_rk, 0.0_rk]
        periodic_Xc(4,:) = [ 0.0_rk,-1.0_rk, 0.0_rk]
        periodic_Xc(5,:) = periodic_Xc(1,:)
        periodic_Xc(6,:) = periodic_Xc(2,:)
        periodic_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        call periodic_curve%set(&
            knot           = periodic_knot,&
            Xc             = periodic_Xc,&
            Wc             = periodic_Wc,&
            degree         = 2,&
            wrap_parameters = .true.)

        periodic_surface = extrude(periodic_curve, [0.0_rk,0.0_rk,1.0_rk])
        call periodic_curve%err%print()
        call periodic_surface%err%print()

        call ut%test(ti)%check(&
            name     = "curve extrusion periodic directions",&
            res      = [periodic_surface%get_parameter_wrapping(1),periodic_surface%get_parameter_wrapping(2)],&
            expected = [.true.,.false.],&
            msg      = "Curve extrusion periodic directions are incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0038


    subroutine forcad_geometry_0039(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_curve) :: periodic_curve
        type(nurbs_surface) :: periodic_surface
        type(nurbs_volume) :: periodic_volume
        real(rk) :: periodic_knot(9), periodic_Xc(6,3), periodic_Wc(6)

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        periodic_Xc(1,:) = [ 1.0_rk, 0.0_rk, 0.0_rk]
        periodic_Xc(2,:) = [ 0.0_rk, 1.0_rk, 0.0_rk]
        periodic_Xc(3,:) = [-1.0_rk, 0.0_rk, 0.0_rk]
        periodic_Xc(4,:) = [ 0.0_rk,-1.0_rk, 0.0_rk]
        periodic_Xc(5,:) = periodic_Xc(1,:)
        periodic_Xc(6,:) = periodic_Xc(2,:)
        periodic_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        call periodic_curve%set(&
            knot           = periodic_knot,&
            Xc             = periodic_Xc,&
            Wc             = periodic_Wc,&
            degree         = 2,&
            wrap_parameters = .true.)

        periodic_surface = extrude(periodic_curve, [0.0_rk,0.0_rk,1.0_rk])

        periodic_volume = extrude(periodic_surface, [0.0_rk,0.0_rk,1.0_rk])
        call periodic_curve%err%print()
        call periodic_surface%err%print()
        call periodic_volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface extrusion periodic directions",&
            res      = [periodic_volume%get_parameter_wrapping(1),periodic_volume%get_parameter_wrapping(2),&
                periodic_volume%get_parameter_wrapping(3)],&
            expected = [.true.,.false.,.false.],&
            msg      = "Surface extrusion must preserve source periodic directions only.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0039


    subroutine forcad_geometry_0040(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        type(nurbs_curve) :: periodic_curve
        type(nurbs_surface) :: periodic_surface
        type(nurbs_volume) :: periodic_volume
        real(rk) :: periodic_knot(9), periodic_Xc(6,3), periodic_Wc(6)

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        periodic_Xc(1,:) = [ 1.0_rk, 0.0_rk, 0.0_rk]
        periodic_Xc(2,:) = [ 0.0_rk, 1.0_rk, 0.0_rk]
        periodic_Xc(3,:) = [-1.0_rk, 0.0_rk, 0.0_rk]
        periodic_Xc(4,:) = [ 0.0_rk,-1.0_rk, 0.0_rk]
        periodic_Xc(5,:) = periodic_Xc(1,:)
        periodic_Xc(6,:) = periodic_Xc(2,:)
        periodic_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        call periodic_curve%set(&
            knot           = periodic_knot,&
            Xc             = periodic_Xc,&
            Wc             = periodic_Wc,&
            degree         = 2,&
            wrap_parameters = .true.)

        periodic_surface = extrude(periodic_curve, [0.0_rk,0.0_rk,1.0_rk])

        periodic_volume = extrude(periodic_surface, [0.0_rk,0.0_rk,1.0_rk])

        periodic_surface = revolve(&
            curve          = periodic_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = 0.5_rk*geometry_pi)
        call periodic_curve%err%print()
        call periodic_surface%err%print()
        call periodic_volume%err%print()

        call ut%test(ti)%check(&
            name     = "curve revolution periodic directions",&
            res      = [periodic_surface%get_parameter_wrapping(1),periodic_surface%get_parameter_wrapping(2)],&
            expected = [.true.,.false.],&
            msg      = "Curve revolution periodic directions are incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0040


    subroutine forcad_geometry_0041(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        type(nurbs_curve) :: periodic_curve
        type(nurbs_surface) :: periodic_surface
        type(nurbs_volume) :: periodic_volume
        real(rk) :: periodic_knot(9), periodic_Xc(6,3), periodic_Wc(6)

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        periodic_Xc(1,:) = [ 1.0_rk, 0.0_rk, 0.0_rk]
        periodic_Xc(2,:) = [ 0.0_rk, 1.0_rk, 0.0_rk]
        periodic_Xc(3,:) = [-1.0_rk, 0.0_rk, 0.0_rk]
        periodic_Xc(4,:) = [ 0.0_rk,-1.0_rk, 0.0_rk]
        periodic_Xc(5,:) = periodic_Xc(1,:)
        periodic_Xc(6,:) = periodic_Xc(2,:)
        periodic_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        call periodic_curve%set(&
            knot           = periodic_knot,&
            Xc             = periodic_Xc,&
            Wc             = periodic_Wc,&
            degree         = 2,&
            wrap_parameters = .true.)

        periodic_surface = extrude(periodic_curve, [0.0_rk,0.0_rk,1.0_rk])

        periodic_volume = extrude(periodic_surface, [0.0_rk,0.0_rk,1.0_rk])

        periodic_surface = revolve(&
            curve          = periodic_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        periodic_surface = sweep(periodic_curve, periodic_curve)
        call periodic_curve%err%print()
        call periodic_surface%err%print()
        call periodic_volume%err%print()

        call ut%test(ti)%check(&
            name     = "curve sweep periodic directions",&
            res      = [periodic_surface%get_parameter_wrapping(1),periodic_surface%get_parameter_wrapping(2)],&
            expected = [.true.,.true.],&
            msg      = "A sweep must preserve periodic profile and spine directions.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0041


    subroutine forcad_geometry_0042(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        type(nurbs_curve) :: periodic_curve
        type(nurbs_surface) :: periodic_surface
        type(nurbs_volume) :: periodic_volume
        real(rk) :: periodic_knot(9), periodic_Xc(6,3), periodic_Wc(6)

        periodic_knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk]
        periodic_Xc(1,:) = [ 1.0_rk, 0.0_rk, 0.0_rk]
        periodic_Xc(2,:) = [ 0.0_rk, 1.0_rk, 0.0_rk]
        periodic_Xc(3,:) = [-1.0_rk, 0.0_rk, 0.0_rk]
        periodic_Xc(4,:) = [ 0.0_rk,-1.0_rk, 0.0_rk]
        periodic_Xc(5,:) = periodic_Xc(1,:)
        periodic_Xc(6,:) = periodic_Xc(2,:)
        periodic_Wc = [1.0_rk, 1.4_rk, 0.8_rk, 1.2_rk, 1.0_rk, 1.4_rk]
        call periodic_curve%set(&
            knot           = periodic_knot,&
            Xc             = periodic_Xc,&
            Wc             = periodic_Wc,&
            degree         = 2,&
            wrap_parameters = .true.)

        periodic_surface = extrude(periodic_curve, [0.0_rk,0.0_rk,1.0_rk])

        periodic_volume = extrude(periodic_surface, [0.0_rk,0.0_rk,1.0_rk])

        periodic_surface = revolve(&
            curve          = periodic_curve,&
            axis_point     = [0.0_rk,0.0_rk,0.0_rk],&
            axis_direction = [0.0_rk,0.0_rk,1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        periodic_surface = sweep(periodic_curve, periodic_curve)

        periodic_surface = loft([periodic_curve,periodic_curve], degree=1)
        call periodic_curve%err%print()
        call periodic_surface%err%print()
        call periodic_volume%err%print()

        call ut%test(ti)%check(&
            name     = "curve loft periodic directions",&
            res      = [periodic_surface%get_parameter_wrapping(1),periodic_surface%get_parameter_wrapping(2)],&
            expected = [.true.,.false.],&
            msg      = "Curve loft periodic directions are incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0042


    subroutine forcad_geometry_0043(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()

        call ut%test(ti)%check(&
            name     = "polynomial curve extrusion",&
            res      = result_surface%err%ok .and. .not. result_surface%is_rational(),&
            expected = .true.,&
            msg      = "Polynomial curve extrusion is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0043


    subroutine forcad_geometry_0044(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "polynomial surface extrusion",&
            res      = result_volume%err%ok .and. .not. result_volume%is_rational(),&
            expected = .true.,&
            msg      = "Polynomial surface extrusion is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0044


    subroutine forcad_geometry_0045(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "invalid curve extrusion source",&
            res      = result_surface%err%ok,&
            expected = .false.,&
            msg      = "Extrusion must reject a curve carrying an active diagnostic.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0045


    subroutine forcad_geometry_0046(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "invalid surface extrusion source",&
            res      = result_volume%err%ok,&
            expected = .false.,&
            msg      = "Extrusion must reject a surface carrying an active diagnostic.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0046


    subroutine forcad_geometry_0047(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "polynomial curve revolution",&
            res      = result_surface%err%ok .and. result_surface%is_rational(),&
            expected = .true.,&
            msg      = "Polynomial curve revolution is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0047


    subroutine forcad_geometry_0048(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "polynomial surface revolution",&
            res      = result_volume%err%ok .and. result_volume%is_rational(),&
            expected = .true.,&
            msg      = "Polynomial surface revolution is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0048


    subroutine forcad_geometry_0049(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "invalid curve revolution source",&
            res      = result_surface%err%ok,&
            expected = .false.,&
            msg      = "Revolution must reject a curve carrying an active diagnostic.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0049


    subroutine forcad_geometry_0050(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "invalid surface revolution source",&
            res      = result_volume%err%ok,&
            expected = .false.,&
            msg      = "Revolution must reject a surface carrying an active diagnostic.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0050


    subroutine forcad_geometry_0051(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "polynomial curve sweep",&
            res      = result_surface%err%ok .and. .not. result_surface%is_rational(),&
            expected = .true.,&
            msg      = "Polynomial curve sweep is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0051


    subroutine forcad_geometry_0052(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)

        result_volume = sweep(polynomial_surface, axis_curve)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "polynomial surface sweep",&
            res      = result_volume%err%ok .and. .not. result_volume%is_rational(),&
            expected = .true.,&
            msg      = "Polynomial surface sweep is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0052


    subroutine forcad_geometry_0053(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)

        result_volume = sweep(polynomial_surface, axis_curve)

        result_surface = sweep(invalid_curve, axis_curve)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "invalid curve sweep source",&
            res      = result_surface%err%ok,&
            expected = .false.,&
            msg      = "Sweep must reject a curve profile carrying an active diagnostic.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0053


    subroutine forcad_geometry_0054(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)

        result_volume = sweep(polynomial_surface, axis_curve)

        result_surface = sweep(invalid_curve, axis_curve)

        result_volume = sweep(invalid_surface, axis_curve)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "invalid surface sweep source",&
            res      = result_volume%err%ok,&
            expected = .false.,&
            msg      = "Invalid surface sweep source is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0054


    subroutine forcad_geometry_0055(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)

        result_volume = sweep(polynomial_surface, axis_curve)

        result_surface = sweep(invalid_curve, axis_curve)

        result_volume = sweep(invalid_surface, axis_curve)

        result_volume = sweep(&
            profile = polynomial_surface,&
            spine   = axis_curve,&
            origin  = [0.0_rk, 0.0_rk])
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface sweep origin diagnostic",&
            res      = result_volume%err%ok,&
            expected = .false.,&
            msg      = "Surface sweep origin diagnostic is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0055


    subroutine forcad_geometry_0056(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)

        result_volume = sweep(polynomial_surface, axis_curve)

        result_surface = sweep(invalid_curve, axis_curve)

        result_volume = sweep(invalid_surface, axis_curve)

        result_volume = sweep(&
            profile = polynomial_surface,&
            spine   = axis_curve,&
            origin  = [0.0_rk, 0.0_rk])

        result_surface = loft([axis_curve], degree=1)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "curve loft section count",&
            res      = result_surface%err%ok,&
            expected = .false.,&
            msg      = "A curve loft must require at least two sections.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0056


    subroutine forcad_geometry_0057(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)

        result_volume = sweep(polynomial_surface, axis_curve)

        result_surface = sweep(invalid_curve, axis_curve)

        result_volume = sweep(invalid_surface, axis_curve)

        result_volume = sweep(&
            profile = polynomial_surface,&
            spine   = axis_curve,&
            origin  = [0.0_rk, 0.0_rk])

        result_surface = loft([axis_curve], degree=1)

        result_surface = loft(&
            sections   = curve_sections,&
            parameters = [0.0_rk, 1.0_rk],&
            degree     = 1)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "curve loft parameter count",&
            res      = result_surface%err%ok,&
            expected = .false.,&
            msg      = "A curve loft must require one parameter per section.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0057


    subroutine forcad_geometry_0058(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)

        result_volume = sweep(polynomial_surface, axis_curve)

        result_surface = sweep(invalid_curve, axis_curve)

        result_volume = sweep(invalid_surface, axis_curve)

        result_volume = sweep(&
            profile = polynomial_surface,&
            spine   = axis_curve,&
            origin  = [0.0_rk, 0.0_rk])

        result_surface = loft([axis_curve], degree=1)

        result_surface = loft(&
            sections   = curve_sections,&
            parameters = [0.0_rk, 1.0_rk],&
            degree     = 1)

        result_surface = loft(curve_sections, degree=size(curve_sections))
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "curve loft degree range",&
            res      = result_surface%err%ok,&
            expected = .false.,&
            msg      = "Curve loft degree range is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0058


    subroutine forcad_geometry_0059(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)

        result_volume = sweep(polynomial_surface, axis_curve)

        result_surface = sweep(invalid_curve, axis_curve)

        result_volume = sweep(invalid_surface, axis_curve)

        result_volume = sweep(&
            profile = polynomial_surface,&
            spine   = axis_curve,&
            origin  = [0.0_rk, 0.0_rk])

        result_surface = loft([axis_curve], degree=1)

        result_surface = loft(&
            sections   = curve_sections,&
            parameters = [0.0_rk, 1.0_rk],&
            degree     = 1)

        result_surface = loft(curve_sections, degree=size(curve_sections))

        result_surface = loft([axis_curve, axis_curve], degree=1)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "polynomial curve loft",&
            res      = result_surface%err%ok .and. .not. result_surface%is_rational(),&
            expected = .true.,&
            msg      = "Polynomial curve loft is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0059


    subroutine forcad_geometry_0060(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)

        result_volume = sweep(polynomial_surface, axis_curve)

        result_surface = sweep(invalid_curve, axis_curve)

        result_volume = sweep(invalid_surface, axis_curve)

        result_volume = sweep(&
            profile = polynomial_surface,&
            spine   = axis_curve,&
            origin  = [0.0_rk, 0.0_rk])

        result_surface = loft([axis_curve], degree=1)

        result_surface = loft(&
            sections   = curve_sections,&
            parameters = [0.0_rk, 1.0_rk],&
            degree     = 1)

        result_surface = loft(curve_sections, degree=size(curve_sections))

        result_surface = loft([axis_curve, axis_curve], degree=1)

        result_volume = loft([polynomial_surface], degree=1)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface loft section count",&
            res      = result_volume%err%ok,&
            expected = .false.,&
            msg      = "A surface loft must require at least two sections.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0060


    subroutine forcad_geometry_0061(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)

        result_volume = sweep(polynomial_surface, axis_curve)

        result_surface = sweep(invalid_curve, axis_curve)

        result_volume = sweep(invalid_surface, axis_curve)

        result_volume = sweep(&
            profile = polynomial_surface,&
            spine   = axis_curve,&
            origin  = [0.0_rk, 0.0_rk])

        result_surface = loft([axis_curve], degree=1)

        result_surface = loft(&
            sections   = curve_sections,&
            parameters = [0.0_rk, 1.0_rk],&
            degree     = 1)

        result_surface = loft(curve_sections, degree=size(curve_sections))

        result_surface = loft([axis_curve, axis_curve], degree=1)

        result_volume = loft([polynomial_surface], degree=1)

        result_volume = loft(&
            sections   = surface_sections,&
            parameters = [0.0_rk, 1.0_rk],&
            degree     = 1)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface loft parameter count",&
            res      = result_volume%err%ok,&
            expected = .false.,&
            msg      = "A surface loft must require one parameter per section.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0061


    subroutine forcad_geometry_0062(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)

        result_volume = sweep(polynomial_surface, axis_curve)

        result_surface = sweep(invalid_curve, axis_curve)

        result_volume = sweep(invalid_surface, axis_curve)

        result_volume = sweep(&
            profile = polynomial_surface,&
            spine   = axis_curve,&
            origin  = [0.0_rk, 0.0_rk])

        result_surface = loft([axis_curve], degree=1)

        result_surface = loft(&
            sections   = curve_sections,&
            parameters = [0.0_rk, 1.0_rk],&
            degree     = 1)

        result_surface = loft(curve_sections, degree=size(curve_sections))

        result_surface = loft([axis_curve, axis_curve], degree=1)

        result_volume = loft([polynomial_surface], degree=1)

        result_volume = loft(&
            sections   = surface_sections,&
            parameters = [0.0_rk, 1.0_rk],&
            degree     = 1)

        result_volume = loft(surface_sections, degree=size(surface_sections))
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "surface loft degree range",&
            res      = result_volume%err%ok,&
            expected = .false.,&
            msg      = "Surface loft degree range is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0062


    subroutine forcad_geometry_0063(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: geometry_pi = acos(-1.0_rk)
        real(rk), parameter :: geometry_linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: geometry_axis_Xc(2,3) = reshape([&
            2.0_rk,2.0_rk,0.0_rk,0.0_rk,-1.0_rk,1.0_rk], [2,3])
        type(nurbs_curve) :: axis_curve, invalid_curve, curve_sections(4)
        type(nurbs_surface) :: polynomial_surface, result_surface
        type(nurbs_surface) :: invalid_surface, surface_sections(3)
        type(nurbs_volume) :: result_volume

        call axis_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            degree = 1)
        call polynomial_surface%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_curve%set(&
            knot   = geometry_linear_knot,&
            Xc     = geometry_axis_Xc,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        invalid_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])
        curve_sections = axis_curve
        surface_sections = polynomial_surface

        result_surface = extrude(axis_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(polynomial_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = extrude(invalid_curve, [0.0_rk, 1.0_rk, 0.0_rk])

        result_volume = extrude(invalid_surface, [0.0_rk, 0.0_rk, 1.0_rk])

        result_surface = revolve(&
            curve          = axis_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_volume = revolve(&
            surface        = polynomial_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = 0.5_rk*geometry_pi)

        result_surface = revolve(&
            curve          = invalid_curve,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_volume = revolve(&
            surface        = invalid_surface,&
            axis_point     = [0.0_rk, 0.0_rk, 0.0_rk],&
            axis_direction = [0.0_rk, 0.0_rk, 1.0_rk],&
            angle          = geometry_pi)

        result_surface = sweep(axis_curve, axis_curve)

        result_volume = sweep(polynomial_surface, axis_curve)

        result_surface = sweep(invalid_curve, axis_curve)

        result_volume = sweep(invalid_surface, axis_curve)

        result_volume = sweep(&
            profile = polynomial_surface,&
            spine   = axis_curve,&
            origin  = [0.0_rk, 0.0_rk])

        result_surface = loft([axis_curve], degree=1)

        result_surface = loft(&
            sections   = curve_sections,&
            parameters = [0.0_rk, 1.0_rk],&
            degree     = 1)

        result_surface = loft(curve_sections, degree=size(curve_sections))

        result_surface = loft([axis_curve, axis_curve], degree=1)

        result_volume = loft([polynomial_surface], degree=1)

        result_volume = loft(&
            sections   = surface_sections,&
            parameters = [0.0_rk, 1.0_rk],&
            degree     = 1)

        result_volume = loft(surface_sections, degree=size(surface_sections))

        result_volume = loft([polynomial_surface, polynomial_surface], degree=1)
        call axis_curve%err%print()
        call invalid_curve%err%print()
        call curve_sections(1)%err%print()
        call curve_sections(2)%err%print()
        call curve_sections(3)%err%print()
        call curve_sections(4)%err%print()
        call polynomial_surface%err%print()
        call result_surface%err%print()
        call invalid_surface%err%print()
        call surface_sections(1)%err%print()
        call surface_sections(2)%err%print()
        call surface_sections(3)%err%print()
        call result_volume%err%print()

        call ut%test(ti)%check(&
            name     = "polynomial surface loft",&
            res      = result_volume%err%ok .and. .not. result_volume%is_rational(),&
            expected = .true.,&
            msg      = "Polynomial surface loft is incorrect.",&
            group    = "forcad_geometry")
        ti = ti + 1

    end subroutine forcad_geometry_0063


    subroutine run_geometry_tests(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call forcad_geometry_0001(ut, ti)
        call forcad_geometry_0002(ut, ti)
        call forcad_geometry_0003(ut, ti)
        call forcad_geometry_0004(ut, ti)
        call forcad_geometry_0005(ut, ti)
        call forcad_geometry_0006(ut, ti)
        call forcad_geometry_0007(ut, ti)
        call forcad_geometry_0008(ut, ti)
        call forcad_geometry_0009(ut, ti)
        call forcad_geometry_0010(ut, ti)
        call forcad_geometry_0011(ut, ti)
        call forcad_geometry_0012(ut, ti)
        call forcad_geometry_0013(ut, ti)
        call forcad_geometry_0014(ut, ti)
        call forcad_geometry_0015(ut, ti)
        call forcad_geometry_0016(ut, ti)
        call forcad_geometry_0017(ut, ti)
        call forcad_geometry_0018(ut, ti)
        call forcad_geometry_0019(ut, ti)
        call forcad_geometry_0020(ut, ti)
        call forcad_geometry_0021(ut, ti)
        call forcad_geometry_0022(ut, ti)
        call forcad_geometry_0023(ut, ti)
        call forcad_geometry_0024(ut, ti)
        call forcad_geometry_0025(ut, ti)
        call forcad_geometry_0026(ut, ti)
        call forcad_geometry_0027(ut, ti)
        call forcad_geometry_0028(ut, ti)
        call forcad_geometry_0029(ut, ti)
        call forcad_geometry_0030(ut, ti)
        call forcad_geometry_0031(ut, ti)
        call forcad_geometry_0032(ut, ti)
        call forcad_geometry_0033(ut, ti)
        call forcad_geometry_0034(ut, ti)
        call forcad_geometry_0035(ut, ti)
        call forcad_geometry_0036(ut, ti)
        call forcad_geometry_0037(ut, ti)
        call forcad_geometry_0038(ut, ti)
        call forcad_geometry_0039(ut, ti)
        call forcad_geometry_0040(ut, ti)
        call forcad_geometry_0041(ut, ti)
        call forcad_geometry_0042(ut, ti)
        call forcad_geometry_0043(ut, ti)
        call forcad_geometry_0044(ut, ti)
        call forcad_geometry_0045(ut, ti)
        call forcad_geometry_0046(ut, ti)
        call forcad_geometry_0047(ut, ti)
        call forcad_geometry_0048(ut, ti)
        call forcad_geometry_0049(ut, ti)
        call forcad_geometry_0050(ut, ti)
        call forcad_geometry_0051(ut, ti)
        call forcad_geometry_0052(ut, ti)
        call forcad_geometry_0053(ut, ti)
        call forcad_geometry_0054(ut, ti)
        call forcad_geometry_0055(ut, ti)
        call forcad_geometry_0056(ut, ti)
        call forcad_geometry_0057(ut, ti)
        call forcad_geometry_0058(ut, ti)
        call forcad_geometry_0059(ut, ti)
        call forcad_geometry_0060(ut, ti)
        call forcad_geometry_0061(ut, ti)
        call forcad_geometry_0062(ut, ti)
        call forcad_geometry_0063(ut, ti)

    end subroutine run_geometry_tests

end module forcad_test_geometry
