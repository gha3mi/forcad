module forcad_test_multipatch

    use forcad_kinds, only: rk
    use forcad_multipatch, only: multipatch_chain_rule_coefficients, multipatch_composition_coefficients, &
        multipatch_connection, multipatch_compatible_trace_space, multipatch_growth_capacity, &
        multipatch_projective_weight_scale, multipatch_valid_reparameterization
    use forunittest, only: unit_tests

    implicit none

    private
    public :: run_multipatch_tests

contains

    subroutine forcad_multipatch_0001(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn

        call conn%set(1, 2, 3, 4, continuity=2, reverse=.true., swap=.true., flip=[.true.])

        call ut%test(ti)%check(&
            name     = "connection patch a",&
            res      = conn%get_patch_a(),&
            expected = 1,&
            msg      = "patch_a metadata was not stored",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0001


    subroutine forcad_multipatch_0002(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn

        call conn%set(1, 2, 3, 4, continuity=2, reverse=.true., swap=.true., flip=[.true.])

        call ut%test(ti)%check(&
            name     = "connection side a",&
            res      = conn%get_side_a(),&
            expected = 2,&
            msg      = "side_a metadata was not stored",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0002


    subroutine forcad_multipatch_0003(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn

        call conn%set(1, 2, 3, 4, continuity=2, reverse=.true., swap=.true., flip=[.true.])

        call ut%test(ti)%check(&
            name     = "connection patch b",&
            res      = conn%get_patch_b(),&
            expected = 3,&
            msg      = "patch_b metadata was not stored",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0003


    subroutine forcad_multipatch_0004(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn

        call conn%set(1, 2, 3, 4, continuity=2, reverse=.true., swap=.true., flip=[.true.])

        call ut%test(ti)%check(&
            name     = "connection side b",&
            res      = conn%get_side_b(),&
            expected = 4,&
            msg      = "side_b metadata was not stored",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0004


    subroutine forcad_multipatch_0005(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn

        call conn%set(1, 2, 3, 4, continuity=2, reverse=.true., swap=.true., flip=[.true.])

        call ut%test(ti)%check(&
            name     = "connection continuity",&
            res      = conn%get_continuity(),&
            expected = 2,&
            msg      = "continuity metadata was not stored",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0005


    subroutine forcad_multipatch_0006(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn

        call conn%set(1, 2, 3, 4, continuity=2, reverse=.true., swap=.true., flip=[.true.])

        call ut%test(ti)%check(&
            name     = "connection orientation",&
            res      = merge(1, 0, conn%is_reversed() .and. conn%is_swapped() .and. conn%is_flipped(1) .and. .not. &
                conn%is_flipped(2) .and. .not. conn%is_flipped(3)),&
            expected = 1,&
            msg      = "orientation metadata was not stored or bounded correctly",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0006


    subroutine forcad_multipatch_0007(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call ut%test(ti)%check(&
            name     = "growth empty",&
            res      = multipatch_growth_capacity(0),&
            expected = 16,&
            msg      = "empty multipatch capacity growth should start at 16",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0007


    subroutine forcad_multipatch_0008(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call ut%test(ti)%check(&
            name     = "growth populated",&
            res      = multipatch_growth_capacity(17),&
            expected = 34,&
            msg      = "nonempty multipatch capacity growth should double",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0008


    subroutine forcad_multipatch_0009(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call ut%test(ti)%check(&
            name     = "trace compatible",&
            res      = multipatch_compatible_trace_space([0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 3.0_rk], 4, 1, [10.0_rk,&
                10.0_rk, 11.0_rk, 12.0_rk, 13.0_rk, 13.0_rk], 4, 1, .false., 1.0e-12_rk),&
            expected = .true.,&
            msg      = "affine-equivalent trace knot vectors should be compatible",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0009


    subroutine forcad_multipatch_0010(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call ut%test(ti)%check(&
            name     = "trace reverse compatible",&
            res      = multipatch_compatible_trace_space([0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 3.0_rk], 4, 1, [0.0_rk, 0.0_rk,&
                1.0_rk, 2.0_rk, 3.0_rk, 3.0_rk], 4, 1, .true., 1.0e-12_rk),&
            expected = .true.,&
            msg      = "Reversed trace spaces must be compatible.",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0010


    subroutine forcad_multipatch_0011(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: weight_scale
        logical :: ok

        call multipatch_projective_weight_scale([2.0_rk, 4.0_rk, 6.0_rk], [1.0_rk, 2.0_rk, 3.0_rk], &
            1.0e-12_rk, weight_scale, ok)

        call ut%test(ti)%check(&
            name     = "weight scale compatible",&
            res      = merge(1, 0, ok .and. abs(weight_scale - 2.0_rk) <= 1.0e-12_rk),&
            expected = 1,&
            msg      = "projective weights should report the common scale",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0011


    subroutine forcad_multipatch_0012(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: weight_scale
        logical :: ok

        call multipatch_projective_weight_scale([2.0_rk, 4.0_rk, 6.0_rk], [1.0_rk, 2.0_rk, 3.0_rk], &
            1.0e-12_rk, weight_scale, ok)
        call multipatch_projective_weight_scale([2.0_rk, 4.0_rk, 7.0_rk], [1.0_rk, 2.0_rk, 3.0_rk], &
            1.0e-12_rk, weight_scale, ok)

        call ut%test(ti)%check(&
            name     = "weight scale incompatible",&
            res      = ok,&
            expected = .false.,&
            msg      = "non-projective weights should be rejected",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0012


    subroutine forcad_multipatch_0013(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn, conn_copy

        call conn%set(1, 2, 3, 4, continuity=2, reverse=.true., swap=.true., flip=[.true.])
        conn_copy = conn
        call conn%set(9, 1, 10, 2, continuity=0)

        call ut%test(ti)%check(&
            name     = "connection assignment",&
            res      = merge(1, 0, conn_copy%get_patch_a() == 1 .and. conn_copy%get_patch_b() == 3 .and. &
                conn_copy%get_continuity() == 2 .and. conn_copy%is_reversed() .and. conn_copy%is_swapped() .and. &
                conn_copy%is_flipped(1)),&
            expected = 1,&
            msg      = "Connection assignment is incorrect.",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0013


    subroutine forcad_multipatch_0014(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn
        real(rk), allocatable :: reparameterization(:)

        call conn%set(&
            patch_a            = 1,&
            side_a             = 2,&
            patch_b            = 2,&
            side_b             = 1,&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        reparameterization = conn%get_reparameterization()

        call ut%test(ti)%check(&
            name     = "G3 connection metadata",&
            res      = conn%is_geometric(),&
            expected = .true.,&
            msg      = "geometric continuity was not stored",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0014


    subroutine forcad_multipatch_0015(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn
        real(rk), allocatable :: reparameterization(:)

        call conn%set(&
            patch_a            = 1,&
            side_a             = 2,&
            patch_b            = 2,&
            side_b             = 1,&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        reparameterization = conn%get_reparameterization()

        call ut%test(ti)%check(&
            name     = "G3 reparameterization metadata",&
            res      = reparameterization,&
            expected = [0.5_rk, -0.25_rk, 0.1875_rk],&
            tol      = 1.0e-14_rk,&
            msg      = "geometric transition jet was not stored",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0015


    subroutine forcad_multipatch_0016(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn
        real(rk), allocatable :: reparameterization(:)

        call conn%set(&
            patch_a            = 1,&
            side_a             = 2,&
            patch_b            = 2,&
            side_b             = 1,&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        reparameterization = conn%get_reparameterization()

        call ut%test(ti)%check(&
            name     = "G3 reparameterization validity",&
            res      = multipatch_valid_reparameterization(conn, 1),&
            expected = .true.,&
            msg      = "a finite orientation-preserving G3 transition jet must be valid",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0016


    subroutine forcad_multipatch_0017(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn
        real(rk), allocatable :: chain(:,:), reparameterization(:)

        call conn%set(&
            patch_a            = 1,&
            side_a             = 2,&
            patch_b            = 2,&
            side_b             = 1,&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        reparameterization = conn%get_reparameterization()

        call multipatch_chain_rule_coefficients(conn, 1, chain)

        call ut%test(ti)%check(&
            name     = "G3 Faa di Bruno coefficients",&
            res      = [chain(0,0), chain(1,1), chain(2,1), chain(2,2), chain(3,1), chain(3,2), chain(3,3)],&
            expected = [1.0_rk, 0.5_rk, -0.25_rk, 0.25_rk, 0.1875_rk, -0.375_rk, 0.125_rk],&
            tol      = 1.0e-14_rk,&
            msg      = "arbitrary-order chain-rule coefficients are incorrect",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0017


    subroutine forcad_multipatch_0018(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn
        real(rk), allocatable :: chain(:,:)

        call conn%set(1, 2, 2, 2, continuity=3, geometric=.false.)
        call multipatch_chain_rule_coefficients(conn, -1, chain)

        call ut%test(ti)%check(&
            name     = "C3 signed affine coefficients",&
            res      = [chain(0,0), chain(1,1), chain(2,2), chain(3,3)],&
            expected = [1.0_rk, -1.0_rk, 1.0_rk, -1.0_rk],&
            tol      = 1.0e-14_rk,&
            msg      = "C3 orientation signs are incorrect",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0018


    subroutine forcad_multipatch_0019(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn
        real(rk), allocatable :: reparameterization(:)

        call conn%set(1, 2, 2, 1, continuity=0)
        reparameterization = conn%get_reparameterization()

        call ut%test(ti)%check(&
            name     = "empty transition jet metadata",&
            res      = size(reparameterization),&
            expected = 0,&
            msg      = "Empty transition jet metadata is incorrect.",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0019


    subroutine forcad_multipatch_0020(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn

        call conn%set(1, 2, 2, 1, continuity=-1, geometric=.true.)

        call ut%test(ti)%check(&
            name     = "negative geometric continuity",&
            res      = multipatch_valid_reparameterization(conn, 1),&
            expected = .false.,&
            msg      = "A geometric transition cannot use negative continuity.",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0020


    subroutine forcad_multipatch_0021(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(multipatch_connection) :: conn

        call conn%set(1, 2, 2, 1, continuity=0, geometric=.true.)

        call ut%test(ti)%check(&
            name     = "G0 empty transition jet",&
            res      = multipatch_valid_reparameterization(conn, 1),&
            expected = .true.,&
            msg      = "G0 empty transition jet is incorrect.",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0021


    subroutine forcad_multipatch_0022(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: transition(2,3)
        real(rk), allocatable :: coefficient(:,:,:,:)

        transition(1,:) = [2.0_rk,5.0_rk,11.0_rk]
        transition(2,:) = [3.0_rk,7.0_rk,13.0_rk]
        call multipatch_composition_coefficients(transition, coefficient)

        call ut%test(ti)%check(&
            name     = "multivariate composition coefficients",&
            res      = [&
                coefficient(2,1,0,0),&
                coefficient(2,0,1,0),&
                coefficient(2,1,1,0),&
                coefficient(3,1,1,0)],&
            expected = [5.0_rk,7.0_rk,12.0_rk,87.0_rk],&
            tol      = 1.0e-14_rk,&
            msg      = "Multivariate Faa di Bruno coefficients are incorrect.",&
            group    = "forcad_multipatch")
        ti = ti + 1

    end subroutine forcad_multipatch_0022


    subroutine run_multipatch_tests(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call forcad_multipatch_0001(ut, ti)
        call forcad_multipatch_0002(ut, ti)
        call forcad_multipatch_0003(ut, ti)
        call forcad_multipatch_0004(ut, ti)
        call forcad_multipatch_0005(ut, ti)
        call forcad_multipatch_0006(ut, ti)
        call forcad_multipatch_0007(ut, ti)
        call forcad_multipatch_0008(ut, ti)
        call forcad_multipatch_0009(ut, ti)
        call forcad_multipatch_0010(ut, ti)
        call forcad_multipatch_0011(ut, ti)
        call forcad_multipatch_0012(ut, ti)
        call forcad_multipatch_0013(ut, ti)
        call forcad_multipatch_0014(ut, ti)
        call forcad_multipatch_0015(ut, ti)
        call forcad_multipatch_0016(ut, ti)
        call forcad_multipatch_0017(ut, ti)
        call forcad_multipatch_0018(ut, ti)
        call forcad_multipatch_0019(ut, ti)
        call forcad_multipatch_0020(ut, ti)
        call forcad_multipatch_0021(ut, ti)
        call forcad_multipatch_0022(ut, ti)

    end subroutine run_multipatch_tests

end module forcad_test_multipatch
