module forcad_test_kinds

    use forcad_kinds, only: rk
    use forunittest, only: unit_tests

    implicit none

    private
    public :: run_kinds_tests

contains

    subroutine forcad_kinds_0001(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call ut%test(ti)%check(&
            name     = "rk positive",&
            res      = merge(1, 0, rk > 0),&
            expected = 1,&
            msg      = "forcad_kinds must expose a valid real kind",&
            group    = "forcad_kinds")
        ti = ti + 1

    end subroutine forcad_kinds_0001


    subroutine forcad_kinds_0002(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call ut%test(ti)%check(&
            name     = "rk precision",&
            res      = merge(1, 0, precision(1.0_rk) >= 6),&
            expected = 1,&
            msg      = "rk precision must support at least single precision",&
            group    = "forcad_kinds")
        ti = ti + 1

    end subroutine forcad_kinds_0002


    subroutine forcad_kinds_0003(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call ut%test(ti)%check(&
            name     = "rk epsilon",&
            res      = merge(1, 0, epsilon(1.0_rk) > 0.0_rk .and. epsilon(1.0_rk) < 1.0_rk),&
            expected = 1,&
            msg      = "rk epsilon must be finite and less than one",&
            group    = "forcad_kinds")
        ti = ti + 1

    end subroutine forcad_kinds_0003


    subroutine forcad_kinds_0004(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call ut%test(ti)%check(&
            name     = "rk range",&
            res      = merge(1, 0, range(1.0_rk) > 0),&
            expected = 1,&
            msg      = "rk range must be positive",&
            group    = "forcad_kinds")
        ti = ti + 1

    end subroutine forcad_kinds_0004


    subroutine run_kinds_tests(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call forcad_kinds_0001(ut, ti)
        call forcad_kinds_0002(ut, ti)
        call forcad_kinds_0003(ut, ti)
        call forcad_kinds_0004(ut, ti)

    end subroutine run_kinds_tests

end module forcad_test_kinds
