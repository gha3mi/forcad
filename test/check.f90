program check

    use forcad_test_geometry, only: run_geometry_tests
    use forcad_test_kinds, only: run_kinds_tests
    use forcad_test_multipatch, only: run_multipatch_tests
    use forcad_test_multipatch_curve, only: run_multipatch_curve_tests
    use forcad_test_multipatch_surface, only: run_multipatch_surface_tests
    use forcad_test_multipatch_volume, only: run_multipatch_volume_tests
    use forcad_test_nurbs_curve, only: run_nurbs_curve_tests
    use forcad_test_nurbs_surface, only: run_nurbs_surface_tests
    use forcad_test_nurbs_volume, only: run_nurbs_volume_tests
    use forcad_test_utils, only: run_utils_tests
    use forunittest, only: unit_tests

    implicit none

    integer, parameter :: ntests = 1381

    type(unit_tests) :: ut
    integer :: ti

    call ut%initialize(n=ntests)
    ti = 1

    call run_geometry_tests(ut, ti)
    call run_kinds_tests(ut, ti)
    call run_multipatch_tests(ut, ti)
    call run_multipatch_curve_tests(ut, ti)
    call run_multipatch_surface_tests(ut, ti)
    call run_multipatch_volume_tests(ut, ti)
    call run_nurbs_curve_tests(ut, ti)
    call run_nurbs_surface_tests(ut, ti)
    call run_nurbs_volume_tests(ut, ti)
    call run_utils_tests(ut, ti)

    if (ti /= ntests + 1) error stop "check: incorrect test count"
    call ut%summary(verbose=3, required_score=100.0, stop_fail=.true.)

end program check
