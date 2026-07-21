module forcad_test_utils

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_positive_inf, ieee_quiet_nan, &
        ieee_value
    use, intrinsic :: iso_fortran_env, only: file_storage_size, int8
    use forcad_kinds, only: rk
    use forcad_utils, only: basis_bspline, basis_bspline_der, basis_bspline_2der, basis_bspline_der_order, &
        basis_bspline_der_all_active, basis_bernstein, compute_multiplicity, ndgrid, dyad, kron, unique, findspan, &
        compute_knot_vector, insert_knot_A_5_1, elevate_degree_A_5_9, hexahedron_Xc, tetragon_Xc, elemConn_C0, &
        elemConn_Cn, rotation, det, inv, gauss_leg, export_vtk_legacy, solve, solve_spd_banded, sparse_left_matmul, &
        repelem, linspace, eye, kron_eye, valid_knot_vector, knot_start, knot_end, knot_tolerance, active_knots, &
        active_span_count, active_knot_multiplicity, fill_uniform, boundary_index, boundary_layer_index, &
        infer_knot_shape, map_parameter, disjoint_set_union, disjoint_set_map, remove_knots_A_5_8, periodic_topology
    use forunittest, only: unit_tests

    implicit none

    private
    public :: run_utils_tests

contains

    subroutine forcad_utils_0001(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B_ref(4)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call ut%test(ti)%check(&
            name     = "basis_bspline",&
            res      = B4,&
            expected = B_ref,&
            msg      = "Partition of unity and shape check",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0001


    subroutine forcad_utils_0002(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B_ref(4)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call ut%test(ti)%check(&
            name     = "basis_sum",&
            res      = sum(B4),&
            expected = 1.0_rk,&
            msg      = "Partition of unity",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0002


    subroutine forcad_utils_0003(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: B_ref(4)
        real(rk) :: dB_ref(4)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call ut%test(ti)%check(&
            name     = "basis_bspline_der",&
            res      = dB,&
            expected = dB_ref,&
            msg      = "1st derivative shape check",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0003


    subroutine forcad_utils_0004(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: B_ref(4)
        real(rk) :: dB_ref(4)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call ut%test(ti)%check(&
            name     = "basis_bspline_der_B",&
            res      = dB,&
            expected = dB_ref,&
            msg      = "1st derivative alternate",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0004


    subroutine forcad_utils_0005(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call ut%test(ti)%check(&
            name     = "basis_bspline_2der_A",&
            res      = d2B,&
            expected = d2B_ref,&
            msg      = "2nd derivative shape check A",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0005


    subroutine forcad_utils_0006(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call ut%test(ti)%check(&
            name     = "basis_bspline_2der_B",&
            res      = d2B,&
            expected = d2B_ref,&
            msg      = "2nd derivative shape check B",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0006


    subroutine forcad_utils_0007(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        call ut%test(ti)%check(&
            name     = "basis_bspline_2der_C",&
            res      = d2B,&
            expected = d2B_ref,&
            msg      = "2nd derivative shape check C",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0007


    subroutine forcad_utils_0008(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        call ut%test(ti)%check(&
            name     = "basis_out_of_bounds_low",&
            res      = sum(B4),&
            expected = 0.0_rk,&
            msg      = "Out-of-bounds left",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0008


    subroutine forcad_utils_0009(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        call ut%test(ti)%check(&
            name     = "basis_out_of_bounds_high",&
            res      = sum(B4),&
            expected = 0.0_rk,&
            msg      = "Out-of-bounds right",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0009


    subroutine forcad_utils_0010(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        call ut%test(ti)%check(&
            name     = "basis_bernstein_sum",&
            res      = sum(basis_bernstein(0.3_rk, 3)),&
            expected = 1.0_rk,&
            msg      = "Bernstein basis partition of unity",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0010


    subroutine forcad_utils_0011(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        call ut%test(ti)%check(&
            name     = "compute_multiplicity1",&
            res      = compute_multiplicity([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]),&
            expected = [2, 2],&
            msg      = "Multiplicity vector",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0011


    subroutine forcad_utils_0012(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        call ut%test(ti)%check(&
            name     = "compute_multiplicity2",&
            res      = compute_multiplicity([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk], 1.0_rk),&
            expected = 3,&
            msg      = "Multiplicity at value",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0012


    subroutine forcad_utils_0013(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        call ut%test(ti)%check(&
            name     = "ndgrid2_shape",&
            res      = shape(Xt2),&
            expected = [4,2],&
            msg      = "ndgrid2 shape",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0013


    subroutine forcad_utils_0014(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        call ut%test(ti)%check(&
            name     = "ndgrid3_value",&
            res      = Xt3(:,3),&
            expected = [100.0_rk, 100.0_rk, 100.0_rk, 100.0_rk],&
            msg      = "ndgrid3 constant Z",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0014


    subroutine forcad_utils_0015(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        call ut%test(ti)%check(&
            name     = "dyad_check",&
            res      = M,&
            expected = reshape([3.0_rk, 6.0_rk, 4.0_rk, 8.0_rk, 5.0_rk, 10.0_rk], [2,3]),&
            msg      = "Outer product a .dyad. b",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0015


    subroutine forcad_utils_0016(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        call ut%test(ti)%check(&
            name     = "kron_vector",&
            res      = w,&
            expected = [3.0_rk, 4.0_rk, 6.0_rk, 8.0_rk],&
            msg      = "u .kron. v",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0016


    subroutine forcad_utils_0017(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        call ut%test(ti)%check(&
            name     = "kron_matrix",&
            res      = Bk,&
            expected = reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk, 3.0_rk, 4.0_rk, 6.0_rk, 8.0_rk], [4,2]),&
            msg      = "u .kron. A",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0017


    subroutine forcad_utils_0018(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        call ut%test(ti)%check(&
            name     = "unique_integer",&
            res      = unique([1,2,2,3,1,4]),&
            expected = [1,2,3,4],&
            msg      = "Unique integers",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0018


    subroutine forcad_utils_0019(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        call ut%test(ti)%check(&
            name     = "unique_real",&
            res      = unique([1.0_rk, 1.0_rk + 1e-16_rk, 2.0_rk, 1.0_rk, 3.0_rk]),&
            expected = [1.0_rk, 2.0_rk, 3.0_rk],&
            msg      = "Unique real with tolerance",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0019


    subroutine forcad_utils_0020(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        call ut%test(ti)%check(&
            name     = "findspan_middle",&
            res      = findspan(4, 2, 0.5_rk, knot),&
            expected = 3,&
            msg      = "Find span index at Xt=0.5",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0020


    subroutine forcad_utils_0021(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        call ut%test(ti)%check(&
            name     = "compute_knot_vector",&
            res      = compute_knot_vector([0.0_rk, 1.0_rk, 2.0_rk], 2, [-1, 1, -1]),&
            expected = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk],&
            msg      = "Knot vector construction",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0021


    subroutine forcad_utils_0022(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        call ut%test(ti)%check(&
            name     = "repelem",&
            res      = repelem([1.0_rk, 2.0_rk, 3.0_rk], [2, 1, 3]),&
            expected = [1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 3.0_rk, 3.0_rk],&
            msg      = "Repeat vector elements",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0022


    subroutine forcad_utils_0023(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        call ut%test(ti)%check(&
            name     = "hexahedron_Xc_shape",&
            res      = shape(hexahedron_Xc([1.0_rk, 1.0_rk, 1.0_rk], [2,2,2])),&
            expected = [8,3],&
            msg      = "3D grid shape",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0023


    subroutine forcad_utils_0024(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        call ut%test(ti)%check(&
            name     = "tetragon_Xc_shape",&
            res      = shape(tetragon_Xc([1.0_rk, 1.0_rk], [2,2])),&
            expected = [4,3],&
            msg      = "2D grid shape",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0024


    subroutine forcad_utils_0025(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        call ut%test(ti)%check(&
            name     = "rotation_identity",&
            res      = rotation(0.0_rk, 0.0_rk, 0.0_rk),&
            expected = reshape([1.0_rk,0.0_rk,0.0_rk, 0.0_rk,1.0_rk,0.0_rk, 0.0_rk,0.0_rk,1.0_rk], [3,3]),&
            msg      = "Rotation(0,0,0) = I",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0025


    subroutine forcad_utils_0026(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        call ut%test(ti)%check(&
            name     = "det_2x2",&
            res      = det(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])),&
            expected = -2.0_rk,&
            msg      = "Determinant 2x2",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0026


    subroutine forcad_utils_0027(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        call ut%test(ti)%check(&
            name     = "inv_2x2",&
            res      = inv(reshape([4.0_rk, 7.0_rk, 2.0_rk, 6.0_rk], [2,2])),&
            expected = reshape([0.6_rk, -0.7_rk, -0.2_rk, 0.4_rk], [2,2]),&
            msg      = "Inverse of 2x2",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0027


    subroutine forcad_utils_0028(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        call ut%test(ti)%check(&
            name     = "solve_linear_system",&
            res      = R,&
            expected = R_expected,&
            msg      = "Solving A.X = B",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0028


    subroutine forcad_utils_0029(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        call ut%test(ti)%check(&
            name     = "insert_knot_A_5_1_nq",&
            res      = nq,&
            expected = 3,&
            msg      = "Inserted knot, new number of control points",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0029


    subroutine forcad_utils_0030(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call ut%test(ti)%check(&
            name     = "elevate_degree_nc",&
            res      = nc,&
            expected = 3,&
            msg      = "New number of control points after elevation",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0030


    subroutine forcad_utils_0031(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call ut%test(ti)%check(&
            name     = "gauss_legendre_1D_pts",&
            res      = size(vec),&
            expected = 3,&
            msg      = "Number of Gauss points (1D)",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0031


    subroutine forcad_utils_0032(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call ut%test(ti)%check(&
            name     = "gauss_legendre_2D_shape",&
            res      = shape(Xksi),&
            expected = [9,2],&
            msg      = "Gauss points shape (2D)",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0032


    subroutine forcad_utils_0033(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call ut%test(ti)%check(&
            name     = "gauss_legendre_3D_size",&
            res      = size(Xksi,1),&
            expected = 8,&
            msg      = "Number of Gauss points (3D)",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0033


    subroutine forcad_utils_0034(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call ut%test(ti)%check(&
            name     = "binary legacy VTK export",&
            res      = .true.,&
            expected = .true.,&
            msg      = "Binary legacy VTK export must complete successfully.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0034


    subroutine forcad_utils_0035(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        call ut%test(ti)%check(&
            name     = "binary and ASCII legacy VTK export",&
            res      = .true.,&
            expected = .true.,&
            msg      = "Binary and ASCII legacy VTK exports must complete successfully.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0035


    subroutine forcad_utils_0036(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        call ut%test(ti)%check(&
            name     = "linspace_uniform",&
            res      = linspace(0.0_rk, 1.0_rk, 5),&
            expected = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk],&
            msg      = "Uniform linspace from 0 to 1",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0036


    subroutine forcad_utils_0037(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        call ut%test(ti)%check(&
            name     = "kron3_product",&
            res      = out,&
            expected = [12.0_rk, 15.0_rk, 24.0_rk, 30.0_rk],&
            msg      = "kron3 result values",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0037


    subroutine forcad_utils_0038(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call ut%test(ti)%check(&
            name     = "elemConn_C0_L",&
            res      = conn1D,&
            expected = reshape([1,3,2,4,3,5], [2,3]),&
            msg      = "Linear C0 connectivity",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0038


    subroutine forcad_utils_0039(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        call ut%test(ti)%check(&
            name     = "elemConn_Cn_L",&
            res      = conn1D,&
            expected = reshape([1,2,2,3,3,4], [2,3]),&
            msg      = "Linear Cn connectivity",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0039


    subroutine forcad_utils_0040(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call ut%test(ti)%check(&
            name     = "elemConn_C0_S",&
            res      = shape(conn2D),&
            expected = [4,9],&
            msg      = "Surface C0 connectivity shape",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0040


    subroutine forcad_utils_0041(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        call ut%test(ti)%check(&
            name     = "elemConn_Cn_S",&
            res      = shape(conn2D),&
            expected = [4,9],&
            msg      = "Surface Cn connectivity shape",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0041


    subroutine forcad_utils_0042(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call ut%test(ti)%check(&
            name     = "elemConn_C0_V",&
            res      = shape(conn3D),&
            expected = [8,27],&
            msg      = "Volume C0 connectivity shape",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0042


    subroutine forcad_utils_0043(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        call ut%test(ti)%check(&
            name     = "elemConn_Cn_V",&
            res      = shape(conn3D),&
            expected = [8,27],&
            msg      = "Volume Cn connectivity shape",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0043


    subroutine forcad_utils_0044(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        call ut%test(ti)%check(&
            name     = "inv_3x3",&
            res      = matmul(A2, A_inv),&
            expected = eye(3),&
            msg      = "A . inv(A) = I for 3x3",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0044


    subroutine forcad_utils_0045(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        call ut%test(ti)%check(&
            name     = "inv_rectangular_3x2",&
            res      = matmul(A_inv, A2),&
            expected = eye(2),&
            msg      = "inv(A) . A = I for 3x2",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0045


    subroutine forcad_utils_0046(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        call ut%test(ti)%check(&
            name     = "inv_rectangular_2x3",&
            res      = matmul(A2, A_inv),&
            expected = eye(2),&
            msg      = "A . inv(A) = I for 2x3",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0046


    subroutine forcad_utils_0047(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        call ut%test(ti)%check(&
            name     = "inv_identity",&
            res      = A_inv,&
            expected = A2,&
            msg      = "inv(I) = I",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0047


    subroutine forcad_utils_0048(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        call ut%test(ti)%check(&
            name     = "inv_rectangular_2x3_proj",&
            res      = matmul(A_inv, A2),&
            expected = transpose(matmul(A_inv, A2)),&
            msg      = "inv(A) . A is symmetric projection (2x3)",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0048


    subroutine forcad_utils_0049(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        call ut%test(ti)%check(&
            name     = "kron_eye_block_diag",&
            res      = R,&
            expected = reshape([ 1.0_rk,0.0_rk,2.0_rk,0.0_rk, 0.0_rk,1.0_rk,0.0_rk,2.0_rk, 3.0_rk,0.0_rk,4.0_rk,0.0_rk,&
                0.0_rk,3.0_rk,0.0_rk,4.0_rk], [4,4]),&
            msg      = "Kronecker with identity blocks",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0049


    subroutine forcad_utils_0050(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        call ut%test(ti)%check(&
            name     = "kron_t1_t1_values",&
            res      = w2,&
            expected = [-5.0_rk, 3.0_rk, 0.0_rk, 0.0_rk, 10.0_rk, -6.0_rk],&
            msg      = "kron(u,v) concatenates u(i)*v blocks",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0050


    subroutine forcad_utils_0051(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        call ut%test(ti)%check(&
            name     = "kron_t1_t1_size",&
            res      = size(w2),&
            expected = size(u2)*size(v2),&
            msg      = "Length of kron(u,v)",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0051


    subroutine forcad_utils_0052(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        call ut%test(ti)%check(&
            name     = "kron_t1_t1_noncommutative",&
            res      = maxval(abs(w2 - kron(v2, u2))) <= 16.0_rk*epsilon(1.0_rk),&
            expected = .false.,&
            msg      = "kron(u,v) /= kron(v,u) for vectors",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0052


    subroutine forcad_utils_0053(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        call ut%test(ti)%check(&
            name     = "kron3_values_mixed",&
            res      = out,&
            expected = [0.0_rk, 6.0_rk, 0.0_rk, -12.0_rk],&
            msg      = "ordering: (u1*v1*w1, u1*v2*w1, u2*v1*w1, u2*v2*w1)",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0053


    subroutine forcad_utils_0054(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        call ut%test(ti)%check(&
            name     = "kron3_size",&
            res      = size(out),&
            expected = size(u2) * size(v2) * size(K3),&
            msg      = "length of kron3 output",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0054


    subroutine forcad_utils_0055(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        call ut%test(ti)%check(&
            name     = "kron3_associativity_vec",&
            res      = all(abs(out - A) <= epsilon(0.0_rk)),&
            expected = .true.,&
            msg      = "kron3 equals kron(u, kron(v,w))",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0055


    subroutine forcad_utils_0056(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        call ut%test(ti)%check(&
            name     = "kron3_values_ordering",&
            res      = out,&
            expected = [14.0_rk, 0.0_rk, -28.0_rk, -7.0_rk, 0.0_rk, 14.0_rk],&
            msg      = "iterate i (u), then j (v), then k (w)",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0056


    subroutine forcad_utils_0057(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        call ut%test(ti)%check(&
            name     = "kron3_noncommutative",&
            res      = maxval(abs(kron(u2, v2, K3) - kron(v2, u2, K3))) <= 16.0_rk*epsilon(1.0_rk),&
            expected = .false.,&
            msg      = "kron3(u,v,w) /= kron3(v,u,w)",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0057


    subroutine forcad_utils_0058(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call ut%test(ti)%check(&
            name     = "valid_knot_vector_nonopen",&
            res      = valid_knot_vector(knot, nc, degree),&
            expected = .true.,&
            msg      = "Unclamped uniform knot vector should be valid",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0058


    subroutine forcad_utils_0059(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call ut%test(ti)%check(&
            name     = "active_knots_nonopen",&
            res      = active_knots(knot, nc, degree),&
            expected = [2.0_rk, 3.0_rk, 4.0_rk],&
            msg      = "Active knots should exclude inactive end spans",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0059


    subroutine forcad_utils_0060(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call ut%test(ti)%check(&
            name     = "active_knot_multiplicity_nonopen",&
            res      = active_knot_multiplicity(knot, nc, degree),&
            expected = [1, 1, 1],&
            msg      = "Active_knot_multiplicity_nonopen is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0060


    subroutine forcad_utils_0061(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call ut%test(ti)%check(&
            name     = "active_knot_bounds_nonopen",&
            res      = [knot_start(knot, nc, degree), knot_end(knot, nc, degree)],&
            expected = [2.0_rk, 4.0_rk],&
            msg      = "Active parameter domain should use knot(p+1):knot(nc+1)",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0061


    subroutine forcad_utils_0062(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        call ut%test(ti)%check(&
            name     = "basis_nonopen_sum",&
            res      = sum(B4),&
            expected = 1.0_rk,&
            msg      = "Unclamped basis should partition unity inside active domain",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0062


    subroutine forcad_utils_0063(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call ut%test(ti)%check(&
            name     = "basis_nonopen_inactive_span",&
            res      = sum(B4),&
            expected = 0.0_rk,&
            msg      = "Unclamped basis should be zero outside active domain",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0063


    subroutine forcad_utils_0064(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call ut%test(ti)%check(&
            name     = "compute_multiplicity_large_scaled",&
            res      = compute_multiplicity([1.0e6_rk, 1.0e6_rk + 8.0_rk*epsilon(1.0_rk)*1.0e6_rk, 2.0e6_rk, 2.0e6_rk]),&
            expected = [1, 1, 2],&
            msg      = "Distinct representable knots must remain topologically distinct",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0064


    subroutine forcad_utils_0065(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call ut%test(ti)%check(&
            name     = "compute_multiplicity_value_large_scaled",&
            res      = compute_multiplicity([1.0e6_rk, 1.0e6_rk + 8.0_rk*epsilon(1.0_rk)*1.0e6_rk, 2.0e6_rk, 2.0e6_rk], 1.0e6_rk),&
            expected = 1,&
            msg      = "Multiplicity at a value must use exact knot topology",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0065


    subroutine forcad_utils_0066(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call ut%test(ti)%check(&
            name     = "unique_real_large_scaled",&
            res      = unique([1.0e6_rk, 1.0e6_rk + 8.0_rk*epsilon(1.0_rk)*1.0e6_rk, 2.0e6_rk]),&
            expected = [1.0e6_rk, 2.0e6_rk],&
            msg      = "Unique real tolerance should scale with value magnitude",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0066


    subroutine forcad_utils_0067(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call ut%test(ti)%check(&
            name     = "fill_uniform_range",&
            res      = Xu,&
            expected = [-1.0_rk, -0.5_rk, 0.0_rk, 0.5_rk, 1.0_rk],&
            msg      = "Uniform fill should include both end points",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0067


    subroutine forcad_utils_0068(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call ut%test(ti)%check(&
            name     = "fill_uniform_single",&
            res      = Xu_single,&
            expected = [7.0_rk],&
            msg      = "Single-value fill should return the lower bound",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0068


    subroutine forcad_utils_0069(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:)
        integer, allocatable :: conn2D(:,:)
        integer, allocatable :: conn3D(:,:)
        integer, allocatable :: idx1(:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call ut%test(ti)%check(&
            name     = "boundary_index_2d_umax",&
            res      = idx1,&
            expected = [3, 6],&
            msg      = "2D u_max boundary indices",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0069


    subroutine forcad_utils_0070(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call ut%test(ti)%check(&
            name     = "boundary_index_3d_wmax",&
            res      = idx2,&
            expected = [5, 6, 7, 8],&
            msg      = "3D w_max boundary indices",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0070


    subroutine forcad_utils_0071(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call ut%test(ti)%check(&
            name     = "valid_knot_vector_rejects_decreasing",&
            res      = valid_knot_vector([0.0_rk, 1.0_rk, 0.5_rk, 1.0_rk], 2, 1),&
            expected = .false.,&
            msg      = "Knot vectors must be nondecreasing",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0071


    subroutine forcad_utils_0072(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer :: degree2(2)
        integer :: nc2(2)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call ut%test(ti)%check(&
            name     = "infer_knot_shape_2d",&
            res      = merge(1, 0, ok .and. all(degree2 == [1, 1]) .and. all(nc2 == [2, 2])),&
            expected = 1,&
            msg      = "2D knot shape inference failed",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0072


    subroutine forcad_utils_0073(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer :: degree2(2)
        integer :: degree3(3)
        integer :: nc2(2)
        integer :: nc3(3)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 8, degree3, nc3, ok)

        call ut%test(ti)%check(&
            name     = "infer_knot_shape_3d",&
            res      = merge(1, 0, ok .and. all(degree3 == [1, 1, 1]) .and. all(nc3 == [2, 2, 2])),&
            expected = 1,&
            msg      = "3D knot shape inference failed",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0073


    subroutine forcad_utils_0074(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer, allocatable :: dsu_map(:)
        integer :: degree2(2), degree3(3), nc2(2), nc3(3), parent(6)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 8, degree3, nc3, ok)

        parent = [(k, k = 1, 6)]
        call disjoint_set_union(parent, 1, 3)
        call disjoint_set_union(parent, 2, 4)
        call disjoint_set_union(parent, 3, 5)
        call disjoint_set_map(parent, dsu_map)

        call ut%test(ti)%check(&
            name     = "disjoint_set_map",&
            res      = dsu_map,&
            expected = [1, 2, 1, 2, 1, 3],&
            msg      = "Disjoint-set compression/map failed",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0074


    subroutine forcad_utils_0075(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer, allocatable :: dsu_map(:)
        integer :: degree2(2), degree3(3), nc2(2), nc3(3), parent(6)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 8, degree3, nc3, ok)

        parent = [(k, k = 1, 6)]
        call disjoint_set_union(parent, 1, 3)
        call disjoint_set_union(parent, 2, 4)
        call disjoint_set_union(parent, 3, 5)
        call disjoint_set_map(parent, dsu_map)

        call boundary_index([3, 2], 4, .true., idx1)

        call ut%test(ti)%check(&
            name     = "boundary_index_2d_reverse",&
            res      = idx1,&
            expected = [6, 5, 4],&
            msg      = "2D reversed boundary indices",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0075


    subroutine forcad_utils_0076(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer, allocatable :: dsu_map(:)
        integer :: degree2(2), degree3(3), nc2(2), nc3(3), parent(6)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 8, degree3, nc3, ok)

        parent = [(k, k = 1, 6)]
        call disjoint_set_union(parent, 1, 3)
        call disjoint_set_union(parent, 2, 4)
        call disjoint_set_union(parent, 3, 5)
        call disjoint_set_map(parent, dsu_map)

        call boundary_index([3, 2], 4, .true., idx1)

        call boundary_index([2, 3, 2], 2, .true., .true., [.true., .false.], idx2)

        call ut%test(ti)%check(&
            name     = "boundary_index_3d_swap_flip_reverse",&
            res      = idx2,&
            expected = [8, 2, 10, 4, 12, 6],&
            msg      = "3D swapped/flipped/reversed boundary indices",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0076


    subroutine forcad_utils_0077(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer, allocatable :: dsu_map(:)
        integer :: degree2(2), degree3(3), nc2(2), nc3(3), parent(6)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 8, degree3, nc3, ok)

        parent = [(k, k = 1, 6)]
        call disjoint_set_union(parent, 1, 3)
        call disjoint_set_union(parent, 2, 4)
        call disjoint_set_union(parent, 3, 5)
        call disjoint_set_map(parent, dsu_map)

        call boundary_index([3, 2], 4, .true., idx1)

        call boundary_index([2, 3, 2], 2, .true., .true., [.true., .false.], idx2)

        call ut%test(ti)%check(&
            name     = "eye_identity",&
            res      = eye(3),&
            expected = reshape([1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk], [3, 3]),&
            msg      = "Identity matrix construction failed",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0077


    subroutine forcad_utils_0078(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer, allocatable :: dsu_map(:)
        integer :: degree2(2), degree3(3), nc2(2), nc3(3), parent(6)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 8, degree3, nc3, ok)

        parent = [(k, k = 1, 6)]
        call disjoint_set_union(parent, 1, 3)
        call disjoint_set_union(parent, 2, 4)
        call disjoint_set_union(parent, 3, 5)
        call disjoint_set_map(parent, dsu_map)

        call boundary_index([3, 2], 4, .true., idx1)

        call boundary_index([2, 3, 2], 2, .true., .true., [.true., .false.], idx2)

        A_band(0,:) = [2.0_rk, 2.0_rk, 2.0_rk]
        A_band(1,:) = [-1.0_rk, -1.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call ut%test(ti)%check(&
            name     = "solve_spd_banded_tridiagonal",&
            res      = Xband(:,1),&
            expected = [1.0_rk, 1.0_rk, 1.0_rk],&
            msg      = "Banded SPD solver failed on tridiagonal system",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0078


    subroutine forcad_utils_0079(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer, allocatable :: dsu_map(:)
        integer :: degree2(2), degree3(3), nc2(2), nc3(3), parent(6)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 8, degree3, nc3, ok)

        parent = [(k, k = 1, 6)]
        call disjoint_set_union(parent, 1, 3)
        call disjoint_set_union(parent, 2, 4)
        call disjoint_set_union(parent, 3, 5)
        call disjoint_set_map(parent, dsu_map)

        call boundary_index([3, 2], 4, .true., idx1)

        call boundary_index([2, 3, 2], 2, .true., .true., [.true., .false.], idx2)

        A_band(0,:) = [2.0_rk, 2.0_rk, 2.0_rk]
        A_band(1,:) = [-1.0_rk, -1.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call ut%test(ti)%check(&
            name     = "valid_knot_vector_rejects_overfull_run",&
            res      = valid_knot_vector([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk], 2, 1),&
            expected = .false.,&
            msg      = "Knot multiplicity larger than degree+1 should be invalid",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0079


    subroutine forcad_utils_0080(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer, allocatable :: dsu_map(:)
        integer :: degree2(2), degree3(3), nc2(2), nc3(3), parent(6)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 8, degree3, nc3, ok)

        parent = [(k, k = 1, 6)]
        call disjoint_set_union(parent, 1, 3)
        call disjoint_set_union(parent, 2, 4)
        call disjoint_set_union(parent, 3, 5)
        call disjoint_set_map(parent, dsu_map)

        call boundary_index([3, 2], 4, .true., idx1)

        call boundary_index([2, 3, 2], 2, .true., .true., [.true., .false.], idx2)

        A_band(0,:) = [2.0_rk, 2.0_rk, 2.0_rk]
        A_band(1,:) = [-1.0_rk, -1.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call ut%test(ti)%check(&
            name     = "knot_start_end_clamped",&
            res      = [knot_start([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 2, 1), knot_end([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 2, 1)],&
            expected = [0.0_rk, 1.0_rk],&
            msg      = "Clamped active knot bounds failed",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0080


    subroutine forcad_utils_0081(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer, allocatable :: dsu_map(:)
        integer :: degree2(2), degree3(3), nc2(2), nc3(3), parent(6)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 8, degree3, nc3, ok)

        parent = [(k, k = 1, 6)]
        call disjoint_set_union(parent, 1, 3)
        call disjoint_set_union(parent, 2, 4)
        call disjoint_set_union(parent, 3, 5)
        call disjoint_set_map(parent, dsu_map)

        call boundary_index([3, 2], 4, .true., idx1)

        call boundary_index([2, 3, 2], 2, .true., .true., [.true., .false.], idx2)

        A_band(0,:) = [2.0_rk, 2.0_rk, 2.0_rk]
        A_band(1,:) = [-1.0_rk, -1.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call ut%test(ti)%check(&
            name     = "active_span_count_nonopen",&
            res      = active_span_count(knot, nc, degree),&
            expected = 2,&
            msg      = "Active_span_count_nonopen is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0081


    subroutine forcad_utils_0082(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer, allocatable :: dsu_map(:)
        integer :: degree2(2), degree3(3), nc2(2), nc3(3), parent(6)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 8, degree3, nc3, ok)

        parent = [(k, k = 1, 6)]
        call disjoint_set_union(parent, 1, 3)
        call disjoint_set_union(parent, 2, 4)
        call disjoint_set_union(parent, 3, 5)
        call disjoint_set_map(parent, dsu_map)

        call boundary_index([3, 2], 4, .true., idx1)

        call boundary_index([2, 3, 2], 2, .true., .true., [.true., .false.], idx2)

        A_band(0,:) = [2.0_rk, 2.0_rk, 2.0_rk]
        A_band(1,:) = [-1.0_rk, -1.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call ut%test(ti)%check(&
            name     = "active_knot_multiplicity_clamped",&
            res      = active_knot_multiplicity([0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk], 5, 2),&
            expected = [3, 2, 3],&
            msg      = "Active multiplicity should count sorted knot runs in one pass",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0082


    subroutine forcad_utils_0083(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer, allocatable :: dsu_map(:)
        integer :: degree2(2), degree3(3), nc2(2), nc3(3), parent(6)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 8, degree3, nc3, ok)

        parent = [(k, k = 1, 6)]
        call disjoint_set_union(parent, 1, 3)
        call disjoint_set_union(parent, 2, 4)
        call disjoint_set_union(parent, 3, 5)
        call disjoint_set_map(parent, dsu_map)

        call boundary_index([3, 2], 4, .true., idx1)

        call boundary_index([2, 3, 2], 2, .true., .true., [.true., .false.], idx2)

        A_band(0,:) = [2.0_rk, 2.0_rk, 2.0_rk]
        A_band(1,:) = [-1.0_rk, -1.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call ut%test(ti)%check(&
            name     = "repelem_rejects_negative_count",&
            res      = size(repelem([1.0_rk], [-1])),&
            expected = 0,&
            msg      = "Negative repetition counts must not form an invalid result shape",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0083


    subroutine forcad_utils_0084(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: utils_vtk_file = "vtk/test_output.vtk"
        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: dB(4)
        real(rk) :: d2B(4)
        real(rk) :: A4(2,2)
        real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
        real(rk) :: u(2), v(2), w(4)
        real(rk), allocatable :: u2(:), v2(:), w2(:)
        real(rk) :: A2x2(2,2), Bk(4,2)
        real(rk), allocatable :: A(:), vec(:), M(:,:)
        real(rk), allocatable :: X1(:), X2(:), X3(:)
        real(rk) :: Xu(5), Xu_single(1)
        real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
        real(rk), allocatable :: R(:,:), R_expected(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
        integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:), idx1(:), idx2(:)
        integer, allocatable :: dsu_map(:)
        integer :: degree2(2), degree3(3), nc2(2), nc3(3), parent(6)
        integer :: nq, p, rr, s, k
        real(rk), allocatable :: A2(:,:), A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)
        logical :: ok

        degree = 2
        nc     = 4
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        Xt     = 0.5_rk

        B4 = basis_bspline(Xt, knot, nc, degree)
        B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
        dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

        call basis_bspline_der(Xt, knot, nc, degree, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
        d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

        call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

        call basis_bspline_2der(Xt, knot, nc, degree, d2B)

        Xt = -0.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        Xt = 1.1_rk
        B4 = basis_bspline(Xt, knot, nc, degree)

        X1 = [1.0_rk, 2.0_rk]
        X2 = [10.0_rk, 20.0_rk]
        X3 = [100.0_rk]

        call ndgrid(X1, X2, Xt2)
        call ndgrid(X1, X2, X3, Xt3)

        A = [1.0_rk, 2.0_rk]
        vec = [3.0_rk, 4.0_rk, 5.0_rk]
        M = dyad(A, vec)

        u = [1.0_rk, 2.0_rk]
        v = [3.0_rk, 4.0_rk]
        w = kron(u, v)

        A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
        Bk   = kron(u, A2x2)

        A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
        R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
        R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

        p = 2
        rr = 1
        s  = 0
        k  = 2
        allocate(knot_in(0:5))
        knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(0:2,1:2))
        Pw(0,:) = [0.0_rk, 0.0_rk]
        Pw(1,:) = [0.5_rk, 0.5_rk]
        Pw(2,:) = [1.0_rk, 1.0_rk]

        call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

        deallocate(knot_in, knot_out, Pw, Qw)
        knot_in = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        allocate(Pw(1:2,1:2))
        Pw(1,:) = [0.0_rk, 0.0_rk]
        Pw(2,:) = [0.5_rk, 0.5_rk]
        call elevate_degree_A_5_9(t=1, knot=knot_in, degree=1, Xcw=Pw, nc_new=nc, knot_new=knot_out, Xcw_new=Qw)

        call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="binary")

        call export_vtk_legacy(filename=utils_vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
            elemConn=reshape([1,2], [1,2]), vtkCellType=3, encoding="ascii")

        K1 = [1.0_rk, 2.0_rk]
        K2 = [3.0_rk]
        K3 = [4.0_rk, 5.0_rk]
        out = kron(K1, K2, K3)

        conn1D = elemConn_C0(5, 2)

        call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

        conn2D = elemConn_C0(5, 5, 2, 2)

        call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], conn2D)

        conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

        call elemConn_Cn(5, 5, 5, 2, 2, 2, &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [0.0_rk, 0.5_rk, 1.0_rk], &
            [3,1,3], [3,1,3], [3,1,3], conn3D)

        allocate(A2(3,3))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
            0.0_rk, 1.0_rk, 4.0_rk, &
            5.0_rk, 6.0_rk, 0.0_rk], [3,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(3,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(2,3))
        A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
        A_inv = inv(A2)

        deallocate(A2, A_inv)
        allocate(A2(4,4))
        A2 = eye(4)
        A_inv = inv(A2)

        if (allocated(A2)) deallocate(A2)
        allocate(A2(2,2))
        A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk],[2,2])
        R = kron_eye(A2, 2)

        u2 = [-1.0_rk, 0.0_rk, 2.0_rk]
        v2 = [ 5.0_rk, -3.0_rk]
        w2 = kron(u2, v2)

        u2 = [-1.0_rk, 2.0_rk]
        v2 = [0.0_rk, 3.0_rk]
        K3 = [-2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]
        out = kron(u2, v2, K3)
        w2  = kron(v2, K3)
        A   = kron(u2, w2)

        u2 = [2.0_rk, -1.0_rk]
        v2 = [7.0_rk]
        K3 = [1.0_rk, 0.0_rk, -2.0_rk]
        out = kron(u2, v2, K3)

        u2 = [1.0_rk, 2.0_rk]
        v2 = [3.0_rk, 4.0_rk]
        K3 = [5.0_rk]

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        B4 = basis_bspline(2.5_rk, knot, nc, degree)

        B4 = basis_bspline(1.5_rk, knot, nc, degree)

        call fill_uniform(-1.0_rk, 1.0_rk, Xu)

        call fill_uniform(7.0_rk, 9.0_rk, Xu_single)

        call boundary_index([3, 2], 2, .false., idx1)

        call boundary_index([2, 2, 2], 6, .false., .false., [.false., .false.], idx2)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            4, degree2, nc2, ok)

        call infer_knot_shape([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], &
            [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 8, degree3, nc3, ok)

        parent = [(k, k = 1, 6)]
        call disjoint_set_union(parent, 1, 3)
        call disjoint_set_union(parent, 2, 4)
        call disjoint_set_union(parent, 3, 5)
        call disjoint_set_map(parent, dsu_map)

        call boundary_index([3, 2], 4, .true., idx1)

        call boundary_index([2, 3, 2], 2, .true., .true., [.true., .false.], idx2)

        A_band(0,:) = [2.0_rk, 2.0_rk, 2.0_rk]
        A_band(1,:) = [-1.0_rk, -1.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 0.0_rk, 1.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call ut%test(ti)%check(&
            name     = "compute_knot_vector_rejects_negative_repetition",&
            res      = size(compute_knot_vector([0.0_rk, 1.0_rk], 1, [2, 2])),&
            expected = 0,&
            msg      = "Compute_knot_vector_rejects_negative_repetition is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0084


    subroutine forcad_utils_0085(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: degree2(2)
        integer :: nc2(2)
        logical :: ok
        real(rk), allocatable :: empty(:)

        allocate(empty(0))
        call infer_knot_shape(empty, [0.0_rk, 0.0_rk], 0, degree2, nc2, ok)

        call ut%test(ti)%check(&
            name     = "infer_knot_shape_2d_rejects_empty_knot",&
            res      = ok,&
            expected = .false.,&
            msg      = "Infer_knot_shape_2d_rejects_empty_knot is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0085


    subroutine forcad_utils_0086(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer :: degree2(2)
        integer :: degree3(3)
        integer :: nc2(2)
        integer :: nc3(3)
        logical :: ok
        real(rk), allocatable :: empty(:)

        allocate(empty(0))
        call infer_knot_shape(empty, [0.0_rk, 0.0_rk], 0, degree2, nc2, ok)

        call infer_knot_shape(empty, [0.0_rk, 0.0_rk], [0.0_rk, 0.0_rk], 0, degree3, nc3, ok)

        call ut%test(ti)%check(&
            name     = "infer_knot_shape_3d_rejects_empty_knot",&
            res      = ok,&
            expected = .false.,&
            msg      = "Infer_knot_shape_3d_rejects_empty_knot is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0086


    subroutine forcad_utils_0087(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call ut%test(ti)%check(&
            name     = "valid_knot_vector_rejects_nan",&
            res      = valid_knot_vector([0.0_rk, ieee_value(1.0_rk, ieee_quiet_nan), 1.0_rk, 1.0_rk], 2, 1),&
            expected = .false.,&
            msg      = "Valid_knot_vector_rejects_nan is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0087


    subroutine forcad_utils_0088(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call ut%test(ti)%check(&
            name     = "valid_knot_vector_rejects_infinity",&
            res      = valid_knot_vector([0.0_rk, ieee_value(1.0_rk, ieee_positive_inf), 1.0_rk, 1.0_rk], 2, 1),&
            expected = .false.,&
            msg      = "Valid_knot_vector_rejects_infinity is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0088


    subroutine forcad_utils_0089(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        call ut%test(ti)%check(&
            name     = "basis_derivative_order_above_degree",&
            res      = B4,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk],&
            msg      = "Basis_derivative_order_above_degree is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0089


    subroutine forcad_utils_0090(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        call ut%test(ti)%check(&
            name     = "valid_knot_vector_tiny_scale",&
            res      = valid_knot_vector(knot, nc, degree),&
            expected = .true.,&
            msg      = "Valid_knot_vector_tiny_scale is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0090


    subroutine forcad_utils_0091(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        call ut%test(ti)%check(&
            name     = "basis_tiny_scale_partition",&
            res      = sum(B4),&
            expected = 1.0_rk,&
            tol      = 1.0e-12_rk,&
            msg      = "Basis_tiny_scale_partition is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0091


    subroutine forcad_utils_0092(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        call ut%test(ti)%check(&
            name     = "basis_rejects_short_knot_vector",&
            res      = B4,&
            expected = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk],&
            msg      = "Basis_rejects_short_knot_vector is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0092


    subroutine forcad_utils_0093(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        call ut%test(ti)%check(&
            name     = "findspan_rejects_short_knot_vector",&
            res      = findspan(nc - 1, degree, 0.5_rk, [0.0_rk, 0.0_rk, 1.0_rk]),&
            expected = 0,&
            msg      = "Findspan_rejects_short_knot_vector is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0093


    subroutine forcad_utils_0094(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        call ut%test(ti)%check(&
            name     = "findspan_rejects_short_degree_zero_knot_vector",&
            res      = findspan(1, 0, 0.5_rk, [0.0_rk, 1.0_rk]),&
            expected = 0,&
            msg      = "Degree-zero span search still needs knot(n+2) to be present",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0094


    subroutine forcad_utils_0095(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        call ut%test(ti)%check(&
            name     = "basis_rejects_degree_not_less_than_nc",&
            res      = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], 2, 2),&
            expected = [0.0_rk, 0.0_rk],&
            msg      = "B-spline basis should reject invalid spaces with degree >= nc",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0095


    subroutine forcad_utils_0096(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        integer, allocatable :: conn1D(:,:)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call ut%test(ti)%check(&
            name     = "elemConn_Cn_nonopen_actual_multiplicity",&
            res      = conn1D,&
            expected = reshape([1,2,2,3,3,4], [2,3]),&
            msg      = "ElemConn_Cn_nonopen_actual_multiplicity is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0096


    subroutine forcad_utils_0097(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        integer, allocatable :: conn1D(:,:)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call ut%test(ti)%check(&
            name     = "linspace_rejects_nonpositive_count",&
            res      = size(linspace(0.0_rk, 1.0_rk, 0)),&
            expected = 0,&
            msg      = "Linspace_rejects_nonpositive_count is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0097


    subroutine forcad_utils_0098(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        integer, allocatable :: conn1D(:,:)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call ut%test(ti)%check(&
            name     = "elemConn_C0_rejects_invalid_degree",&
            res      = shape(elemConn_C0(4, 0)),&
            expected = [0, 0],&
            msg      = "ElemConn_C0_rejects_invalid_degree is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0098


    subroutine forcad_utils_0099(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        integer, allocatable :: conn1D(:,:)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call ut%test(ti)%check(&
            name     = "elemConn_Cn_rejects_invalid_node_count",&
            res      = shape(conn1D),&
            expected = [0, 0],&
            msg      = "ElemConn_Cn_rejects_invalid_node_count is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0099


    subroutine forcad_utils_0100(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk), allocatable :: R(:,:)
        integer, allocatable :: conn1D(:,:)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call ut%test(ti)%check(&
            name     = "sparse_left_matmul_values",&
            res      = R,&
            expected = reshape([22.0_rk, 10.0_rk, 34.0_rk, 16.0_rk], [2, 2]),&
            msg      = "Sparse_left_matmul_values is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0100


    subroutine forcad_utils_0101(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk), allocatable :: R(:,:)
        integer, allocatable :: conn1D(:,:)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        call ut%test(ti)%check(&
            name     = "sparse_left_matmul_rejects_shape_mismatch",&
            res      = shape(R),&
            expected = [0, 0],&
            msg      = "Sparse_left_matmul_rejects_shape_mismatch is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0101


    subroutine forcad_utils_0102(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        integer, allocatable :: conn1D(:,:)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call ut%test(ti)%check(&
            name     = "basis_repeated_active_start_sum",&
            res      = sum(B5),&
            expected = 1.0_rk,&
            msg      = "Basis_repeated_active_start_sum is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0102


    subroutine forcad_utils_0103(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        integer, allocatable :: conn1D(:,:)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        call ut%test(ti)%check(&
            name     = "basis_repeated_active_start_derivative_sum",&
            res      = sum(dB5),&
            expected = 0.0_rk,&
            msg      = "Basis_repeated_active_start_derivative_sum is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0103


    subroutine forcad_utils_0104(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        integer, allocatable :: conn1D(:,:)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        call ut%test(ti)%check(&
            name     = "det_unsupported_shape_returns_nan",&
            res      = merge(1, 0, ieee_is_nan(det(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3])))),&
            expected = 1,&
            msg      = "Unsupported determinant shapes should not terminate",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0104


    subroutine forcad_utils_0105(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        call ut%test(ti)%check(&
            name     = "inv_singular_returns_empty",&
            res      = shape(A_inv),&
            expected = [0, 0],&
            msg      = "Inv_singular_returns_empty is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0105


    subroutine forcad_utils_0106(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        call ut%test(ti)%check(&
            name     = "solve_nonsquare_returns_empty",&
            res      = shape(A_inv),&
            expected = [0, 0],&
            msg      = "Invalid solve shapes should not terminate",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0106


    subroutine forcad_utils_0107(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        call ut%test(ti)%check(&
            name     = "solve_singular_returns_empty",&
            res      = shape(A_inv),&
            expected = [0, 0],&
            msg      = "Solve_singular_returns_empty is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0107


    subroutine forcad_utils_0108(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        call ut%test(ti)%check(&
            name     = "solve_spd_banded_mismatch_returns_empty",&
            res      = shape(Xband),&
            expected = [0, 0],&
            msg      = "Invalid banded solve shapes should not terminate",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0108


    subroutine forcad_utils_0109(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        A_band(0,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        A_band(1,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 2.0_rk, 3.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call ut%test(ti)%check(&
            name     = "solve_spd_banded_non_spd_returns_empty",&
            res      = shape(Xband),&
            expected = [0, 0],&
            msg      = "Solve_spd_banded_non_spd_returns_empty is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0109


    subroutine forcad_utils_0110(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: X1(:)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        real(rk), allocatable :: Wksi(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        A_band(0,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        A_band(1,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 2.0_rk, 3.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call gauss_leg([1.0_rk, 0.0_rk], 1, X1, Wksi)

        call ut%test(ti)%check(&
            name     = "gauss_invalid_interval_returns_nan_points",&
            res      = merge(1, 0, ieee_is_nan(X1(1))),&
            expected = 1,&
            msg      = "Invalid quadrature intervals should not terminate",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0110


    subroutine forcad_utils_0111(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: X1(:)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        real(rk), allocatable :: Wksi(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        A_band(0,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        A_band(1,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 2.0_rk, 3.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call gauss_leg([1.0_rk, 0.0_rk], 1, X1, Wksi)

        call ut%test(ti)%check(&
            name     = "gauss_invalid_interval_returns_nan_weights",&
            res      = merge(1, 0, ieee_is_nan(Wksi(1))),&
            expected = 1,&
            msg      = "Invalid quadrature intervals should set NaN weights",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0111


    subroutine forcad_utils_0112(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: X1(:)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        real(rk), allocatable :: Wksi(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        A_band(0,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        A_band(1,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 2.0_rk, 3.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call gauss_leg([1.0_rk, 0.0_rk], 1, X1, Wksi)

        call export_vtk_legacy("vtk/test_invalid_export.vtk", reshape([1.0_rk, 2.0_rk, 3.0_rk], [1, 3]), &
            reshape([1, 1, 1], [1, 3]), 12, encoding="invalid")

        call ut%test(ti)%check(&
            name     = "export_vtk_invalid_encoding_returns",&
            res      = 1,&
            expected = 1,&
            msg      = "Invalid VTK export inputs should not terminate",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0112


    subroutine forcad_utils_0113(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: X1(:)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        real(rk), allocatable :: Wksi(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        A_band(0,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        A_band(1,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 2.0_rk, 3.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call gauss_leg([1.0_rk, 0.0_rk], 1, X1, Wksi)

        call export_vtk_legacy("vtk/test_invalid_export.vtk", reshape([1.0_rk, 2.0_rk, 3.0_rk], [1, 3]), &
            reshape([1, 1, 1], [1, 3]), 12, encoding="invalid")

        knot_in = [&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk + 16.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk)]

        call ut%test(ti)%check(&
            name     = "valid_knot_vector_translated_tight",&
            res      = valid_knot_vector(knot_in, 4, 2),&
            expected = .true.,&
            msg      = "Translation must not collapse representable knot spans",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0113


    subroutine forcad_utils_0114(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: X1(:)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        real(rk), allocatable :: Wksi(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        A_band(0,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        A_band(1,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 2.0_rk, 3.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call gauss_leg([1.0_rk, 0.0_rk], 1, X1, Wksi)

        call export_vtk_legacy("vtk/test_invalid_export.vtk", reshape([1.0_rk, 2.0_rk, 3.0_rk], [1, 3]), &
            reshape([1, 1, 1], [1, 3]), 12, encoding="invalid")

        knot_in = [&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk + 16.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk)]

        call ut%test(ti)%check(&
            name     = "active_span_count_translated_tight",&
            res      = active_span_count(knot_in, 4, 2),&
            expected = 2,&
            msg      = "Both representable translated spans must remain active",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0114


    subroutine forcad_utils_0115(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: X1(:)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        real(rk), allocatable :: Wksi(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        A_band(0,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        A_band(1,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 2.0_rk, 3.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call gauss_leg([1.0_rk, 0.0_rk], 1, X1, Wksi)

        call export_vtk_legacy("vtk/test_invalid_export.vtk", reshape([1.0_rk, 2.0_rk, 3.0_rk], [1, 3]), &
            reshape([1, 1, 1], [1, 3]), 12, encoding="invalid")

        knot_in = [&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk + 16.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk)]

        B4 = basis_bspline(1.0e6_rk + 8.0_rk*spacing(1.0e6_rk), knot_in, 4, 2)

        call ut%test(ti)%check(&
            name     = "basis_sum_translated_tight",&
            res      = sum(B4),&
            expected = 1.0_rk,&
            tol      = 64.0_rk*epsilon(1.0_rk),&
            msg      = "Tight translated basis must retain partition of unity",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0115


    subroutine forcad_utils_0116(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: X1(:)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        real(rk), allocatable :: Wksi(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)
        real(rk) :: ders_all(0:6,0:5), derivative_sums(0:6), cubic_derivatives(0:6), cubic_coeff(0:5)
        integer :: first_active
        integer :: derivative_order
        integer :: point

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        A_band(0,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        A_band(1,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 2.0_rk, 3.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call gauss_leg([1.0_rk, 0.0_rk], 1, X1, Wksi)

        call export_vtk_legacy("vtk/test_invalid_export.vtk", reshape([1.0_rk, 2.0_rk, 3.0_rk], [1, 3]), &
            reshape([1, 1, 1], [1, 3]), 12, encoding="invalid")

        knot_in = [&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk + 16.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk)]

        B4 = basis_bspline(1.0e6_rk + 8.0_rk*spacing(1.0e6_rk), knot_in, 4, 2)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(12))
        knot_in(1:6) = 0.0_rk
        knot_in(7:12) = 1.0_rk
        call basis_bspline_der_all_active(&
            Xt     = 0.3_rk,&
            knot   = knot_in,&
            nc     = 6,&
            degree = 5,&
            nder   = 6,&
            first  = first_active,&
            ders   = ders_all)

        cubic_coeff = [0.0_rk, 0.0_rk, 0.0_rk, 0.1_rk, 0.4_rk, 1.0_rk]
        derivative_sums = 0.0_rk
        cubic_derivatives = 0.0_rk
        do derivative_order = 0, 6
            do point = 0, 5
                derivative_sums(derivative_order) = derivative_sums(derivative_order) + ders_all(derivative_order,point)
                cubic_derivatives(derivative_order) = cubic_derivatives(derivative_order) + &
                    ders_all(derivative_order,point)*cubic_coeff(point)
            end do
        end do

        call ut%test(ti)%check(&
            name     = "basis_all_derivatives_partition",&
            res      = derivative_sums,&
            expected = [1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 1024.0_rk*epsilon(1.0_rk),&
            msg      = "Basis_all_derivatives_partition is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0116


    subroutine forcad_utils_0117(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: X1(:)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        real(rk), allocatable :: Wksi(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)
        real(rk) :: ders_all(0:6,0:5), derivative_sums(0:6), cubic_derivatives(0:6), cubic_coeff(0:5)
        integer :: first_active
        integer :: derivative_order
        integer :: point

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        A_band(0,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        A_band(1,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 2.0_rk, 3.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call gauss_leg([1.0_rk, 0.0_rk], 1, X1, Wksi)

        call export_vtk_legacy("vtk/test_invalid_export.vtk", reshape([1.0_rk, 2.0_rk, 3.0_rk], [1, 3]), &
            reshape([1, 1, 1], [1, 3]), 12, encoding="invalid")

        knot_in = [&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk + 16.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk)]

        B4 = basis_bspline(1.0e6_rk + 8.0_rk*spacing(1.0e6_rk), knot_in, 4, 2)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(12))
        knot_in(1:6) = 0.0_rk
        knot_in(7:12) = 1.0_rk
        call basis_bspline_der_all_active(&
            Xt     = 0.3_rk,&
            knot   = knot_in,&
            nc     = 6,&
            degree = 5,&
            nder   = 6,&
            first  = first_active,&
            ders   = ders_all)

        cubic_coeff = [0.0_rk, 0.0_rk, 0.0_rk, 0.1_rk, 0.4_rk, 1.0_rk]
        derivative_sums = 0.0_rk
        cubic_derivatives = 0.0_rk
        do derivative_order = 0, 6
            do point = 0, 5
                derivative_sums(derivative_order) = derivative_sums(derivative_order) + ders_all(derivative_order,point)
                cubic_derivatives(derivative_order) = cubic_derivatives(derivative_order) + &
                    ders_all(derivative_order,point)*cubic_coeff(point)
            end do
        end do

        call ut%test(ti)%check(&
            name     = "basis_all_derivatives_cubic_reproduction",&
            res      = cubic_derivatives,&
            expected = [0.027_rk, 0.27_rk, 1.8_rk, 6.0_rk, 0.0_rk, 0.0_rk, 0.0_rk],&
            tol      = 512.0_rk*epsilon(1.0_rk),&
            msg      = "Basis_all_derivatives_cubic_reproduction is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0117


    subroutine forcad_utils_0118(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: X1(:)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:)
        real(rk), allocatable :: Wksi(:)
        integer, allocatable :: conn1D(:,:)
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)
        real(rk) :: ders_all(0:6,0:5), derivative_sums(0:6), cubic_derivatives(0:6), cubic_coeff(0:5)
        integer :: first_active
        integer :: derivative_order
        integer :: point

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        A_band(0,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        A_band(1,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 2.0_rk, 3.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call gauss_leg([1.0_rk, 0.0_rk], 1, X1, Wksi)

        call export_vtk_legacy("vtk/test_invalid_export.vtk", reshape([1.0_rk, 2.0_rk, 3.0_rk], [1, 3]), &
            reshape([1, 1, 1], [1, 3]), 12, encoding="invalid")

        knot_in = [&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk + 16.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk)]

        B4 = basis_bspline(1.0e6_rk + 8.0_rk*spacing(1.0e6_rk), knot_in, 4, 2)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(12))
        knot_in(1:6) = 0.0_rk
        knot_in(7:12) = 1.0_rk
        call basis_bspline_der_all_active(&
            Xt     = 0.3_rk,&
            knot   = knot_in,&
            nc     = 6,&
            degree = 5,&
            nder   = 6,&
            first  = first_active,&
            ders   = ders_all)

        cubic_coeff = [0.0_rk, 0.0_rk, 0.0_rk, 0.1_rk, 0.4_rk, 1.0_rk]
        derivative_sums = 0.0_rk
        cubic_derivatives = 0.0_rk
        do derivative_order = 0, 6
            do point = 0, 5
                derivative_sums(derivative_order) = derivative_sums(derivative_order) + ders_all(derivative_order,point)
                cubic_derivatives(derivative_order) = cubic_derivatives(derivative_order) + &
                    ders_all(derivative_order,point)*cubic_coeff(point)
            end do
        end do

        call ut%test(ti)%check(&
            name     = "basis_all_derivatives_first_active",&
            res      = first_active,&
            expected = 1,&
            msg      = "Open Bezier basis must start at the first control point",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0118


    subroutine forcad_utils_0119(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: X1(:)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Wksi(:)
        integer, allocatable :: conn1D(:,:)
        integer :: p
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)
        real(rk) :: ders_all(0:6,0:5), derivative_sums(0:6), cubic_derivatives(0:6), cubic_coeff(0:5)
        real(rk) :: weight
        integer :: first_active
        integer :: derivative_order
        integer :: point
        logical :: ok

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        A_band(0,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        A_band(1,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 2.0_rk, 3.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call gauss_leg([1.0_rk, 0.0_rk], 1, X1, Wksi)

        call export_vtk_legacy("vtk/test_invalid_export.vtk", reshape([1.0_rk, 2.0_rk, 3.0_rk], [1, 3]), &
            reshape([1, 1, 1], [1, 3]), 12, encoding="invalid")

        knot_in = [&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk + 16.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk)]

        B4 = basis_bspline(1.0e6_rk + 8.0_rk*spacing(1.0e6_rk), knot_in, 4, 2)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(12))
        knot_in(1:6) = 0.0_rk
        knot_in(7:12) = 1.0_rk
        call basis_bspline_der_all_active(&
            Xt     = 0.3_rk,&
            knot   = knot_in,&
            nc     = 6,&
            degree = 5,&
            nder   = 6,&
            first  = first_active,&
            ders   = ders_all)

        cubic_coeff = [0.0_rk, 0.0_rk, 0.0_rk, 0.1_rk, 0.4_rk, 1.0_rk]
        derivative_sums = 0.0_rk
        cubic_derivatives = 0.0_rk
        do derivative_order = 0, 6
            do point = 0, 5
                derivative_sums(derivative_order) = derivative_sums(derivative_order) + ders_all(derivative_order,point)
                cubic_derivatives(derivative_order) = cubic_derivatives(derivative_order) + &
                    ders_all(derivative_order,point)*cubic_coeff(point)
            end do
        end do

        if (allocated(knot_in)) deallocate(knot_in)
        if (allocated(knot_out)) deallocate(knot_out)
        if (allocated(Pw)) deallocate(Pw)
        if (allocated(Qw)) deallocate(Qw)
        p = 180
        allocate(knot_in(2*(p+1)), Pw(p+1,3))
        knot_in(1:p+1) = 0.0_rk
        knot_in(p+2:) = 1.0_rk
        do point = 1, p+1
            Xt = real(point-1,rk)/real(p,rk)
            weight = 1.0_rk + 0.2_rk*Xt
            Pw(point,1) = weight*(Xt*Xt + 0.1_rk*Xt)
            Pw(point,2) = weight*(Xt - 0.25_rk*Xt*Xt)
            Pw(point,3) = weight
        end do

        call elevate_degree_A_5_9(&
            t        = 3,&
            knot     = knot_in,&
            degree   = p,&
            Xcw      = Pw,&
            nc_new   = nc,&
            knot_new = knot_out,&
            Xcw_new  = Qw,&
            success  = ok)

        call ut%test(ti)%check(&
            name     = "elevate_degree_high_degree_shape",&
            res      = ok .and. nc == p + 4 .and. size(knot_out) == 2*(p+4) .and. all(ieee_is_finite(Qw)),&
            expected = .true.,&
            msg      = "Elevate_degree_high_degree_shape is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0119


    subroutine forcad_utils_0120(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: Xt
        integer  :: nc, degree
        real(rk) :: knot(7)
        real(rk) :: B4(4)
        real(rk) :: B5(5)
        real(rk) :: dB5(5)
        real(rk), allocatable :: X1(:)
        real(rk), allocatable :: R(:,:)
        real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:)
        real(rk), allocatable :: Wksi(:)
        integer, allocatable :: conn1D(:,:)
        integer :: p
        real(rk), allocatable :: A_inv(:,:)
        real(rk), allocatable :: Xband(:,:)
        real(rk) :: A_band(0:1,3), B_rhs(3,1)
        real(rk) :: ders_all(0:6,0:5), derivative_sums(0:6), cubic_derivatives(0:6), cubic_coeff(0:5)
        real(rk) :: sample_parameters(5), old_point(3), new_point(3), elevation_error, weight
        real(rk), allocatable :: B_old(:), B_new(:)
        integer :: first_active, derivative_order, point, sample, component
        logical :: ok

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]

        call basis_bspline_der_order(2.5_rk, knot, nc, degree, degree + 1, B4)

        degree = 2
        nc = 4
        knot = 1.0e-20_rk*[0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk]

        B4 = basis_bspline(0.5e-20_rk, knot, nc, degree)

        B4 = basis_bspline(0.5_rk, [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], nc, degree)

        degree = 2
        nc = 4
        knot = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk]
        call elemConn_Cn(nc, degree, active_knots(knot, nc, degree), active_knot_multiplicity(knot, nc, degree), conn1D)

        call elemConn_Cn(0, 2, [0.0_rk, 1.0_rk], [1], conn1D)

        call sparse_left_matmul(reshape([1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, 0.0_rk], [2, 3]), &
            reshape([4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 9.0_rk], [3, 2]), R)

        call sparse_left_matmul(reshape([1.0_rk, 2.0_rk], [1, 2]), reshape([1.0_rk], [1, 1]), R)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(8))
        knot_in = [0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk]

        B5 = basis_bspline(1.0_rk, knot_in, 5, 2)

        call basis_bspline_der_order(1.0_rk, knot_in, 5, 2, 1, dB5)

        A_inv = inv(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [2, 3]), &
            reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_inv = solve(reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk], [2, 2]), reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        Xband = solve_spd_banded(A_band, reshape([1.0_rk, 2.0_rk], [2, 1]))

        A_band = 0.0_rk
        A_band(0,:) = [0.0_rk, 0.0_rk, 1.0_rk]
        A_band(1,:) = [1.0_rk, 0.0_rk, 0.0_rk]
        B_rhs(:,1) = [1.0_rk, 2.0_rk, 3.0_rk]
        Xband = solve_spd_banded(A_band, B_rhs)

        call gauss_leg([1.0_rk, 0.0_rk], 1, X1, Wksi)

        call export_vtk_legacy("vtk/test_invalid_export.vtk", reshape([1.0_rk, 2.0_rk, 3.0_rk], [1, 3]), &
            reshape([1, 1, 1], [1, 3]), 12, encoding="invalid")

        knot_in = [&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk,&
            1.0e6_rk + 16.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk),&
            1.0e6_rk + 32.0_rk*spacing(1.0e6_rk)]

        B4 = basis_bspline(1.0e6_rk + 8.0_rk*spacing(1.0e6_rk), knot_in, 4, 2)

        if (allocated(knot_in)) deallocate(knot_in)
        allocate(knot_in(12))
        knot_in(1:6) = 0.0_rk
        knot_in(7:12) = 1.0_rk
        call basis_bspline_der_all_active(&
            Xt     = 0.3_rk,&
            knot   = knot_in,&
            nc     = 6,&
            degree = 5,&
            nder   = 6,&
            first  = first_active,&
            ders   = ders_all)

        cubic_coeff = [0.0_rk, 0.0_rk, 0.0_rk, 0.1_rk, 0.4_rk, 1.0_rk]
        derivative_sums = 0.0_rk
        cubic_derivatives = 0.0_rk
        do derivative_order = 0, 6
            do point = 0, 5
                derivative_sums(derivative_order) = derivative_sums(derivative_order) + ders_all(derivative_order,point)
                cubic_derivatives(derivative_order) = cubic_derivatives(derivative_order) + &
                    ders_all(derivative_order,point)*cubic_coeff(point)
            end do
        end do

        if (allocated(knot_in)) deallocate(knot_in)
        if (allocated(knot_out)) deallocate(knot_out)
        if (allocated(Pw)) deallocate(Pw)
        if (allocated(Qw)) deallocate(Qw)
        p = 180
        allocate(knot_in(2*(p+1)), Pw(p+1,3))
        knot_in(1:p+1) = 0.0_rk
        knot_in(p+2:) = 1.0_rk
        do point = 1, p+1
            Xt = real(point-1,rk)/real(p,rk)
            weight = 1.0_rk + 0.2_rk*Xt
            Pw(point,1) = weight*(Xt*Xt + 0.1_rk*Xt)
            Pw(point,2) = weight*(Xt - 0.25_rk*Xt*Xt)
            Pw(point,3) = weight
        end do

        call elevate_degree_A_5_9(&
            t        = 3,&
            knot     = knot_in,&
            degree   = p,&
            Xcw      = Pw,&
            nc_new   = nc,&
            knot_new = knot_out,&
            Xcw_new  = Qw,&
            success  = ok)

        allocate(B_old(p+1), B_new(nc))
        sample_parameters = [0.0_rk, 0.1_rk, 0.37_rk, 0.83_rk, 1.0_rk]
        elevation_error = 0.0_rk
        do sample = 1, size(sample_parameters)
            B_old = basis_bspline(sample_parameters(sample), knot_in, p+1, p)
            B_new = basis_bspline(sample_parameters(sample), knot_out, nc, p+3)
            old_point = 0.0_rk
            new_point = 0.0_rk
            do component = 1, 3
                do point = 1, p+1
                    old_point(component) = old_point(component) + B_old(point)*Pw(point,component)
                end do
                do point = 1, nc
                    new_point(component) = new_point(component) + B_new(point)*Qw(point,component)
                end do
            end do
            elevation_error = max(elevation_error, maxval(abs(new_point-old_point))/max(1.0_rk,maxval(abs(old_point))))
        end do

        call ut%test(ti)%check(&
            name     = "elevate_degree_high_degree_geometry",&
            res      = elevation_error,&
            expected = 0.0_rk,&
            tol      = 4096.0_rk*epsilon(1.0_rk),&
            msg      = "High-degree elevation must preserve homogeneous geometry",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0120


    subroutine forcad_utils_0121(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NFAMILY = 11
        character(len=32), parameter :: family_name(NFAMILY) = [character(len=32) :: &
            "uniform unclamped",&
            "nonuniform unclamped A",&
            "uniform clamped",&
            "nonuniform clamped",&
            "uniform left clamped",&
            "uniform right clamped",&
            "nonuniform left clamped",&
            "nonuniform right clamped",&
            "uniform partially clamped",&
            "nonuniform partially clamped",&
            "nonuniform unclamped B"]
        real(rk), parameter :: knot_family(9,NFAMILY) = reshape([&
            0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk,&
            -1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 2.5_rk, 2.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 3.5_rk, 3.5_rk,&
            0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 6.0_rk,&
            0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 3.1_rk, 3.1_rk, 4.4_rk,&
            0.0_rk, 0.2_rk, 0.7_rk, 1.6_rk, 2.2_rk, 3.8_rk, 4.0_rk, 5.2_rk, 6.1_rk],&
            [9,NFAMILY])
        real(rk) :: family_knot(9)
        integer :: family
        do family = 1, NFAMILY
            family_knot = knot_family(:,family)

            call ut%test(ti)%check(&
                name     = trim(family_name(family))//": valid",&
                res      = valid_knot_vector(family_knot, 6, 2),&
                expected = .true.,&
                msg      = "Each knot-vector family must satisfy the validity rules.",&
                group    = "forcad_utils")
            ti = ti + 1
        end do

    end subroutine forcad_utils_0121


    subroutine forcad_utils_0122(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NFAMILY = 11
        character(len=32), parameter :: family_name(NFAMILY) = [character(len=32) :: &
            "uniform unclamped",&
            "nonuniform unclamped A",&
            "uniform clamped",&
            "nonuniform clamped",&
            "uniform left clamped",&
            "uniform right clamped",&
            "nonuniform left clamped",&
            "nonuniform right clamped",&
            "uniform partially clamped",&
            "nonuniform partially clamped",&
            "nonuniform unclamped B"]
        real(rk), parameter :: knot_family(9,NFAMILY) = reshape([&
            0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk,&
            -1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 2.5_rk, 2.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 3.5_rk, 3.5_rk,&
            0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 6.0_rk,&
            0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 3.1_rk, 3.1_rk, 4.4_rk,&
            0.0_rk, 0.2_rk, 0.7_rk, 1.6_rk, 2.2_rk, 3.8_rk, 4.0_rk, 5.2_rk, 6.1_rk],&
            [9,NFAMILY])
        real(rk) :: family_knot(9)
        real(rk), allocatable :: family_active(:)
        integer, allocatable :: family_multiplicity(:)
        integer :: family
        do family = 1, NFAMILY
            family_knot = knot_family(:,family)

            family_active = active_knots(family_knot, 6, 2)
            family_multiplicity = active_knot_multiplicity(family_knot, 6, 2)

            call ut%test(ti)%check(&
                name     = trim(family_name(family))//": active domain",&
                res      = [family_active(1), family_active(size(family_active))],&
                expected = [family_knot(3), family_knot(7)],&
                tol      = 0.0_rk,&
                msg      = "The active domain must be [U_p,U_(n+1)]",&
                group    = "forcad_utils")
            ti = ti + 1
        end do

    end subroutine forcad_utils_0122


    subroutine forcad_utils_0123(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NFAMILY = 11
        character(len=32), parameter :: family_name(NFAMILY) = [character(len=32) :: &
            "uniform unclamped",&
            "nonuniform unclamped A",&
            "uniform clamped",&
            "nonuniform clamped",&
            "uniform left clamped",&
            "uniform right clamped",&
            "nonuniform left clamped",&
            "nonuniform right clamped",&
            "uniform partially clamped",&
            "nonuniform partially clamped",&
            "nonuniform unclamped B"]
        real(rk), parameter :: knot_family(9,NFAMILY) = reshape([&
            0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk,&
            -1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 2.5_rk, 2.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 3.5_rk, 3.5_rk,&
            0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 6.0_rk,&
            0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 3.1_rk, 3.1_rk, 4.4_rk,&
            0.0_rk, 0.2_rk, 0.7_rk, 1.6_rk, 2.2_rk, 3.8_rk, 4.0_rk, 5.2_rk, 6.1_rk],&
            [9,NFAMILY])
        real(rk) :: family_knot(9)
        real(rk), allocatable :: family_active(:)
        integer, allocatable :: family_multiplicity(:)
        integer :: family
        do family = 1, NFAMILY
            family_knot = knot_family(:,family)

            family_active = active_knots(family_knot, 6, 2)
            family_multiplicity = active_knot_multiplicity(family_knot, 6, 2)

            call ut%test(ti)%check(&
                name     = trim(family_name(family))//": active spans",&
                res      = size(family_active),&
                expected = active_span_count(family_knot, 6, 2) + 1,&
                msg      = "Distinct active knots must delimit every nonzero active span",&
                group    = "forcad_utils")
            ti = ti + 1
        end do

    end subroutine forcad_utils_0123


    subroutine forcad_utils_0124(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NFAMILY = 11
        character(len=32), parameter :: family_name(NFAMILY) = [character(len=32) :: &
            "uniform unclamped",&
            "nonuniform unclamped A",&
            "uniform clamped",&
            "nonuniform clamped",&
            "uniform left clamped",&
            "uniform right clamped",&
            "nonuniform left clamped",&
            "nonuniform right clamped",&
            "uniform partially clamped",&
            "nonuniform partially clamped",&
            "nonuniform unclamped B"]
        real(rk), parameter :: knot_family(9,NFAMILY) = reshape([&
            0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk,&
            -1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 2.5_rk, 2.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 3.5_rk, 3.5_rk,&
            0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 6.0_rk,&
            0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 3.1_rk, 3.1_rk, 4.4_rk,&
            0.0_rk, 0.2_rk, 0.7_rk, 1.6_rk, 2.2_rk, 3.8_rk, 4.0_rk, 5.2_rk, 6.1_rk],&
            [9,NFAMILY])
        real(rk) :: family_knot(9)
        real(rk), allocatable :: family_active(:)
        integer, allocatable :: family_multiplicity(:)
        integer :: family
        do family = 1, NFAMILY
            family_knot = knot_family(:,family)

            family_active = active_knots(family_knot, 6, 2)
            family_multiplicity = active_knot_multiplicity(family_knot, 6, 2)

            call ut%test(ti)%check(&
                name     = trim(family_name(family))//": multiplicities",&
                res      = size(family_multiplicity),&
                expected = size(family_active),&
                msg      = "Every active knot must have one reported multiplicity",&
                group    = "forcad_utils")
            ti = ti + 1
        end do

    end subroutine forcad_utils_0124


    subroutine forcad_utils_0125(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NFAMILY = 11
        real(rk), parameter :: FAMILY_TOL = 4096.0_rk*epsilon(1.0_rk)
        character(len=32), parameter :: family_name(NFAMILY) = [character(len=32) :: &
            "uniform unclamped",&
            "nonuniform unclamped A",&
            "uniform clamped",&
            "nonuniform clamped",&
            "uniform left clamped",&
            "uniform right clamped",&
            "nonuniform left clamped",&
            "nonuniform right clamped",&
            "uniform partially clamped",&
            "nonuniform partially clamped",&
            "nonuniform unclamped B"]
        real(rk), parameter :: knot_family(9,NFAMILY) = reshape([&
            0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk,&
            -1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 2.5_rk, 2.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 3.5_rk, 3.5_rk,&
            0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 6.0_rk,&
            0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 3.1_rk, 3.1_rk, 4.4_rk,&
            0.0_rk, 0.2_rk, 0.7_rk, 1.6_rk, 2.2_rk, 3.8_rk, 4.0_rk, 5.2_rk, 6.1_rk],&
            [9,NFAMILY])
        real(rk) :: family_knot(9), family_parameter, family_derivative(6)
        real(rk), allocatable :: family_active(:), family_basis(:)
        integer, allocatable :: family_multiplicity(:)
        integer :: family
        do family = 1, NFAMILY
            family_knot = knot_family(:,family)

            family_active = active_knots(family_knot, 6, 2)
            family_multiplicity = active_knot_multiplicity(family_knot, 6, 2)

            family_parameter = 0.5_rk*(family_active(1) + family_active(2))
            family_basis = basis_bspline(family_parameter, family_knot, 6, 2)
            call basis_bspline_der_order(&
                Xt     = family_parameter,&
                knot   = family_knot,&
                nc     = 6,&
                degree = 2,&
                nder   = 1,&
                B      = family_derivative)

            call ut%test(ti)%check(&
                name     = trim(family_name(family))//": partition",&
                res      = sum(family_basis),&
                expected = 1.0_rk,&
                tol      = FAMILY_TOL,&
                msg      = "Basis functions must form a partition of unity in an active span",&
                group    = "forcad_utils")
            ti = ti + 1
        end do

    end subroutine forcad_utils_0125


    subroutine forcad_utils_0126(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NFAMILY = 11
        real(rk), parameter :: FAMILY_TOL = 4096.0_rk*epsilon(1.0_rk)
        character(len=32), parameter :: family_name(NFAMILY) = [character(len=32) :: &
            "uniform unclamped",&
            "nonuniform unclamped A",&
            "uniform clamped",&
            "nonuniform clamped",&
            "uniform left clamped",&
            "uniform right clamped",&
            "nonuniform left clamped",&
            "nonuniform right clamped",&
            "uniform partially clamped",&
            "nonuniform partially clamped",&
            "nonuniform unclamped B"]
        real(rk), parameter :: knot_family(9,NFAMILY) = reshape([&
            0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk,&
            -1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 2.5_rk, 2.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 3.5_rk, 3.5_rk,&
            0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 6.0_rk,&
            0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 3.1_rk, 3.1_rk, 4.4_rk,&
            0.0_rk, 0.2_rk, 0.7_rk, 1.6_rk, 2.2_rk, 3.8_rk, 4.0_rk, 5.2_rk, 6.1_rk],&
            [9,NFAMILY])
        real(rk) :: family_knot(9), family_parameter, family_derivative(6)
        real(rk), allocatable :: family_active(:), family_basis(:)
        integer, allocatable :: family_multiplicity(:)
        integer :: family
        do family = 1, NFAMILY
            family_knot = knot_family(:,family)

            family_active = active_knots(family_knot, 6, 2)
            family_multiplicity = active_knot_multiplicity(family_knot, 6, 2)

            family_parameter = 0.5_rk*(family_active(1) + family_active(2))
            family_basis = basis_bspline(family_parameter, family_knot, 6, 2)
            call basis_bspline_der_order(&
                Xt     = family_parameter,&
                knot   = family_knot,&
                nc     = 6,&
                degree = 2,&
                nder   = 1,&
                B      = family_derivative)

            call ut%test(ti)%check(&
                name     = trim(family_name(family))//": derivative sum",&
                res      = sum(family_derivative),&
                expected = 0.0_rk,&
                tol      = FAMILY_TOL,&
                msg      = "First basis derivatives must sum to zero in an active span",&
                group    = "forcad_utils")
            ti = ti + 1
        end do

    end subroutine forcad_utils_0126


    subroutine forcad_utils_0127(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NFAMILY = 11
        real(rk), parameter :: knot_family(9,NFAMILY) = reshape([&
            0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk,&
            -1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 2.5_rk, 2.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 3.5_rk, 3.5_rk,&
            0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 6.0_rk,&
            0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 3.1_rk, 3.1_rk, 4.4_rk,&
            0.0_rk, 0.2_rk, 0.7_rk, 1.6_rk, 2.2_rk, 3.8_rk, 4.0_rk, 5.2_rk, 6.1_rk],&
            [9,NFAMILY])
        real(rk) :: family_knot(9), family_parameter, family_derivative(6)
        real(rk) :: map_knot(8), map_Pw(5,2)
        real(rk), allocatable :: family_active(:), family_basis(:)
        real(rk), allocatable :: map_knot_new(:), map_Pw_new(:,:), insertion_map(:,:)
        integer, allocatable :: family_multiplicity(:)
        integer :: family, map_last
        do family = 1, NFAMILY
            family_knot = knot_family(:,family)

            family_active = active_knots(family_knot, 6, 2)
            family_multiplicity = active_knot_multiplicity(family_knot, 6, 2)

            family_parameter = 0.5_rk*(family_active(1) + family_active(2))
            family_basis = basis_bspline(family_parameter, family_knot, 6, 2)
            call basis_bspline_der_order(&
                Xt     = family_parameter,&
                knot   = family_knot,&
                nc     = 6,&
                degree = 2,&
                nder   = 1,&
                B      = family_derivative)
        end do

        map_knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        map_Pw = reshape([&
            0.0_rk, 0.0_rk,&
            0.2_rk, 0.5_rk,&
            0.5_rk, 0.2_rk,&
            0.8_rk, 0.7_rk,&
            1.0_rk, 1.0_rk], shape(map_Pw), order=[2,1])
        call insert_knot_A_5_1(&
            p  = 2,&
            UP = map_knot,&
            Pw = map_Pw,&
            u  = 0.5_rk,&
            k  = 4,&
            s  = 2,&
            r  = 1,&
            nq = map_last,&
            UQ = map_knot_new,&
            Qw = map_Pw_new,&
            T  = insertion_map)

        call ut%test(ti)%check(&
            name     = "full-break insertion last index",&
            res      = map_last,&
            expected = 5,&
            msg      = "Full-break insertion last index is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0127


    subroutine forcad_utils_0128(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, parameter :: NFAMILY = 11
        real(rk), parameter :: knot_family(9,NFAMILY) = reshape([&
            0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 8.0_rk,&
            -1.4_rk, 0.0_rk, 0.3_rk, 1.1_rk, 2.6_rk, 3.0_rk, 4.4_rk, 4.7_rk, 5.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 2.5_rk, 2.5_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk,&
            0.0_rk, 0.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 4.1_rk, 5.0_rk,&
            -2.0_rk,-1.0_rk, 0.0_rk, 0.2_rk, 1.4_rk, 2.0_rk, 3.5_rk, 3.5_rk, 3.5_rk,&
            0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 6.0_rk,&
            0.0_rk, 0.0_rk, 0.4_rk, 1.3_rk, 1.7_rk, 2.5_rk, 3.1_rk, 3.1_rk, 4.4_rk,&
            0.0_rk, 0.2_rk, 0.7_rk, 1.6_rk, 2.2_rk, 3.8_rk, 4.0_rk, 5.2_rk, 6.1_rk],&
            [9,NFAMILY])
        real(rk) :: family_knot(9), family_parameter, family_derivative(6)
        real(rk) :: map_knot(8), map_Pw(5,2)
        real(rk), allocatable :: family_active(:), family_basis(:)
        real(rk), allocatable :: map_knot_new(:), map_Pw_new(:,:), insertion_map(:,:)
        integer, allocatable :: family_multiplicity(:)
        integer :: family, map_last
        do family = 1, NFAMILY
            family_knot = knot_family(:,family)

            family_active = active_knots(family_knot, 6, 2)
            family_multiplicity = active_knot_multiplicity(family_knot, 6, 2)

            family_parameter = 0.5_rk*(family_active(1) + family_active(2))
            family_basis = basis_bspline(family_parameter, family_knot, 6, 2)
            call basis_bspline_der_order(&
                Xt     = family_parameter,&
                knot   = family_knot,&
                nc     = 6,&
                degree = 2,&
                nder   = 1,&
                B      = family_derivative)
        end do

        map_knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        map_Pw = reshape([&
            0.0_rk, 0.0_rk,&
            0.2_rk, 0.5_rk,&
            0.5_rk, 0.2_rk,&
            0.8_rk, 0.7_rk,&
            1.0_rk, 1.0_rk], shape(map_Pw), order=[2,1])
        call insert_knot_A_5_1(&
            p  = 2,&
            UP = map_knot,&
            Pw = map_Pw,&
            u  = 0.5_rk,&
            k  = 4,&
            s  = 2,&
            r  = 1,&
            nq = map_last,&
            UQ = map_knot_new,&
            Qw = map_Pw_new,&
            T  = insertion_map)

        call ut%test(ti)%check(&
            name     = "full-break insertion transformation",&
            res      = map_Pw_new,&
            expected = matmul(insertion_map, map_Pw),&
            tol      = 16.0_rk*epsilon(1.0_rk),&
            msg      = "Full-break insertion transformation is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0128


    subroutine forcad_utils_0129(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call ut%test(ti)%check(&
            name     = "boundary_layer_index_1d_left",&
            res      = idx1,&
            expected = [3],&
            msg      = "The second layer from the left boundary must use index three",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0129


    subroutine forcad_utils_0130(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call ut%test(ti)%check(&
            name     = "boundary_layer_index_1d_right",&
            res      = idx1,&
            expected = [6],&
            msg      = "The first layer from the right boundary must use index six",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0130


    subroutine forcad_utils_0131(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call ut%test(ti)%check(&
            name     = "boundary_layer_index_1d_invalid_side",&
            res      = idx1,&
            expected = [0],&
            msg      = "Boundary_layer_index_1d_invalid_side is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0131


    subroutine forcad_utils_0132(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: utils_periodic_knot_test(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        integer, allocatable :: idx1(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call ut%test(ti)%check(&
            name     = "map_parameter_wrapping_interior",&
            res      = map_parameter(3.5_rk, utils_periodic_knot_test, 6, 2, .true.),&
            expected = 3.5_rk,&
            tol      = 0.0_rk,&
            msg      = "An interior periodic parameter must remain unchanged",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0132


    subroutine forcad_utils_0133(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: utils_periodic_knot_test(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        integer, allocatable :: idx1(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call ut%test(ti)%check(&
            name     = "map_parameter_wrapping_endpoint",&
            res      = map_parameter(6.0_rk, utils_periodic_knot_test, 6, 2, .true.),&
            expected = 2.0_rk,&
            tol      = 0.0_rk,&
            msg      = "The periodic upper endpoint must map to the lower endpoint",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0133


    subroutine forcad_utils_0134(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: utils_periodic_knot_test(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        integer, allocatable :: idx1(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call ut%test(ti)%check(&
            name     = "map_parameter_wrapping_above_domain",&
            res      = map_parameter(10.5_rk, utils_periodic_knot_test, 6, 2, .true.),&
            expected = 2.5_rk,&
            tol      = 0.0_rk,&
            msg      = "Map_parameter_wrapping_above_domain is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0134


    subroutine forcad_utils_0135(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: utils_periodic_knot_test(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        integer, allocatable :: idx1(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call ut%test(ti)%check(&
            name     = "map_parameter_wrapping_below_domain",&
            res      = map_parameter(-1.5_rk, utils_periodic_knot_test, 6, 2, .true.),&
            expected = 2.5_rk,&
            tol      = 0.0_rk,&
            msg      = "Map_parameter_wrapping_below_domain is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0135


    subroutine forcad_utils_0136(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: utils_periodic_knot_test(9) = [&
            0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk,8.0_rk]
        integer, allocatable :: idx1(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call ut%test(ti)%check(&
            name     = "map_parameter_without_wrapping",&
            res      = map_parameter(10.5_rk, utils_periodic_knot_test, 6, 2, .false.),&
            expected = 10.5_rk,&
            tol      = 0.0_rk,&
            msg      = "Map_parameter_without_wrapping is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0136


    subroutine forcad_utils_0137(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call boundary_index([3, 2], 3, .false., idx1)

        call ut%test(ti)%check(&
            name     = "boundary_index_2d_vmin",&
            res      = idx1,&
            expected = [1, 2, 3],&
            msg      = "Boundary_index_2d_vmin is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0137


    subroutine forcad_utils_0138(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)
        integer, allocatable :: idx2(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call boundary_index([3, 2], 3, .false., idx1)

        call boundary_index([2, 3, 2], 3, .false., .false., [.false., .false.], idx2)

        call ut%test(ti)%check(&
            name     = "boundary_index_3d_vmin",&
            res      = idx2,&
            expected = [1, 2, 7, 8],&
            msg      = "The minimum-v volume face must preserve u-fast control ordering.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0138


    subroutine forcad_utils_0139(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)
        integer, allocatable :: idx2(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call boundary_index([3, 2], 3, .false., idx1)

        call boundary_index([2, 3, 2], 3, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 5, .false., .false., [.false., .false.], idx2)

        call ut%test(ti)%check(&
            name     = "boundary_index_3d_wmin",&
            res      = idx2,&
            expected = [1, 2, 3, 4, 5, 6],&
            msg      = "The minimum-w volume face must contain the first control layer.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0139


    subroutine forcad_utils_0140(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)
        integer, allocatable :: idx2(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call boundary_index([3, 2], 3, .false., idx1)

        call boundary_index([2, 3, 2], 3, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 5, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 7, .false., .false., [.false., .false.], idx2)

        call ut%test(ti)%check(&
            name     = "boundary_index_3d_invalid_side",&
            res      = size(idx2),&
            expected = 0,&
            msg      = "An invalid volume side must return an empty boundary index set.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0140


    subroutine forcad_utils_0141(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)
        integer, allocatable :: idx2(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call boundary_index([3, 2], 3, .false., idx1)

        call boundary_index([2, 3, 2], 3, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 5, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 7, .false., .false., [.false., .false.], idx2)

        call boundary_layer_index([3, 4], 3, 1, .false., idx1)

        call ut%test(ti)%check(&
            name     = "boundary_layer_index_2d_vmin",&
            res      = idx1,&
            expected = [4, 5, 6],&
            msg      = "Boundary_layer_index_2d_vmin is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0141


    subroutine forcad_utils_0142(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)
        integer, allocatable :: idx2(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call boundary_index([3, 2], 3, .false., idx1)

        call boundary_index([2, 3, 2], 3, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 5, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 7, .false., .false., [.false., .false.], idx2)

        call boundary_layer_index([3, 4], 3, 1, .false., idx1)

        call boundary_layer_index([3, 4], 4, 1, .true., idx1)

        call ut%test(ti)%check(&
            name     = "boundary_layer_index_2d_vmax_reverse",&
            res      = idx1,&
            expected = [9, 8, 7],&
            msg      = "Boundary_layer_index_2d_vmax_reverse is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0142


    subroutine forcad_utils_0143(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)
        integer, allocatable :: idx2(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call boundary_index([3, 2], 3, .false., idx1)

        call boundary_index([2, 3, 2], 3, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 5, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 7, .false., .false., [.false., .false.], idx2)

        call boundary_layer_index([3, 4], 3, 1, .false., idx1)

        call boundary_layer_index([3, 4], 4, 1, .true., idx1)

        call boundary_layer_index([2, 4, 2], 3, 1, .false., .false., [.false., .false.], idx2)

        call ut%test(ti)%check(&
            name     = "boundary_layer_index_3d_vmin",&
            res      = idx2,&
            expected = [3, 4, 11, 12],&
            msg      = "The first interior minimum-v volume layer has incorrect indices.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0143


    subroutine forcad_utils_0144(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)
        integer, allocatable :: idx2(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call boundary_index([3, 2], 3, .false., idx1)

        call boundary_index([2, 3, 2], 3, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 5, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 7, .false., .false., [.false., .false.], idx2)

        call boundary_layer_index([3, 4], 3, 1, .false., idx1)

        call boundary_layer_index([3, 4], 4, 1, .true., idx1)

        call boundary_layer_index([2, 4, 2], 3, 1, .false., .false., [.false., .false.], idx2)

        call boundary_layer_index([2, 3, 3], 5, 1, .false., .false., [.false., .false.], idx2)

        call ut%test(ti)%check(&
            name     = "boundary_layer_index_3d_wmin",&
            res      = idx2,&
            expected = [7, 8, 9, 10, 11, 12],&
            msg      = "The first interior minimum-w volume layer has incorrect indices.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0144


    subroutine forcad_utils_0145(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)
        integer, allocatable :: idx2(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call boundary_index([3, 2], 3, .false., idx1)

        call boundary_index([2, 3, 2], 3, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 5, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 7, .false., .false., [.false., .false.], idx2)

        call boundary_layer_index([3, 4], 3, 1, .false., idx1)

        call boundary_layer_index([3, 4], 4, 1, .true., idx1)

        call boundary_layer_index([2, 4, 2], 3, 1, .false., .false., [.false., .false.], idx2)

        call boundary_layer_index([2, 3, 3], 5, 1, .false., .false., [.false., .false.], idx2)

        call boundary_layer_index([2, 3, 3], 7, 1, .false., .false., [.false., .false.], idx2)

        call ut%test(ti)%check(&
            name     = "boundary_layer_index_3d_invalid_side",&
            res      = size(idx2),&
            expected = 0,&
            msg      = "Boundary_layer_index_3d_invalid_side is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0145


    subroutine forcad_utils_0146(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        integer, allocatable :: idx1(:)
        integer, allocatable :: idx2(:)

        call boundary_layer_index(7, 1, 2, .false., idx1)

        call boundary_layer_index(7, 2, 1, .true., idx1)

        call boundary_layer_index(7, 3, 0, .false., idx1)

        call boundary_index([3, 2], 3, .false., idx1)

        call boundary_index([2, 3, 2], 3, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 5, .false., .false., [.false., .false.], idx2)

        call boundary_index([2, 3, 2], 7, .false., .false., [.false., .false.], idx2)

        call boundary_layer_index([3, 4], 3, 1, .false., idx1)

        call boundary_layer_index([3, 4], 4, 1, .true., idx1)

        call boundary_layer_index([2, 4, 2], 3, 1, .false., .false., [.false., .false.], idx2)

        call boundary_layer_index([2, 3, 3], 5, 1, .false., .false., [.false., .false.], idx2)

        call boundary_layer_index([2, 3, 3], 7, 1, .false., .false., [.false., .false.], idx2)

        call ut%test(ti)%check(&
            name     = "knot_tolerance_reversed_bounds",&
            res      = knot_tolerance([0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk], 5, 1) > 0.0_rk,&
            expected = .true.,&
            msg      = "Knot tolerance must normalize reversed index bounds.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0146


    subroutine forcad_utils_0147(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), allocatable :: empty_knot(:)

        allocate(empty_knot(0))

        call ut%test(ti)%check(&
            name     = "knot_tolerance_empty_vector",&
            res      = knot_tolerance(empty_knot, 1, 1),&
            expected = 0.0_rk,&
            tol      = 0.0_rk,&
            msg      = "An empty knot vector must have zero comparison tolerance.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0147


    subroutine forcad_utils_0148(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: identity_knot(0:3), identity_controls(2,2)
        real(rk), allocatable :: refined_knot(:), refined_controls(:,:), refinement_map(:,:)
        integer :: last_control

        identity_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        identity_controls = reshape([0.0_rk, 1.0_rk, 0.0_rk, 2.0_rk], [2,2])
        call insert_knot_A_5_1(&
            p  = 1,&
            UP = identity_knot,&
            Pw = identity_controls,&
            u  = 0.5_rk,&
            k  = 1,&
            s  = 0,&
            r  = 0,&
            nq = last_control,&
            UQ = refined_knot,&
            Qw = refined_controls,&
            T  = refinement_map)

        call ut%test(ti)%check(&
            name     = "zero knot insertion identity",&
            res      = last_control == 1 .and. maxval(abs(refined_knot-identity_knot)) <= tiny(1.0_rk) .and. &
                maxval(abs(refined_controls-identity_controls)) <= tiny(1.0_rk) .and. maxval(abs(refinement_map-eye(2))) <= &
                tiny(1.0_rk),&
            expected = .true.,&
            msg      = "Zero knot insertion identity is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0148


    subroutine forcad_utils_0149(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: identity_knot(0:3), identity_controls(2,2)
        real(rk), allocatable :: refined_knot(:), refined_controls(:,:), refinement_map(:,:)
        integer :: last_control, elevated_nc
        logical :: success

        identity_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        identity_controls = reshape([0.0_rk, 1.0_rk, 0.0_rk, 2.0_rk], [2,2])
        call insert_knot_A_5_1(&
            p  = 1,&
            UP = identity_knot,&
            Pw = identity_controls,&
            u  = 0.5_rk,&
            k  = 1,&
            s  = 0,&
            r  = 0,&
            nq = last_control,&
            UQ = refined_knot,&
            Qw = refined_controls,&
            T  = refinement_map)

        call elevate_degree_A_5_9(&
            t        = 0,&
            knot     = identity_knot,&
            degree   = 1,&
            Xcw      = identity_controls,&
            nc_new   = elevated_nc,&
            knot_new = refined_knot,&
            Xcw_new  = refined_controls,&
            Tmap     = refinement_map,&
            success  = success)

        call ut%test(ti)%check(&
            name     = "zero degree elevation identity",&
            res      = success .and. elevated_nc == 2 .and. maxval(abs(refined_knot-identity_knot)) <= tiny(1.0_rk) .and. &
                maxval(abs(refined_controls-identity_controls)) <= tiny(1.0_rk) .and. maxval(abs(refinement_map-eye(2))) <= &
                tiny(1.0_rk),&
            expected = .true.,&
            msg      = "Zero degree elevation identity is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0149


    subroutine forcad_utils_0150(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: identity_knot(0:3), identity_controls(2,2)
        real(rk), allocatable :: refined_knot(:), refined_controls(:,:), refinement_map(:,:)
        integer :: last_control, elevated_nc
        logical :: success

        identity_knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        identity_controls = reshape([0.0_rk, 1.0_rk, 0.0_rk, 2.0_rk], [2,2])
        call insert_knot_A_5_1(&
            p  = 1,&
            UP = identity_knot,&
            Pw = identity_controls,&
            u  = 0.5_rk,&
            k  = 1,&
            s  = 0,&
            r  = 0,&
            nq = last_control,&
            UQ = refined_knot,&
            Qw = refined_controls,&
            T  = refinement_map)

        call elevate_degree_A_5_9(&
            t        = 0,&
            knot     = identity_knot,&
            degree   = 1,&
            Xcw      = identity_controls,&
            nc_new   = elevated_nc,&
            knot_new = refined_knot,&
            Xcw_new  = refined_controls,&
            Tmap     = refinement_map,&
            success  = success)

        call elevate_degree_A_5_9(&
            t        = -1,&
            knot     = identity_knot,&
            degree   = 1,&
            Xcw      = identity_controls,&
            nc_new   = elevated_nc,&
            knot_new = refined_knot,&
            Xcw_new  = refined_controls,&
            Tmap     = refinement_map,&
            success  = success)

        call ut%test(ti)%check(&
            name     = "negative degree elevation diagnostic",&
            res      = .not. success .and. elevated_nc == 0 .and. size(refined_knot) == 0 .and. size(refined_controls,1) == 0 &
                .and. size(refinement_map) == 0,&
            expected = .true.,&
            msg      = "Negative degree elevation must return defined empty outputs.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0150


    subroutine forcad_utils_0151(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        character(len=*), parameter :: vtk_file = "vtk/test_binary_endian.vtk"
        character(len=*), parameter :: points_header = &
            "# vtk DataFile Version 2.0"//achar(10)//&
            "Generated by ForCAD"//achar(10)//&
            "BINARY"//achar(10)//&
            "DATASET UNSTRUCTURED_GRID"//achar(10)//&
            "POINTS 2 double"//achar(10)
        character(len=*), parameter :: cells_header = achar(10)//"CELLS 1 3"//achar(10)
        character(len=*), parameter :: types_header = achar(10)//"CELL_TYPES 1"//achar(10)
        character(len=*), parameter :: field_header = &
            achar(10)//"POINT_DATA 2"//achar(10)//&
            "SCALARS temperature double"//achar(10)//&
            "LOOKUP_TABLE default"//achar(10)
        integer, parameter :: expected_points(48) = [&
            63, 240, 0, 0, 0, 0, 0, 0,&
            64,   0, 0, 0, 0, 0, 0, 0,&
            64,   8, 0, 0, 0, 0, 0, 0,&
            64,  16, 0, 0, 0, 0, 0, 0,&
            64,  20, 0, 0, 0, 0, 0, 0,&
            64,  24, 0, 0, 0, 0, 0, 0]
        integer, parameter :: expected_cells(12) = [&
            0, 0, 0, 2,&
            0, 0, 0, 0,&
            0, 0, 0, 1]
        integer, parameter :: expected_type(4) = [0, 0, 0, 3]
        integer, parameter :: expected_field(16) = [&
            64, 28, 0, 0, 0, 0, 0, 0,&
            64, 32, 0, 0, 0, 0, 0, 0]
        character(len=len(points_header)) :: actual_points_header
        character(len=len(cells_header)) :: actual_cells_header
        character(len=len(types_header)) :: actual_types_header
        character(len=len(field_header)) :: actual_field_header
        character :: actual_line_feed
        character(len=11) :: field_names(1)
        integer(int8) :: actual_points(size(expected_points))
        integer(int8) :: actual_cells(size(expected_cells))
        integer(int8) :: actual_type(size(expected_type))
        integer(int8) :: actual_field(size(expected_field))
        real(rk) :: points(2,3), point_data(2,1)
        integer :: bit, file_unit, i, io_status
        integer :: elem_conn(1,2)
        logical :: exists, valid

        points(1,1) = 1.0_rk
        points(1,2) = 2.0_rk
        points(1,3) = 3.0_rk
        points(2,1) = 4.0_rk
        points(2,2) = 5.0_rk
        points(2,3) = 6.0_rk
        elem_conn(1,1) = 1
        elem_conn(1,2) = 2
        point_data(1,1) = 7.0_rk
        point_data(2,1) = 8.0_rk
        field_names(1) = "temperature"

        open(&
            newunit = file_unit,&
            file    = vtk_file,&
            access  = "stream",&
            form    = "unformatted",&
            status  = "old",&
            iostat  = io_status)
        if (io_status == 0) close(file_unit, status="delete")

        call export_vtk_legacy(&
            filename    = vtk_file,&
            points      = points,&
            elemConn    = elem_conn,&
            vtkCellType = 3,&
            point_data  = point_data,&
            field_names = field_names,&
            encoding    = "BINARY")

        actual_points_header = ""
        actual_cells_header = ""
        actual_types_header = ""
        actual_field_header = ""
        actual_line_feed = achar(0)
        actual_points = 0_int8
        actual_cells = 0_int8
        actual_type = 0_int8
        actual_field = 0_int8
        io_status = 1
        inquire(file=vtk_file, exist=exists)
        valid = exists .and. file_storage_size == 8 .and. bit_size(0_int8) == 8
        if (valid) then
            open(&
                newunit = file_unit,&
                file    = vtk_file,&
                access  = "stream",&
                form    = "unformatted",&
                action  = "read",&
                status  = "old",&
                iostat  = io_status)
            if (io_status == 0) then
                read(file_unit, iostat=io_status) actual_points_header
                if (io_status == 0) read(file_unit, iostat=io_status) actual_points
                if (io_status == 0) read(file_unit, iostat=io_status) actual_cells_header
                if (io_status == 0) read(file_unit, iostat=io_status) actual_cells
                if (io_status == 0) read(file_unit, iostat=io_status) actual_types_header
                if (io_status == 0) read(file_unit, iostat=io_status) actual_type
                if (io_status == 0) read(file_unit, iostat=io_status) actual_field_header
                if (io_status == 0) read(file_unit, iostat=io_status) actual_field
                if (io_status == 0) read(file_unit, iostat=io_status) actual_line_feed
                close(file_unit)
            end if
        end if

        valid = valid .and. io_status == 0 .and. actual_points_header == points_header .and. &
            actual_cells_header == cells_header .and. actual_types_header == types_header .and. &
            actual_field_header == field_header .and. actual_line_feed == achar(10)
        do i = 1, size(expected_points)
            do bit = 0, 7
                valid = valid .and. (btest(actual_points(i),bit) .eqv. btest(expected_points(i),bit))
            end do
        end do
        do i = 1, size(expected_cells)
            do bit = 0, 7
                valid = valid .and. (btest(actual_cells(i),bit) .eqv. btest(expected_cells(i),bit))
            end do
        end do
        do i = 1, size(expected_type)
            do bit = 0, 7
                valid = valid .and. (btest(actual_type(i),bit) .eqv. btest(expected_type(i),bit))
            end do
        end do
        do i = 1, size(expected_field)
            do bit = 0, 7
                valid = valid .and. (btest(actual_field(i),bit) .eqv. btest(expected_field(i),bit))
            end do
        end do

        call ut%test(ti)%check(&
            name     = "binary VTK big-endian representation",&
            res      = valid,&
            expected = .true.,&
            msg      = "Binary VTK big-endian representation is incorrect.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0151


    subroutine forcad_utils_0152(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(8) = [&
            0.0_rk,0.0_rk,0.0_rk,0.5_rk,0.5_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Pw(5,1) = reshape([&
            0.0_rk,0.0_rk,1.0_rk,0.0_rk,0.0_rk], [5,1])
        real(rk), allocatable :: knot_new(:), Pw_new(:,:)
        logical :: preserved
        integer :: t

        call remove_knots_A_5_8(&
            p        = 2,&
            knot     = knot,&
            Pw       = Pw,&
            u        = 0.5_rk,&
            r        = 5,&
            s        = 2,&
            num      = 1,&
            t        = t,&
            knot_new = knot_new,&
            Pw_new   = Pw_new)

        preserved = t == 0 .and. size(knot_new) == size(knot) .and. &
            size(Pw_new,1) == size(Pw,1) .and. size(Pw_new,2) == size(Pw,2)
        if (preserved) preserved = maxval(abs(knot_new-knot)) <= epsilon(1.0_rk) .and. &
            maxval(abs(Pw_new-Pw)) <= epsilon(1.0_rk)

        call ut%test(ti)%check(&
            name     = "A5.8 rejects geometry-changing removal",&
            res      = preserved,&
            expected = .true.,&
            msg      = "A5.8 accepted a removal that changes the curve.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0152


    subroutine forcad_utils_0153(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(9) = [&
            -2.2_rk,-1.6_rk,0.0_rk,0.7_rk,1.8_rk,2.4_rk,4.0_rk,4.7_rk,5.8_rk]
        real(rk), parameter :: Xc(6,2) = reshape([&
            1.0_rk,0.0_rk,-1.0_rk,0.0_rk,1.0_rk,0.0_rk,&
            0.0_rk,1.0_rk,0.0_rk,-1.0_rk,0.0_rk,1.0_rk], [6,2])
        real(rk), parameter :: Wc(6) = [1.0_rk,1.4_rk,0.8_rk,1.2_rk,1.0_rk,1.4_rk]

        call ut%test(ti)%check(&
            name     = "verified periodic topology",&
            res      = periodic_topology(knot, 2, Xc, Wc),&
            expected = .true.,&
            msg      = "A complete cyclic representation was not recognized.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0153


    subroutine forcad_utils_0154(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk) :: knot(9)
        real(rk), parameter :: Xc(6,2) = reshape([&
            1.0_rk,0.0_rk,-1.0_rk,0.0_rk,1.0_rk,0.0_rk,&
            0.0_rk,1.0_rk,0.0_rk,-1.0_rk,0.0_rk,1.0_rk], [6,2])

        knot = [-2.2_rk,-1.6_rk,0.0_rk,0.7_rk,1.8_rk,2.4_rk,4.0_rk,4.7_rk,5.9_rk]

        call ut%test(ti)%check(&
            name     = "periodic knot extension rejection",&
            res      = periodic_topology(knot, 2, Xc),&
            expected = .false.,&
            msg      = "Inconsistent cyclic knot spacing was accepted.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0154


    subroutine forcad_utils_0155(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(9) = [&
            -2.2_rk,-1.6_rk,0.0_rk,0.7_rk,1.8_rk,2.4_rk,4.0_rk,4.7_rk,5.8_rk]
        real(rk) :: Xc(6,2)

        Xc = reshape([&
            1.0_rk,0.0_rk,-1.0_rk,0.0_rk,1.0_rk,0.0_rk,&
            0.0_rk,1.0_rk,0.0_rk,-1.0_rk,0.0_rk,1.25_rk], [6,2])

        call ut%test(ti)%check(&
            name     = "periodic control closure rejection",&
            res      = periodic_topology(knot, 2, Xc),&
            expected = .false.,&
            msg      = "Mismatched repeated controls were accepted.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0155


    subroutine forcad_utils_0156(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(9) = [&
            -2.2_rk,-1.6_rk,0.0_rk,0.7_rk,1.8_rk,2.4_rk,4.0_rk,4.7_rk,5.8_rk]
        real(rk), parameter :: Xc(6,2) = reshape([&
            1.0_rk,0.0_rk,-1.0_rk,0.0_rk,1.0_rk,0.0_rk,&
            0.0_rk,1.0_rk,0.0_rk,-1.0_rk,0.0_rk,1.0_rk], [6,2])
        real(rk) :: Wc(6)

        Wc = [1.0_rk,1.4_rk,0.8_rk,1.2_rk,1.0_rk,1.5_rk]

        call ut%test(ti)%check(&
            name     = "periodic weight closure rejection",&
            res      = periodic_topology(knot, 2, Xc, Wc),&
            expected = .false.,&
            msg      = "Mismatched repeated weights were accepted.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0156


    subroutine forcad_utils_0157(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: A(2,2) = reshape([&
            1.0e-20_rk,0.0_rk,1.0e-20_rk,1.0e-20_rk], [2,2])
        real(rk), parameter :: B(2,1) = reshape([&
            2.0e-20_rk,1.0e-20_rk], [2,1])
        real(rk), parameter :: expected(2,1) = reshape([&
            1.0_rk,1.0_rk], [2,1])
        real(rk), allocatable :: X(:,:)

        X = solve(A, B)

        call ut%test(ti)%check(&
            name     = "small nonsymmetric dense solve",&
            res      = X,&
            expected = expected,&
            tol      = 64.0_rk*epsilon(1.0_rk),&
            msg      = "A small nonsymmetric system must use pivoted LU.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0157


    subroutine forcad_utils_0158(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: A(2,2) = reshape([&
            2.0e-20_rk,1.0e-20_rk,1.0e-20_rk,3.0e-20_rk], [2,2])
        real(rk), parameter :: B(2,1) = reshape([&
            0.0_rk,-5.0e-20_rk], [2,1])
        real(rk), parameter :: expected(2,1) = reshape([&
            1.0_rk,-2.0_rk], [2,1])
        real(rk), allocatable :: X(:,:)

        X = solve(A, B)

        call ut%test(ti)%check(&
            name     = "small SPD dense solve",&
            res      = X,&
            expected = expected,&
            tol      = 64.0_rk*epsilon(1.0_rk),&
            msg      = "Uniform scaling must preserve the SPD solution.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0158


    subroutine forcad_utils_0159(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: A(2,2) = reshape([&
            1.0_rk,0.0_rk,1.0_rk,1.0e-20_rk], [2,2])
        real(rk), parameter :: B(2,1) = reshape([&
            2.0_rk,1.0e-20_rk], [2,1])
        real(rk), parameter :: expected(2,1) = reshape([&
            1.0_rk,1.0_rk], [2,1])
        real(rk), allocatable :: X(:,:)

        X = solve(A, B)

        call ut%test(ti)%check(&
            name     = "scale-separated LU pivots",&
            res      = X,&
            expected = expected,&
            tol      = 64.0_rk*epsilon(1.0_rk),&
            msg      = "A nonzero LU pivot must not be rejected by scale.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0159


    subroutine forcad_utils_0160(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: A(2,2) = reshape([&
            1.0e-8_rk,0.0_rk,0.0_rk,1.0e-8_rk], [2,2])
        real(rk), parameter :: expected(2,2) = reshape([&
            1.0e8_rk,0.0_rk,0.0_rk,1.0e8_rk], [2,2])
        real(rk), allocatable :: A_inv(:,:)

        A_inv = inv(A)

        call ut%test(ti)%check(&
            name     = "small-scale inverse",&
            res      = A_inv,&
            expected = expected,&
            tol      = 64.0_rk*epsilon(1.0_rk)*1.0e8_rk,&
            msg      = "A well-conditioned inverse must be scale invariant.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0160


    subroutine forcad_utils_0161(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: A(3,3) = reshape([&
            1.0e120_rk,0.0_rk,0.0_rk,&
            0.0_rk,1.0e120_rk,0.0_rk,&
            0.0_rk,0.0_rk,1.0e120_rk], [3,3])
        real(rk), parameter :: expected(3,3) = reshape([&
            1.0e-120_rk,0.0_rk,0.0_rk,&
            0.0_rk,1.0e-120_rk,0.0_rk,&
            0.0_rk,0.0_rk,1.0e-120_rk], [3,3])
        real(rk), allocatable :: A_inv(:,:)

        A_inv = inv(A)

        call ut%test(ti)%check(&
            name     = "large-scale inverse",&
            res      = A_inv,&
            expected = expected,&
            tol      = 64.0_rk*epsilon(1.0_rk)*1.0e-120_rk,&
            msg      = "Inverse scaling must avoid determinant overflow.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0161


    subroutine forcad_utils_0162(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: A(2,2) = reshape([&
            1.0_rk,0.0_rk,0.0_rk,epsilon(1.0_rk)], [2,2])
        real(rk), allocatable :: A_inv(:,:)

        A_inv = inv(A)

        call ut%test(ti)%check(&
            name     = "ill-conditioned inverse rejection",&
            res      = shape(A_inv),&
            expected = [0,0],&
            msg      = "A numerically singular inverse must be rejected.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0162


    subroutine forcad_utils_0163(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: A_band(0:1,3) = reshape([&
            4.0e-20_rk,1.0e-20_rk,&
            3.0e-20_rk,1.0e-20_rk,&
            2.0e-20_rk,0.0_rk], [2,3])
        real(rk), parameter :: B(3,1) = reshape([&
            6.0e-20_rk,1.0e-19_rk,8.0e-20_rk], [3,1])
        real(rk), parameter :: expected(3,1) = reshape([&
            1.0_rk,2.0_rk,3.0_rk], [3,1])
        real(rk), allocatable :: X(:,:)

        X = solve_spd_banded(A_band, B)

        call ut%test(ti)%check(&
            name     = "small-scale banded SPD solve",&
            res      = X,&
            expected = expected,&
            tol      = 128.0_rk*epsilon(1.0_rk),&
            msg      = "Banded Cholesky must be invariant under scaling.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0163


    subroutine forcad_utils_0164(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: A_band(0:1,2) = reshape([&
            1.0e-20_rk,1.0e-20_rk,&
            1.0e-20_rk,0.0_rk], [2,2])
        real(rk), parameter :: B(2,1) = reshape([&
            2.0e-20_rk,2.0e-20_rk], [2,1])
        real(rk), allocatable :: X(:,:)

        X = solve_spd_banded(A_band, B)

        call ut%test(ti)%check(&
            name     = "singular banded SPD rejection",&
            res      = shape(X),&
            expected = [0,0],&
            msg      = "A rank-deficient band matrix must be rejected.",&
            group    = "forcad_utils")
        ti = ti + 1

    end subroutine forcad_utils_0164


    subroutine run_utils_tests(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call forcad_utils_0001(ut, ti)
        call forcad_utils_0002(ut, ti)
        call forcad_utils_0003(ut, ti)
        call forcad_utils_0004(ut, ti)
        call forcad_utils_0005(ut, ti)
        call forcad_utils_0006(ut, ti)
        call forcad_utils_0007(ut, ti)
        call forcad_utils_0008(ut, ti)
        call forcad_utils_0009(ut, ti)
        call forcad_utils_0010(ut, ti)
        call forcad_utils_0011(ut, ti)
        call forcad_utils_0012(ut, ti)
        call forcad_utils_0013(ut, ti)
        call forcad_utils_0014(ut, ti)
        call forcad_utils_0015(ut, ti)
        call forcad_utils_0016(ut, ti)
        call forcad_utils_0017(ut, ti)
        call forcad_utils_0018(ut, ti)
        call forcad_utils_0019(ut, ti)
        call forcad_utils_0020(ut, ti)
        call forcad_utils_0021(ut, ti)
        call forcad_utils_0022(ut, ti)
        call forcad_utils_0023(ut, ti)
        call forcad_utils_0024(ut, ti)
        call forcad_utils_0025(ut, ti)
        call forcad_utils_0026(ut, ti)
        call forcad_utils_0027(ut, ti)
        call forcad_utils_0028(ut, ti)
        call forcad_utils_0029(ut, ti)
        call forcad_utils_0030(ut, ti)
        call forcad_utils_0031(ut, ti)
        call forcad_utils_0032(ut, ti)
        call forcad_utils_0033(ut, ti)
        call forcad_utils_0034(ut, ti)
        call forcad_utils_0035(ut, ti)
        call forcad_utils_0036(ut, ti)
        call forcad_utils_0037(ut, ti)
        call forcad_utils_0038(ut, ti)
        call forcad_utils_0039(ut, ti)
        call forcad_utils_0040(ut, ti)
        call forcad_utils_0041(ut, ti)
        call forcad_utils_0042(ut, ti)
        call forcad_utils_0043(ut, ti)
        call forcad_utils_0044(ut, ti)
        call forcad_utils_0045(ut, ti)
        call forcad_utils_0046(ut, ti)
        call forcad_utils_0047(ut, ti)
        call forcad_utils_0048(ut, ti)
        call forcad_utils_0049(ut, ti)
        call forcad_utils_0050(ut, ti)
        call forcad_utils_0051(ut, ti)
        call forcad_utils_0052(ut, ti)
        call forcad_utils_0053(ut, ti)
        call forcad_utils_0054(ut, ti)
        call forcad_utils_0055(ut, ti)
        call forcad_utils_0056(ut, ti)
        call forcad_utils_0057(ut, ti)
        call forcad_utils_0058(ut, ti)
        call forcad_utils_0059(ut, ti)
        call forcad_utils_0060(ut, ti)
        call forcad_utils_0061(ut, ti)
        call forcad_utils_0062(ut, ti)
        call forcad_utils_0063(ut, ti)
        call forcad_utils_0064(ut, ti)
        call forcad_utils_0065(ut, ti)
        call forcad_utils_0066(ut, ti)
        call forcad_utils_0067(ut, ti)
        call forcad_utils_0068(ut, ti)
        call forcad_utils_0069(ut, ti)
        call forcad_utils_0070(ut, ti)
        call forcad_utils_0071(ut, ti)
        call forcad_utils_0072(ut, ti)
        call forcad_utils_0073(ut, ti)
        call forcad_utils_0074(ut, ti)
        call forcad_utils_0075(ut, ti)
        call forcad_utils_0076(ut, ti)
        call forcad_utils_0077(ut, ti)
        call forcad_utils_0078(ut, ti)
        call forcad_utils_0079(ut, ti)
        call forcad_utils_0080(ut, ti)
        call forcad_utils_0081(ut, ti)
        call forcad_utils_0082(ut, ti)
        call forcad_utils_0083(ut, ti)
        call forcad_utils_0084(ut, ti)
        call forcad_utils_0085(ut, ti)
        call forcad_utils_0086(ut, ti)
        call forcad_utils_0087(ut, ti)
        call forcad_utils_0088(ut, ti)
        call forcad_utils_0089(ut, ti)
        call forcad_utils_0090(ut, ti)
        call forcad_utils_0091(ut, ti)
        call forcad_utils_0092(ut, ti)
        call forcad_utils_0093(ut, ti)
        call forcad_utils_0094(ut, ti)
        call forcad_utils_0095(ut, ti)
        call forcad_utils_0096(ut, ti)
        call forcad_utils_0097(ut, ti)
        call forcad_utils_0098(ut, ti)
        call forcad_utils_0099(ut, ti)
        call forcad_utils_0100(ut, ti)
        call forcad_utils_0101(ut, ti)
        call forcad_utils_0102(ut, ti)
        call forcad_utils_0103(ut, ti)
        call forcad_utils_0104(ut, ti)
        call forcad_utils_0105(ut, ti)
        call forcad_utils_0106(ut, ti)
        call forcad_utils_0107(ut, ti)
        call forcad_utils_0108(ut, ti)
        call forcad_utils_0109(ut, ti)
        call forcad_utils_0110(ut, ti)
        call forcad_utils_0111(ut, ti)
        call forcad_utils_0112(ut, ti)
        call forcad_utils_0113(ut, ti)
        call forcad_utils_0114(ut, ti)
        call forcad_utils_0115(ut, ti)
        call forcad_utils_0116(ut, ti)
        call forcad_utils_0117(ut, ti)
        call forcad_utils_0118(ut, ti)
        call forcad_utils_0119(ut, ti)
        call forcad_utils_0120(ut, ti)
        call forcad_utils_0121(ut, ti)
        call forcad_utils_0122(ut, ti)
        call forcad_utils_0123(ut, ti)
        call forcad_utils_0124(ut, ti)
        call forcad_utils_0125(ut, ti)
        call forcad_utils_0126(ut, ti)
        call forcad_utils_0127(ut, ti)
        call forcad_utils_0128(ut, ti)
        call forcad_utils_0129(ut, ti)
        call forcad_utils_0130(ut, ti)
        call forcad_utils_0131(ut, ti)
        call forcad_utils_0132(ut, ti)
        call forcad_utils_0133(ut, ti)
        call forcad_utils_0134(ut, ti)
        call forcad_utils_0135(ut, ti)
        call forcad_utils_0136(ut, ti)
        call forcad_utils_0137(ut, ti)
        call forcad_utils_0138(ut, ti)
        call forcad_utils_0139(ut, ti)
        call forcad_utils_0140(ut, ti)
        call forcad_utils_0141(ut, ti)
        call forcad_utils_0142(ut, ti)
        call forcad_utils_0143(ut, ti)
        call forcad_utils_0144(ut, ti)
        call forcad_utils_0145(ut, ti)
        call forcad_utils_0146(ut, ti)
        call forcad_utils_0147(ut, ti)
        call forcad_utils_0148(ut, ti)
        call forcad_utils_0149(ut, ti)
        call forcad_utils_0150(ut, ti)
        call forcad_utils_0151(ut, ti)
        call forcad_utils_0152(ut, ti)
        call forcad_utils_0153(ut, ti)
        call forcad_utils_0154(ut, ti)
        call forcad_utils_0155(ut, ti)
        call forcad_utils_0156(ut, ti)
        call forcad_utils_0157(ut, ti)
        call forcad_utils_0158(ut, ti)
        call forcad_utils_0159(ut, ti)
        call forcad_utils_0160(ut, ti)
        call forcad_utils_0161(ut, ti)
        call forcad_utils_0162(ut, ti)
        call forcad_utils_0163(ut, ti)
        call forcad_utils_0164(ut, ti)

    end subroutine run_utils_tests

end module forcad_test_utils
