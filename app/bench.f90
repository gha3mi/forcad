program benchmark_bspline

    use bspline_kinds, only: rk
    use bspline_basis
    use bspline_utils, only: print_results, write_results_csv
    use fortime, only: timer
    use iso_fortran_env, only: compiler_version, compiler_options

    implicit none

    !! Setup parameters for the benchmark
    !------------------------------------------------------------------
    ! Names of the methods being benchmarked
    character(*), parameter :: method_names(*) = ["b1", "b2", "b3"]
    ! Minimum degree of B-spline basis functions
    integer, parameter :: degree_min = 1
    ! Maximum degree of B-spline basis functions
    integer, parameter :: degree_max = 5
    ! Minimum number of control points
    integer, parameter :: nc_min = 0
    ! Maximum number of control points
    integer, parameter :: nc_max = 1e6
    ! Step size for control points
    integer, parameter :: nc_step = nc_max/100
    ! Number of repetitions for benchmarking
    integer, parameter :: reps = 20

    integer, parameter :: nmethods = size(method_names)
    integer :: num_knots, i, method, nc, degree
    real(rk), allocatable :: knot(:), B(:)
    type(timer) :: t((nc_max-nc_min)/nc_step+1, degree_max-degree_min+1, nmethods)
    real(rk) :: Xt

    call random_number(Xt) ! Generate a random point in [0, 1]

    print '(a)', '============================================='
    print '(a)', ' B-spline Basis Function Benchmark Report'
    print '(a)', '============================================='
    print '(a,i12)', ' Degree of B-spline Basis     : ', degree
    print '(a,i12)', ' Number of Repetitions        : ', reps
    print '(a,i12)', ' Number of Methods Benchmarked: ', nmethods
    print '(a)', '---------------------------------------------'
    print '(a)', ' Compiler Information:'
    print '(a)', '---------------------------------------------'
    print '(a)', trim(compiler_version())
    print '(a)', trim(compiler_options())
    print '(a)', '============================================='

    do degree = degree_min, degree_max
        do nc = nc_min, nc_max, nc_step
            print'(a,g0, ",",1x ,a,g0)', "degree=", degree, "nc=", nc
            num_knots = nc+degree+1
            allocate (knot(num_knots), source=[(real(i-1, rk)/(num_knots-1), i=1, num_knots)])
            allocate (B(nc))
            do method = 1, nmethods
                B = 0.0_rk

                ! Start the timer for the current method
                call t((nc-nc_min)/nc_step+1, degree, method)%timer_start()
                do i = 1, reps
                    select case (method)
                    case (1); B = bspline_basis1(Xt, knot, nc, degree)
                    case (2); B = bspline_basis2(Xt, knot, nc, degree)
                    case (3); B = bspline_basis3(Xt, knot, nc, degree)
                    end select
                    ! to prevent loop invariant optimization
                    write (*, '(a)', advance='no') "."
                end do
                ! Stop the timer for the current method
                call t((nc-nc_min)/nc_step+1, degree, method)%timer_stop(nloops=reps, message=method_names(method), print=.true.)
            end do
            deallocate (knot, B)
        end do
    end do

    ! Print the timing results
    call print_results(degree_min, degree_max, nc_min, nc_max, nc_step, reps, nmethods, method_names, t)
    call write_results_csv(degree_min, degree_max, nc_min, nc_max, nc_step, reps, nmethods, method_names, t)
end program
