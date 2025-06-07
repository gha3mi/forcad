module bspline_utils

    use bspline_kinds, only: rk
    use fortime, only: timer

    implicit none

    private
    public :: print_results

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    subroutine print_results(degree_min, degree_max, nc_min, nc_max, nc_step, reps, nmethods, method_names, t)
        integer, intent(in) :: degree_min, degree_max, nc_min, nc_max, nc_step, reps, nmethods
        character(*), dimension(nmethods), intent(in) :: method_names
        integer :: degree, nc, method
        type(timer), intent(in) :: t((nc_max-nc_min)/nc_step+1, degree_max-degree_min+1, nmethods)

        print '(a)', '---------------------------------------------'
        print '(a)', ' Timing Results:'
        print '(a)', '---------------------------------------------'
        print '(a)', ' Degree | Control Points | Method | Time (s)'
        print '(a)', '---------------------------------------------'
        do degree = degree_min, degree_max
            do nc = nc_min, nc_max, nc_step
                do method = 1, nmethods
                    write (*, '(i3, 2x, i12, 2x, a10, 2x, f10.8)') &
                        degree, nc, trim(method_names(method)), t((nc-nc_min)/nc_step+1, degree, method)%elapsed_time
                end do
            end do
        end do
    end subroutine print_results
    !===============================================================================

end module
