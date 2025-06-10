module bspline_basis

    use bspline_kinds, only: rk

    implicit none

    private
    public :: &
        bspline_basis1, &
        bspline_basis2, &
        bspline_basis3, &
        bspline_basis4

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function bspline_basis1(Xt, knot, nc, degree) result(B)
        integer, intent(in) :: degree
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc
        real(rk), intent(in) :: Xt
        real(rk) :: B(nc), B_curr(nc)
        integer :: i, p
        real(rk) :: knot_i, knot_i1, knot_ip, knot_ip1, knot_last

        B = 0.0_rk
        B_curr = 0.0_rk

        ! Degree 0 initialization
        do concurrent(i=1:nc)
            knot_i = knot(i)
            knot_i1 = knot(i+1)
            knot_last = knot(size(knot))

            if ((Xt >= knot_i .and. Xt < knot_i1) .or. (Xt == knot_last .and. Xt == knot_i1)) B(i) = 1.0_rk
        end do

        ! Recursion for higher degrees
        do p = 1, degree
            B_curr = 0.0_rk
            do concurrent(i=1:nc)
                knot_i = knot(i)
                knot_i1 = knot(i+1)
                knot_ip = knot(i+p)
                knot_ip1 = knot(i+p+1)

                if (knot_ip /= knot_i) B_curr(i) = (Xt-knot_i)/(knot_ip-knot_i)*B(i)
                if (i < nc .and. knot_ip1 /= knot_i1) B_curr(i) = B_curr(i)+(knot_ip1-Xt)/(knot_ip1-knot_i1)*B(i+1)
            end do
            B = B_curr
        end do
    end function
    !===============================================================================

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function bspline_basis2(Xt, knot, nc, degree) result(B)
        integer, intent(in) :: degree
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc
        real(rk), intent(in) :: Xt
        real(rk) :: temp, knot_i, knot_i1, knot_ip, knot_ip1
        real(rk) :: N(nc, 0:degree)
        integer :: i, p
        real(rk) :: B(nc)

        temp = abs(Xt-knot(size(knot)))
        N = 0.0_rk

        do p = 0, degree
            do concurrent(i=1:nc)
                knot_i = knot(i)
                knot_i1 = knot(i+1)
                knot_ip = knot(i+p)
                knot_ip1 = knot(i+p+1)

                if (temp /= tiny(0.0_rk) .and. Xt >= knot_i .and. Xt <= knot_i1) N(i, 0) = 1.0_rk
                if (knot_ip /= knot_i) N(i, p) = (Xt-knot_i)/(knot_ip-knot_i)*N(i, p-1)
                if (i < nc .and. knot_ip1 /= knot_i1) N(i, p) = N(i, p)+(knot_ip1-Xt)/(knot_ip1-knot_i1)*N(i+1, p-1)
            end do
        end do
        B = N(:, degree)
    end function
    !===============================================================================

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function bspline_basis3(Xt, knot, nc, degree) result(B)
        integer, intent(in) :: degree
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc
        real(rk), intent(in) :: Xt
        real(rk) :: knot_i, knot_i1, knot_ip, knot_ip1, knot_last
        real(rk) :: N(nc, 0:degree)
        integer :: i, p
        real(rk) :: B(nc)

        N = 0.0_rk

        ! Initialize basis functions for degree 0
        do concurrent(i=1:nc)
            knot_i = knot(i)
            knot_i1 = knot(i+1)
            knot_last = knot(size(knot))

            if ((Xt >= knot_i .and. Xt < knot_i1) .or. (Xt == knot_last .and. Xt == knot_i1)) N(i, 0) = 1.0_rk
        end do

        do p = 1, degree
            do concurrent(i=1:nc)
                knot_i = knot(i)
                knot_i1 = knot(i+1)
                knot_ip = knot(i+p)
                knot_ip1 = knot(i+p+1)

                if (knot_ip /= knot_i) N(i, p) = (Xt-knot_i)/(knot_ip-knot_i)*N(i, p-1)
                if (i < nc .and. knot_ip1 /= knot_i1) N(i, p) = N(i, p)+(knot_ip1-Xt)/(knot_ip1-knot_i1)*N(i+1, p-1)
            end do
        end do
        B = N(:, degree)
    end function
    !===============================================================================

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function bspline_basis4(Xt, knot, nc, degree) result(B)
        integer, intent(in) :: degree
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc
        real(rk), intent(in) :: Xt
        real(rk) :: B(nc)
        integer :: span, j, r, low, mid, high, nk
        real(rk) :: left(degree), right(degree)
        real(rk) :: N(0:degree)
        real(rk) :: saved, temp1, temp2, l1, r1
        integer :: i

        if (nc == 0) then
            B = 0.0_rk
            return
        end if

        nk = size(knot)
        B = 0.0_rk

        ! Check if Xt is outside the domain
        if (Xt < knot(1) .or. Xt > knot(nk)) error stop "Xt is outside the knot range"

        ! Find knot span
        if (Xt == knot(nk)) then
            span = nk-1
        else
            low = 1
            high = nk
            do while (low < high)
                mid = (low+high)/2
                if (Xt < knot(mid)) then
                    high = mid
                else
                    low = mid+1
                end if
            end do
            span = low-1
        end if

        ! Degree 0 case (TODO: check this)
        if (degree == 0) then
            if (span >= 1 .and. span <= nc) B(span) = 1.0_rk
            return
        end if

        ! Precompute differences for recurrence
        do concurrent(j=1:degree)
            left(j) = Xt-knot(span-j+1)
            right(j) = knot(span+j)-Xt
        end do

        ! Recurrence for higher degrees
        N(0) = 1.0_rk
        do j = 1, degree
            saved = 0.0_rk
            do r = 0, j-1
                l1 = left(j-r)
                r1 = right(r+1)
                temp2 = r1+l1
                if (abs(temp2) <= tiny(1.0_rk)) then
                    temp1 = 0.0_rk
                else
                    temp1 = N(r)/temp2
                end if
                N(r) = saved+r1*temp1
                saved = l1*temp1
            end do
            N(j) = saved
        end do

        do concurrent(i=max(1, span-degree):min(nc, span))
            B(i) = N(i-(span-degree))
        end do
    end function
    !===============================================================================

end module
