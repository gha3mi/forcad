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
        real(rk) :: saved, temp
        integer :: i, index_start

        if (nc == 0) then
            B = 0.0_rk
            return
        end if

        B = 0.0_rk
        nk = size(knot)

        if (Xt < knot(1) .or. Xt > knot(nk)) then
            B = 0.0_rk
            return
        end if

        ! Find span
        if (Xt == knot(nk)) then
            span = nk-degree-1
        else
            low = degree+1
            high = nk-degree
            do while (low <= high)
                mid = (low+high)/2
                if (Xt >= knot(mid) .and. Xt < knot(mid+1)) then
                    span = mid
                    exit
                else if (Xt < knot(mid)) then
                    high = mid-1
                else
                    low = mid+1
                end if
            end do
        end if

        ! Cox-de Boor recursion
        N = 0.0_rk
        N(0) = 1.0_rk
        do j = 1, degree
            left(j) = Xt-knot(span+1-j)
            right(j) = knot(span+j)-Xt
            saved = 0.0_rk
            do r = 0, j-1
                temp = N(r)/(right(r+1)+left(j-r))
                N(r) = saved+right(r+1)*temp
                saved = left(j-r)*temp
            end do
            N(j) = saved
        end do

        index_start = span-degree
        do i = 0, degree
            if (index_start+i >= 1 .and. index_start+i <= nc) then
                B(index_start+i) = N(i)
            end if
        end do
    end function
    !===============================================================================

end module
