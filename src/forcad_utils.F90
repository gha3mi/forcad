!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> This module contains parameters, functions and subroutines that are used in the library.
module forcad_utils

    use forcad_kinds, only: rk

    implicit none

    private
    public basis_bernstein, basis_bspline, elemConn_C0, kron, ndgrid, compute_multiplicity, compute_knot_vector, &
        basis_bspline_der, insert_knot_A_5_1, findspan, elevate_degree_A_5_9, hexahedron_Xc, tetragon_Xc, remove_knots_A_5_8, &
        elemConn_Cn, unique, rotation, basis_bspline_2der, det, inv, dyad, gauss_leg, export_vtk_legacy, solve, &
        repelem, linspace, eye

    !===============================================================================
    interface elemConn_C0
        module procedure cmp_elemConn_C0_L
        module procedure cmp_elemConn_C0_S
        module procedure cmp_elemConn_C0_V
    end interface
    !===============================================================================


    !===============================================================================
    interface elemConn_Cn
        module procedure cmp_elemConn_Cn_L
        module procedure cmp_elemConn_Cn_S
        module procedure cmp_elemConn_Cn_V
    end interface
    !===============================================================================


    !===============================================================================
    interface ndgrid
        module procedure ndgrid2
        module procedure ndgrid3
    end interface
    !===============================================================================


    !===============================================================================
    interface compute_multiplicity
        module procedure compute_multiplicity1
        module procedure compute_multiplicity2
    end interface
    !===============================================================================


    !===============================================================================
    interface unique
        module procedure unique_integer
        module procedure unique_real
    end interface
    !===============================================================================


    !===============================================================================
    interface dyad
        module procedure dyad_t1_t1
    end interface
    !===============================================================================


    !===============================================================================
    interface gauss_leg
        module procedure gauss_legendre_1D
        module procedure gauss_legendre_2D
        module procedure gauss_legendre_3D
    end interface
    !===============================================================================


    !===============================================================================
    interface kron
        module procedure kron_t1_t1
        module procedure kron_t1_t2
        module procedure kron3
    end interface
    !===============================================================================


    !===============================================================================
    interface basis_bspline_der
        module procedure basis_bspline_der_A
        module procedure basis_bspline_der_B
    end interface
    !===============================================================================


    !===============================================================================
    interface basis_bspline_2der
        module procedure basis_bspline_2der_A
        module procedure basis_bspline_2der_B
        module procedure basis_bspline_2der_C
    end interface
    !===============================================================================

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function basis_bspline(Xt, knot, nc, degree) result(B)
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
        ! if (Xt == knot(nk)) then
        if (abs(Xt - knot(nk)) < 2.0_rk*epsilon(0.0_rk)) then
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


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine basis_bspline_der_A(Xt, knot, nc, degree, dB, B)
        integer, intent(in) :: degree
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc
        real(rk), intent(in) :: Xt
        real(rk), intent(out) :: dB(nc)
        real(rk), intent(out) :: B(nc)
        integer :: i, p
        real(rk) :: Xth_i, Xth_i1, Xth_ip, Xth_ip1, Xth_last
        real(rk) :: B_curr(nc)
        real(rk) :: dB_curr(nc)

        B = 0.0_rk
        dB = 0.0_rk

        ! Degree 0 initialization
        do concurrent(i = 1:nc)
            Xth_i    = knot(i)
            Xth_i1   = knot(i + 1)
            Xth_last = knot(size(knot))

            ! if ((Xt >= Xth_i .and. Xt < Xth_i1) .or. (Xt == Xth_last .and. Xt == Xth_i1)) then
            if ((Xt >= Xth_i .and. Xt < Xth_i1) .or. (abs(Xt - Xth_last) < 2.0_rk*epsilon(0.0_rk) .and. abs(Xt - Xth_i1) < 2.0_rk*epsilon(0.0_rk))) then
                B(i) = 1.0_rk
                dB(i) = 0.0_rk
            end if
        end do

        ! Recursion for higher degrees
        do p = 1, degree
            B_curr = 0.0_rk
            dB_curr = 0.0_rk
            do concurrent(i = 1:nc)
                Xth_i   = knot(i)
                Xth_i1  = knot(i+1)
                Xth_ip  = knot(i+p)
                Xth_ip1 = knot(i+p+1)

                ! if (Xth_ip /= Xth_i) then
                if (abs(Xth_ip - Xth_i) > 2.0_rk*epsilon(0.0_rk)) then
                    B_curr(i) = (Xt-Xth_i)/(Xth_ip-Xth_i)*B(i)
                    dB_curr(i) = B(i)/(Xth_ip-Xth_i)+(Xt-Xth_i)/(Xth_ip-Xth_i)*dB(i)
                end if
                ! if (i < nc .and. Xth_ip1 /= Xth_i1) then
                if (i < nc .and. abs(Xth_ip1 - Xth_i1) > 2.0_rk*epsilon(0.0_rk)) then
                    B_curr(i) = B_curr(i)+(Xth_ip1-Xt)/(Xth_ip1-Xth_i1)*B(i+1)
                    dB_curr(i) = dB_curr(i)-B(i+1)/(Xth_ip1-Xth_i1)+ &
                                (Xth_ip1-Xt)/(Xth_ip1-Xth_i1)*dB(i+1)
                end if
            end do
            B = B_curr
            dB = dB_curr
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine basis_bspline_der_B(Xt, knot, nc, degree, dB)
        integer, intent(in) :: degree
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc
        real(rk), intent(in) :: Xt
        real(rk), intent(out) :: dB(nc)
        real(rk) :: B(nc)
        integer :: i, p
        real(rk) :: Xth_i, Xth_i1, Xth_ip, Xth_ip1, Xth_last
        real(rk) :: B_curr(nc)
        real(rk) :: dB_curr(nc)

        B = 0.0_rk
        dB = 0.0_rk

        ! Degree 0 initialization
        do concurrent(i = 1:nc)
            Xth_i    = knot(i)
            Xth_i1   = knot(i + 1)
            Xth_last = knot(size(knot))

            ! if ((Xt >= Xth_i .and. Xt < Xth_i1) .or. (Xt == Xth_last .and. Xt == Xth_i1)) then
            if ((Xt >= Xth_i .and. Xt < Xth_i1) .or. (abs(Xt - Xth_last) < 2.0_rk*epsilon(0.0_rk) .and. abs(Xt - Xth_i1) < 2.0_rk*epsilon(0.0_rk))) then
                B(i) = 1.0_rk
                dB(i) = 0.0_rk
            end if
        end do

        ! Recursion for higher degrees
        do p = 1, degree
            B_curr = 0.0_rk
            dB_curr = 0.0_rk
            do concurrent(i = 1:nc)
                Xth_i   = knot(i)
                Xth_i1  = knot(i+1)
                Xth_ip  = knot(i+p)
                Xth_ip1 = knot(i+p+1)

                ! if (Xth_ip /= Xth_i) then
                if (abs(Xth_ip - Xth_i) > 2.0_rk*epsilon(0.0_rk)) then
                    B_curr(i) = (Xt-Xth_i)/(Xth_ip-Xth_i)*B(i)
                    dB_curr(i) = B(i)/(Xth_ip-Xth_i)+(Xt-Xth_i)/(Xth_ip-Xth_i)*dB(i)
                end if
                ! if (i < nc .and. Xth_ip1 /= Xth_i1) then
                if (i < nc .and. abs(Xth_ip1 - Xth_i1) > 2.0_rk*epsilon(0.0_rk)) then
                    B_curr(i) = B_curr(i)+(Xth_ip1-Xt)/(Xth_ip1-Xth_i1)*B(i+1)
                    dB_curr(i) = dB_curr(i)-B(i+1)/(Xth_ip1-Xth_i1)+ &
                                (Xth_ip1-Xt)/(Xth_ip1-Xth_i1)*dB(i+1)
                end if
            end do
            B = B_curr
            dB = dB_curr
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine basis_bspline_2der_A(Xt, knot, nc, degree, d2B, dB, B)
        integer, intent(in) :: degree
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc
        real(rk), intent(in) :: Xt
        real(rk), intent(out) :: d2B(nc)
        real(rk), intent(out) :: dB(nc)
        real(rk), intent(out) :: B(nc)
        integer :: i, p
        real(rk) :: Xth_i, Xth_i1, Xth_ip, Xth_ip1, Xth_last
        real(rk) :: B_curr(nc)
        real(rk) :: dB_curr(nc)
        real(rk) :: d2B_curr(nc)

        B = 0.0_rk
        dB = 0.0_rk
        d2B = 0.0_rk

        ! Degree 0 initialization
        do concurrent(i = 1:nc)
            Xth_i = knot(i)
            Xth_i1 = knot(i + 1)
            Xth_last = knot(size(knot))

            ! if ((Xt >= Xth_i .and. Xt < Xth_i1) .or. (Xt == Xth_last .and. Xt == Xth_i1)) then
            if ((Xt >= Xth_i .and. Xt < Xth_i1) .or. (abs(Xt - Xth_last) < 2.0_rk*epsilon(0.0_rk) .and. abs(Xt - Xth_i1) < 2.0_rk*epsilon(0.0_rk))) then
                B(i) = 1.0_rk
                dB(i) = 0.0_rk
                d2B(i) = 0.0_rk
            end if
        end do

        ! Recursion for higher degrees
        do p = 1, degree
            B_curr = 0.0_rk
            dB_curr = 0.0_rk
            d2B_curr = 0.0_rk
            do concurrent(i = 1:nc)
                Xth_i = knot(i)
                Xth_i1 = knot(i + 1)
                Xth_ip = knot(i + p)
                Xth_ip1 = knot(i + p + 1)

                ! if (Xth_ip /= Xth_i) then
                if (abs(Xth_ip - Xth_i) > 2.0_rk*epsilon(0.0_rk)) then
                    B_curr(i) = (Xt - Xth_i)/(Xth_ip - Xth_i)*B(i)
                    dB_curr(i) = B(i)/(Xth_ip - Xth_i) + (Xt - Xth_i)/(Xth_ip - Xth_i)*dB(i)
                    d2B_curr(i) = (2*dB(i) + (Xt - Xth_i)*d2B(i))/(Xth_ip - Xth_i)
                end if

                ! if (i < nc .and. Xth_ip1 /= Xth_i1) then
                if (i < nc .and. abs(Xth_ip1 - Xth_i1) > 2.0_rk*epsilon(0.0_rk)) then
                    B_curr(i) = B_curr(i) + (Xth_ip1 - Xt)/(Xth_ip1 - Xth_i1)*B(i + 1)
                    dB_curr(i) = dB_curr(i) - B(i + 1)/(Xth_ip1 - Xth_i1) + (Xth_ip1 - Xt)/(Xth_ip1 - Xth_i1)*dB(i + 1)
                    d2B_curr(i) = d2B_curr(i) - (2*dB(i + 1) - (Xth_ip1 - Xt)*d2B(i + 1))/(Xth_ip1 - Xth_i1)
                end if
            end do
            B = B_curr
            dB = dB_curr
            d2B = d2B_curr
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine basis_bspline_2der_B(Xt, knot, nc, degree, d2B, dB)
        integer, intent(in) :: degree
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc
        real(rk), intent(in) :: Xt
        real(rk), intent(out) :: d2B(nc)
        real(rk), intent(out) :: dB(nc)
        real(rk) :: B(nc)
        integer :: i, p
        real(rk) :: Xth_i, Xth_i1, Xth_ip, Xth_ip1, Xth_last
        real(rk) :: B_curr(nc)
        real(rk) :: dB_curr(nc)
        real(rk) :: d2B_curr(nc)

        B = 0.0_rk
        dB = 0.0_rk
        d2B = 0.0_rk

        ! Degree 0 initialization
        do concurrent(i = 1:nc)
            Xth_i = knot(i)
            Xth_i1 = knot(i + 1)
            Xth_last = knot(size(knot))

            ! if ((Xt >= Xth_i .and. Xt < Xth_i1) .or. (Xt == Xth_last .and. Xt == Xth_i1)) then
            if ((Xt >= Xth_i .and. Xt < Xth_i1) .or. (abs(Xt - Xth_last) < 2.0_rk*epsilon(0.0_rk) .and. abs(Xt - Xth_i1) < 2.0_rk*epsilon(0.0_rk))) then
                B(i) = 1.0_rk
                dB(i) = 0.0_rk
                d2B(i) = 0.0_rk
            end if
        end do

        ! Recursion for higher degrees
        do p = 1, degree
            B_curr = 0.0_rk
            dB_curr = 0.0_rk
            d2B_curr = 0.0_rk
            do concurrent(i = 1:nc)
                Xth_i = knot(i)
                Xth_i1 = knot(i + 1)
                Xth_ip = knot(i + p)
                Xth_ip1 = knot(i + p + 1)

                ! if (Xth_ip /= Xth_i) then
                if (abs(Xth_ip - Xth_i) > 2.0_rk*epsilon(0.0_rk)) then
                    B_curr(i) = (Xt - Xth_i)/(Xth_ip - Xth_i)*B(i)
                    dB_curr(i) = B(i)/(Xth_ip - Xth_i) + (Xt - Xth_i)/(Xth_ip - Xth_i)*dB(i)
                    d2B_curr(i) = (2*dB(i) + (Xt - Xth_i)*d2B(i))/(Xth_ip - Xth_i)
                end if

                ! if (i < nc .and. Xth_ip1 /= Xth_i1) then
                if (i < nc .and. abs(Xth_ip1 - Xth_i1) > 2.0_rk*epsilon(0.0_rk)) then
                    B_curr(i) = B_curr(i) + (Xth_ip1 - Xt)/(Xth_ip1 - Xth_i1)*B(i + 1)
                    dB_curr(i) = dB_curr(i) - B(i + 1)/(Xth_ip1 - Xth_i1) + (Xth_ip1 - Xt)/(Xth_ip1 - Xth_i1)*dB(i + 1)
                    d2B_curr(i) = d2B_curr(i) - (2*dB(i + 1) - (Xth_ip1 - Xt)*d2B(i + 1))/(Xth_ip1 - Xth_i1)
                end if
            end do
            B = B_curr
            dB = dB_curr
            d2B = d2B_curr
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine basis_bspline_2der_C(Xt, knot, nc, degree, d2B)
        integer, intent(in) :: degree
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc
        real(rk), intent(in) :: Xt
        real(rk), intent(out) :: d2B(nc)
        real(rk) :: dB(nc)
        real(rk) :: B(nc)
        integer :: i, p
        real(rk) :: Xth_i, Xth_i1, Xth_ip, Xth_ip1, Xth_last
        real(rk) :: B_curr(nc)
        real(rk) :: dB_curr(nc)
        real(rk) :: d2B_curr(nc)

        B = 0.0_rk
        dB = 0.0_rk
        d2B = 0.0_rk

        ! Degree 0 initialization
        do concurrent(i = 1:nc)
            Xth_i = knot(i)
            Xth_i1 = knot(i + 1)
            Xth_last = knot(size(knot))

            ! if ((Xt >= Xth_i .and. Xt < Xth_i1) .or. (Xt == Xth_last .and. Xt == Xth_i1)) then
            if ((Xt >= Xth_i .and. Xt < Xth_i1) .or. (abs(Xt - Xth_last) < 2.0_rk*epsilon(0.0_rk) .and. abs(Xt - Xth_i1) < 2.0_rk*epsilon(0.0_rk))) then
            B(i) = 1.0_rk
                dB(i) = 0.0_rk
                d2B(i) = 0.0_rk
            end if
        end do

        ! Recursion for higher degrees
        do p = 1, degree
            B_curr = 0.0_rk
            dB_curr = 0.0_rk
            d2B_curr = 0.0_rk
            do concurrent(i = 1:nc)
                Xth_i = knot(i)
                Xth_i1 = knot(i + 1)
                Xth_ip = knot(i + p)
                Xth_ip1 = knot(i + p + 1)

                ! if (Xth_ip /= Xth_i) then
                if (abs(Xth_ip - Xth_i) > 2.0_rk*epsilon(0.0_rk)) then
                    B_curr(i) = (Xt - Xth_i)/(Xth_ip - Xth_i)*B(i)
                    dB_curr(i) = B(i)/(Xth_ip - Xth_i) + (Xt - Xth_i)/(Xth_ip - Xth_i)*dB(i)
                    d2B_curr(i) = (2*dB(i) + (Xt - Xth_i)*d2B(i))/(Xth_ip - Xth_i)
                end if

                ! if (i < nc .and. Xth_ip1 /= Xth_i1) then
                if (i < nc .and. abs(Xth_ip1 - Xth_i1) > 2.0_rk*epsilon(0.0_rk)) then
                    B_curr(i) = B_curr(i) + (Xth_ip1 - Xt)/(Xth_ip1 - Xth_i1)*B(i + 1)
                    dB_curr(i) = dB_curr(i) - B(i + 1)/(Xth_ip1 - Xth_i1) + (Xth_ip1 - Xt)/(Xth_ip1 - Xth_i1)*dB(i + 1)
                    d2B_curr(i) = d2B_curr(i) - (2*dB(i + 1) - (Xth_ip1 - Xt)*d2B(i + 1))/(Xth_ip1 - Xth_i1)
                end if
            end do
            B = B_curr
            dB = dB_curr
            d2B = d2B_curr
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function basis_bernstein(Xt, nc) result(B)
        real(rk), intent(in) :: Xt
        integer, intent(in) :: nc
        real(rk), allocatable :: B(:)
        integer  :: p, degree

        degree = nc - 1

        allocate(B(nc), source=0.0_rk)

        do concurrent (p = 0:degree)
            B(p+1) = gamma(real(nc, kind=rk))/(gamma(real(p+1, kind=rk))*gamma(real(nc-p, kind=rk)))
            ! if (Xt == 0.0_rk .and. p == 0) then
            if (abs(Xt) < 2.0_rk*epsilon(0.0_rk) .and. p == 0) then
                B(p+1) = B(p+1)*(1.0_rk-Xt)**(degree-p)
            ! else if (Xt == 0.0_rk .and. degree-p == 0) then
            else if (abs(Xt) < 2.0_rk*epsilon(0.0_rk) .and. degree-p == 0) then
                B(p+1) = B(p+1)*(Xt**p)
            else
                B(p+1) = B(p+1)*(Xt**p)*(1.0_rk-Xt)**(degree-p)
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function kron_t1_t1(u,v) result(w)
        real(rk), intent(in), contiguous :: u(:), v(:)
        real(rk) :: w(size(u)*size(v))
        integer :: i, j, n

        n = size(v)

        do concurrent(i = 1:size(u), j = 1:n)
            w((i-1)*n + j) = u(i)*v(j)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function kron_t1_t2(u,A) result(B)
        real(rk), intent(in), contiguous :: u(:)
        real(rk), intent(in), contiguous :: A(:,:)
        real(rk) :: B(size(u)*size(A,1), size(A,2))
        integer :: i, j, k, m, n, r, c

        m = size(u)
        r = size(A, 1)
        c = size(A, 2)

        do concurrent (i=1:m, j=1:r, k=1:c)
            B((i-1)*r + j, k) = u(i) * A(j, k)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function kron3(u, v, w) result(out)
        real(rk), intent(in), contiguous :: u(:), v(:), w(:)
        real(rk) :: out(size(u)*size(v)*size(w))
        integer :: i, j, k, nv, nw

        nv = size(v)
        nw = size(w)

        do concurrent(i = 1:size(u), j = 1:nv, k = 1:nw)
            out(((i-1)*nv + j -1)*nw + k) = u(i) * v(j) * w(k)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine ndgrid2(X_dir1,X_dir2, Xt)
        real(rk), intent(in), contiguous :: X_dir1(:), X_dir2(:)
        real(rk), allocatable, intent(out) :: Xt(:,:)
        integer :: s1, s2, i, j

        s1 = size(X_dir1)
        s2 = size(X_dir2)
        allocate(Xt(s1*s2,2))
        do concurrent (j = 1:s2, i = 1: s1)
            Xt((j - 1) * s1 + i,1) = X_dir1(i)
            Xt((j - 1) * s1 + i,2) = X_dir2(j)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine ndgrid3(X_dir1,X_dir2,X_dir3, Xt)
        real(rk), intent(in), contiguous :: X_dir1(:), X_dir2(:), X_dir3(:)
        real(rk), allocatable, intent(out) :: Xt(:,:)
        integer :: s1, s2, s3, i, j, k

        s1 = size(X_dir1)
        s2 = size(X_dir2)
        s3 = size(X_dir3)
        allocate(Xt(s1*s2*s3,3))
        do concurrent (k = 1:s3, j = 1:s2, i = 1: s1)
            Xt(((k - 1) * s2 + (j - 1)) * s1 + i,1) = X_dir1(i)
            Xt(((k - 1) * s2 + (j - 1)) * s1 + i,2) = X_dir2(j)
            Xt(((k - 1) * s2 + (j - 1)) * s1 + i,3) = X_dir3(k)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function repelem(a, b) result(c)
        real(rk), intent(in), contiguous :: a(:)
        integer, intent(in), contiguous :: b(:)
        real(rk) :: c(sum(b))
        integer :: i, l, n

        l = 0
        do i = 1, size(a)
            n = b(i)
            c(l+1:l+n) = a(i)
            l = l + n
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function linspace(a, b, n) result(x)
        real(rk), intent(in) :: a, b
        integer, intent(in) :: n
        real(rk), allocatable :: x(:)
        integer :: i

        if (n < 1) error stop "linspace: n must be ≥ 1"
        allocate(x(n))

        if (n == 1) then
            x(1) = a
        else
            do concurrent(i = 1:n)
                x(i) = a + (i - 1) * (b - a) / real(n - 1, rk)
            end do
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elemConn_C0_L(nnode,p) result(elemConn)
        integer, intent(in) :: nnode
        integer, intent(in) :: p
        integer, allocatable :: elemConn(:,:)
        integer :: i
        integer, allocatable :: nodes(:)

        if (mod(nnode-1,p) /= 0) error stop 'cmp_elemConn_C0_L: nnode-1 must be divisible by p'

        allocate(elemConn( (nnode-1) / p ,p+1))
        nodes = [(i, i=1,nnode)]
        do concurrent (i = 1:nnode-p:p)
            elemConn((i-1)/p+1, :) = reshape(nodes(i:i+p), [p + 1])
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elemConn_C0_S(nnode1, nnode2, p1, p2) result(elemConn)
        integer, intent(in) :: nnode1, nnode2, p1, p2
        integer, allocatable :: elemConn(:,:)
        integer :: i, j, nelem1, nelem2, nnel
        integer, allocatable :: nodes(:,:)

        if (mod(nnode1 - 1, p1) /= 0) error stop 'cmp_elemConn_C0_S: nnode1-1 must be divisible by p1'
        if (mod(nnode2 - 1, p2) /= 0) error stop 'cmp_elemConn_C0_S: nnode2-1 must be divisible by p2'

        nelem1 = (nnode1 - 1) / p1
        nelem2 = (nnode2 - 1) / p2
        nnel   = (p1 + 1) * (p2 + 1)

        allocate(elemConn(nelem1 * nelem2, nnel))
        nodes = reshape([(i, i = 1, nnode1 * nnode2)], [nnode1, nnode2])

#if defined(__NVCOMPILER)
        do i = 1,nnode1 - p1,p1
            do j = 1,nnode2 - p2,p2
            elemConn(((j-1)/p2)*nelem1+(i-1)/p1+1, :) = reshape(nodes(i:i+p1, j:j+p2), [nnel])
            end do
        end do
#else
        do concurrent (j = 1:nnode2 - p2:p2, i = 1:nnode1 - p1:p1)
            elemConn(((j-1)/p2)*nelem1+(i-1)/p1+1, :) = reshape(nodes(i:i+p1, j:j+p2), [nnel])
        end do
#endif
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elemConn_C0_V(nnode1,nnode2,nnode3,p1,p2,p3) result(elemConn)
        integer, intent(in) :: nnode1, nnode2, nnode3, p1, p2, p3
        integer, allocatable :: elemConn(:,:)
        integer :: i, j, k, nnel, nelem1, nelem2, nelem3
        integer, allocatable :: nodes(:,:,:)

        if (mod(nnode1-1,p1) /= 0) error stop 'cmp_elemConn_C0_V: nnode1-1 must be divisible by p1'
        if (mod(nnode2-1,p2) /= 0) error stop 'cmp_elemConn_C0_V: nnode2-1 must be divisible by p2'
        if (mod(nnode3-1,p3) /= 0) error stop 'cmp_elemConn_C0_V: nnode3-1 must be divisible by p3'

        nelem1 = (nnode1 - 1) / p1
        nelem2 = (nnode2 - 1) / p2
        nelem3 = (nnode3 - 1) / p3

        nnel = (p1 + 1) * (p2 + 1) * (p3 + 1)
        allocate(elemConn(nelem1 * nelem2 * nelem3, nnel))
        nodes = reshape([(i, i=1, nnode1 * nnode2 * nnode3)], [nnode1, nnode2, nnode3])
#if defined(__NVCOMPILER)
        do k = 1,nnode3 - p3,p3
            do j = 1,nnode2 - p2,p2
                do i = 1,nnode1 - p1,p1
                    elemConn(((k-1)/p3)*(nelem1*nelem2)+((j-1)/p2)*nelem1+((i-1)/p1)+1, :) = reshape(nodes(i:i+p1, j:j+p2, k:k+p3), [nnel])
                end do
            end do
        end do
#else
        do concurrent (k = 1:nnode3-p3:p3, j = 1:nnode2-p2:p2, i = 1:nnode1-p1:p1)
            elemConn(((k-1)/p3)*(nelem1*nelem2)+((j-1)/p2)*nelem1+((i-1)/p1)+1, :) = reshape(nodes(i:i+p1, j:j+p2, k:k+p3), [nnel])
        end do
#endif
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine cmp_elemConn_Cn_L(nnode, p, Xth, vecKnot_mul, elemConn)
        integer, intent(in) :: p, nnode
        integer, intent(in), contiguous :: vecKnot_mul(:)
        real(rk), intent(in), contiguous :: Xth(:)
        integer, allocatable, intent(out) :: elemConn(:,:)
        integer, allocatable :: nodes(:)
        integer :: i, m, nnel, nelem

        nnel  = p + 1
        nelem = size(Xth) - 1

        allocate(nodes(nnode))
        nodes = [(i, i=1, nnode)]

        allocate(elemConn(nelem, nnel))

        do concurrent (i = 1:nelem)
            m = -p + sum(vecKnot_mul(1:i))
            elemConn(i,:) = nodes(m:m+p)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine cmp_elemConn_Cn_S(nnode1, nnode2, p1, p2, &
        Xth1, Xth2, vecKnot_mul1, vecKnot_mul2, elemConn)
        integer, intent(in) :: p1, p2, nnode1, nnode2
        integer, intent(in), contiguous :: vecKnot_mul1(:), vecKnot_mul2(:)
        real(rk), intent(in), contiguous :: Xth1(:), Xth2(:)
        integer, allocatable, intent(out) :: elemConn(:,:)
        integer, allocatable :: nodes(:,:), nodes_vec(:)
        integer :: nnd_total, i, j, l, m, n, nnel1, nnel2, nelem1, nelem2, nelem

        nnel1 = p1 + 1
        nnel2 = p2 + 1
        nelem1 = size(Xth1) - 1
        nelem2 = size(Xth2) - 1
        nelem  = nelem1 * nelem2

        nnd_total = nnode1 * nnode2
        allocate(nodes_vec(nnd_total))
        nodes_vec = [(i, i=1, nnd_total)]
        nodes = reshape(nodes_vec, [nnode1, nnode2])

        allocate(elemConn(nelem, nnel1 * nnel2))

#if defined(__NVCOMPILER)
        do j = 1, nelem2
            do i = 1, nelem1
                m = -p1 + sum(vecKnot_mul1(1:i))
                n = -p2 + sum(vecKnot_mul2(1:j))
                l = (j - 1) * nelem1 + i
                elemConn(l,:) = reshape(nodes(m:m+p1, n:n+p2), [nnel1 * nnel2])
            end do
        end do
#else
        do concurrent (j = 1:nelem2, i = 1:nelem1)
            m = -p1 + sum(vecKnot_mul1(1:i))
            n = -p2 + sum(vecKnot_mul2(1:j))
            l = (j - 1) * nelem1 + i
            elemConn(l,:) = reshape(nodes(m:m+p1, n:n+p2), [nnel1 * nnel2])
        end do
#endif
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine cmp_elemConn_Cn_V(nnode1,nnode2,nnode3,p1,p2,p3,&
        Xth1,Xth2,Xth3,vecKnot_mul1,vecKnot_mul2,vecKnot_mul3, elemConn)
        integer, intent(in) :: p1, p2, p3, nnode1, nnode2, nnode3
        integer, intent(in), contiguous :: vecKnot_mul1(:), vecKnot_mul2(:), vecKnot_mul3(:)
        real(rk), intent(in), contiguous :: Xth1(:), Xth2(:), Xth3(:)
        integer, allocatable, intent(out) :: elemConn(:,:)
        integer, allocatable :: nodes(:,:,:), nodes_vec(:)
        integer :: nnd_total, i, j, k, l, nnel1, nnel2, nnel3, nnel, m, n, o, nelem1, nelem2, nelem3, nelem

        nnel1 = p1 + 1
        nnel2 = p2 + 1
        nnel3 = p3 + 1

        nnd_total = nnode1*nnode2*nnode3
        allocate(nodes_vec(nnd_total))

        nodes_vec = [(i, i=1, nnd_total)]
        nodes = reshape(nodes_vec,[nnode1,nnode2,nnode3])

        nelem1 = size(Xth1) - 1
        nelem2 = size(Xth2) - 1
        nelem3 = size(Xth3) - 1

        nelem = nelem1*nelem2*nelem3

        nnel = nnel1*nnel2*nnel3
        allocate(elemConn(nelem,nnel))
#if defined(__NVCOMPILER)
        do k = 1, nelem3
            do j = 1, nelem2
                do i = 1, nelem1
                    o = -p3 + sum(vecKnot_mul3(1:k))
                    n = -p2 + sum(vecKnot_mul2(1:j))
                    m = -p1 + sum(vecKnot_mul1(1:i))
                    l = (k - 1) * nelem1 * nelem2 + (j - 1) * nelem1 + i
                    elemConn(l, :) = reshape( nodes(m:m+p1, n:n+p2, o:o+p3), [nnel] )
                end do
            end do
        end do
#else
        do concurrent (k = 1:nelem3, j = 1:nelem2, i = 1:nelem1)
            o = -p3 + sum(vecKnot_mul3(1:k))
            n = -p2 + sum(vecKnot_mul2(1:j))
            m = -p1 + sum(vecKnot_mul1(1:i))
            l = (k - 1) * nelem1 * nelem2 + (j - 1) * nelem1 + i
            elemConn(l, :) = reshape( nodes(m:m+p1, n:n+p2, o:o+p3), [nnel] )
        end do
#endif
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_multiplicity1(knot) result(multiplicity)
        real(rk), intent(in), contiguous :: knot(:)
        integer, allocatable :: multiplicity(:)
        integer :: i, count

        count = 1
        do i = 2, size(knot)
            ! if (knot(i) /= knot(i-1)) count = count + 1
            if (abs(knot(i) - knot(i-1)) > 2.0_rk*epsilon(0.0_rk)) count = count + 1
        end do

        allocate(multiplicity(count))

        multiplicity(1) = 1
        count = 1

        do i = 2, size(knot)
            ! if (knot(i) /= knot(i-1)) then
            if (abs(knot(i) - knot(i-1)) > 2.0_rk*epsilon(0.0_rk)) then
                count = count + 1
                multiplicity(count) = 1
            else
                multiplicity(count) = multiplicity(count) + 1
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_multiplicity2(knot, Xth) result(multiplicity)
        real(rk), intent(in), contiguous :: knot(:)
        real(rk), intent(in) :: Xth
        integer :: multiplicity
        integer :: i, count, size_knot

        size_knot = size(knot)
        multiplicity = 0
        i = 1
        do while (i <= size_knot)
            ! if (knot(i) == Xth) then
            if (abs(knot(i) - Xth) < 2.0_rk*epsilon(0.0_rk)) then
                count = 1
                ! do while (i + count <= size_knot .and. knot(i + count) == Xth)
                do while ((i + count) <= size_knot)
                    if (abs(knot(i + count) - Xth) >= 2.0_rk*epsilon(0.0_rk)) exit
                    count = count + 1
                end do
                if (count > multiplicity) then
                    multiplicity = count
                end if
                i = i + count
            else
                i = i + 1
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_knot_vector(Xth_dir, degree, continuity) result(knot)
        real(rk), intent(in), contiguous :: Xth_dir(:)
        integer, intent(in) :: degree
        integer, intent(in), contiguous :: continuity(:)
        real(rk), allocatable :: knot(:)

        knot = repelem(Xth_dir, (degree - continuity))
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine insert_knot_A_5_1(p, UP, Pw, u, k, s, r, nq, UQ, Qw)
        integer, intent(in) :: p, k, s, r
        real(rk), intent(in), contiguous :: UP(0:), Pw(0:,:)
        real(rk), intent(in) :: u
        real(rk), allocatable, intent(out) :: UQ(:), Qw(:,:)
        integer, intent(out) :: nq
        integer :: i, j, L, mp, d, np
        real(rk), allocatable :: Rw(:,:)
        real(rk) :: alpha

        d = size(Pw, 2)
        np  = size(Pw, 1) - 1
        mp  = np + p + 1
        nq  = np + r

        allocate(UQ(0:mp+r))
        allocate(Qw(0:nq,1:d))
        allocate(Rw(0:p ,1:d))

        UQ(0:k) = UP(0:k)
        UQ(k+1:k+r) = u
        UQ(k+1+r:mp+r) = UP(k+1:mp)
        Qw(0:k-p,:) = Pw(0:k-p,:)
        Qw(k-s+r:np+r,:) = Pw(k-s:np,:)
        Rw(0:p-s,:) = Pw(k-p:k-s,:)
        do j = 1, r
            L = k-p+j
            do i = 0, p-j-s
                alpha = (u - UP(L+i)) / (UP(i+k+1) - UP(L+i))
                Rw(i,:) = alpha*Rw(i+1,:) + (1.0_rk - alpha) * Rw(i,:)
            end do
            Qw(L,:) = Rw(0,:)
            Qw(k+r-j-s,:) = Rw(p-j-s,:)
        end do
        Qw(L+1:k-s-1,:) = Rw(1:k-s-1-L,:)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function findspan(n,degree,Xth,knot) result(s)
        integer, intent(in) :: n, degree
        real(rk), intent(in) :: Xth
        real(rk), intent(in), contiguous :: knot(:)
        integer :: s
        integer :: low, high, mid
        ! if (Xth == knot(n+2)) then
        if (abs(Xth - knot(n+2)) < 2.0_rk*epsilon(0.0_rk)) then
            s = n
            return
        end if
        low = degree
        high = n + 1
        mid = (low + high) / 2
        do while (Xth < knot(mid+1) .or. Xth >= knot(mid+2))
            if (Xth < knot(mid+1)) then
                high = mid
            else
                low = mid
            end if
            mid = (low + high) / 2
        end do
        s = mid
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine elevate_degree_A_5_9(t, knot, degree, Xcw, nc_new, knot_new, Xcw_new)
        integer, intent(in) :: t
        real(rk), intent(in), contiguous :: Xcw(:,:), knot(:)
        integer, intent(in) :: degree
        integer, intent(out) :: nc_new
        real(rk), allocatable, intent(out) :: Xcw_new(:,:), knot_new(:)
        real(rk), allocatable :: bezalfs(:,:), bpts(:,:), ebpts(:,:), Nextbpts(:,:), alfs(:)
        real(rk) :: iinv, alpha1, alpha2, Xth1, Xth2, numer, den
        integer :: n, lbz, rbz, sv, tr, kj, first, knoti, last, alpha3, d, nc
        integer :: i, j, q, s, m, ph, ph2, mpi, mh, r, a, b, Xcwi, oldr, mul
        integer, allocatable :: mlp(:)

        nc = size(Xcw,1)
        d = size(Xcw,2)
        mlp = compute_multiplicity(knot)
        mlp = mlp + t
        nc_new = sum(mlp) - (mlp(1)-1) - 1
        allocate(Xcw_new(nc_new,d), source=0.0_rk)
        allocate(bezalfs(degree+1,degree+t+1), source=0.0_rk)
        allocate(bpts(degree+1,d), source=0.0_rk)
        allocate(ebpts(degree+t+1,d), source=0.0_rk)
        allocate(Nextbpts(degree+1,d), source=0.0_rk)
        allocate(alfs(degree), source=0.0_rk)
        n = nc - 1
        m = n + degree + 1
        ph = degree + t
        ph2 = ph / 2
        bezalfs(1,1) = 1.0_rk
        bezalfs(degree+1,ph+1) = 1.0_rk
        do i = 1,ph2
            iinv = 1.0_rk/bincoeff(ph,i)
            mpi = min(degree,i)
            do j = max(0,i-t),mpi
                bezalfs(j+1,i+1) = iinv*bincoeff(degree,j)*bincoeff(t,i-j)
            end do
        end do
        do i = ph2+1,ph-1
            mpi = min(degree,i)
            do j = max(0,i-t),mpi
                bezalfs(j+1,i+1) = bezalfs(degree-j+1,ph-i+1)
            end do
        end do
        mh = ph
        knoti = ph+1
        r = -1
        a = degree
        b = degree+1
        Xcwi = 1
        Xth1 = knot(1)
        Xcw_new(1,:) = Xcw(1,:)
        allocate(knot_new(sum(mlp)), source=0.0_rk)
        knot_new(1:ph+1) = Xth1
        do i = 0,degree
            bpts(i+1,:) = Xcw(i+1,:)
        end do
        do while (b<m)
            i = b
            ! do while (b<m .and. knot(b+1) == knot(b+2))
            do while (b<m .and. abs(knot(b+1)-knot(b+2)) < 2.0_rk*epsilon(0.0_rk))
                b = b + 1
                if (b+2 > size(knot)) then
                    exit
                end if
            end do
            mul = b - i + 1
            mh = mh + mul + t
            Xth2 = knot(b+1)
            oldr = r
            r = degree - mul
            if (oldr > 0) then
                lbz = (oldr+2)/2
            else
                lbz = 1
            end if
            if (r > 0) then
                rbz = ph - (r+1)/2
            else
                rbz = ph
            end if
            if (r>0) then
                numer = Xth2 - Xth1
                do q = degree,mul+1,-1
                    alfs(q-mul) = numer / (knot(a+q+1)-Xth1)
                end do
                do j = 1,r
                    sv = r - j
                    s = mul + j
                    do q = degree,s,-1
                        bpts(q+1,:) = (1.0_rk-alfs(q-s+1))*bpts(q,:) + alfs(q-s+1)*bpts(q+1,:)
                    end do
                    Nextbpts(sv+1,:) = bpts(degree+1,:)
                end do
            end if
            do i = lbz,ph
                ebpts(i+1,:) = 0.0_rk
                mpi = min(degree,i)
                do j = max(0,i-t),mpi
                    ebpts(i+1,:) = bezalfs(j+1,i+1)*bpts(j+1,:) + ebpts(i+1,:)
                end do
            end do
            if (oldr > 1) then
                first = knoti - 2
                last = knoti
                den = Xth2 - Xth1
                alpha3 = floor((Xth2-knot(knoti)) / den)
                do tr = 1,oldr-1
                    i = first
                    j = last
                    kj = j - knoti + 1
                    do while (j-i > tr)
                        if (i < Xcwi) then
                            alpha1 = (Xth2-knot(i+1))/(Xth1-knot(i+1))
                            Xcw_new(i+1,:) = (1-alpha1)*Xcw_new(i,:) + alpha1*Xcw_new(i+1,:)
                        end if
                        if (j >= lbz) then
                            if (j-tr <= knoti-ph+oldr) then
                                alpha2 = (Xth2-knot_new(j-tr+1)) / den
                                ebpts(kj+1,:) = alpha2*ebpts(kj+1,:) + (1-alpha2)*ebpts(kj+2,:)
                            else
                                ebpts(kj+1,:) = (1-alpha3)*ebpts(:,kj+2) + alpha3*ebpts(kj+1,:)
                            end if
                        end if
                        i = i + 1
                        j = j - 1
                        kj = kj - 1
                    end do
                    first = first - 1
                    last = last + 1
                end do
            end if
            if (a /= degree) then
                do i = 0,ph-oldr-1
                    knot_new(knoti+1) = Xth1
                    knoti = knoti + 1
                end do
            end if
            do j = lbz,rbz
                Xcw_new(Xcwi+1,:) = ebpts(j+1,:)
                Xcwi = Xcwi + 1
            end do
            if (b<m) then
                do j = 0,r-1
                    bpts(j+1,:) = Nextbpts(j+1,:)
                end do
                do j = r,degree
                    bpts(j+1,:) = Xcw(b-degree+j+1,:)
                end do
                a = b
                b = b+1
                Xth1 = Xth2
            else
                do i = 0,ph
                    knot_new(knoti+i+1) = Xth2
                end do
            end if
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure elemental function bincoeff(n,k) result(b)
        integer, intent(in) :: n, k
        real(rk) :: b
        b = floor(0.5_rk+exp(factln(n)-factln(k)-factln(n-k)))
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure elemental function factln(n) result(f)
        integer, intent(in) :: n
        real(rk) :: f
        if (n <= 1) then
            f = 0.0_rk
            return
        end if
        f = log(gamma(real(n+1,rk)))
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function hexahedron_Xc(L, nc) result(Xc)
        real(rk), intent(in) :: L(3)
        integer, intent(in) :: nc(3)
        real(rk), allocatable :: Xc(:,:)
        real(rk) :: dx, dy, dz
        integer :: i, j, k, nci

        dx = L(1) / real(nc(1)-1, rk)
        dy = L(2) / real(nc(2)-1, rk)
        dz = L(3) / real(nc(3)-1, rk)

        allocate(Xc(nc(1) * nc(2) * nc(3), 3))
        do concurrent (k = 0:nc(3)-1, j = 0:nc(2)-1, i = 0:nc(1)-1)
            Xc(i + j * nc(1) + k * nc(1) * nc(2) + 1, 1) = real(i, rk) * dx
            Xc(i + j * nc(1) + k * nc(1) * nc(2) + 1, 2) = real(j, rk) * dy
            Xc(i + j * nc(1) + k * nc(1) * nc(2) + 1, 3) = real(k, rk) * dz
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function tetragon_Xc(L, nc) result(Xc)
        real(rk), intent(in) :: L(2)
        integer, intent(in) :: nc(2)
        real(rk), allocatable :: Xc(:,:)
        real(rk) :: dx, dy
        integer :: i, j, nci

        dx = L(1) / real(nc(1)-1, rk)
        dy = L(2) / real(nc(2)-1, rk)

        allocate(Xc(nc(1) * nc(2), 3))
        do concurrent (j = 0:nc(2)-1, i = 0:nc(1)-1)
            Xc(i + j * nc(1) + 1, 1) = real(i, rk) * dx
            Xc(i + j * nc(1) + 1, 2) = real(j, rk) * dy
            Xc(i + j * nc(1) + 1, 3) = 0.0_rk
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine remove_knots_A_5_8(p,knot,Pw,u,r,s,num,t, knot_new, Pw_new)
        real(rk), intent(in) :: u
        integer, intent(in) :: p, r, s, num
        real(rk), intent(in), contiguous :: knot(:)
        real(rk), intent(in), contiguous :: Pw(:,:)
        real(rk), allocatable, intent(out) :: knot_new(:)
        real(rk), allocatable, intent(out) :: Pw_new(:,:)
        real(rk), allocatable :: Pw_copy(:,:), knot_copy(:)
        integer, intent(out) :: t
        real(rk) :: tol, alfi, alfj
        real(rk), allocatable :: temp(:,:)
        integer :: i, j, ii, jj, remflag, off, first, last, ord, fout, m, k, n, nc, d, tt

        d  = size(Pw,2)
        nc = size(Pw,1)
        n = nc
        m = n+p+1
        ord = p+1
        fout = (2*r-s-p)/2
        last = r-s
        first = r-p

        Pw_copy = Pw
        knot_copy = knot

        ! TODO:
        tol = 1.0e-6_rk * minval(Pw(:,d))/(1.0_rk + maxval(sqrt(sum(Pw**2, 2))))

        allocate(temp(2*p+1,d), source=0.0_rk)
        t = 0
        do tt = 0, num-1
            off = first-1
            temp(1,:) = Pw_copy(off,:)
            temp(last+1-off+1,:) = Pw_copy(last+1,:)
            i = first
            j = last
            ii = 1
            jj = last-off
            remflag = 0
            do while (j-i>t)
                alfi = (u-knot_copy(i))/(knot_copy(i+ord+t)-knot_copy(i))
                alfj = (u-knot_copy(j-t))/(knot_copy(j+ord)-knot_copy(j-t))
                temp(ii+1,:) = (Pw_copy(i,:)-(1.0_rk-alfi)*temp(ii-1+1,:))/alfi
                temp(jj+1,:) = (Pw_copy(j,:)-alfj*temp(jj+1+1,:))/(1.0_rk-alfj)
                i = i+1
                ii = ii+1
                j = j-1
                jj = jj-1
            end do
            if (j-i<=t) then
                if (norm2(temp(ii-1+1,:) - temp(jj+1+1,:))<=tol) then
                    remflag = 1
                else
                    alfi = (u-knot_copy(i))/(knot_copy(i+ord+t)-knot_copy(i))
                    if (norm2(Pw_copy(i,:) - (alfi*temp(ii+t+1+1,:)+(1.0_rk-alfi)*temp(ii-1+1,:)))<=tol) then
                        remflag = 1
                    end if
                end if
            end if
            if (remflag == 0) then
                exit
            else
                i = first
                j = last
                do while(j-i>t)
                    Pw_copy(i,:) = temp(i-off+1,:)
                    Pw_copy(j,:) = temp(j-off+1,:)
                    i = i+1
                    j = j-1
                end do
            end if
            first = first-1
            last = last+1
            t=t+1
        end do
        if (t==0) then
            return
        end if
        do k = r+1, m
            knot_copy(k-t) = knot_copy(k)
        end do
        j = fout
        i = j
        do k = 1, t-1
            if (mod(k,2)==1) then
                i = i+1
            else
                j = j-1
            end if
        end do
        do k = i+1, n
            Pw_copy(j,:) = Pw_copy(k,:)
            j = j+1
        end do
        knot_new = knot_copy(1:size(knot_copy)-t)
        Pw_new = Pw_copy(1:size(Pw_copy,1)-t,:)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function unique_integer(vec) result(output)
        integer, dimension(:), intent(in), contiguous :: vec
        integer, dimension(:), allocatable :: output
        integer :: i, j, k
        allocate(output(0))
        do i = 1, size(vec)
            k = 0
            do j = 1, size(output)
                if (vec(i) == output(j)) then
                    k = k + 1
                    exit
                end if
            end do
            if (k == 0) then
                output = [output, vec(i)]
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function unique_real(vec) result(output)
        real(rk), dimension(:), intent(in), contiguous :: vec
        real(rk), dimension(:), allocatable :: output
        integer :: i, j, k
        allocate(output(0))
        do i = 1, size(vec)
            k = 0
            do j = 1, size(output)
                ! if (vec(i) == output(j)) then
                if (abs(vec(i) - output(j)) < 2.0_rk*epsilon(0.0_rk)) then
                    k = k + 1
                    exit
                end if
            end do
            if (k == 0) then
                output = [output, vec(i)]
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function rotation(alpha, beta, theta) result(R)
        real(rk), intent(in) :: alpha, beta, theta
        real(rk), dimension(3,3) :: R

        R(1,1) = cosd(beta)*cosd(theta)
        R(2,1) = cosd(beta)*sind(theta)
        R(3,1) = -sind(beta)
        R(1,2) = sind(alpha)*sind(beta)*cosd(theta) - cosd(alpha)*sind(theta)
        R(2,2) = sind(alpha)*sind(beta)*sind(theta) + cosd(alpha)*cosd(theta)
        R(3,2) = sind(alpha)*cosd(beta)
        R(1,3) = cosd(alpha)*sind(beta)*cosd(theta) + sind(alpha)*sind(theta)
        R(2,3) = cosd(alpha)*sind(beta)*sind(theta) - sind(alpha)*cosd(theta)
        R(3,3) = cosd(alpha)*cosd(beta)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function det(A) result(detA)
        real(rk), intent(in) :: A(:,:)
        real(rk) :: detA

        if (size(A,1) == size(A,2)) then
            select case(size(A,1))
              case(2)
                detA = A(1,1)*A(2,2) - A(1,2)*A(2,1)
              case(3)
                detA = &
                    + A(1,1)*( A(2,2)*A(3,3) - A(2,3)*A(3,2) )&
                    - A(1,2)*( A(2,1)*A(3,3) - A(2,3)*A(3,1) )&
                    + A(1,3)*( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
            end select
        elseif (size(A,1) == 3 .and. size(A,2) == 2) then
            detA = &
                + A(1,1) * ( A(2,2) * 1.0_rk - A(3,2) * 1.0_rk )&
                - A(1,2) * ( A(2,1) * 1.0_rk - A(3,1) * 1.0_rk )&
                + 1.0_rk * ( A(2,1) * A(3,2) - A(3,1) * A(2,2) )
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    recursive pure function inv(A) result(A_inv)
        real(rk), intent(in) :: A(:,:)
        real(rk), allocatable :: A_inv(:,:)

        if (size(A,1) == size(A,2)) then
            select case(size(A,1))
              case(2)
                allocate(A_inv(size(A,1),size(A,2)))
                A_inv(1,1) =  A(2,2)
                A_inv(1,2) = -A(1,2)
                A_inv(2,1) = -A(2,1)
                A_inv(2,2) =  A(1,1)
                A_inv = A_inv/det(A)
              case(3)
                allocate(A_inv(size(A,1),size(A,2)))
                A_inv(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
                A_inv(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
                A_inv(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
                A_inv(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
                A_inv(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
                A_inv(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
                A_inv(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
                A_inv(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
                A_inv(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)
                A_inv = A_inv/det(A)
              case default
                A_inv = solve(A,eye(size(A,1)))
            end select
        elseif (size(A,1)>size(A,2)) then
            allocate(A_inv(size(A,2),size(A,1)))
            A_inv = transpose(A)
            A_inv = matmul(inv(matmul(A_inv, A)), A_inv)
        elseif (size(A,1)<size(A,2)) then
            allocate(A_inv(size(A,2),size(A,1)))
            A_inv = transpose(A)
            A_inv = matmul(A_inv, inv(matmul(A, A_inv)))
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function eye(n) result(I)
        integer, intent(in) :: n
        real(rk), allocatable :: I(:,:)

        ! local variables
        integer :: k

        allocate(I(n,n), source=0.0_rk)
        do concurrent (k = 1: n)
            I(k, k) = 1.0_rk
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function dyad_t1_t1(a, b) result(c)
        real(rk), intent(in), contiguous :: a(:)
        real(rk), intent(in), contiguous :: b(:)
        real(rk), allocatable :: c(:,:)
        integer :: i

        allocate(c(size(a), size(b)))
        do concurrent(i = 1:size(c,1))
            c(i, :) = a(i) * b(:)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine gauss_legendre_1D(interval, degree, Xksi, Wksi)
        real(rk), intent(in) :: interval(2)
        integer, intent(in) :: degree
        real(rk), allocatable, intent(out) :: Xksi(:), Wksi(:)

        allocate(Xksi(degree+1), Wksi(degree+1))

        call gauss_legendre(Xksi, Wksi, interval)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine gauss_legendre_2D(interval1, interval2, degree, Xksi, Wksi)
        real(rk), intent(in) :: interval1(2), interval2(2)
        integer, intent(in) :: degree(2)
        real(rk), allocatable, intent(out) :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: Xksi1(:), Wksi1(:), Xksi2(:), Wksi2(:)

        allocate(Xksi1(degree(1)+1), Wksi1(degree(1)+1))
        allocate(Xksi2(degree(2)+1), Wksi2(degree(2)+1))

        call gauss_legendre(Xksi1, Wksi1, interval1)
        call gauss_legendre(Xksi2, Wksi2, interval2)

        call ndgrid(Xksi1, Xksi2, Xksi)
        Wksi = kron(Wksi1, Wksi2)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine gauss_legendre_3D(interval1, interval2, interval3, degree, Xksi, Wksi)
        real(rk), intent(in) :: interval1(2), interval2(2), interval3(2)
        integer, intent(in) :: degree(3)
        real(rk), allocatable, intent(out) :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: Xksi1(:), Wksi1(:), Xksi2(:), Wksi2(:), Xksi3(:), Wksi3(:)

        allocate(Xksi1(degree(1)+1), Wksi1(degree(1)+1))
        allocate(Xksi2(degree(2)+1), Wksi2(degree(2)+1))
        allocate(Xksi3(degree(3)+1), Wksi3(degree(3)+1))

        call gauss_legendre(Xksi1, Wksi1, interval1)
        call gauss_legendre(Xksi2, Wksi2, interval2)
        call gauss_legendre(Xksi3, Wksi3, interval3)

        call ndgrid(Xksi1, Xksi2, Xksi3, Xksi)
        Wksi = kron(kron(Wksi3, Wksi2), Wksi1)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine gauss_legendre(x, w, interval)
        real(rk), intent(out) :: x(:), w(:)
        real(rk), intent(in) :: interval(2)
        real(rk) :: xi, delta, p_next, dp_next, p_prev, p_curr, dp_prev, dp_curr, midpoint, half_length
        integer :: i, j, k, n
        real(rk), parameter :: pi = acos(-1.0_rk)
        real(rk), parameter :: tol = 4.0_rk*epsilon(1.0_rk)
        integer, parameter :: maxit = 100
        logical :: converged

        if (interval(1) >= interval(2)) error stop "gauss_legendre: Invalid interval, interval(1) must be less than interval(2)"
        n = size(x)
        ! Gauss-Legendre points are symmetric, only compute half
        do concurrent (i = 1:(n+1)/2)
            ! Initial guess (Chebyshev approximation)
            xi = -cos(pi * (i-0.25_rk)/(n+0.5_rk))
            ! Newton iteration
            j = 0
            converged = .false.
            do while (.not. converged .and. j < maxit)
                j = j + 1
                ! Compute Legendre polynomial and derivative via recurrence
                p_prev = 1.0_rk        ! P_0(xi)
                p_curr = xi            ! P_1(xi)
                dp_prev = 0.0_rk       ! P_0d(xi)
                dp_curr = 1.0_rk       ! P_1d(xi)
                do k = 2, n
                    p_next  = ((2*k-1)*xi*p_curr-(k-1)*p_prev)/k
                    dp_next = ((2*k-1)*(xi*dp_curr+p_curr)-(k-1)*dp_prev)/k
                    p_prev  = p_curr
                    p_curr  = p_next
                    dp_prev = dp_curr
                    dp_curr = dp_next
                end do
                ! Newton correction
                delta = -p_curr / dp_curr
                xi = xi + delta
                ! Check for convergence
                converged = (abs(delta) <= tol*abs(xi))
            end do
#if defined(__NVCOMPILER)
            ! if (.not. converged) error stop "gauss_legendre: Newton iteration did not converge"
#else
            if (.not. converged) error stop "gauss_legendre: Newton iteration did not converge"
#endif
            ! Store symmetric nodes and weights
            x(i)     = xi
            x(n+1-i) = -xi
            w(i)     = 2.0_rk/((1.0_rk-xi**2)*dp_curr**2)
            w(n+1-i) = w(i)
        end do
        ! Transform from [-1,1] to [interval(1), interval(2)]
        midpoint    =0.5_rk*(interval(1)+interval(2))
        half_length =0.5_rk*(interval(2)-interval(1))
        x = midpoint+half_length*x
        w = half_length*w
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_vtk_legacy(filename, points, elemConn, vtkCellType, point_data, field_names, encoding)
        character(len=*), intent(in) :: filename
        real(rk), intent(in) :: points(:,:)
        integer, intent(in) :: elemConn(:,:)
        integer, intent(in) :: vtkCellType
        real(rk), intent(in), optional :: point_data(:,:)     ! [npoints, nfields]
        character(len=*), intent(in), optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding

        integer :: i, j, ne, np, nn, n, nunit
        character(len=6) :: encoding_
        integer, parameter :: dp = kind(1.0d0)

        ne = size(elemConn, 1)
        nn = size(elemConn, 2)
        np = size(points, 1)
        n = ne*(nn+1)

        if (present(encoding)) then
            select case (trim(encoding))
                case ('ascii')
                    encoding_ = 'ASCII'
                case ('binary')
                    encoding_ = 'BINARY'
                case default
                    error stop 'Invalid encoding type. Use "ASCII" or "BINARY".'
            end select
        else
            encoding_ = 'BINARY'
        end if


        if (trim(encoding_) == 'ASCII') then
            open(newunit=nunit, file=filename, action='write')
            write(nunit,'(a)') '# vtk DataFile Version 2.0'
            write(nunit,'(a)') 'Generated by ForCAD'
            write(nunit,'(a)') 'ASCII'
            write(nunit,'(a)') 'DATASET UNSTRUCTURED_GRID'

            write(nunit,'(a,1x,g0,1x,a)') 'POINTS', np, 'double'
            if (size(points,2) == 2) then
                write(nunit,'(g0,1x,g0,1x,g0)') (points(i,1), points(i,2), 0.0_rk , i = 1, np)
            elseif (size(points,2) == 3) then
                write(nunit,'(g0,1x,g0,1x,g0)') (points(i,1), points(i,2), points(i,3) , i = 1, np)
            else
                error stop 'Invalid dimension for points.'
            end if

            write(nunit,'(a,1x,g0,1x,g0)') 'CELLS', ne, n
            select case (nn)
            case (2)
                write(nunit,'(g0,1x,g0,1x,g0)')&
                    (2, elemConn(i,1)-1,elemConn(i,2)-1, i = 1, ne)
            case (4)
                write(nunit,'(g0,1x,g0,1x,g0,1x,g0)')&
                    (4, elemConn(i,1)-1,elemConn(i,2)-1,elemConn(i,4)-1,elemConn(i,3)-1, i = 1, ne)
            case (8)
                write(nunit,'(g0,1x,g0,1x,g0,1x,g0,1x,g0,1x,g0,1x,g0,1x,g0,1x,g0)')&
                    (8, elemConn(i,1)-1,elemConn(i,2)-1,elemConn(i,4)-1,elemConn(i,3)-1,&
                    elemConn(i,5)-1,elemConn(i,6)-1,elemConn(i,8)-1,elemConn(i,7)-1, i = 1, ne)
            case default
                error stop 'Invalid number of nodes per element.'
            end select

            write(nunit,'(a,1x,g0)') 'CELL_TYPES', ne
            write(nunit,'(g0)') (vtkCellType , i = 1, ne)

            if (present(point_data) .and. present(field_names)) then
                write(nunit, '(a,1x,g0)') 'POINT_DATA', size(point_data,1)
                do i = 1, size(point_data,2)
                    write(nunit, '(a,1x,a,1x,a)') 'SCALARS', trim(field_names(i)), 'double'
                    write(nunit, '(a)') 'LOOKUP_TABLE default'
                    write(nunit, '(g0)') (point_data(j,i), j = 1, size(point_data,1))
                end do
            end if

            close(nunit)
        end if


        if (trim(encoding_) == 'BINARY') then
            open(newunit=nunit, file=filename, form='formatted', action='write')
            write(nunit,'(a)') '# vtk DataFile Version 2.0'
            write(nunit,'(a)') 'Generated by ForCAD'
            write(nunit,'(a)') 'BINARY'
            write(nunit,'(a)') 'DATASET UNSTRUCTURED_GRID'
            close(nunit)

            open(newunit=nunit, file=filename, form='formatted', action='write', position='append')
            write(nunit,'(a,1x,g0,1x,a)') 'POINTS', np, 'double'
            close(nunit)
            open(newunit=nunit, file=filename, position='append', access="stream", form="unformatted",&
                action="write", convert="big_endian", status="unknown")
            if (size(points,2) == 2) then
                write(nunit) (real(points(i,1),dp), real(points(i,2),dp), real(0.0_rk,dp) , i = 1, np)
            elseif (size(points,2) == 3) then
                write(nunit) (real(points(i,1),dp), real(points(i,2),dp), real(points(i,3),dp) , i = 1, np)
            else
                error stop 'Invalid dimension for points.'
            end if
            close(nunit)

            open(newunit=nunit, file=filename, form='formatted', action='write', position='append')
            write(nunit,'(a,1x,g0,1x,g0)') 'CELLS', ne, n
            close(nunit)
            open(newunit=nunit, file=filename, position='append', access="stream", form="unformatted",&
                action="write", convert="big_endian", status="unknown")
            select case (nn)
            case (2)
                write(nunit)&
                    (2, elemConn(i,1)-1,elemConn(i,2)-1, i = 1, ne)
            case (4)
                write(nunit)&
                    (4, elemConn(i,1)-1,elemConn(i,2)-1,elemConn(i,4)-1,elemConn(i,3)-1, i = 1, ne)
            case (8)
                write(nunit)&
                    (8, elemConn(i,1)-1,elemConn(i,2)-1,elemConn(i,4)-1,elemConn(i,3)-1,&
                    elemConn(i,5)-1,elemConn(i,6)-1,elemConn(i,8)-1,elemConn(i,7)-1, i = 1, ne)
            case default
                error stop 'Invalid number of nodes per element.'
            end select
            close(nunit)

            open(newunit=nunit, file=filename, form='formatted', action='write', position='append')
            write(nunit,'(a,1x,g0)') 'CELL_TYPES', ne
            close(nunit)
            open(newunit=nunit, file=filename, position='append', access="stream", form="unformatted",&
                action="write", convert="big_endian", status="unknown")
            write(nunit) (vtkCellType, i=1, ne)
            close(nunit)

            if (present(point_data) .and. present(field_names)) then
            open(newunit=nunit, file=filename, form='formatted', action='write', position='append')
            write(nunit, '(a,1x,g0)') 'POINT_DATA', size(point_data,1)
            close(nunit)

            do i = 1, size(point_data,2)
                open(newunit=nunit, file=filename, form='formatted', action='write', position='append')
                write(nunit, '(a,1x,a,1x,a)') 'SCALARS', trim(field_names(i)), 'double'
                write(nunit, '(a)') 'LOOKUP_TABLE default'
                close(nunit)

                open(newunit=nunit, file=filename, position='append', access='stream', form='unformatted', &
                    action='write', convert='big_endian', status='unknown')
                write(nunit) (real(point_data(j,i),dp), j=1,size(point_data,1))
                close(nunit)
            end do
        end if

        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function solve(A, B) result(X)
        real(rk), intent(in), contiguous :: A(:,:), B(:,:)
        real(rk), allocatable :: X(:,:)

        integer :: n, m, i, j, k, p
        real(rk), allocatable :: L(:,:), Y(:,:)
        real(rk) :: sum

        p = size(A,1)
        n = size(A,2)
        m = size(B,2)

        if (p /= size(B,1)) error stop "solve: A and B row mismatch"

        allocate(L(n,n), Y(n,m), X(n,m), source=0.0_rk)

        do i = 1, n
            do j = 1, i
                sum = A(i,j)
                do k = 1, j-1
                    sum = sum - L(i,k) * L(j,k)
                end do
                if (i == j) then
                    if (sum <= 0.0_rk) error stop "solve: Matrix not positive definite"
                    L(i,j) = sqrt(sum)
                else
                    L(i,j) = sum / L(j,j)
                end if
            end do
        end do

        ! Forward substitution: L·Y = AtB
        do j = 1,m
            do i = 1, n
                sum = B(i,j)
                do k = 1, i-1
                    sum = sum - L(i,k) * Y(k,j)
                end do
                Y(i,j) = sum / L(i,i)
            end do
        end do

        ! Backward substitution: Lᵗ·X = Y
        do j = 1,m
            do i = n, 1, -1
                sum = Y(i,j)
                do k = i+1, n
                    sum = sum - L(k,i) * X(k,j)
                end do
                X(i,j) = sum / L(i,i)
            end do
        end do
    end function
    !===============================================================================

end module forcad_utils
