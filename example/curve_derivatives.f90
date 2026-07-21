!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Evaluate first, second, arbitrary-order, and active curve derivatives.
program curve_derivatives

    use forcad, only: rk, nurbs_curve

    implicit none

    real(rk), parameter :: center(3) = [0.0_rk, 0.0_rk, 0.0_rk]
    real(rk), parameter :: Xt = 0.37_rk

    type(nurbs_curve) :: curve
    integer :: first_active
    real(rk), allocatable :: basis(:), first_derivative(:), second_derivative(:)
    real(rk), allocatable :: third_derivative(:), first_derivative_active(:)

    call curve%set_circle(&
        center = center,&
        radius = 2.0_rk)
    call curve%derivative(&
        Xt   = Xt,&
        dTgc = first_derivative,&
        Tgc  = basis)
    call curve%derivative2(&
        Xt     = Xt,&
        d2Tgc  = second_derivative)
    call curve%derivative_order(&
        Xt    = Xt,&
        order = 3,&
        Dgc   = third_derivative)
    call curve%derivative_order_active(&
        Xt           = Xt,&
        order        = 1,&
        first_active = first_active,&
        Dgc          = first_derivative_active)

    write(*,"(a,l1)") "Valid curve             : ", curve%err%ok
    write(*,"(a,i0)") "First active control     : ", first_active
    write(*,"(a,es12.4)") "Basis partition error    : ", abs(sum(basis)-1.0_rk)
    write(*,"(a,es12.4)") "First derivative sum     : ", abs(sum(first_derivative))
    write(*,"(a,es12.4)") "Second derivative sum    : ", abs(sum(second_derivative))
    write(*,"(a,es12.4)") "Third derivative sum     : ", abs(sum(third_derivative))
    write(*,"(a,es12.4)") "Active derivative sum    : ", abs(sum(first_derivative_active))
    write(*,"(a,*(es12.4,1x))") "Active derivative values : ", first_derivative_active

    call curve%finalize()

end program curve_derivatives
