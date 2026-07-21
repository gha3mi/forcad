!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Evaluate pure, mixed, and locally active volume derivatives.
program volume_derivatives

    use forcad, only: rk, nurbs_volume

    implicit none

    integer, parameter :: mixed_order(3) = [1, 1, 1]
    real(rk), parameter :: center(3) = [0.0_rk, 0.0_rk, 0.0_rk]
    real(rk), parameter :: Xt(3) = [0.37_rk, 0.43_rk, 0.61_rk]

    type(nurbs_volume) :: volume
    integer :: first_active(3)
    real(rk), allocatable :: basis(:), first_derivative(:,:), second_derivative(:,:)
    real(rk), allocatable :: mixed_derivative(:), mixed_derivative_active(:)

    call volume%set_ring(&
        center  = center,&
        radius1 = 1.0_rk,&
        radius2 = 2.0_rk,&
        length  = 1.0_rk)
    call volume%derivative(&
        Xt   = Xt,&
        dTgc = first_derivative,&
        Tgc  = basis)
    call volume%derivative2(&
        Xt     = Xt,&
        d2Tgc  = second_derivative)
    call volume%derivative_order(&
        Xt    = Xt,&
        order = mixed_order,&
        Dgc   = mixed_derivative)
    call volume%derivative_order_active(&
        Xt           = Xt,&
        order        = mixed_order,&
        first_active = first_active,&
        Dgc          = mixed_derivative_active)

    write(*,"(a,l1)") "Valid volume            : ", volume%err%ok
    write(*,"(a,3(i0,1x))") "First active controls    : ", first_active
    write(*,"(a,es12.4)") "Basis partition error    : ", abs(sum(basis)-1.0_rk)
    write(*,"(a,es12.4)") "First derivative sum     : ", maxval(abs(sum(first_derivative,dim=1)))
    write(*,"(a,es12.4)") "Pure second derivative   : ", maxval(abs(sum(second_derivative,dim=1)))
    write(*,"(a,es12.4)") "Mixed derivative sum     : ", abs(sum(mixed_derivative))
    write(*,"(a,es12.4)") "Active mixed sum         : ", abs(sum(mixed_derivative_active))
    write(*,"(a,*(es12.4,1x))") "Active mixed values      : ", mixed_derivative_active

    call volume%finalize()

end program volume_derivatives
