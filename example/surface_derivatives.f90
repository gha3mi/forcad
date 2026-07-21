!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Evaluate pure, mixed, and locally active surface derivatives.
program surface_derivatives

    use forcad, only: rk, nurbs_surface

    implicit none

    integer, parameter :: mixed_order(2) = [1, 1]
    real(rk), parameter :: center(3) = [0.0_rk, 0.0_rk, 0.0_rk]
    real(rk), parameter :: Xt(2) = [0.37_rk, 0.43_rk]

    type(nurbs_surface) :: surface
    integer :: first_active(2)
    real(rk), allocatable :: basis(:), first_derivative(:,:), second_derivative(:,:)
    real(rk), allocatable :: mixed_derivative(:), mixed_derivative_active(:)

    call surface%set_ring(&
        center  = center,&
        radius1 = 1.0_rk,&
        radius2 = 2.0_rk)
    call surface%derivative(&
        Xt   = Xt,&
        dTgc = first_derivative,&
        Tgc  = basis)
    call surface%derivative2(&
        Xt     = Xt,&
        d2Tgc  = second_derivative)
    call surface%derivative_order(&
        Xt    = Xt,&
        order = mixed_order,&
        Dgc   = mixed_derivative)
    call surface%derivative_order_active(&
        Xt           = Xt,&
        order        = mixed_order,&
        first_active = first_active,&
        Dgc          = mixed_derivative_active)

    write(*,"(a,l1)") "Valid surface           : ", surface%err%ok
    write(*,"(a,2(i0,1x))") "First active controls    : ", first_active
    write(*,"(a,es12.4)") "Basis partition error    : ", abs(sum(basis)-1.0_rk)
    write(*,"(a,es12.4)") "First derivative sum     : ", maxval(abs(sum(first_derivative,dim=1)))
    write(*,"(a,es12.4)") "Pure second derivative   : ", maxval(abs(sum(second_derivative,dim=1)))
    write(*,"(a,es12.4)") "Mixed derivative sum     : ", abs(sum(mixed_derivative))
    write(*,"(a,es12.4)") "Active mixed sum         : ", abs(sum(mixed_derivative_active))
    write(*,"(a,*(es12.4,1x))") "Active mixed values      : ", mixed_derivative_active

    call surface%finalize()

end program surface_derivatives
