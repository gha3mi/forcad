!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Evaluate dense and locally active rational curve basis functions.
program curve_basis

    use forcad, only: rk, nurbs_curve

    implicit none

    real(rk), parameter :: center(3) = [0.0_rk, 0.0_rk, 0.0_rk]
    real(rk), parameter :: Xt = 0.37_rk

    type(nurbs_curve) :: curve
    integer :: first_active
    real(rk), allocatable :: basis(:), basis_active(:)

    call curve%set_circle(&
        center = center,&
        radius = 2.0_rk)
    call curve%basis(&
        Xt  = Xt,&
        Tgc = basis)
    call curve%derivative_order_active(&
        Xt           = Xt,&
        order        = 0,&
        first_active = first_active,&
        Dgc          = basis_active)

    write(*,"(a,l1)") "Valid curve          : ", curve%err%ok
    write(*,"(a,i0)") "Dense basis count     : ", size(basis)
    write(*,"(a,i0)") "Active basis count    : ", size(basis_active)
    write(*,"(a,i0)") "First active control  : ", first_active
    write(*,"(a,es12.4)") "Dense partition error : ", abs(sum(basis)-1.0_rk)
    write(*,"(a,es12.4)") "Active partition error: ", abs(sum(basis_active)-1.0_rk)
    write(*,"(a,*(es12.4,1x))") "Active basis values   : ", basis_active

    call curve%finalize()

end program curve_basis
