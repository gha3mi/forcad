!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Evaluate dense and locally active rational surface basis functions.
program surface_basis

    use forcad, only: rk, nurbs_surface

    implicit none

    integer, parameter :: derivative_order(2) = [0, 0]
    real(rk), parameter :: center(3) = [0.0_rk, 0.0_rk, 0.0_rk]
    real(rk), parameter :: Xt(2) = [0.37_rk, 0.43_rk]

    type(nurbs_surface) :: surface
    integer :: first_active(2)
    real(rk), allocatable :: basis(:), basis_active(:)

    call surface%set_ring(&
        center  = center,&
        radius1 = 1.0_rk,&
        radius2 = 2.0_rk)
    call surface%basis(&
        Xt  = Xt,&
        Tgc = basis)
    call surface%derivative_order_active(&
        Xt           = Xt,&
        order        = derivative_order,&
        first_active = first_active,&
        Dgc          = basis_active)

    write(*,"(a,l1)") "Valid surface        : ", surface%err%ok
    write(*,"(a,i0)") "Dense basis count     : ", size(basis)
    write(*,"(a,i0)") "Active basis count    : ", size(basis_active)
    write(*,"(a,2(i0,1x))") "First active controls : ", first_active
    write(*,"(a,es12.4)") "Dense partition error : ", abs(sum(basis)-1.0_rk)
    write(*,"(a,es12.4)") "Active partition error: ", abs(sum(basis_active)-1.0_rk)
    write(*,"(a,*(es12.4,1x))") "Active basis values   : ", basis_active

    call surface%finalize()

end program surface_basis
