!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Evaluate an element-local surface ansatz at one tensor quadrature point.
program surface_ansatz

    use forcad, only: rk, nurbs_surface

    implicit none

    integer, parameter :: ngauss(2) = [4, 4]
    real(rk), parameter :: center(3) = [0.0_rk, 0.0_rk, 0.0_rk]

    type(nurbs_surface) :: surface
    integer :: element_id
    integer, allocatable :: elem(:,:)
    real(rk) :: dA
    real(rk), allocatable :: basis(:), gradient(:,:), hessian(:,:,:)

    call surface%set_ring(&
        center  = center,&
        radius1 = 1.0_rk,&
        radius2 = 2.0_rk)
    elem = surface%cmp_elem()
    element_id = max(1, (size(elem,1)+1)/2)
    call surface%ansatz(&
        ie          = element_id,&
        ig          = 6,&
        Tgc         = basis,&
        dTgc_dXg    = gradient,&
        d2Tgc_dXg2  = hessian,&
        dA          = dA,&
        ngauss      = ngauss)

    write(*,"(a,l1)") "Valid surface           : ", surface%err%ok
    write(*,"(a,i0)") "Element                  : ", element_id
    write(*,"(a,i0)") "Local basis count        : ", size(basis)
    write(*,"(a,es12.4)") "Basis partition error    : ", abs(sum(basis)-1.0_rk)
    write(*,"(a,es12.4)") "Gradient partition error : ", maxval(abs(sum(gradient,dim=1)))
    write(*,"(a,es12.4)") "Hessian partition error  : ", maxval(abs(sum(hessian,dim=1)))
    write(*,"(a,es12.4)") "Weighted area measure    : ", dA

    call surface%finalize()

end program surface_ansatz
