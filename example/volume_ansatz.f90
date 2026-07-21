!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Evaluate an element-local volume ansatz at one tensor quadrature point.
program volume_ansatz

    use forcad, only: rk, nurbs_volume

    implicit none

    integer, parameter :: ngauss(3) = [3, 3, 3]
    real(rk), parameter :: center(3) = [0.0_rk, 0.0_rk, 0.0_rk]

    type(nurbs_volume) :: volume
    integer :: element_id
    integer, allocatable :: elem(:,:)
    real(rk) :: dV
    real(rk), allocatable :: basis(:), gradient(:,:), hessian(:,:,:)

    call volume%set_ring(&
        center  = center,&
        radius1 = 1.0_rk,&
        radius2 = 2.0_rk,&
        length  = 1.0_rk)
    elem = volume%cmp_elem()
    element_id = max(1, (size(elem,1)+1)/2)
    call volume%ansatz(&
        ie          = element_id,&
        ig          = 14,&
        Tgc         = basis,&
        dTgc_dXg    = gradient,&
        d2Tgc_dXg2  = hessian,&
        dV          = dV,&
        ngauss      = ngauss)

    write(*,"(a,l1)") "Valid volume            : ", volume%err%ok
    write(*,"(a,i0)") "Element                  : ", element_id
    write(*,"(a,i0)") "Local basis count        : ", size(basis)
    write(*,"(a,es12.4)") "Basis partition error    : ", abs(sum(basis)-1.0_rk)
    write(*,"(a,es12.4)") "Gradient partition error : ", maxval(abs(sum(gradient,dim=1)))
    write(*,"(a,es12.4)") "Hessian partition error  : ", maxval(abs(sum(hessian,dim=1)))
    write(*,"(a,es12.4)") "Weighted volume measure  : ", dV

    call volume%finalize()

end program volume_ansatz
