!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Evaluate an element-local curve ansatz at one quadrature point.
program curve_ansatz

    use forcad, only: rk, nurbs_curve

    implicit none

    integer, parameter :: ngauss = 5
    real(rk), parameter :: center(3) = [0.0_rk, 0.0_rk, 0.0_rk]

    type(nurbs_curve) :: curve
    integer :: element_id
    integer, allocatable :: elem(:,:)
    real(rk) :: dL
    real(rk), allocatable :: basis(:), gradient(:,:), hessian(:)

    call curve%set_circle(&
        center = center,&
        radius = 2.0_rk)
    elem = curve%cmp_elem()
    element_id = max(1, (size(elem,1)+1)/2)
    call curve%ansatz(&
        ie          = element_id,&
        ig          = 3,&
        Tgc         = basis,&
        dTgc_dXg    = gradient,&
        d2Tgc_dXg2  = hessian,&
        dL          = dL,&
        ngauss      = ngauss)

    write(*,"(a,l1)") "Valid curve             : ", curve%err%ok
    write(*,"(a,i0)") "Element                  : ", element_id
    write(*,"(a,i0)") "Local basis count        : ", size(basis)
    write(*,"(a,es12.4)") "Basis partition error    : ", abs(sum(basis)-1.0_rk)
    write(*,"(a,es12.4)") "Gradient partition error : ", maxval(abs(sum(gradient,dim=1)))
    write(*,"(a,es12.4)") "Hessian partition error  : ", abs(sum(hessian))
    write(*,"(a,es12.4)") "Weighted line measure    : ", dL

    call curve%finalize()

end program curve_ansatz
