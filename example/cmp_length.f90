program compute_length

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: shape
    real(rk) :: length
    real(rk) :: Xc(2,3)

    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [2.0_rk, 0.0_rk, 0.0_rk]

    call shape%set(&
        knot=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
        Xc=Xc)

    call shape%cmp_length(length)
    print*, length
end program
