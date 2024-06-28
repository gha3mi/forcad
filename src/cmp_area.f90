program compute_area

    use forcad

    implicit none

    type(nurbs_surface) :: shape
    real(rk) :: area
    real(rk) :: Xc(4,3)

    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [2.0_rk, 0.0_rk, 0.0_rk]
    Xc(3,:) = [0.0_rk, 2.0_rk, 0.0_rk]
    Xc(4,:) = [2.0_rk, 2.0_rk, 0.0_rk]

    call shape%set(&
        knot1=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
        knot2=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
        Xc=Xc)

    call shape%cmp_area(area)
    print*, area
end program
