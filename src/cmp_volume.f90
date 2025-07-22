program compute_volume

    use forcad, only: rk, nurbs_volume

    implicit none

    type(nurbs_volume) :: shape
    real(rk) :: volume
    real(rk) :: Xc(8,3)

    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [2.0_rk, 0.0_rk, 0.0_rk]
    Xc(3,:) = [0.0_rk, 2.0_rk, 0.0_rk]
    Xc(4,:) = [2.0_rk, 2.0_rk, 0.0_rk]
    Xc(5,:) = [0.0_rk, 0.0_rk, 2.0_rk]
    Xc(6,:) = [2.0_rk, 0.0_rk, 2.0_rk]
    Xc(7,:) = [0.0_rk, 2.0_rk, 2.0_rk]
    Xc(8,:) = [2.0_rk, 2.0_rk, 2.0_rk]

    call shape%set(&
        knot1=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
        knot2=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
        knot3=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
        Xc=Xc)

    call shape%cmp_volume(volume)
    print*, volume
end program
