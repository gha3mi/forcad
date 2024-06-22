program fdm_test_curve

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: curve              !! Declare a NURBS curve object
    real(rk), allocatable :: Xc(:,:), Wc(:) !! Arrays for control points and weights
    real(rk) :: knot(6)                     !! Array for knot vector
    real(rk) :: Xtp, tol, Xt, Xtm
    real(rk), allocatable :: Tgc(:), dTgc(:), Tgcp(:), dTgcp(:), Tgcm(:), dTgcm(:), CFD(:), BFD(:), FFD(:), d2Tgc(:), d2Tgcp(:), d2Tgcm(:)

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS curve
    !-----------------------------------------------------------------------------

    !> Define control points for the NURBS curve
    allocate(Xc(3, 3))
    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [0.0_rk, 5.0_rk, 0.0_rk]
    Xc(3,:) = [5.0_rk, 5.0_rk, 0.0_rk]

    !> Define weights for the control points (optional)
    allocate(Wc(3))
    Wc = [1.0_rk, 1.1_rk, 1.0_rk]

    !> Define knot vector
    knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    !> Set knot vector, control points, and weights for the NURBS curve object.
    !> Wc is optional
    call curve%set(knot, Xc, Wc)

    !-----------------------------------------------------------------------------
    ! Creating the NURBS curve
    !-----------------------------------------------------------------------------

    !> Generate the NURBS curve with a resolution of 20
    call curve%create(res = 20)

    !-----------------------------------------------------------------------------
    ! Finite Difference Derivative Test
    !-----------------------------------------------------------------------------

    tol = 1.0e-6_rk
    Xt = 0.5_rk
    call curve%derivative2(Xt=Xt, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
    Xtm = Xt - tol
    call curve%derivative2(Xt=Xtm, d2Tgc=d2Tgcm, dTgc=dTgcm, Tgc=Tgcm)
    Xtp = Xt + tol
    call curve%derivative2(Xt=Xtp, d2Tgc=d2Tgcp, dTgc=dTgcp, Tgc=Tgcp)
    print *, 'Tolerance:', tol
    print *, 'Error BFD dTgc:  ', norm2((Tgc  - Tgcm)/tol          - dTgc)
    print *, 'Error CFD dTgc:  ', norm2((Tgcp - Tgcm)/(2.0_rk*tol) - dTgc)
    print *, 'Error FFD dTgc:  ', norm2((Tgcp - Tgc )/tol          - dTgc)
    print *, 'Error BFD d2Tgc: ', norm2((dTgc  - dTgcm)/tol          - d2Tgc)
    print *, 'Error CFD d2Tgc: ', norm2((dTgcp - dTgcm)/(2.0_rk*tol) - d2Tgc)
    print *, 'Error FFD d2Tgc: ', norm2((dTgcp - dTgc )/tol          - d2Tgc)

    !> Finalize the NURBS curve object
    call curve%finalize()
    deallocate(Xc, Wc)









    !-----------------------------------------------------------------------------
    ! Setting up the NURBS curve
    !-----------------------------------------------------------------------------

    !> Define control points for the NURBS curve
    allocate(Xc(3, 3))
    Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    Xc(2,:) = [0.0_rk, 5.0_rk, 0.0_rk]
    Xc(3,:) = [5.0_rk, 5.0_rk, 0.0_rk]

    !> Define knot vector
    knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

    !> Set knot vector, control points for the NURBS curve object.
    !> Wc is optional
    call curve%set(knot, Xc)

    !-----------------------------------------------------------------------------
    ! Creating the NURBS curve
    !-----------------------------------------------------------------------------

    !> Generate the NURBS curve with a resolution of 20
    call curve%create(res = 20)

    !-----------------------------------------------------------------------------
    ! Finite Difference Derivative Test
    !-----------------------------------------------------------------------------

    tol = 1.0e-6_rk
    Xt = 0.5_rk
    call curve%derivative2(Xt=Xt, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
    Xtm = Xt - tol
    call curve%derivative2(Xt=Xtm, d2Tgc=d2Tgcm, dTgc=dTgcm, Tgc=Tgcm)
    Xtp = Xt + tol
    call curve%derivative2(Xt=Xtp, d2Tgc=d2Tgcp, dTgc=dTgcp, Tgc=Tgcp)
    print *, 'Tolerance:', tol
    print *, 'Error BFD dTgc:  ', norm2((Tgc  - Tgcm)/tol          - dTgc)
    print *, 'Error CFD dTgc:  ', norm2((Tgcp - Tgcm)/(2.0_rk*tol) - dTgc)
    print *, 'Error FFD dTgc:  ', norm2((Tgcp - Tgc )/tol          - dTgc)
    print *, 'Error BFD d2Tgc: ', norm2((dTgc  - dTgcm)/tol          - d2Tgc)
    print *, 'Error CFD d2Tgc: ', norm2((dTgcp - dTgcm)/(2.0_rk*tol) - d2Tgc)
    print *, 'Error FFD d2Tgc: ', norm2((dTgcp - dTgc )/tol          - d2Tgc)

    !> Finalize the NURBS curve object
    call curve%finalize()
    deallocate(Xc)
end program
