program fdm_test_volume

    use forcad, only: rk, nurbs_volume

    implicit none

    type(nurbs_volume) :: volume            !! Declare a NURBS volume object
    real(rk), allocatable :: Wc(:)          !! Weights for the control points
    real(rk) :: Xt(3), tol, Xtm(3), Xtp(3)
    real(rk), allocatable :: Tgc(:), dTgc(:,:), Tgcp(:), dTgcp(:,:), Tgcm(:), dTgcm(:,:), d2Tgc(:,:), d2Tgcp(:,:), d2Tgcm(:,:)
    real(rk), allocatable :: CFD(:,:), BFD(:,:), FFD(:,:), CFD2(:,:), BFD2(:,:), FFD2(:,:)
    integer :: i

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS volume
    !-----------------------------------------------------------------------------
    allocate(Wc(64))
    Wc = 1.0_rk
    Wc(10) = 0.5_rk
    call volume%set_hexahedron(L=[2.0_rk, 4.0_rk, 8.0_rk], nc=[4,4,4], Wc=Wc)

    !-----------------------------------------------------------------------------
    ! Creating the NURBS volume
    !-----------------------------------------------------------------------------

    !> Generate the NURBS volume with a resolution of 20
    call volume%create(20, 20, 20)

    !-----------------------------------------------------------------------------
    ! Finite Difference Derivative Test
    !-----------------------------------------------------------------------------

    allocate(CFD(64,3), BFD(64,3), FFD(64,3))
    allocate(CFD2(3*64,3), BFD2(3*64,3), FFD2(3*64,3))

    tol = 1.0e-6_rk

    Xt(1) = 0.5_rk
    Xt(2) = 0.3_rk
    Xt(3) = 0.7_rk
    call volume%derivative2(Xt=Xt, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
    do i = 1, 3
        Xtm = Xt
        Xtm(i) = Xt(i) - tol
        call volume%derivative2(Xt=Xtm, d2Tgc=d2Tgcm, dTgc=dTgcm, Tgc=Tgcm)

        Xtp = Xt
        Xtp(i) = Xt(i) + tol
        call volume%derivative2(Xt=Xtp, d2Tgc=d2Tgcp, dTgc=dTgcp, Tgc=Tgcp)

        BFD(:,i) = (Tgc  - Tgcm)/tol
        CFD(:,i) = (Tgcp - Tgcm)/(2.0_rk*tol)
        FFD(:,i) = (Tgcp - Tgc )/tol

        BFD2(:,i) = reshape((dTgc  - dTgcm)/tol,          shape=[3*64])
        CFD2(:,i) = reshape((dTgcp - dTgcm)/(2.0_rk*tol), shape=[3*64])
        FFD2(:,i) = reshape((dTgcp - dTgc )/tol,          shape=[3*64])
    end do

    print *, 'Tolerance:', tol
    print *, 'Error BFD dTgc:  ', norm2(BFD - dTgc)
    print *, 'Error CFD dTgc:  ', norm2(CFD - dTgc)
    print *, 'Error FFD dTgc:  ', norm2(FFD - dTgc)
    print *, 'Error BFD d2Tgc: ', norm2(BFD2 - d2Tgc)
    print *, 'Error CFD d2Tgc: ', norm2(CFD2 - d2Tgc)
    print *, 'Error FFD d2Tgc: ', norm2(FFD2 - d2Tgc)

    !> Finalize the NURBS volume object
    call volume%finalize()
    deallocate(CFD, BFD, FFD, CFD2, BFD2, FFD2)
    deallocate(Wc)


















    !-----------------------------------------------------------------------------
    ! Setting up the NURBS volume
    !-----------------------------------------------------------------------------
    call volume%set_hexahedron(L=[2.0_rk, 4.0_rk, 8.0_rk], nc=[4,4,4])

    !-----------------------------------------------------------------------------
    ! Creating the NURBS volume
    !-----------------------------------------------------------------------------

    !> Generate the NURBS volume with a resolution of 20
    call volume%create(20, 20, 20)

    !-----------------------------------------------------------------------------
    ! Finite Difference Derivative Test
    !-----------------------------------------------------------------------------

    allocate(CFD(64,3), BFD(64,3), FFD(64,3))
    allocate(CFD2(3*64,3), BFD2(3*64,3), FFD2(3*64,3))

    tol = 1.0e-6_rk

    Xt(1) = 0.5_rk
    Xt(2) = 0.3_rk
    Xt(3) = 0.7_rk
    call volume%derivative2(Xt=Xt, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
    do i = 1, 3
        Xtm = Xt
        Xtm(i) = Xt(i) - tol
        call volume%derivative2(Xt=Xtm, d2Tgc=d2Tgcm, dTgc=dTgcm, Tgc=Tgcm)

        Xtp = Xt
        Xtp(i) = Xt(i) + tol
        call volume%derivative2(Xt=Xtp, d2Tgc=d2Tgcp, dTgc=dTgcp, Tgc=Tgcp)

        BFD(:,i) = (Tgc  - Tgcm)/tol
        CFD(:,i) = (Tgcp - Tgcm)/(2.0_rk*tol)
        FFD(:,i) = (Tgcp - Tgc )/tol

        BFD2(:,i) = reshape((dTgc  - dTgcm)/tol,          shape=[3*64])
        CFD2(:,i) = reshape((dTgcp - dTgcm)/(2.0_rk*tol), shape=[3*64])
        FFD2(:,i) = reshape((dTgcp - dTgc )/tol,          shape=[3*64])
    end do

    print *, 'Tolerance:', tol
    print *, 'Error BFD dTgc:  ', norm2(BFD - dTgc)
    print *, 'Error CFD dTgc:  ', norm2(CFD - dTgc)
    print *, 'Error FFD dTgc:  ', norm2(FFD - dTgc)
    print *, 'Error BFD d2Tgc: ', norm2(BFD2 - d2Tgc)
    print *, 'Error CFD d2Tgc: ', norm2(CFD2 - d2Tgc)
    print *, 'Error FFD d2Tgc: ', norm2(FFD2 - d2Tgc)

    !> Finalize the NURBS volume object
    call volume%finalize()
    deallocate(CFD, BFD, FFD, CFD2, BFD2, FFD2)

end program
