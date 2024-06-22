program fdm_test_surface

    use forcad, only: rk, nurbs_surface

    implicit none

    type(nurbs_surface) :: surface              !! Declare a NURBS surface object
    real(rk), allocatable :: Xc(:,:), Wc(:) !! Arrays for control points and weights
    real(rk) :: Xtp(2), tol, Xt(2), Xtm(2)
    real(rk), allocatable :: Tgc(:), dTgc(:,:), Tgcp(:), dTgcp(:,:), Tgcm(:), dTgcm(:,:), d2Tgc(:,:), d2Tgcp(:,:), d2Tgcm(:,:)
    real(rk), allocatable :: CFD(:,:), BFD(:,:), FFD(:,:), CFD2(:,:), BFD2(:,:), FFD2(:,:)
    integer :: i

    !-----------------------------------------------------------------------------
    ! Setting up the NURBS surface
    !-----------------------------------------------------------------------------
    Wc = [1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, 1.0_rk, 0.5_rk, 1.0_rk,&
        0.2_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.2_rk]
    call surface%set_tetragon(L=[5.0_rk, 8.0_rk], nc=[4,4], Wc=Wc)

    !-----------------------------------------------------------------------------
    ! Creating the NURBS surface
    !-----------------------------------------------------------------------------

    !> Generate the NURBS surface with a resolution of 20
    call surface%create(20, 20)

    !-----------------------------------------------------------------------------
    ! Finite Difference Derivative Test
    !-----------------------------------------------------------------------------

    allocate(CFD(16,2), BFD(16,2), FFD(16,2))
    allocate(CFD2(2*16,2), BFD2(2*16,2), FFD2(2*16,2))

    tol = 1.0e-6_rk

    Xt(1) = 0.5_rk
    Xt(2) = 0.3_rk
    call surface%derivative2(Xt=Xt, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)

    do i = 1, 2
        Xtm = Xt
        Xtm(i) = Xt(i) - tol
        call surface%derivative2(Xt=Xtm, d2Tgc=d2Tgcm, dTgc=dTgcm, Tgc=Tgcm)

        Xtp = Xt
        Xtp(i) = Xt(i) + tol
        call surface%derivative2(Xt=Xtp, d2Tgc=d2Tgcp, dTgc=dTgcp, Tgc=Tgcp)

        BFD(:,i) = (Tgc  - Tgcm)/tol
        CFD(:,i) = (Tgcp - Tgcm)/(2.0_rk*tol)
        FFD(:,i) = (Tgcp - Tgc )/tol

        BFD2(:,i) = reshape((dTgc  - dTgcm)/tol,          shape=[2*16])
        CFD2(:,i) = reshape((dTgcp - dTgcm)/(2.0_rk*tol), shape=[2*16])
        FFD2(:,i) = reshape((dTgcp - dTgc )/tol,          shape=[2*16])
    end do

    print *, 'Tolerance:', tol
    print *, 'Error BFD dTgc:  ', norm2(BFD - dTgc)
    print *, 'Error CFD dTgc:  ', norm2(CFD - dTgc)
    print *, 'Error FFD dTgc:  ', norm2(FFD - dTgc)
    print *, 'Error BFD d2Tgc: ', norm2(BFD2 - d2Tgc)
    print *, 'Error CFD d2Tgc: ', norm2(CFD2 - d2Tgc)
    print *, 'Error FFD d2Tgc: ', norm2(FFD2 - d2Tgc)

    !> Finalize the NURBS surface object
    call surface%finalize()
    deallocate(CFD, BFD, FFD, CFD2, BFD2, FFD2)



















    !-----------------------------------------------------------------------------
    ! Setting up the NURBS surface
    !-----------------------------------------------------------------------------
    call surface%set_tetragon(L=[5.0_rk, 8.0_rk], nc=[4,4])

    !-----------------------------------------------------------------------------
    ! Creating the NURBS surface
    !-----------------------------------------------------------------------------

    !> Generate the NURBS surface with a resolution of 20
    call surface%create(20, 20)

    !-----------------------------------------------------------------------------
    ! Finite Difference Derivative Test
    !-----------------------------------------------------------------------------

    allocate(CFD(16,2), BFD(16,2), FFD(16,2))
    allocate(CFD2(2*16,2), BFD2(2*16,2), FFD2(2*16,2))

    tol = 1.0e-6_rk

    Xt(1) = 0.5_rk
    Xt(2) = 0.3_rk
    call surface%derivative2(Xt=Xt, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)

    do i = 1, 2
        Xtm = Xt
        Xtm(i) = Xt(i) - tol
        call surface%derivative2(Xt=Xtm, d2Tgc=d2Tgcm, dTgc=dTgcm, Tgc=Tgcm)

        Xtp = Xt
        Xtp(i) = Xt(i) + tol
        call surface%derivative2(Xt=Xtp, d2Tgc=d2Tgcp, dTgc=dTgcp, Tgc=Tgcp)

        BFD(:,i) = (Tgc  - Tgcm)/tol
        CFD(:,i) = (Tgcp - Tgcm)/(2.0_rk*tol)
        FFD(:,i) = (Tgcp - Tgc )/tol

        BFD2(:,i) = reshape((dTgc  - dTgcm)/tol,          shape=[2*16])
        CFD2(:,i) = reshape((dTgcp - dTgcm)/(2.0_rk*tol), shape=[2*16])
        FFD2(:,i) = reshape((dTgcp - dTgc )/tol,          shape=[2*16])
    end do

    print *, 'Tolerance:', tol
    print *, 'Error BFD dTgc:  ', norm2(BFD - dTgc)
    print *, 'Error CFD dTgc:  ', norm2(CFD - dTgc)
    print *, 'Error FFD dTgc:  ', norm2(FFD - dTgc)
    print *, 'Error BFD d2Tgc: ', norm2(BFD2 - d2Tgc)
    print *, 'Error CFD d2Tgc: ', norm2(CFD2 - d2Tgc)
    print *, 'Error FFD d2Tgc: ', norm2(FFD2 - d2Tgc)

    !> Finalize the NURBS surface object
    call surface%finalize()
    deallocate(CFD, BFD, FFD, CFD2, BFD2, FFD2)

end program
