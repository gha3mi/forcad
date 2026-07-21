!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program surface_fit_nurbs

   use forcad, only: rk, nurbs_surface, ndgrid

   implicit none

   type(nurbs_surface) :: nrb
   integer :: n(2), ndata, i
   real(rk), parameter :: pi = acos(-1.0_rk)
   real(rk), allocatable :: Xdata(:,:)
   real(rk), allocatable :: Xt1(:), Xt2(:), Xt(:,:)
   real(rk), allocatable :: Xg_eval(:,:)
   real(rk) :: err1, err2, err3, rms

   n = [60, 30]

   ! create parametric grid points
   allocate(Xt1(n(1)), Xt2(n(2)))
   do concurrent (i = 1: n(1))
      Xt1(i) = real(i-1, rk)/real(n(1)-1, rk)
   end do
   do concurrent (i = 1: n(2))
      Xt2(i) = real(i-1, rk)/real(n(2)-1, rk)
   end do
   call ndgrid(Xt1, Xt2, Xt)

   ! data points to be fitted
   ndata = n(1)*n(2)
   allocate(Xdata(ndata, 3))
   do i = 1, ndata
      Xdata(i,1) = 0.0_rk + (1.0_rk + 0.35_rk*cos(2.0_rk*pi*Xt(i,2))) * cos(2.0_rk*pi*Xt(i,1))
      Xdata(i,2) = 0.0_rk + (1.0_rk + 0.35_rk*cos(2.0_rk*pi*Xt(i,2))) * sin(2.0_rk*pi*Xt(i,1))
      Xdata(i,3) = 0.0_rk +  0.35_rk*sin(2.0_rk*pi*Xt(i,2))
   end do

   ! set up NURBS surface
   call nrb%set(&
      degree      = [4, 4],&
      Xth_dir1    = [0.0_rk, 0.1_rk, 0.2_rk, 0.3_rk, 0.4_rk, 0.5_rk, 0.6_rk, 0.7_rk, 0.8_rk, 0.9_rk, 1.0_rk],&
      Xth_dir2    = [0.0_rk, 0.2_rk, 0.4_rk, 0.6_rk, 0.8_rk, 1.0_rk],&
      continuity1 = [ -1   ,   1   ,   1   ,   1   ,   1   ,   1   ,   1   ,   1   ,   1   ,   1   ,  -1   ],&
      continuity2 = [ -1   ,   1   ,   1   ,   1   ,   1   ,  -1   ])

   write(*,"(a)") "========================================"
   write(*,"(a)") "NURBS Surface Configuration"
   write(*,"(a)") "----------------------------------------"
   write(*,"(a,2(i0,a))") "Degrees    : ", nrb%get_degree(1), ", ", nrb%get_degree(2)
   write(*,"(a,2(i0,a))") "Control pts: ", nrb%get_nc(1), " x ", nrb%get_nc(2)
   write(*,"(a,2(i0,a))") "Data grid  : ", n(1), " x ", n(2)

   write(*,"(a)") "----------------------------------------"
   write(*,"(a)") "Continuity"
   write(*,"(a,*(i3,1x))") "  dir1:", nrb%get_continuity(1)
   write(*,"(a,*(i3,1x))") "  dir2:", nrb%get_continuity(2)

   write(*,"(a)") "----------------------------------------"
   write(*,"(a)") "Knot vectors"
   write(*,"(a,*(f5.2,1x))") "  dir1:", nrb%get_knot(1)
   write(*,"(a,*(f5.2,1x))") "  dir2:", nrb%get_knot(2)
   write(*,"(a)") "========================================"

   write(*,"(a)") "Fitting least squares surface..."
   call nrb%lsq_fit_nurbs(&
      Xt        = Xt,     &
      Xdata     = Xdata,  &
      ndata     = n,      &
      maxit     = 100,    &
      tol       = sqrt(epsilon(0.0_rk)), &
      lambda_xc = sqrt(epsilon(0.0_rk)), &
      reg_logw  = sqrt(epsilon(0.0_rk)) )
   write(*,"(a)") "Fitting complete."

   ! create NURBS surface
   call nrb%create(Xt1=Xt1, Xt2=Xt2)
   Xg_eval = nrb%get_Xg()

   ! Compute errors
   err1 = norm2(Xg_eval(:,1) - Xdata(:,1)) / max( norm2(Xdata(:,1)), epsilon(0.0_rk) )
   err2 = norm2(Xg_eval(:,2) - Xdata(:,2)) / max( norm2(Xdata(:,2)), epsilon(0.0_rk) )
   err3 = norm2(Xg_eval(:,3) - Xdata(:,3)) / max( norm2(Xdata(:,3)), epsilon(0.0_rk) )
   rms  = sqrt((err1**2 + err2**2 + err3**2)/3.0_rk)

   write(*,"(a)") "========================================"
   write(*,"(a)") "Fitting Error Report"
   write(*,"(a)") "----------------------------------------"
   write(*,"(a,e13.6)") "Rel. error (dir1):", err1
   write(*,"(a,e13.6)") "Rel. error (dir2):", err2
   write(*,"(a,e13.6)") "Rel. error (dir3):", err3
   write(*,"(a,e13.6)") "Total RMS error  :", rms
   write(*,"(a)") "========================================"

   ! Export results
   call nrb%export_Xc("vtk/surface_fit_nurbs_Xc.vtk")
   call nrb%export_Xg("vtk/surface_fit_nurbs_Xg.vtk")
   call nrb%export_Xth_in_Xg("vtk/surface_fit_nurbs_Xth.vtk", res=20)
   call nrb%show("vtk/surface_fit_nurbs_Xc.vtk", "vtk/surface_fit_nurbs_Xg.vtk", "vtk/surface_fit_nurbs_Xth.vtk")

end program surface_fit_nurbs
