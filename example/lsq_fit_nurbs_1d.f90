program lsq_fit_nurbs_1d

   use forcad, only: rk, nurbs_curve

   implicit none

   type(nurbs_curve) :: nrb
   integer  :: n, i
   real(rk), parameter :: pi = acos(-1.0_rk)
   real(rk), allocatable :: Xdata(:,:)
   real(rk), allocatable :: Xt(:)
   real(rk), allocatable :: Xg_eval(:,:)
   real(rk) :: err1, err2, err3, rms

   n = 100

   ! create parametric grid points
   allocate(Xt(n))
   do concurrent (i = 1: n)
      Xt(i) = real(i-1, rk) / real(n-1, rk)
   end do

   ! data points to be fitted
   allocate(Xdata(n, 3))
   do concurrent (i = 1: n)
      Xdata(i,1) = 0.0_rk + 1.0_rk * cos(2.0_rk*pi * Xt(i) )
      Xdata(i,2) = 0.0_rk + 1.0_rk * sin(2.0_rk*pi * Xt(i) )
      Xdata(i,3) = 0.0_rk
   end do

   call nrb%set(&
      degree     = 7,&
      Xth_dir    = [0.0_rk, 0.1_rk, 0.2_rk, 0.3_rk, 0.4_rk, 0.5_rk, 0.6_rk, 0.7_rk, 0.8_rk, 0.9_rk, 1.0_rk],&
      continuity = [ -1   ,   2   ,   2   ,   2   ,   2   ,   2   ,   2   ,   2   ,   2   ,   2   ,  -1   ])

   print'(a)', "========================================"
   print'(a)', "NURBS Curve Configuration"
   print'(a)', "----------------------------------------"
   print'(a,i0,a)', "Degrees    : ", nrb%get_degree()
   print'(a,i0,a)', "Control pts: ", nrb%get_nc()
   print'(a,i0,a)', "Data grid  : ", n

   print'(a)', "----------------------------------------"
   print'(a)', "Continuity"
   print'(a,*(i3,1x))', "  dir1:", nrb%get_continuity()

   print'(a)', "----------------------------------------"
   print'(a)', "Knot vectors"
   print'(a,*(f5.2,1x))', "  dir1:", nrb%get_knot()
   print'(a)', "========================================"

   print'(a)', "Fitting least squares curve..."
   call nrb%lsq_fit_nurbs(&
      Xt        = Xt,     &
      Xdata     = Xdata,  &
      ndata     = n,      &
      maxit     = 100,    &
      tol       = sqrt(epsilon(0.0_rk)), &
      lambda_xc = sqrt(epsilon(0.0_rk)), &
      reg_logw  = sqrt(epsilon(0.0_rk)) )
   print'(a)', "Fitting complete."

   ! create NURBS curve
   call nrb%create(Xt=Xt)
   Xg_eval = nrb%get_Xg()

   ! Compute errors
   err1 = norm2(Xg_eval(:,1) - Xdata(:,1)) / max( norm2(Xdata(:,1)), epsilon(0.0_rk) )
   err2 = norm2(Xg_eval(:,2) - Xdata(:,2)) / max( norm2(Xdata(:,2)), epsilon(0.0_rk) )
   err3 = norm2(Xg_eval(:,3) - Xdata(:,3)) / max( norm2(Xdata(:,3)), epsilon(0.0_rk) )
   rms  = sqrt((err1**2 + err2**2 + err3**2)/3.0_rk)

   print'(a)', "========================================"
   print'(a)', "Fitting Error Report"
   print'(a)', "----------------------------------------"
   print'(a,e13.6)', "Rel. error (dir1):", err1
   print'(a,e13.6)', "Rel. error (dir2):", err2
   print'(a,e13.6)', "Rel. error (dir3):", err3
   print'(a,e13.6)', "Total RMS error  :", rms
   print'(a)', "========================================"

   ! Export results
   call nrb%export_Xc("vtk/lsq_fit_bspline_1d_Xc.vtk")
   call nrb%export_Xg("vtk/lsq_fit_bspline_1d_Xg.vtk")
   call nrb%show("vtk/lsq_fit_bspline_1d_Xc.vtk", "vtk/lsq_fit_bspline_1d_Xg.vtk")

end program
