program lsq_fit_bspline_1d

   use forcad, only: rk, nurbs_curve

   implicit none

   type(nurbs_curve) :: bsp
   integer :: n, i
   real(rk), parameter :: pi = acos(-1.0_rk)
   real(rk), allocatable :: Xdata(:,:)
   real(rk), allocatable :: Xt(:)
   real(rk), allocatable :: Xg_eval(:,:)
   real(rk) :: err1, err2, err3, rms

   n = 42

   ! create parametric grid points
   allocate(Xt(n))
   do concurrent (i = 1: n)
      Xt(i) = real(i-1, rk)/real(n-1, rk)
   end do

   ! data points to be fitted
   allocate(Xdata(n, 3))
   do i = 1, n
      Xdata(i,1) = Xt(i)
      Xdata(i,2) = 0.3_rk * sin(4.0_rk * pi * Xt(i))
      Xdata(i,3) = 0.3_rk * cos(4.0_rk * pi * Xt(i))
   end do

   ! set up B-Spline curve
   ! Xth_dir(1) = minval(Xt), Xth_dir(2) = maxval(Xt)
   call bsp%set(&
      degree     = 5,&
      Xth_dir    = [0.0_rk, 0.1_rk, 0.2_rk, 0.3_rk, 0.4_rk, 0.5_rk, 0.6_rk, 0.7_rk, 0.8_rk, 0.9_rk, 1.0_rk],&
      continuity = [ -1   ,   1   ,   1   ,   1   ,   1   ,   1   ,   1   ,   1   ,   1   ,   1   ,  -1   ])

   print'(a)', "========================================"
   print'(a)', "B-Spline Curve Configuration"
   print'(a)', "----------------------------------------"
   print'(a,i0,a)', "Degrees    : ", bsp%get_degree()
   print'(a,i0,a)', "Control pts: ", bsp%get_nc()
   print'(a,i0,a)', "Data grid  : ", n

   print'(a)', "----------------------------------------"
   print'(a)', "Continuity"
   print'(a,*(i3,1x))', "  dir1:", bsp%get_continuity()

   print'(a)', "----------------------------------------"
   print'(a)', "Knot vectors"
   print'(a,*(f5.2,1x))', "  dir1:", bsp%get_knot()
   print'(a)', "========================================"

   print'(a)', "Fitting least squares curve..."
   call bsp%lsq_fit_bspline(Xt, Xdata, n)
   print'(a)', "Fitting complete."

   ! create B-Spline curve
   call bsp%create(Xt=Xt)
   Xg_eval = bsp%get_Xg()

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
   call bsp%export_Xc("vtk/lsq_fit_bspline_1d_Xc.vtk")
   call bsp%export_Xg("vtk/lsq_fit_bspline_1d_Xg.vtk")
   call bsp%show("vtk/lsq_fit_bspline_1d_Xc.vtk", "vtk/lsq_fit_bspline_1d_Xg.vtk")

end program
