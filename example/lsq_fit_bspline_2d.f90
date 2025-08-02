program lsq_fit_bspline_2d

   use forcad, only: rk, nurbs_surface
   use forcad_utils, only: ndgrid

   implicit none

   type(nurbs_surface) :: bsp
   integer :: n(2), ndata, i
   real(rk), parameter :: pi = acos(-1.0_rk)
   real(rk), allocatable :: Xdata(:,:)
   real(rk), allocatable :: Xt1(:), Xt2(:), Xt(:,:)
   real(rk), allocatable :: Xg_eval(:,:)
   real(rk) :: err1, err2, err3, rms

   n = [14,14]

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
      Xdata(i,1) = Xt(i,1)
      Xdata(i,2) = Xt(i,2)
      Xdata(i,3) = 0.1_rk * sin(2.0_rk*pi*Xt(i,1)) * cos(2.0_rk*pi*Xt(i,2))
   end do

   ! set up B-Spline surface
   ! Xth_dir1(1) = minval(Xt1), Xth_dir1(2) = maxval(Xt1)
   ! Xth_dir2(1) = minval(Xt2), Xth_dir2(2) = maxval(Xt2)
   call bsp%set(&
      degree      = [4, 4],&
      Xth_dir1    = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk],&
      Xth_dir2    = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk],&
      continuity1 = [ -1   ,   1    ,   1   ,   1    ,  -1   ],&
      continuity2 = [ -1   ,   1    ,   1   ,   1    ,  -1   ])

   print'(a)', "========================================"
   print'(a)', "B-Spline Surface Configuration"
   print'(a)', "----------------------------------------"
   print'(a,2(i0,a))', "Degrees    : ", bsp%get_degree(1), ", ", bsp%get_degree(2)
   print'(a,2(i0,a))', "Control pts: ", bsp%get_nc(1), " x ", bsp%get_nc(2)
   print'(a,2(i0,a))', "Data grid  : ", n(1), " x ", n(2)

   print'(a)', "----------------------------------------"
   print'(a)', "Continuity"
   print'(a,*(i3,1x))', "  dir1:", bsp%get_continuity(1)
   print'(a,*(i3,1x))', "  dir2:", bsp%get_continuity(2)

   print'(a)', "----------------------------------------"
   print'(a)', "Knot vectors"
   print'(a,*(f5.2,1x))', "  dir1:", bsp%get_knot(1)
   print'(a,*(f5.2,1x))', "  dir2:", bsp%get_knot(2)
   print'(a)', "========================================"

   print'(a)', "Fitting least squares surface..."
   call bsp%lsq_fit_bspline(Xt, Xdata, n)
   print'(a)', "Fitting complete."

   ! create B-Spline surface
   call bsp%create(n(1), n(2))
   Xg_eval = bsp%get_Xg()

   ! Compute errors
   err1 = norm2(Xg_eval(:,1) - Xdata(:,1)) / norm2(Xdata(:,1))
   err2 = norm2(Xg_eval(:,2) - Xdata(:,2)) / norm2(Xdata(:,2))
   err3 = norm2(Xg_eval(:,3) - Xdata(:,3)) / norm2(Xdata(:,3))
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
   call bsp%export_Xc("vtk/lsq_fit_bspline_2d_Xc.vtk")
   call bsp%export_Xg("vtk/lsq_fit_bspline_2d_Xg.vtk")
   call bsp%show("vtk/lsq_fit_bspline_2d_Xc.vtk", "vtk/lsq_fit_bspline_2d_Xg.vtk")

end program
