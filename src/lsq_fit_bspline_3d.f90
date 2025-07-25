program lsq_fit_bspline_3d

   use forcad_kinds, only: rk
   use forcad_utils, only: ndgrid
   use forcad_nurbs_volume, only: nurbs_volume

   implicit none
   type(nurbs_volume) :: bsp
   integer :: n(3), ndata, i
   real(rk), parameter :: pi = acos(-1.0_rk)
   real(rk), allocatable :: Xdata(:,:)
   real(rk), allocatable :: Xt1(:), Xt2(:), Xt3(:), Xt(:,:)
   real(rk), allocatable :: Xg_eval(:,:)
   real(rk) :: err1, err2, err3, rms

   n = [10,10,10]

   ! create parametric grid points
   allocate(Xt1(n(1)), Xt2(n(2)), Xt3(n(3)))
   do concurrent (i = 1:n(1))
      Xt1(i) = real(i - 1, rk) / real(n(1) - 1, rk)
   end do
   do concurrent (i = 1:n(2))
      Xt2(i) = real(i - 1, rk) / real(n(2) - 1, rk)
   end do
   do concurrent (i = 1:n(3))
      Xt3(i) = real(i - 1, rk) / real(n(3) - 1, rk)
   end do
   call ndgrid(Xt1, Xt2, Xt3, Xt)

   ! data points to be fitted
   ndata = n(1) * n(2) * n(3)
   allocate(Xdata(ndata, 3))
   do i = 1, ndata
      Xdata(i,1) = Xt(i,1) + 0.1_rk * sin(2.0_rk * pi * Xt(i,2))
      Xdata(i,2) = Xt(i,2) + 0.1_rk * sin(2.0_rk * pi * Xt(i,3))
      Xdata(i,3) = Xt(i,3) + 0.1_rk * sin(2.0_rk * pi * Xt(i,1))
   end do

   ! set up B-Spline volume
   call bsp%set(&
      degree      = [3, 3, 3],&
      Xth_dir1    = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk],&
      Xth_dir2    = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk],&
      Xth_dir3    = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk],&
      continuity1 = [ -1   ,   1    ,   1   ,   1    ,  -1   ],&
      continuity2 = [ -1   ,   1    ,   1   ,   1    ,  -1   ],&
      continuity3 = [ -1   ,   1    ,   1   ,   1    ,  -1   ])

   print'(a)', "========================================"
   print'(a)', "B-Spline Volume Configuration"
   print'(a)', "----------------------------------------"
   print'(a,3(i0,a))', "Degrees    : ", bsp%get_degree(1), ", ", bsp%get_degree(2), ", ", bsp%get_degree(3)
   print'(a,3(i0,a))', "Control pts: ", bsp%get_nc(1), " x ", bsp%get_nc(2), " x ", bsp%get_nc(3)
   print'(a,3(i0,a))', "Data grid  : ", n(1), " x ", n(2), " x ", n(3)

   print'(a)', "----------------------------------------"
   print'(a)', "Continuity"
   print'(a,*(i3,1x))', "  dir1:", bsp%get_continuity(1)
   print'(a,*(i3,1x))', "  dir2:", bsp%get_continuity(2)
   print'(a,*(i3,1x))', "  dir3:", bsp%get_continuity(3)

   print'(a)', "----------------------------------------"
   print'(a)', "Knot vectors"
   print'(a,*(f5.2,1x))', "  dir1:", bsp%get_knot(1)
   print'(a,*(f5.2,1x))', "  dir2:", bsp%get_knot(2)
   print'(a,*(f5.2,1x))', "  dir3:", bsp%get_knot(3)
   print'(a)', "========================================"

   print'(a)', "Fitting least squares volume..."
   call bsp%lsq_fit_bspline(Xt, Xdata, n)
   print'(a)', "Fitting complete."

   ! create B-Spline volume
   ! call bsp%create(n(1), n(2), n(3))
   call bsp%create(Xt1=Xt1, Xt2=Xt2, Xt3=Xt3)
   Xg_eval = bsp%get_Xg()

   ! Compute relative errors in each direction
   err1 = norm2(Xg_eval(:,1) - Xdata(:,1)) / norm2(Xdata(:,1))
   err2 = norm2(Xg_eval(:,2) - Xdata(:,2)) / norm2(Xdata(:,2))
   err3 = norm2(Xg_eval(:,3) - Xdata(:,3)) / norm2(Xdata(:,3))
   rms  = sqrt((err1**2 + err2**2 + err3**2) / 3.0_rk)

   ! Report
   print'(a)', "========================================"
   print'(a)', "Fitting Error Report"
   print'(a)', "----------------------------------------"
   print'(a,e13.6)', "Rel. error (dir1):", err1
   print'(a,e13.6)', "Rel. error (dir2):", err2
   print'(a,e13.6)', "Rel. error (dir3):", err3
   print'(a,e13.6)', "Total RMS error  :", rms
   print'(a)', "========================================"

   ! Export results and visualize
   call bsp%export_Xc("vtk/lsq_fit_bspline_3d_Xc.vtk")
   call bsp%export_Xg("vtk/lsq_fit_bspline_3d_Xg.vtk")
   call bsp%show("vtk/lsq_fit_bspline_3d_Xc.vtk", "vtk/lsq_fit_bspline_3d_Xg.vtk")

end program
