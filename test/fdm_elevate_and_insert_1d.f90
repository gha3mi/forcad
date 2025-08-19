program fdm_elevate_and_insert_1d

   use forcad,       only: rk, nurbs_curve
   use forcad_utils, only: linspace, kron_eye
   use fortime,      only: timer

   implicit none

   type(nurbs_curve) :: sh0, shr, shfd
   real(rk), allocatable :: Xc(:,:), Wc(:)
   real(rk), allocatable :: Xc0(:,:), Xp(:,:), Xm(:,:)
   real(rk), allocatable :: Xcp_vec(:), Xcm_vec(:)
   real(rk), allocatable :: knot(:)
   real(rk), allocatable :: S1(:,:), S2(:,:)
   real(rk), allocatable :: Bs(:,:), B(:,:), Bfd(:,:)
   real(rk), allocatable :: u(:)
   integer,  allocatable :: r(:)
   real(rk) :: rel_err
   integer  :: dim, nc0, ndof_old, ndof_new, i, d, idx
   type(timer) :: t

   real(rk), parameter :: tol  = 1e-5_rk      !! tolerance of finite differences
   integer,  parameter :: tdeg = 4            !! degrees to elevate
   integer,  parameter :: nins = 10           !! number of knots to insert

   !> set control points
   Xc = generate_Xc(5.0_rk)

   !> set weights
   allocate(Wc(size(Xc,1)), source=1.0_rk)
   Wc(2) = 0.5_rk

   !> set knot vector
   knot = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

   !> set NURBS curve
   call shr%set(knot, Xc, Wc)

   !> deallocate temporary arrays
   deallocate(Xc, Wc, knot)

   !> copy initial NURBS curve (before refinement)
   sh0 = shr

   !> get initial control points, dimension, number of control points and degrees of freedom
   Xc0      = sh0%get_Xc()
   dim      = size(Xc0,2)
   nc0      = size(Xc0,1)
   ndof_old = nc0*dim

   !> elevate degree and get sensitivities (Bs: in compact form, memory efficient)
   call t%timer_start()
   call shr%elevate_degree(tdeg, Bs=S1)
   call t%timer_stop()

   !> set knot vectors to insert
   u = linspace(0.0_rk, 1.0_rk, nins+2)
   u = u(2:nins+1)
   !> multiplicities of knots to insert
   allocate(r(size(u)), source=2)

   !> insert knots and get sensitivities (Bs: in compact form, memory efficient)
   call t%timer_start()
   call shr%insert_knots(u, r, Bs=S2)
   call t%timer_stop()

   !> compute global sensitivities (dXc_old/dXc_new)
   call t%timer_start()
   Bs = matmul(S2, S1)
   B  = kron_eye(Bs, dim)
   call t%timer_stop()

   !> get new degrees of freedom (after refinement)
   ndof_new = size(shr%get_Xc(),1) * dim

   !> start finite difference computations
   allocate(Xp(nc0,dim), Xm(nc0,dim))
   allocate(Xcp_vec(ndof_new), Xcm_vec(ndof_new))
   allocate(Bfd(ndof_new, ndof_old))

   do idx = 1, ndof_old
      Xp = Xc0
      Xm = Xc0
      i = (idx-1)/dim + 1
      d = mod(idx-1, dim) + 1
      Xp(i,d) = Xp(i,d) + tol
      Xm(i,d) = Xm(i,d) - tol

      call shfd%set(sh0%get_knot(), Xp, sh0%get_Wc())
      call shfd%elevate_degree(tdeg)
      call shfd%insert_knots(u, r)
      Xcp_vec = reshape(transpose(shfd%get_Xc()), [ndof_new])

      call shfd%set(sh0%get_knot(), Xm, sh0%get_Wc())
      call shfd%elevate_degree(tdeg)
      call shfd%insert_knots(u, r)
      Xcm_vec = reshape(transpose(shfd%get_Xc()), [ndof_new])

      Bfd(:,idx) = (Xcp_vec - Xcm_vec) * (0.5_rk/tol)
   end do

   !> compute relative error between finite difference and analytical sensitivities
   rel_err = norm2(Bfd - B) / norm2(Bfd)

   print '(a)',        '--- CENTRAL FDM vs Analytic (elevate + insert) ---'
   print '(a,i0)',     '  ndof_old     = ', ndof_old
   print '(a,i0)',     '  ndof_new     = ', ndof_new
   print '(a,1pe12.4)','  ||B||_2      = ', norm2(B)
   print '(a,1pe12.4)','  ||Bfd||_2    = ', norm2(Bfd)
   print '(a,1pe12.4)','  ||B-Bfd||_2  = ', norm2(Bfd - B)
   print '(a,1pe12.4)','  rel l2 error = ', rel_err

   !> finalize
   call shr%finalize()
   call sh0%finalize()
   call shfd%finalize()

contains

   pure function generate_Xc(L) result(control_points)
      real(rk), intent(in) :: L
      real(rk), allocatable :: control_points(:,:)
      real(rk) :: L2
      L2 = L / 2.0_rk
      allocate(control_points(2,3))
      control_points(1,:) = [-L2, 0.0_rk, 0.0_rk]
      control_points(2,:) = [ L2, 0.0_rk, 0.0_rk]
   end function

end program
