program fdm_elevate_and_insert_3d

   use forcad,       only: rk, nurbs_volume
   use forcad_utils, only: linspace, kron_eye
   use fortime,      only: timer

   implicit none

   type(nurbs_volume) :: sh0, shr, shfd
   real(rk), allocatable :: Xc(:,:), Wc(:)
   real(rk), allocatable :: Xc0(:,:), Xp(:,:), Xm(:,:)
   real(rk), allocatable :: Xcp_vec(:), Xcm_vec(:)
   real(rk), allocatable :: knot1(:), knot2(:), knot3(:)
   real(rk), allocatable :: S1(:,:), S2(:,:), S3(:,:)
   real(rk), allocatable :: S4(:,:), S5(:,:), S6(:,:)
   real(rk), allocatable :: Bs(:,:), B(:,:), Bfd(:,:)
   real(rk), allocatable :: u1(:), u2(:), u3(:)
   integer,  allocatable :: r1(:), r2(:), r3(:)
   real(rk) :: rel_err
   integer  :: dim, nc0, ndof_old, ndof_new, i, d, idx
   type(timer) :: t

   real(rk), parameter :: tol = 1e-5_rk             !! tolerance of finite differences
   integer,  parameter :: dg1 = 5, dg2 = 4, dg3 = 3 !! degrees to elevate
   integer,  parameter :: n1 = 8, n2 = 9, n3 = 10   !! number of knots to insert

   !> set control points
   Xc = generate_Xc(5.0_rk)

   !> set weights
   allocate(Wc(size(Xc,1)), source=1.0_rk)
   Wc(2) = 0.5_rk

   !> set knot vectors
   knot1 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
   knot2 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
   knot3 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

   !> set NURBS volume
   call shr%set(knot1, knot2, knot3, Xc, Wc)

   !> deallocate temporary arrays
   deallocate(Xc, Wc, knot1, knot2, knot3)

   !> copy initial NURBS volume (before refinement)
   sh0 = shr

   !> get initial control points, dimension, number of control points and degrees of freedom
   Xc0      = sh0%get_Xc()
   dim      = size(Xc0,2)
   nc0      = size(Xc0,1)
   ndof_old = nc0*dim

   !> elevate degree in three directions and get sensitivities (Bs: in compact form, memory efficient)
   call t%timer_start()
   call shr%elevate_degree(1, dg1, Bs=S1)
   call shr%elevate_degree(2, dg2, Bs=S2)
   call shr%elevate_degree(3, dg3, Bs=S3)
   call t%timer_stop()

   !> set knot vectors to insert
   u1 = linspace(0.0_rk, 1.0_rk, n1+2)
   u1 = u1(2:n1+1)
   u2 = linspace(0.0_rk, 1.0_rk, n2+2)
   u2 = u2(2:n2+1)
   u3 = linspace(0.0_rk, 1.0_rk, n3+2)
   u3 = u3(2:n3+1)
   !> multiplicities of knots to insert
   allocate(r1(size(u1)), source=2)
   allocate(r2(size(u2)), source=1)
   allocate(r3(size(u3)), source=2)

   !> insert knots in three directions and get sensitivities (Bs: in compact form, memory efficient)
   call t%timer_start()
   call shr%insert_knots(1, u1, r1, Bs=S4)
   call shr%insert_knots(2, u2, r2, Bs=S5)
   call shr%insert_knots(3, u3, r3, Bs=S6)
   call t%timer_stop()

   !> compute global sensitivities (dXc_old/dXc_new)
   call t%timer_start()
   Bs = matmul(S6, matmul(S5, matmul(S4, matmul(S3, matmul(S2, S1)))))
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

      call shfd%set(sh0%get_knot(1), sh0%get_knot(2), sh0%get_knot(3), Xp, sh0%get_Wc())
      call shfd%elevate_degree(1, dg1)
      call shfd%elevate_degree(2, dg2)
      call shfd%elevate_degree(3, dg3)
      call shfd%insert_knots(1, u1, r1)
      call shfd%insert_knots(2, u2, r2)
      call shfd%insert_knots(3, u3, r3)
      Xcp_vec = reshape(transpose(shfd%get_Xc()), [ndof_new])

      call shfd%set(sh0%get_knot(1), sh0%get_knot(2), sh0%get_knot(3), Xm, sh0%get_Wc())
      call shfd%elevate_degree(1, dg1)
      call shfd%elevate_degree(2, dg2)
      call shfd%elevate_degree(3, dg3)
      call shfd%insert_knots(1, u1, r1)
      call shfd%insert_knots(2, u2, r2)
      call shfd%insert_knots(3, u3, r3)
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
      allocate(control_points(8, 3))
      control_points(1,:) = [ L2, -L2,  L2]
      control_points(2,:) = [ L2, -L2, -L2]
      control_points(3,:) = [-L2, -L2,  L2]
      control_points(4,:) = [-L2, -L2, -L2]
      control_points(5,:) = [ L2,  L2,  L2]
      control_points(6,:) = [ L2,  L2, -L2]
      control_points(7,:) = [-L2,  L2,  L2]
      control_points(8,:) = [-L2,  L2, -L2]
   end function

end program
