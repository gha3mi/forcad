 !> Solves the 2D Poisson problem using Isogeometric Analysis (IGA).
 !>
 !> This code solves the equation:
 !> \[
 !>   -\Delta u = f \quad \text{in } \Omega, \qquad u = 0 \quad \text{on } \partial \Omega
 !> \]
 !> using a B-spline surface on a rectangular domain.
 !>
 !> The solution is discretized using tensor-product B-spline basis functions
 !> over a structured control net. The global linear system is assembled and
 !> solved using a basic internal Cholesky solver.
 !>
 !> The resulting solution is exported to VTK format for visualization.
 !>
 !> The L2 error norm with respect to the exact solution
 !> is computed and printed.
 !>
 !> @note
 !> This implementation uses B-spline geometry (no rational weights), hence the surface is not a full NURBS.
 !> @endnote
 !>
 !> @warning "Slow solver"
 !> The solver uses an internal Cholesky factorization which is not optimized. For large systems, consider replacing it with a more scalable external solver.
 !> @endwarning
 !>

program poisson_iga_solver_2d
   use forcad, only: rk, nurbs_surface
   use forcad_utils, only: solve
   use fortime, only: timer
   implicit none

   type(nurbs_surface) :: surf                !! NURBS surface object
   real(rk), parameter :: pi = acos(-1.0_rk)  !! Constant \( \pi \)
   integer :: ie                              !! Element index
   integer :: ig                              !! Quadrature (Gauss) point index
   integer :: i                               !! Generic loop index
   integer :: nelem                           !! Number of elements
   integer :: nnelem                          !! Number of local nodes per element
   integer :: nc(2)                           !! Number of control points in each direction
   integer :: nct                             !! Total number of control points
   integer :: res(2)                          !! Visualization resolution
   integer :: dof                             !! Degrees of freedom per control point (1 for scalar field)
   integer :: ndof                            !! Total degrees of freedom
   integer :: m(2)                            !! Mode numbers \( (m_1, m_2) \) for source and exact solution
   integer :: ki(2)                           !! Number of knots to insert in each direction
   integer, allocatable :: elem(:,:)          !! Element connectivity matrix
   integer, allocatable :: elem_e(:)          !! Local connectivity of current element
   integer, allocatable :: dirichlet_id(:)    !! Indices for Dirichlet boundary conditions
   real(rk), allocatable :: K(:,:)            !! Global stiffness matrix \( K \)
   real(rk), allocatable :: b(:)              !! Global right-hand side vector \( b \)
   real(rk), allocatable :: Xc(:,:)           !! Control point coordinates
   real(rk), allocatable :: Xg(:)             !! Physical coordinates at quadrature point
   real(rk), allocatable :: T(:)              !! Basis function values at quadrature point
   real(rk), allocatable :: dT(:,:)           !! Derivatives of basis functions at quadrature point
   real(rk), allocatable :: X(:,:)            !! Global solution vector \( X \)
   real(rk), allocatable :: u_h(:,:)          !! Interpolated solution field \( u(x,y) \) on grid
   real(rk) :: dA                             !! Differential element area \( \text{d}A = J \cdot w \)
   real(rk) :: l2_error                       !! L2 error norm accumulator \( \|u_h - u\|^2 \)
   real(rk) :: L(2)                           !! Domain size \( (L_1, L_2) \in \mathbb{R}^2 \)
   character(len=256) :: filename             !! Filename for VTK export
   type(timer) :: ti                          !! Timer object for performance measurement

   !> Domain size and number of control points
   L  = [1.0_rk, 1.0_rk]
   nc = [10, 10]
   !> Number knots to insert in each direction
   ki = [8, 8]
   !> Mode numbers for the source term and exact solution
   m = [2, 2]
   !> Resolution of the visualization grid
   res = [50, 50]
   !> filename for VTK export
   filename = "vtk/poisson_iga_solver_2d"

   !> Construct the NURBS surface
   !> For simplicity, set_tetragon creates a rectangular surface with uniform knot spacing
   !> For more complex geometries, use surf%set() with knots, continuity,...
   call surf%set_tetragon(L=L, nc=nc)

   !> Insert knots in the first and second directions
   call surf%insert_knots(1, [(real(i,rk)/real(ki(1),rk), i=1,ki(1)-1)], [(1, i=1, ki(1)-1)])
   call surf%insert_knots(2, [(real(i,rk)/real(ki(2),rk), i=1,ki(2)-1)], [(1, i=1, ki(2)-1)])

   !> Extract geometry and mesh structure
   Xc     = surf%get_Xc()
   elem   = surf%cmp_elem()
   nct    = product(surf%get_nc())
   nelem  = size(elem, 1)
   nnelem = size(elem, 2)
   dof    = 1
   ndof   = dof * nct

   print '(a,g0,",",g0)',   "Degree (dir1, dir2)                  : ", surf%get_degree(1), surf%get_degree(2)
   print '(a,g0,"x",g0)',   "Control net size (dir1 x dir2)       : ", nc(1), nc(2)
   print '(a,g0)',          "Total control points                 : ", nct
   print '(a,g0)',          "Number of elements                   : ", nelem
   print '(a,g0)',          "Degrees of freedom (DoFs)            : ", ndof
   print '(a,g0," x ",g0)', "Domain size (L1, L2)                 : ", L(1), L(2)
   print '(a,g0,",",g0)',   "Mode numbers (m1, m2)                : ", m(1), m(2)
   print '(a,g0,",",g0)',   "Knot insertion (dir1, dir2)          : ", ki(1), ki(2)
   print '(a,g0,"x",g0)',   "Visualization resolution (res1, res2): ", res(1), res(2)

   !> Assemble global stiffness matrix and load vector
   allocate(K(nct,nct), b(nct), source=0.0_rk)
   call ti%timer_start()
   !$omp parallel do private(ie, ig, elem_e, T, dT, Xg, dA) shared(K, b)
   do ie = 1, nelem
      elem_e = elem(ie, :)
      do ig = 1, nnelem
         call surf%ansatz(ie, ig, T, dT, dA)
         Xg = matmul(T, Xc(elem_e,:))
         !$omp critical
         b(elem_e)         = b(elem_e)         + T * source_term(Xg, L, m) * dA
         K(elem_e, elem_e) = K(elem_e, elem_e) + matmul(dT, transpose(dT)) * dA
         !$omp end critical
      end do
   end do
   !$omp end parallel do
   call ti%timer_stop(message="Assembly                             : ")

   !> Apply homogeneous Dirichlet boundary conditions
   call ti%timer_start()
   allocate(dirichlet_id(0))
   do i = 1, nct
      if (&
         abs(Xc(i,1)) < 1e-12_rk .or. abs(Xc(i,1)-L(1)) < 1e-12_rk .or. &
         abs(Xc(i,2)) < 1e-12_rk .or. abs(Xc(i,2)-L(2)) < 1e-12_rk) then
         dirichlet_id = [dirichlet_id, i]
      end if
   end do
   K(dirichlet_id, :) = 0.0_rk
   K(:, dirichlet_id) = 0.0_rk
   b(dirichlet_id) = 0.0_rk
   do concurrent (i = 1:size(dirichlet_id))
      K(dirichlet_id(i), dirichlet_id(i)) = 1.0_rk
   end do
   call ti%timer_stop(message="Boundary conditions                  : ")

   !> Solve the linear system K·X = b
   call ti%timer_start()
   X = solve(K, reshape(b, [nct,1]))
   call ti%timer_stop(message="System solution                      : ")

   !> Export solution at control points to VTK
   call surf%export_Xc(filename=trim(filename)//".vtk", point_data=reshape(X, [nct,1]), field_names=["u"])

   !> Interpolate solution and export field
   call surf%create(res1=res(1), res2=res(2))
   call surf%basis(Tgc = u_h)
   u_h = matmul(u_h, reshape(X(:,1), [nct,1]))
   call surf%export_Xg(filename=trim(filename)//"_interp.vtk", point_data=u_h, field_names=["u"])

   !> Compute the L2 error norm
   call ti%timer_start()
   l2_error = 0.0_rk
   do ie = 1, nelem
      elem_e = elem(ie, :)
      do ig = 1, nnelem
         call surf%ansatz(ie, ig, T, dT, dA)
         Xg = matmul(T, Xc(elem_e,:))
         l2_error = l2_error + (dot_product(T, X(elem_e,1)) - exact_solution(Xg, L, m))**2 * dA
      end do
   end do
   call ti%timer_stop(message="L2 error evaluation                  : ")

   print '(a,1pe10.4)', "L2 error norm                        = ", sqrt(l2_error)
   print '(a,a,a,a)', trim(filename)//".vtk", " and ", trim(filename)//"_interp.vtk", " exported"

   call surf%finalize()
   ! deallocate(K, b, Xc, Xg, T, dT, X, u_h)

contains

   !> Computes the source function \( f(x,y) = \sin(m_1 \pi x / L_1) \sin(m_2 \pi y / L_2) \)
   pure function source_term(p, d, n) result(f)
      real(rk), intent(in) :: p(2) !! Coordinates (x, y)
      real(rk), intent(in) :: d(2) !! Domain size (L1, L2)
      integer, intent(in)  :: n(2) !! Mode numbers (m1, m2)
      real(rk) :: f
      f = sin(n(1)*pi*p(1)/d(1)) * sin(n(2)*pi*p(2)/d(2))
   end function

   !> Computes the exact solution corresponding to the source term
   !>
   !> \[
   !> u(x, y) = \frac{1}{\lambda} \sin(m_1 \pi x / L_1) \sin(m_2 \pi y / L_2)
   !> \quad\text{where}\quad
   !> \lambda = \left( \frac{m_1 \pi}{L_1} \right)^2 + \left( \frac{m_2 \pi}{L_2} \right)^2
   !> \]
   pure function exact_solution(p, d, n) result(u)
      real(rk), intent(in) :: p(2) !! Coordinates (x, y)
      real(rk), intent(in) :: d(2) !! Domain size (L1, L2)
      integer, intent(in)  :: n(2) !! Mode numbers (m1, m2)
      real(rk) :: u, lam
      lam = (n(1)*pi/d(1))**2 + (n(2)*pi/d(2))**2
      u = (1.0_rk / lam) * sin(n(1)*pi*p(1)/d(1)) * sin(n(2)*pi*p(2)/d(2))
   end function

end program
