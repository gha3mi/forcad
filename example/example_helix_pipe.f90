!> Example program demonstrating how to sweep a straight, pipe-like NURBS volume
!> into a **cylindrical helix** with a prescribed radius, pitch, and number of turns.
!>
!> The program:
!>   Creates a straight ring extruded along \(z\) (a pipe segment),
!>   Refines it along the axial direction (knot insertion + degree elevation),
!>   Maps the control points onto a helix of radius `rh`, pitch `pitch`, with `nturns` turns,
!>   Exports the resulting NURBS volume to VTK,
!>   Displays the geometry.
program example_helix_pipe

   use forcad, only: rk, nurbs_volume

   implicit none

   type(nurbs_volume) :: ring                                   !! Straight pipe ring.
   type(nurbs_volume) :: shape                                  !! Helical pipe shape.
   real(rk), parameter :: center(3) = [0.0_rk, 0.0_rk, 0.0_rk]  !! Pipe center \((c_x,c_y,c_z)\).
   real(rk), parameter :: r1 = 0.1_rk                           !! Inner radius of the pipe.
   real(rk), parameter :: r2 = 0.2_rk                           !! Outer radius of the pipe.
   real(rk), parameter :: length = 1.0_rk                       !! Length of the straight segment.
   real(rk), parameter :: rh = 1.0_rk                           !! Helix radius.
   real(rk), parameter :: pitch = 1.5_rk                        !! Helix pitch (axial rise per full turn).
   integer,  parameter :: nturns = 3                            !! Number of full turns.

   integer, parameter :: N = 50                                 !! # interior knots to insert along \(z\).
   integer :: i

   ! Initialize a straight pipe ring
   call ring%set_ring(center=center, radius1=r1, radius2=r2, length=length)

   ! Refine along z (dir=3): add N evenly-spaced interior knots (multiplicity 1 keeps C^{p-1}).
   call ring%insert_knots(dir=3, Xth=[(real(i,rk)/real(N+1,rk), i=1,N)], r=[(1, i=1,N)])

   ! elevate degree along z to increase smoothness/flexibility
   call ring%elevate_degree(3,5)

   ! Build the helical shape
   shape = ring
   call build_helix(shape, c=center, rh=rh, p=pitch, n=nturns)

   ! Create the NURBS volume sampling
   call shape%create(50,50,250)

   ! Export the NURBS volume to VTK files
   call shape%export_Xc("vtk/helix_pipe_Xc.vtk")
   call shape%export_Xg("vtk/helix_pipe_Xg.vtk")
   call shape%export_Xth_in_Xg("vtk/helix_pipe_Xth.vtk", res=20)

   ! Show the NURBS volume
   call shape%show("vtk/helix_pipe_Xc.vtk","vtk/helix_pipe_Xg.vtk","vtk/helix_pipe_Xth.vtk")

contains

   !===============================================================================
   !> Map a straight pipe-like NURBS volume onto a **cylindrical helix**.
   !>
   !> Each control point \(\mathbf{X}_c=(x,y,z)\) is first expressed relative to
   !> the input center \(\mathbf{c}=(c_x,c_y,c_z)\). For the axial (k) index,
   !> define a normalized parameter
   !> \[
   !>   s(k) =
   !>   \begin{cases}
   !>     \dfrac{k-1}{n_c^{(z)}-1}, & n_c^{(z)} > 1,\\[4pt]
   !>     0, & \text{otherwise},
   !>   \end{cases}
   !> \]
   !> which runs from 0 (start) to 1 (end) along the volume.
   !>
   !> The helix angle as a function of \(s\) is
   !> \[
   !>   \theta(s) = 2\pi\,n\,s ,
   !> \]
   !> where `n` is the number of turns. For a point at polar coordinates
   !> \(\rho=\sqrt{(x-c_x)^2+(y-c_y)^2}\) and \(\phi=\operatorname{atan2}(y-c_y,x-c_x)\),
   !> the mapping is:
   !> \[
   !>   X = c_x + \bigl(r_h + \rho\cos\phi\bigr)\cos\theta(s),\quad
   !>   Y = c_y + \bigl(r_h + \rho\cos\phi\bigr)\sin\theta(s),\quad
   !>   Z = c_z + \rho\sin\phi + p\,s ,
   !> \]
   !> where \(r_h\) is the helix radius and \(p\) is the pitch (axial rise per turn).
   !>
   !>  Knots are preserved; only the control points are updated.
   pure subroutine build_helix(this, c, rh, p, n)
      type(nurbs_volume), intent(inout) :: this        !! NURBS volume to be helicalized.
      real(rk), intent(in) :: c(3), rh, p              !! Center \(\mathbf{c}\), helix radius \(r_h\), pitch \(p\).
      integer,  intent(in) :: n                        !! Number of turns.
      real(rk), allocatable :: Xc(:,:), X4(:,:,:,:)
      integer :: nc(3), i,j,k, dim
      real(rk) :: s, theta, x0,y0,rho,phi, X,Y,Z
      real(rk), parameter :: pi = acos(-1.0_rk)

      Xc  = this%get_Xc();  dim = size(Xc,2)
      nc  = this%get_nc()
      X4  = reshape(Xc, [nc(1),nc(2),nc(3),dim])

      do k=1,nc(3)
         s = merge( real(k-1,rk)/real(max(1,nc(3)-1),rk), 0._rk, nc(3)>1 )
         theta = 2.0_rk*pi*real(n,rk)*s

         do j=1,nc(2)
            do i=1,nc(1)
               x0  = X4(i,j,k,1) - c(1)
               y0  = X4(i,j,k,2) - c(2)
               rho = sqrt(x0*x0 + y0*y0)
               phi = merge(atan2(y0,x0), 0._rk, rho>0._rk)

               X = (rh + rho*cos(phi))*cos(theta) + c(1)
               Y = (rh + rho*cos(phi))*sin(theta) + c(2)
               Z =  rho*sin(phi) + c(3) + p*s

               X4(i,j,k,1)=X; X4(i,j,k,2)=Y; X4(i,j,k,3)=Z
            end do
         end do
      end do

      Xc = reshape(X4, [product(nc),dim])

      if (this%is_rational()) then
         call this%set(this%get_knot(1), this%get_knot(2), this%get_knot(3), Xc, this%get_Wc())
      else
         call this%set(this%get_knot(1), this%get_knot(2), this%get_knot(3), Xc)
      end if
   end subroutine
   !===============================================================================

end program
