!> Example program demonstrating how to sweep a straight, pipe-like NURBS volume
!> onto a **toroidal pipe** (donut) with an optional **progressive cross-section twist**
!> and a **sinusoidal wobble in the global z-direction** as it travels around the torus.
!>
!> The program:
!>   Creates a straight pipe segment (ring extruded along \(z\)),
!>   Refines it along the axial direction (knot insertion + degree elevation),
!>   Maps the control points onto a torus of major radius `R`,
!>   Applies a cross-section twist of `twist_turns` full rotations over one loop,
!>   Superimposes a vertical wobble \(z_\text{off}(s) = A_z\sin(2\pi\,n_\text{waves}\,s + \varphi)\),
!>   Exports the resulting NURBS volume to VTK and displays it,
!>   Prints the approximate enclosed volume (Gaussian integration).
program example_toroidal_pipe

   use forcad, only: rk, nurbs_volume

   implicit none
   type(nurbs_volume) :: ring  !! straight pipe segment
   type(nurbs_volume) :: shape !! toroidal pipe shape

   ! Straight ring (pipe) parameters
   real(rk), parameter :: center(3)   = [0.0_rk, 0.0_rk, 0.0_rk] !! center of the ring
   real(rk), parameter :: r1          = 0.12_rk                  !! inner radius
   real(rk), parameter :: r2          = 0.20_rk                  !! outer radius
   real(rk), parameter :: length      = 1.0_rk                   !! length of the pipe

   ! Toroid mapping parameters
   real(rk), parameter :: R           = 1.00_rk  !! major radius
   real(rk), parameter :: twist_turns = 1.0_rk   !! cross-section twist (turns over full loop)

   ! z-sine wobble parameters (global Z offset around the torus)
   real(rk), parameter :: Az          = 0.1_rk  !! amplitude of vertical wobble
   integer,  parameter :: nwaves_z    = 6       !! number of sine waves per full loop
   real(rk), parameter :: phase_z     = 0.1_rk  !! phase shift (radians)

   ! Refinement along z (dir=3)
   integer,  parameter :: Nref = 30
   integer :: i

   ! Make a straight pipe ring
   call ring%set_ring(center=center, radius1=r1, radius2=r2, length=length)

   ! Refine along z (keep C^{p-1})
   call ring%insert_knots(dir=3, Xth=[(real(i,rk)/real(Nref+1,rk), i=1,Nref)], r=[(1, i=1,Nref)])
   call ring%elevate_degree(3, 5)

   ! Map onto a torus with twist + z-sine wobble
   shape = ring
   call map_to_torus_sineZ(shape, c=center, R=R, twist_turns=twist_turns, Az=Az, nwaves_z=nwaves_z, phase_z=phase_z)

   ! Export the NURBS volume to VTK files
   call shape%create(30, 30, 120)
   call shape%export_Xc("vtk/example_toroidal_pipe_Xc.vtk")
   call shape%export_Xg("vtk/example_toroidal_pipe_Xg.vtk")
   call shape%export_Xth_in_Xg("vtk/example_toroidal_pipe_Xth.vtk", res=24)

   ! Show the resulting VTK files
   call shape%show("vtk/example_toroidal_pipe_Xc.vtk","vtk/example_toroidal_pipe_Xg.vtk","vtk/example_toroidal_pipe_Xth.vtk")

contains

   !===============================================================================
   !> Map a straight pipe (ring extruded in \(z\)) onto a **toroidal pipe** with:
   !>   • cross-section twist: \(\phi' = \phi + 2\pi\,(\texttt{twist\_turns})\,s\),
   !>   • vertical wobble: \(z_\text{off}(s) = A_z \sin(2\pi\,n_\text{waves}\,s + \varphi)\).
   !>
   !> Parameterization:
   !> For the axial (k) index, define
   !> \[
   !>    s(k) = \begin{cases}
   !>      \dfrac{k-1}{n_c^{(z)}-1}, & n_c^{(z)} > 1,\\[4pt]
   !>      0, & \text{otherwise},
   !>    \end{cases}
   !> \]
   !> and the torus angle (here one full loop since \(\theta=2\pi s\)):
   !> \[
   !>   \theta(s) = 2\pi\,s .
   !> \]
   !>
   !> For each control point \(\mathbf{X}_c=(x,y,z)\), relative to \(\mathbf{c}=(c_x,c_y,c_z)\),
   !> set \(\rho=\sqrt{(x-c_x)^2+(y-c_y)^2}\), \(\phi=\mathrm{atan2}(y-c_y,x-c_x)\), then
   !> \(\phi'=\phi+ 2\pi(\texttt{twist\_turns})\,s\).
   !>
   !> Mapping onto a torus of major radius \(R\):
   !> \[
   !> \begin{aligned}
   !>   X &= c_x + \bigl(R + \rho\cos\phi'\bigr)\cos\theta,\\
   !>   Y &= c_y + \bigl(R + \rho\cos\phi'\bigr)\sin\theta,\\
   !>   Z &= c_z + \rho\sin\phi' + z_\text{off}(s), \qquad
   !>       z_\text{off}(s) = A_z \sin\!\bigl(2\pi\,n_\text{waves}\,s + \varphi\bigr).
   !> \end{aligned}
   !> \]
   !>
   !> Knots/weights are preserved; only the control points are updated.
   pure subroutine map_to_torus_sineZ(this, c, R, twist_turns, Az, nwaves_z, phase_z)
      type(nurbs_volume), intent(inout) :: this
      real(rk), intent(in) :: c(3), R
      integer,  intent(in) :: nwaves_z
      real(rk), intent(in) :: twist_turns, Az, phase_z

      real(rk), allocatable :: Xc(:,:), X4(:,:,:,:)
      integer :: nc(3), i, j, k, dim
      real(rk) :: s, theta, x0, y0, rho, phi, phip, X, Y, Z, z_off
      real(rk), parameter :: pi = acos(-1.0_rk)

      Xc  = this%get_Xc();   dim = size(Xc,2)
      nc  = this%get_nc()
      X4  = reshape(Xc, [nc(1), nc(2), nc(3), dim])

      do k = 1, nc(3)
         s     = merge( real(k-1,rk)/real(max(1,nc(3)-1),rk), 0._rk, nc(3)>1 )
         theta = 2.0_rk*pi*s
         z_off = Az * sin(2.0_rk*pi*real(nwaves_z,rk)*s + phase_z)

         do j = 1, nc(2)
            do i = 1, nc(1)
               x0  = X4(i,j,k,1) - c(1)
               y0  = X4(i,j,k,2) - c(2)
               rho = sqrt(x0*x0 + y0*y0)
               phi = merge(atan2(y0,x0), 0._rk, rho>0._rk)
               phip = phi + 2.0_rk*pi*twist_turns*s

               X = c(1) + (R + rho*cos(phip))*cos(theta)
               Y = c(2) + (R + rho*cos(phip))*sin(theta)
               Z = c(3) +  rho*sin(phip) + z_off

               X4(i,j,k,1) = X
               X4(i,j,k,2) = Y
               X4(i,j,k,3) = Z
            end do
         end do
      end do

      Xc = reshape(X4, [product(nc), dim])
      if (this%is_rational()) then
         call this%set(this%get_knot(1), this%get_knot(2), this%get_knot(3), Xc, this%get_Wc())
      else
         call this%set(this%get_knot(1), this%get_knot(2), this%get_knot(3), Xc)
      end if
   end subroutine
   !===============================================================================

end program
