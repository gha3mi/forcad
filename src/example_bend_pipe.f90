!> Example program demonstrating how to bend a straight pipe-like NURBS volume into a circular arc.
!>
!> The program:
!>   Creates a straight pipe segment (as a ring extruded in \(z\)),
!>   Refines the shape,
!>   Applies bending with different bend angles (90°, 270°, 360°),
!>   Exports the resulting NURBS volumes to VTK,
!>   Displays the geometry.
program example_bend_pipe

   use forcad, only: rk, nurbs_volume

   implicit none
   type(nurbs_volume) :: ring                                    !! Straight pipe ring
   type(nurbs_volume) :: sh                                      !! Bent pipe shape
   real(rk), parameter :: c(3)  = [0.0_rk, 0.0_rk, 0.0_rk]       !! Center of the ring.
   real(rk), parameter :: r1    = 0.3_rk                         !! Inner radius of the pipe.
   real(rk), parameter :: r2    = 0.5_rk                         !! Outer radius of the pipe.
   real(rk), parameter :: l     = 1.0_rk                         !! Length of the straight segment.
   real(rk), parameter :: rb    = 1.5_rk                         !! Bend radius.
   real(rk), parameter :: a(3)  = [90.0_rk, 270.0_rk, 360.0_rk]  !! Bend angles (degrees).

   ! Initialize a straight pipe ring
   call ring%set_ring(c, r1, r2, l)

   ! Refine the pipe ring
   call ring%elevate_degree(3,1)
   call ring%insert_knots(dir=3, Xth=[0.1_rk, 0.2_rk, 0.3_rk, 0.4_rk, 0.5_rk, 0.6_rk, 0.7_rk, 0.8_rk, 0.9_rk], r=[1,1,1,1,1,1,1,1,1])
   call ring%elevate_degree(3,5)

   !===============================================================================
   ! Apply bending with angle 90°

   ! Set the shape
   sh = ring

   ! Build the bend pipe shape
   call bend_pipe(sh, c, l, rb, a(1))

   ! Create the NURBS volume
   call sh%create(30,40,80)

   ! Export the NURBS volume to VTK files
   call sh%export_Xc("vtk/bend_pipe_pipe90_Xc.vtk")
   call sh%export_Xg("vtk/bend_pipe_pipe90_Xg.vtk")
   call sh%export_Xth_in_Xg("vtk/bend_pipe_pipe90_Xth.vtk")

   ! Show the NURBS volume
   call sh%show("vtk/bend_pipe_pipe90_Xc.vtk", "vtk/bend_pipe_pipe90_Xg.vtk", "vtk/bend_pipe_pipe90_Xth.vtk")

   !===============================================================================
   ! Apply bending with angle 270°

   ! Set the shape
   sh = ring

   ! Build the bend pipe shape
   call bend_pipe(sh, c, l, rb, a(2))

   ! Create the NURBS volume
   call sh%create(30,40,80)

   ! Export the NURBS volume to VTK files
   call sh%export_Xc("vtk/bend_pipe_pipe270_Xc.vtk")
   call sh%export_Xg("vtk/bend_pipe_pipe270_Xg.vtk")
   call sh%export_Xth_in_Xg("vtk/bend_pipe_pipe270_Xth.vtk")

   ! Show the NURBS volume
   call sh%show("vtk/bend_pipe_pipe270_Xc.vtk", "vtk/bend_pipe_pipe270_Xg.vtk", "vtk/bend_pipe_pipe270_Xth.vtk")

   !===============================================================================
   ! Apply bending with angle 360°

   ! Set the shape
   sh = ring

   ! Build the bend pipe shape
   call bend_pipe(sh, c, l, rb, a(3))

   ! Create the NURBS volume
   call sh%create(30,40,80)

   ! Export the NURBS volume to VTK files
   call sh%export_Xc("vtk/bend_pipe_pipe360_Xc.vtk")
   call sh%export_Xg("vtk/bend_pipe_pipe360_Xg.vtk")
   call sh%export_Xth_in_Xg("vtk/bend_pipe_pipe360_Xth.vtk")

   ! Show the NURBS volume
   call sh%show("vtk/bend_pipe_pipe360_Xc.vtk", "vtk/bend_pipe_pipe360_Xg.vtk", "vtk/bend_pipe_pipe360_Xth.vtk")

contains

   !===============================================================================
   !> Bend a straight pipe-like NURBS volume into a circular arc.
   !> Each control point \(\mathbf{X}_c=(x,y,z)\) is transformed relative
   !> to the center \(\mathbf{c}=(c_x,c_y,c_z)\).
   !> Let the total bend angle in degrees be \(\alpha^\circ\) (input),
   !> and define the total bend angle in radians
   !> \[
   !>   \alpha_{\max} = \frac{\pi}{180}\,\alpha^\circ .
   !> \]
   !> For each control point:
   !> \[
   !>   \rho = \sqrt{(x-c_x)^2+(y-c_y)^2},\quad
   !>   \phi = \operatorname{atan2}(y-c_y,\,x-c_x),\quad
   !>   Z = z - c_z ,
   !> \]
   !> \[
   !>   \theta(Z) = \alpha_{\max}\,\min\!\bigl(1,\max\!\bigl(0,\tfrac{Z}{L}\bigr)\bigr) .
   !> \]
   !> Mapping:
   !> \[
   !>   x' = c_x + (\rho\cos\phi + R_b)\cos\theta,\quad
   !>   y' = c_y + (\rho\cos\phi + R_b)\sin\theta,\quad
   !>   z' = c_z + \rho\sin\phi .
   !> \]
   !> Knots and weights are preserved; only the control lattice is updated.
   pure subroutine bend_pipe(this, center, length, rbend, angle_deg)
      type(nurbs_volume), intent(inout) :: this !! NURBS volume to be bent.
      real(rk), intent(in) :: center(3)         !! Pipe center coordinates \((c_x,c_y,c_z)\).
      real(rk), intent(in) :: length            !! Length of the straight pipe segment before bending.
      real(rk), intent(in) :: rbend             !! Bend radius \(R_b\), i.e. distance from the bend centerline.
      real(rk), intent(in) :: angle_deg         !! Bend angle in degrees \(\alpha^\circ\).

      real(rk), allocatable :: Xc(:,:), X4(:,:,:,:)
      integer :: nc(3), i,j,k, d
      real(rk) :: x0, y0, z0, theta, cth, sth, rho, phi, ang_tot
      real(rk), parameter :: pi = acos(-1.0_rk)

      Xc  = this%get_Xc()
      d   = size(Xc,2)
      nc  = this%get_nc()
      X4  = reshape(Xc, shape=[nc(1), nc(2), nc(3), d], order=[1,2,3,4])

      ang_tot = angle_deg*pi/180.0_rk

      do k = 1, nc(3)
         z0    = X4(1,1,k,3) - center(3)
         theta = ang_tot * max(0.0_rk, min(1.0_rk, z0/length))
         cth   = cos(theta)
         sth   = sin(theta)

         do j = 1, nc(2)
            do i = 1, nc(1)
               x0  = X4(i,j,k,1) - center(1)
               y0  = X4(i,j,k,2) - center(2)
               rho = sqrt(x0*x0 + y0*y0)
               phi = merge(atan2(y0,x0), 0.0_rk, rho>0.0_rk)
               X4(i,j,k,1) = center(1) + (rho*cos(phi) + rbend) * cth
               X4(i,j,k,2) = center(2) + (rho*cos(phi) + rbend) * sth
               X4(i,j,k,3) = center(3) +  rho*sin(phi)
            end do
         end do
      end do

      Xc = reshape(X4, shape=[product(nc), d], order=[1,2])

      if (this%is_rational()) then
         call this%set(this%get_knot(1), this%get_knot(2), this%get_knot(3), Xc, this%get_Wc())
      else
         call this%set(this%get_knot(1), this%get_knot(2), this%get_knot(3), Xc)
      end if
   end subroutine
   !===============================================================================

end program
