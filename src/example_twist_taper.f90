!> Example program demonstrating how to apply a **progressive twist** and **linear taper**
!> to a hexahedral NURBS volume along its axial (z) direction.
!>
!> The program:
!>   Creates a straight hexahedral block (control points `nc = [7,7,9]` over a box of size `L`),
!>   Copies it to `shape`,
!>   Applies a z-dependent twist up to \(\alpha_{\max}^\circ =\) `twist_deg` and a linear taper to factor `taper`,
!>   Exports the resulting NURBS volume to VTK,
!>   Displays the geometry.
program example_twist_taper

   use forcad, only: rk, nurbs_volume

   implicit none
   type(nurbs_volume) :: sh, hexa
   real(rk), parameter :: L(3) = [1.0_rk, 1.0_rk, 3.0_rk]  !! Domain extents in \(x,y,z\).
   integer,  parameter :: nc(3) = [7, 7, 9]                !! Control point counts per direction.
   real(rk), parameter :: twist_deg = 360.0_rk             !! Total twist angle (degrees) at the top face.
   real(rk), parameter :: taper = 0.1_rk                   !! Target in-plane scale at the top face (0<`taper`≤1).

   ! Initialize a straight hexahedral block
   call hexa%set_hexahedron(L=L, nc=nc)

   ! Work on a copy
   sh = hexa

   ! Build the twisted and tapered shape
   call build_twist_taper(sh, L, nc, twist_deg, taper)

   ! Create the NURBS volume sampling
   call sh%create(30,30,80)

   ! Export the NURBS volume to VTK files
   call sh%export_Xc("vtk/example_twist_taper_Xc.vtk")
   call sh%export_Xg("vtk/example_twist_taper_Xg.vtk")
   call sh%export_Xth_in_Xg("vtk/example_twist_taper_Xth_in_Xg.vtk", res=30)

   ! Show the NURBS volume
   call sh%show("vtk/example_twist_taper_Xc.vtk","vtk/example_twist_taper_Xg.vtk", "vtk/example_twist_taper_Xth_in_Xg.vtk")

contains

   !===============================================================================
   !> Apply a **z-progressive twist** and **linear taper** to a NURBS hexahedron.
   !>
   !> Let the control points be indexed by \(k=1,\dots,n_c^{(z)}\) along \(z\).
   !> Define the normalized axial coordinate
   !> \[
   !>   t(k) = \begin{cases}
   !>     \dfrac{k-1}{n_c^{(z)}-1}, & n_c^{(z)} > 1,\\[4pt]
   !>     0, & \text{otherwise},
   !>   \end{cases}
   !> \]
   !> which varies from 0 at the bottom face to 1 at the top face.
   !>
   !> The total twist angle (in radians) at level \(t\) is
   !> \[
   !>   \theta(t) = \Bigl(\dfrac{\pi}{180}\Bigr)\,\texttt{twist\_deg}\; t ,
   !> \]
   !> i.e. a linear ramp from \(0\) to \(\alpha_{\max} = (\pi/180)\,\texttt{twist\_deg}\).
   !>
   !> The in-plane (x–y) **taper scale** is chosen linear in \(t\):
   !> \[
   !>   s_{xy}(t) = (1-\texttt{taper})(1-t) + \texttt{taper},
   !> \]
   !> so \(s_{xy}(0)=1\) at the bottom and \(s_{xy}(1)=\texttt{taper}\) at the top (shrinking if \(\texttt{taper}<1\)).
   !>
   !> For each control point \(\mathbf{X}_c=(x,y,z)\), first shift to the in-plane
   !> centroid \(\mathbf{c}_{xy}=(c_x,c_y)\) of the box,
   !> apply the scale \(s_{xy}(t)\) and rotation by \(\theta(t)\) about \(\mathbf{c}_{xy}\),
   !> and keep \(z\) unchanged:
   !> \[
   !>   \begin{bmatrix}x'\\y'\end{bmatrix}
   !>   = \begin{bmatrix}c_x\\c_y\end{bmatrix}
   !>     + s_{xy}(t)\,
   !>       \begin{bmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{bmatrix}
   !>       \!\left(\begin{bmatrix}x\\y\end{bmatrix}-\begin{bmatrix}c_x\\c_y\end{bmatrix}\right),
   !>   \qquad
   !>   z' = z .
   !> \]
   !>
   !>  Knots are preserved; only the control points are updated.
   pure subroutine build_twist_taper(this, L, nc, twist_deg, taper)
      type(nurbs_volume), intent(inout) :: this !! Volume to be transformed.
      real(rk), intent(in) :: L(3)              !! Box lengths \((L_x,L_y,L_z)\).
      integer,  intent(in) :: nc(3)             !! Control points sizes.
      real(rk), intent(in) :: twist_deg         !! Total twist at top face (degrees).
      real(rk), intent(in) :: taper             !! In-plane scale at top face (0<`taper`≤1).

      real(rk), allocatable :: Xc(:,:), X4(:,:,:,:)
      real(rk), parameter :: pi = acos(-1.0_rk)
      real(rk) :: t, ang, sxy, ca, sa, x, y, cx, cy
      integer :: i, j, k, d

      Xc  = this%get_Xc()
      d  = size(Xc,2)
      X4 = reshape(Xc, shape=[nc(1), nc(2), nc(3), d], order=[1,2,3,4])

      cx = 0.5_rk*L(1)
      cy = 0.5_rk*L(2)
      do k=1,nc(3)
         t   = merge(real(k-1,rk)/real(max(1,nc(3)-1),rk), 0.0_rk, nc(3)>1)
         ang = (twist_deg*pi/180.0_rk) * t
         sxy = (1.0_rk - taper) * (1.0_rk - t) + taper
         ca  = cos(ang)
         sa  = sin(ang)
         do j=1,nc(2)
            do i=1,nc(1)
               x = (X4(i,j,k,1) - cx)*sxy
               y = (X4(i,j,k,2) - cy)*sxy
               X4(i,j,k,1) =  cx + ca*x - sa*y
               X4(i,j,k,2) =  cy + sa*x + ca*y
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
