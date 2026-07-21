!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
program surface_plate_hole_half

   use forcad, only: rk, nurbs_surface

   implicit none
   type(nurbs_surface) :: plate_hole
   real(rk), allocatable :: Xc(:,:)
   real(rk), allocatable :: Wc(:)

   real(rk), parameter :: radius1 = 2.5_rk
   real(rk), parameter :: radius2 = 3.5_rk
   real(rk), parameter :: length  = 5.0_rk
   real(rk), parameter :: height  = 5.0_rk

   call set_Xc_Wc("ellipse", [radius1, radius2, length, height], Xc, Wc)

   call plate_hole%set(&
      knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk],&
      knot2 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk],&
      Xc    = Xc,&
      Wc    = Wc&
      )

   call plate_hole%create(31, 31)
   call plate_hole%export_Xc("vtk/surface_plate_hole_half_Xc.vtk")
   call plate_hole%export_Xg("vtk/surface_plate_hole_half_Xg.vtk")
   call plate_hole%export_Xth_in_Xg("vtk/surface_plate_hole_half_Xth.vtk")

   call plate_hole%show("vtk/surface_plate_hole_half_Xc.vtk", "vtk/surface_plate_hole_half_Xg.vtk", "vtk/surface_plate_hole_half_Xth.vtk")

contains

   !===============================================================================
   pure subroutine set_Xc_Wc(tp, params, X_c, W_c)
      character(len=*), intent(in) :: tp
      real(rk), intent(in), contiguous :: params(:)
      real(rk), allocatable, intent(out) :: X_c(:,:)
      real(rk), allocatable, intent(out) :: W_c(:)

      real(rk) :: r1, r2, l, h

      select case (tp)
       case("circle")

         r1 = params(1)
         l  = params(3)
         h  = params(4)


         allocate(X_c(21,2))
         X_c(1, :) = [-r1, 0.0_rk]
         X_c(2, :) = [-r1, r1*tand(22.5_rk)]
         X_c(3, :) = [-r1*tand(22.5_rk), r1]
         X_c(4, :) = [0.0_rk, r1]
         X_c(5, :) = [r1*tand(22.5_rk), r1]
         X_c(6, :) = [r1, r1*tand(22.5_rk)]
         X_c(7, :) = [r1, 0.0_rk]
         X_c(8, :) = [-(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(9, :) = [-(r1 + (l-r1)/2.0_rk), (r1 + (h-r1)/2.0_rk)*tand(16.7_rk)]
         X_c(10,:) = [-(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), (r1 + (h-r1)/2.0_rk)]
         X_c(11,:) = [0.0_rk, (r1 + (h-r1)/2.0_rk)]
         X_c(12,:) = [(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), (r1 + (h-r1)/2.0_rk)]
         X_c(13,:) = [(r1 + (l-r1)/2.0_rk), (r1 + (h-r1)/2.0_rk)*tand(16.7_rk)]
         X_c(14,:) = [(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(15,:) = [-l, 0.0_rk]
         X_c(16,:) = [-l, h]
         X_c(17,:) = [-l, h]
         X_c(18,:) = [0.0_rk, h]
         X_c(19,:) = [l, h]
         X_c(20,:) = [l, h]
         X_c(21,:) = [l, 0.0_rk]

         allocate(W_c(21), source=1.0_rk)
         W_c([2,3,5,6]) = (1.0_rk + 1.0_rk/sqrt(2.0_rk))/2.0_rk

       case("ellipse")

         r1 = params(1)
         r2 = params(2)
         l  = params(3)
         h  = params(4)


         allocate(X_c(21,2))
         X_c(1 ,:)  = [-r1, 0.0_rk]
         X_c(2 ,:)  = [-r1, r2*tand(22.5_rk)]
         X_c(3 ,:)  = [-r1*tand(22.5_rk), r2]
         X_c(4 ,:)  = [0.0_rk, r2]
         X_c(5 ,:)  = [r1*tand(22.5_rk), r2]
         X_c(6 ,:)  = [r1, r2*tand(22.5_rk)]
         X_c(7 ,:)  = [r1, 0.0_rk]
         X_c(8 ,:)  = [-(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(9 ,:)  = [-(r1 + (l-r1)/2.0_rk), (r2 + (h-r2)/2.0_rk)*tand(16.7_rk)]
         X_c(10,:)  = [-(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), (r2 + (h-r2)/2.0_rk)]
         X_c(11,:)  = [0.0_rk, (r2 + (h-r2)/2.0_rk)]
         X_c(12,:)  = [(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), (r2 + (h-r2)/2.0_rk)]
         X_c(13,:)  = [(r1 + (l-r1)/2.0_rk), (r2 + (h-r2)/2.0_rk)*tand(16.7_rk)]
         X_c(14,:)  = [(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(15,:)  = [-l, 0.0_rk]
         X_c(16,:)  = [-l, h]
         X_c(17,:)  = [-l, h]
         X_c(18,:)  = [0.0_rk, h]
         X_c(19,:)  = [l, h]
         X_c(20,:)  = [l, h]
         X_c(21,:)  = [l, 0.0_rk]

         allocate(W_c(21), source=1.0_rk)
         W_c([2,3,5,6]) = cosd(22.5_rk)

       case default
      end select
   end subroutine
   !===============================================================================
end program surface_plate_hole_half
