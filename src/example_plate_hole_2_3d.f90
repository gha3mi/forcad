program example_plate_hole_2_3d

   use forcad, only: rk, nurbs_volume

   implicit none
   type(nurbs_volume) :: plate_hole
   real(rk), allocatable :: Xc(:,:)
   real(rk), allocatable :: Wc(:)

   real(rk), parameter :: radius1 = 2.5_rk
   real(rk), parameter :: radius2 = 3.5_rk
   real(rk), parameter :: length  = 5.0_rk
   real(rk), parameter :: height  = 5.0_rk
   real(rk), parameter :: width   = 0.5_rk

   call set_Xc_Wc('ellipse', [radius1, radius2, length, height, width], Xc, Wc)

   call plate_hole%set(&
      knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 4.0_rk],&
      knot2 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk],&
      knot3 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk],&
      Xc    = Xc,&
      Wc    = Wc&
      )

   call plate_hole%create(31, 31, 5)
   call plate_hole%export_Xc("vtk/plate_hole_2_3d_Xc.vtk")
   call plate_hole%export_Xg("vtk/plate_hole_2_3d_Xg.vtk")
   call plate_hole%export_Xth_in_Xg("vtk/plate_hole_2_3d_Xth.vtk")

   call plate_hole%show("vtk/plate_hole_2_3d_Xc.vtk","vtk/plate_hole_2_3d_Xg.vtk","vtk/plate_hole_2_3d_Xth.vtk")

contains

   !===============================================================================
   pure subroutine set_Xc_Wc(tp, params, X_c, W_c)
      character(len=*), intent(in) :: tp
      real(rk), intent(in), contiguous :: params(:)
      real(rk), allocatable, intent(out) :: X_c(:,:)
      real(rk), allocatable, intent(out) :: W_c(:)

      real(rk) :: r1, r2, l, w, h

      select case (tp)
       case('circle')

         r1 = params(1)
         l  = params(3)
         h  = params(4)
         w  = params(5)

         if (r1 < 0.0_rk) error stop 'Radius must be positive'
         if (l  < 0.0_rk) error stop 'Length must be positive'
         if (h  < 0.0_rk) error stop 'Height must be positive'
         if (w  < 0.0_rk) error stop 'Width must be positive'

         allocate(X_c(21*2, 3))
         X_c(1, 1:2) = [-r1, 0.0_rk]
         X_c(2, 1:2) = [-r1, r1*tand(22.5_rk)]
         X_c(3, 1:2) = [-r1*tand(22.5_rk), r1]
         X_c(4, 1:2) = [0.0_rk, r1]
         X_c(5, 1:2) = [r1*tand(22.5_rk), r1]
         X_c(6, 1:2) = [r1, r1*tand(22.5_rk)]
         X_c(7, 1:2) = [r1, 0.0_rk]
         X_c(8, 1:2) = [-(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(9, 1:2) = [-(r1 + (l-r1)/2.0_rk), (r1 + (h-r1)/2.0_rk)*tand(16.7_rk)]
         X_c(10,1:2) = [-(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), (r1 + (h-r1)/2.0_rk)]
         X_c(11,1:2) = [0.0_rk, (r1 + (h-r1)/2.0_rk)]
         X_c(12,1:2) = [(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), (r1 + (h-r1)/2.0_rk)]
         X_c(13,1:2) = [(r1 + (l-r1)/2.0_rk), (r1 + (h-r1)/2.0_rk)*tand(16.7_rk)]
         X_c(14,1:2) = [(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(15,1:2) = [-l, 0.0_rk]
         X_c(16,1:2) = [-l, h]
         X_c(17,1:2) = [-l, h]
         X_c(18,1:2) = [0.0_rk, h]
         X_c(19,1:2) = [l, h]
         X_c(20,1:2) = [l, h]
         X_c(21,1:2) = [l, 0.0_rk]

         X_c(22:42,1:2) = X_c(1:21,1:2)
         X_c(1:21, 3)   = 0.0_rk
         X_c(22:42,3)   = w

         allocate(W_c(21*2), source=1.0_rk)
         W_c([2,3,5,6, 23,24,26,27]) = (1.0_rk + 1.0_rk/sqrt(2.0_rk))/2.0_rk

       case('ellipse')

         r1 = params(1)
         r2 = params(2)
         l  = params(3)
         h  = params(4)
         w  = params(5)

         if (r1 < 0.0_rk) error stop 'Radius1 must be positive'
         if (r2 < 0.0_rk) error stop 'Radius2 must be positive'
         if (l  < 0.0_rk) error stop 'Length must be positive'
         if (h  < 0.0_rk) error stop 'Height must be positive'
         if (w  < 0.0_rk) error stop 'Width must be positive'

         allocate(X_c(21*2, 3))
         X_c(1 ,1:2)  = [-r1, 0.0_rk]
         X_c(2 ,1:2)  = [-r1, r2*tand(22.5_rk)]
         X_c(3 ,1:2)  = [-r1*tand(22.5_rk), r2]
         X_c(4 ,1:2)  = [0.0_rk, r2]
         X_c(5 ,1:2)  = [r1*tand(22.5_rk), r2]
         X_c(6 ,1:2)  = [r1, r2*tand(22.5_rk)]
         X_c(7 ,1:2)  = [r1, 0.0_rk]
         X_c(8 ,1:2)  = [-(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(9 ,1:2)  = [-(r1 + (l-r1)/2.0_rk), (r2 + (h-r2)/2.0_rk)*tand(16.7_rk)]
         X_c(10,1:2)  = [-(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), (r2 + (h-r2)/2.0_rk)]
         X_c(11,1:2)  = [0.0_rk, (r2 + (h-r2)/2.0_rk)]
         X_c(12,1:2)  = [(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), (r2 + (h-r2)/2.0_rk)]
         X_c(13,1:2)  = [(r1 + (l-r1)/2.0_rk), (r2 + (h-r2)/2.0_rk)*tand(16.7_rk)]
         X_c(14,1:2)  = [(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(15,1:2)  = [-l, 0.0_rk]
         X_c(16,1:2)  = [-l, h]
         X_c(17,1:2)  = [-l, h]
         X_c(18,1:2)  = [0.0_rk, h]
         X_c(19,1:2)  = [l, h]
         X_c(20,1:2)  = [l, h]
         X_c(21,1:2)  = [l, 0.0_rk]

         X_c(22:42,1:2) = X_c(1:21,1:2)
         X_c(1:21, 3)   = 0.0_rk
         X_c(22:42,3)   = w

         allocate(W_c(21*2), source=1.0_rk)
         W_c([2,3,5,6, 23,24,26,27]) = cosd(22.5_rk)

       case default
         error stop 'set_Xc_Wc: Invalid type. Valid types are: circle, ellipse'
      end select
   end subroutine
   !===============================================================================

end program
