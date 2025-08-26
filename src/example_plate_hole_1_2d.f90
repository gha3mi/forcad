program example_plate_hole_1_2d

   use forcad, only: rk, nurbs_surface

   implicit none
   type(nurbs_surface) :: plate_hole
   real(rk), allocatable :: Xc(:,:)
   real(rk), allocatable :: Wc(:)

   real(rk), parameter :: radius1 = 2.5_rk
   real(rk), parameter :: radius2 = 3.5_rk
   real(rk), parameter :: length  = 5.0_rk
   real(rk), parameter :: height  = 5.0_rk

   call set_Xc_Wc('ellipse', [radius1, radius2, length, height], Xc, Wc)

   call plate_hole%set(&
      knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 4.0_rk, 5.0_rk, 6.0_rk, 6.0_rk, 7.0_rk, 8.0_rk, 8.0_rk, 8.0_rk],&
      knot2 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk],&
      Xc    = Xc,&
      Wc    = Wc&
      )

   call plate_hole%create(31, 31)
   call plate_hole%export_Xc("vtk/plate_hole_1_2d_Xc.vtk")
   call plate_hole%export_Xg("vtk/plate_hole_1_2d_Xg.vtk")
   call plate_hole%export_Xth_in_Xg("vtk/plate_hole_1_2d_Xth.vtk")

   call plate_hole%show("vtk/plate_hole_1_2d_Xc.vtk","vtk/plate_hole_1_2d_Xg.vtk","vtk/plate_hole_1_2d_Xth.vtk")

contains

   !===============================================================================
   pure subroutine set_Xc_Wc(type, params, X_c, W_c)
      character(len=*), intent(in) :: type
      real(rk), intent(in), contiguous :: params(:)
      real(rk), allocatable, intent(out) :: X_c(:,:)
      real(rk), allocatable, intent(out) :: W_c(:)

      real(rk) :: r1, r2, l, h

      select case (type)
       case('circle')

         r1 = params(1)
         l  = params(3)
         h  = params(4)

         if (r1 < 0.0_rk) error stop 'Radius must be positive'
         if (l  < 0.0_rk) error stop 'Length must be positive'
         if (h  < 0.0_rk) error stop 'Height must be positive'

         allocate(X_c(39,2))
         X_c(1, :) = [-r1, 0.0_rk]
         X_c(2, :) = [-r1, r1*tand(22.5_rk)]
         X_c(3, :) = [-r1*tand(22.5_rk), r1]
         X_c(4, :) = [0.0_rk, r1]
         X_c(5, :) = [r1*tand(22.5_rk), r1]
         X_c(6, :) = [r1, r1*tand(22.5_rk)]
         X_c(7, :) = [r1, 0.0_rk]
         X_c(8, :) = [r1, -r1*tand(22.5_rk)]
         X_c(9, :) = [r1*tand(22.5_rk), -r1]
         X_c(10,:) = [0.0_rk, -r1]
         X_c(11,:) = [-r1*tand(22.5_rk), -r1]
         X_c(12,:) = [-r1, -r1*tand(22.5_rk)]
         X_c(13,:) = [-r1, 0.0_rk]
         X_c(14,:) = [-(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(15,:) = [-(r1 + (l-r1)/2.0_rk), (r1 + (h-r1)/2.0_rk)*tand(16.7_rk)]
         X_c(16,:) = [-(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), (r1 + (h-r1)/2.0_rk)]
         X_c(17,:) = [0.0_rk, (r1 + (h-r1)/2.0_rk)]
         X_c(18,:) = [(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), (r1 + (h-r1)/2.0_rk)]
         X_c(19,:) = [(r1 + (l-r1)/2.0_rk), (r1 + (h-r1)/2.0_rk)*tand(16.7_rk)]
         X_c(20,:) = [(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(21,:) = [(r1 + (l-r1)/2.0_rk), -(r1 + (h-r1)/2.0_rk)*tand(16.7_rk)]
         X_c(22,:) = [(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), -(r1 + (h-r1)/2.0_rk)]
         X_c(23,:) = [0.0_rk, -(r1 + (h-r1)/2.0_rk)]
         X_c(24,:) = [-(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), -(r1 + (h-r1)/2.0_rk)]
         X_c(25,:) = [-(r1 + (l-r1)/2.0_rk), -(r1 + (h-r1)/2.0_rk)*tand(16.7_rk)]
         X_c(26,:) = [-(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(27,:) = [-l, 0.0_rk]
         X_c(28,:) = [-l, h]
         X_c(29,:) = [-l, h]
         X_c(30,:) = [0.0_rk, h]
         X_c(31,:) = [l, h]
         X_c(32,:) = [l, h]
         X_c(33,:) = [l, 0.0_rk]
         X_c(34,:) = [l, -h]
         X_c(35,:) = [l, -h]
         X_c(36,:) = [0.0_rk, -h]
         X_c(37,:) = [-l, -h]
         X_c(38,:) = [-l, -h]
         X_c(39,:) = [-l, 0.0_rk]

         allocate(W_c(39), source=1.0_rk)
         W_c([2,3,5,6,8,9,11,12]) = (1.0_rk + 1.0_rk/sqrt(2.0_rk))/2.0_rk

       case('ellipse')

         r1 = params(1)
         r2 = params(2)
         l  = params(3)
         h  = params(4)

         if (r1 < 0.0_rk) error stop 'Radius1 must be positive'
         if (r2 < 0.0_rk) error stop 'Radius2 must be positive'
         if (l  < 0.0_rk) error stop 'Length must be positive'
         if (h  < 0.0_rk) error stop 'Height must be positive'

         allocate(X_c(39,2))
         X_c(1 ,:) = [-r1, 0.0_rk]
         X_c(2 ,:) = [-r1, r2*tand(22.5_rk)]
         X_c(3 ,:) = [-r1*tand(22.5_rk), r2]
         X_c(4 ,:) = [0.0_rk, r2]
         X_c(5 ,:) = [r1*tand(22.5_rk), r2]
         X_c(6 ,:) = [r1, r2*tand(22.5_rk)]
         X_c(7 ,:) = [r1, 0.0_rk]
         X_c(8 ,:) = [r1, -r2*tand(22.5_rk)]
         X_c(9 ,:) = [r1*tand(22.5_rk), -r2]
         X_c(10,:) = [0.0_rk, -r2]
         X_c(11,:) = [-r1*tand(22.5_rk), -r2]
         X_c(12,:) = [-r1, -r2*tand(22.5_rk)]
         X_c(13,:) = [-r1, 0.0_rk]
         X_c(14,:) = [-(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(15,:) = [-(r1 + (l-r1)/2.0_rk), (r2 + (h-r2)/2.0_rk)*tand(16.7_rk)]
         X_c(16,:) = [-(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), (r2 + (h-r2)/2.0_rk)]
         X_c(17,:) = [0.0_rk, (r2 + (h-r2)/2.0_rk)]
         X_c(18,:) = [(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), (r2 + (h-r2)/2.0_rk)]
         X_c(19,:) = [(r1 + (l-r1)/2.0_rk), (r2 + (h-r2)/2.0_rk)*tand(16.7_rk)]
         X_c(20,:) = [(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(21,:) = [(r1 + (l-r1)/2.0_rk), -(r2 + (h-r2)/2.0_rk)*tand(16.7_rk)]
         X_c(22,:) = [(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), -(r2 + (h-r2)/2.0_rk)]
         X_c(23,:) = [0.0_rk, -(r2 + (h-r2)/2.0_rk)]
         X_c(24,:) = [-(r1 + (l-r1)/2.0_rk)*tand(16.7_rk), -(r2 + (h-r2)/2.0_rk)]
         X_c(25,:) = [-(r1 + (l-r1)/2.0_rk), -(r2 + (h-r2)/2.0_rk)*tand(16.7_rk)]
         X_c(26,:) = [-(r1 + (l-r1)/2.0_rk), 0.0_rk]
         X_c(27,:) = [-l, 0.0_rk]
         X_c(28,:) = [-l, h]
         X_c(29,:) = [-l, h]
         X_c(30,:) = [0.0_rk, h]
         X_c(31,:) = [l, h]
         X_c(32,:) = [l, h]
         X_c(33,:) = [l, 0.0_rk]
         X_c(34,:) = [l, -h]
         X_c(35,:) = [l, -h]
         X_c(36,:) = [0.0_rk, -h]
         X_c(37,:) = [-l, -h]
         X_c(38,:) = [-l, -h]
         X_c(39,:) = [-l, 0.0_rk]

         allocate(W_c(39), source=1.0_rk)
         W_c([2,3,5,6,8,9,11,12]) = cosd(22.5_rk)

       case default
         error stop 'set_Xc_Wc: Invalid type. Valid types are: circle, ellipse'
      end select
   end subroutine
   !===============================================================================

end program
