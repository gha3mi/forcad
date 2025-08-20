program test_nurbs_curve
   use forcad_kinds, only: rk
   use forcad_nurbs_curve, only: nurbs_curve, compute_Tgc, compute_dTgc
   use forunittest, only: unit_tests
   implicit none

   integer, parameter :: NTESTS = 113
   type(unit_tests) :: ut
   real(rk), parameter :: TOL = 1.0e-5_rk
   real(rk), parameter :: PI = acos(-1.0_rk)

   type(nurbs_curve) :: nurbs, bsp
   real(rk), allocatable :: Xc(:,:), Wc(:)
   real(rk), allocatable :: Xg(:,:), Xgb(:,:)
   real(rk) :: knot(6)
   integer, allocatable :: elemConn(:,:)
   real(rk), allocatable :: Tgc(:,:), dTgc(:,:), Tgcb(:,:), dTgcb(:,:), d2Tgc(:,:), d2Tgcb(:,:)
   real(rk), allocatable :: Tgc1(:), dTgc1(:), Tgc1b(:), dTgc1b(:), d2Tgc1(:), d2Tgc1b(:)
   real(rk) :: nearest_Xg(3), nearest_Xt, length, lengthb
   integer :: i, id, ti
   real(rk), allocatable :: Xt(:)
   character(len=*), parameter :: fXc_nurbs = 'vtk/test_nurbs_curve_Xc.vtk'
   character(len=*), parameter :: fXg_nurbs = 'vtk/test_nurbs_curve_Xg.vtk'
   character(len=*), parameter :: fXth_nurbs = 'vtk/test_nurbs_curve_Xth.vtk'
   character(len=*), parameter :: fIgs_nurbs = 'iges/test_nurbs_curve.iges'
   character(len=*), parameter :: fXc_bsp = 'vtk/test_bsp_curve_Xc.vtk'
   character(len=*), parameter :: fXg_bsp = 'vtk/test_bsp_curve_Xg.vtk'
   character(len=*), parameter :: fXth_bsp = 'vtk/test_bsp_curve_Xth.vtk'
   character(len=*), parameter :: fIgs_bsp = 'iges/test_bsp_curve.iges'
   logical :: okXc_nurbs, okXg_nurbs, okXth_nurbs, okIges_nurbs
   logical :: okXc_bsp, okXg_bsp, okXth_bsp, okIges_bsp

   call ut%initialize(NTESTS)
   ti = 1

   ! Initialize curve data
   allocate(Xc(3, 3))
   Xc(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
   Xc(2,:) = [1.0_rk, 0.0_rk, 0.0_rk]
   Xc(3,:) = [2.0_rk, 0.0_rk, 0.0_rk]
   allocate(Wc(3))
   Wc = [1.0_rk, 0.9_rk, 0.8_rk]
   knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

   ! 1) Set and create NURBS and B-spline curves
   call nurbs%set(knot, Xc, Wc)
   call bsp%set(knot, Xc)
   call nurbs%set(degree=nurbs%get_degree(), nc=nurbs%get_nc(), Xc=nurbs%get_Xc(), Wc=nurbs%get_Wc())
   call bsp%set(degree=bsp%get_degree(), nc=bsp%get_nc(), Xc=bsp%get_Xc())
   call nurbs%create(res=23)
   call bsp%create(res=23)


   call ut%test(ti)%check( &
      name="set(): degree==2", &
      res=nurbs%get_degree(), &
      expected=2, &
      msg="NURBS degree not 2", &
      group="setup"); ti=ti+1
   call ut%test(ti)%check( &
      name="set(): degree==2 (B-spline)", &
      res=bsp%get_degree(), &
      expected=2, &
      msg="B-spline degree not 2", &
      group="setup"); ti=ti+1
   call ut%test(ti)%check( &
      name="set(): nc==3", &
      res=nurbs%get_nc(), &
      expected=3, &
      msg="NURBS nc not 3", &
      group="setup"); ti=ti+1
   call ut%test(ti)%check( &
      name="set(): nc==3 (B-spline)", &
      res=bsp%get_nc(), &
      expected=3, &
      msg="B-spline nc not 3", &
      group="setup"); ti=ti+1

   call ut%test(ti)%check( &
      name="set(): knot matches", &
      res=nurbs%get_knot(), &
      expected=knot, &
      tol=TOL, &
      msg="NURBS knot vector mismatch", &
      group="setup"); ti=ti+1

      call ut%test(ti)%check( &
      name="set(): knot matches (B-spline)", &
      res=bsp%get_knot(), &
      expected=knot, &
      tol=TOL, &
      msg="B-spline knot vector mismatch", &
      group="setup"); ti=ti+1

   ! 2) Export tests
   call nurbs%export_Xc(fXc_nurbs)
   call bsp%export_Xc(fXc_bsp)
   call nurbs%export_Xg(fXg_nurbs)
   call bsp%export_Xg(fXg_bsp)
   call nurbs%export_Xth(fXth_nurbs)
   call bsp%export_Xth(fXth_bsp)
   call nurbs%export_iges(fIgs_nurbs)
   call bsp%export_iges(fIgs_bsp)
   inquire(file=fXc_nurbs, exist=okXc_nurbs)
   inquire(file=fXg_nurbs, exist=okXg_nurbs)
   inquire(file=fXth_nurbs, exist=okXth_nurbs)
   inquire(file=fIgs_nurbs, exist=okIges_nurbs)
   inquire(file=fXc_bsp, exist=okXc_bsp)
   inquire(file=fXg_bsp, exist=okXg_bsp)
   inquire(file=fXth_bsp, exist=okXth_bsp)
   inquire(file=fIgs_bsp, exist=okIges_bsp)
   call ut%test(ti)%check( &
      name="export_Xc/Xg/Xth/IGES(): files exist (NURBS)", &
      res=merge(1,0,okXc_nurbs .and. okXg_nurbs .and. okXth_nurbs .and. okIges_nurbs), &
      expected=1, &
      msg="NURBS export did not create all files", &
      group="io"); ti=ti+1
   call ut%test(ti)%check( &
      name="export_Xc/Xg/Xth/IGES(): files exist (B-spline)", &
      res=merge(1,0,okXc_bsp .and. okXg_bsp .and. okXth_bsp .and. okIges_bsp), &
      expected=1, &
      msg="B-spline export did not create all files", &
      group="io"); ti=ti+1

   ! 3) Length tests
   call nurbs%cmp_length(length)
   call bsp%cmp_length(lengthb)
   call ut%test(ti)%check( &
      name="cmp_length(): length == 2 (NURBS)", &
      res=length, &
      expected=2.0_rk, &
      tol=TOL, &
      msg="NURBS arc length incorrect", &
      group="geometry"); ti=ti+1
   call ut%test(ti)%check( &
      name="cmp_length(): length == 2 (B-spline)", &
      res=lengthb, &
      expected=2.0_rk, &
      tol=TOL, &
      msg="B-spline arc length incorrect", &
      group="geometry"); ti=ti+1

   ! 4) Nearest point tests
   call nurbs%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
   call ut%test(ti)%check( &
      name="nearest_point(): point matches [0,0,0] (NURBS)", &
      res=nearest_Xg, &
      expected=[0.0_rk, 0.0_rk, 0.0_rk], &
      tol=TOL, &
      msg="NURBS nearest point mismatch", &
      group="nearest"); ti=ti+1
   call ut%test(ti)%check( &
      name="nearest_point(): param == 0 (NURBS)", &
      res=nearest_Xt, &
      expected=0.0_rk, &
      tol=TOL, &
      msg="NURBS nearest param incorrect", &
      group="nearest"); ti=ti+1
   call ut%test(ti)%check( &
      name="nearest_point(): id == 1 (NURBS)", &
      res=id, &
      expected=1, &
      msg="NURBS nearest id incorrect", &
      group="nearest"); ti=ti+1

   call bsp%nearest_point([0.0_rk, 0.0_rk, 0.5_rk], nearest_Xg, nearest_Xt, id)
   call ut%test(ti)%check( &
      name="nearest_point(): point matches [0,0,0] (B-spline)", &
      res=nearest_Xg, &
      expected=[0.0_rk, 0.0_rk, 0.0_rk], &
      tol=TOL, &
      msg="B-spline nearest point mismatch", &
      group="nearest"); ti=ti+1
   call ut%test(ti)%check( &
      name="nearest_point(): param == 0 (B-spline)", &
      res=nearest_Xt, &
      expected=0.0_rk, &
      tol=TOL, &
      msg="B-spline nearest param incorrect", &
      group="nearest"); ti=ti+1
   call ut%test(ti)%check( &
      name="nearest_point(): id == 1 (B-spline)", &
      res=id, &
      expected=1, &
      msg="B-spline nearest id incorrect", &
      group="nearest"); ti=ti+1

   call nurbs%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1.0e-10_rk, 10, nearest_Xt, nearest_Xg)
   call ut%test(ti)%check( &
      name="nearest_point2(): point matches [0,0,0] (NURBS)", &
      res=nearest_Xg, &
      expected=[0.0_rk, 0.0_rk, 0.0_rk], &
      tol=TOL, &
      msg="NURBS nearest_point2 point mismatch", &
      group="nearest"); ti=ti+1
   call ut%test(ti)%check( &
      name="nearest_point2(): param == 0 (NURBS)", &
      res=nearest_Xt, &
      expected=0.0_rk, &
      tol=TOL, &
      msg="NURBS nearest_point2 param incorrect", &
      group="nearest"); ti=ti+1

   call bsp%nearest_point2([0.0_rk, 0.0_rk, 0.5_rk], 1.0e-10_rk, 10, nearest_Xt, nearest_Xg)
   call ut%test(ti)%check( &
      name="nearest_point2(): point matches [0,0,0] (B-spline)", &
      res=nearest_Xg, &
      expected=[0.0_rk, 0.0_rk, 0.0_rk], &
      tol=TOL, &
      msg="B-spline nearest_point2 point mismatch", &
      group="nearest"); ti=ti+1
   call ut%test(ti)%check( &
      name="nearest_point2(): param == 0 (B-spline)", &
      res=nearest_Xt, &
      expected=0.0_rk, &
      tol=TOL, &
      msg="B-spline nearest_point2 param incorrect", &
      group="nearest"); ti=ti+1

   ! 5) Store initial geometry
   Xg = nurbs%get_Xg()
   Xgb = bsp%get_Xg()

   ! 6) Set with Xth_dir and continuity
   call nurbs%set([0.0_rk, 1.0_rk], 2, [-1, -1], Xc, Wc)
   call bsp%set([0.0_rk, 1.0_rk], 2, [-1, -1], Xc)
   call nurbs%create(res=23)
   call bsp%create(res=23)
   call ut%test(ti)%check( &
      name="set(Xth_dir): geometry preserved (NURBS)", &
      res=nurbs%get_Xg(), &
      expected=Xg, &
      tol=TOL, &
      msg="NURBS geometry changed after set(Xth_dir)", &
      group="setup"); ti=ti+1
   call ut%test(ti)%check( &
      name="set(Xth_dir): geometry preserved (B-spline)", &
      res=bsp%get_Xg(), &
      expected=Xgb, &
      tol=TOL, &
      msg="B-spline geometry changed after set(Xth_dir)", &
      group="setup"); ti=ti+1

   ! 7) Set with Xc and Wc
   call nurbs%set(Xc, Wc)
   call bsp%set(Xc)
   call nurbs%create(res=23)
   call bsp%create(res=23)
   call ut%test(ti)%check( &
      name="set(Xc,Wc): geometry preserved (NURBS)", &
      res=nurbs%get_Xg(), &
      expected=Xg, &
      tol=TOL, &
      msg="NURBS geometry changed after set(Xc,Wc)", &
      group="setup"); ti=ti+1
   call ut%test(ti)%check( &
      name="set(Xc): geometry preserved (B-spline)", &
      res=bsp%get_Xg(), &
      expected=Xgb, &
      tol=TOL, &
      msg="B-spline geometry changed after set(Xc)", &
      group="setup"); ti=ti+1

   ! 8) Create with explicit Xt
   call nurbs%create(Xt=nurbs%get_Xt())
   call bsp%create(Xt=bsp%get_Xt())
   call ut%test(ti)%check( &
      name="create(Xt): geometry preserved (NURBS)", &
      res=nurbs%get_Xg(), &
      expected=Xg, &
      tol=TOL, &
      msg="NURBS geometry changed after create(Xt)", &
      group="sampling"); ti=ti+1
   call ut%test(ti)%check( &
      name="create(Xt): geometry preserved (B-spline)", &
      res=bsp%get_Xg(), &
      expected=Xgb, &
      tol=TOL, &
      msg="B-spline geometry changed after create(Xt)", &
      group="sampling"); ti=ti+1

   ! 9) Export after create(Xt)
   call nurbs%export_Xc(fXc_nurbs)
   call bsp%export_Xc(fXc_bsp)
   call nurbs%export_Xg(fXg_nurbs)
   call bsp%export_Xg(fXg_bsp)
   call nurbs%export_Xth(fXth_nurbs)
   call bsp%export_Xth(fXth_bsp)
   call nurbs%export_iges(fIgs_nurbs)
   call bsp%export_iges(fIgs_bsp)
   inquire(file=fXc_nurbs, exist=okXc_nurbs)
   inquire(file=fXg_nurbs, exist=okXg_nurbs)
   inquire(file=fXth_nurbs, exist=okXth_nurbs)
   inquire(file=fIgs_nurbs, exist=okIges_nurbs)
   inquire(file=fXc_bsp, exist=okXc_bsp)
   inquire(file=fXg_bsp, exist=okXg_bsp)
   inquire(file=fXth_bsp, exist=okXth_bsp)
   inquire(file=fIgs_bsp, exist=okIges_bsp)
   call ut%test(ti)%check( &
      name="export after create(Xt): files exist (NURBS)", &
      res=merge(1,0,okXc_nurbs .and. okXg_nurbs .and. okXth_nurbs .and. okIges_nurbs), &
      expected=1, &
      msg="NURBS export after create(Xt) did not create all files", &
      group="io"); ti=ti+1
   call ut%test(ti)%check( &
      name="export after create(Xt): files exist (B-spline)", &
      res=merge(1,0,okXc_bsp .and. okXg_bsp .and. okXth_bsp .and. okIges_bsp), &
      expected=1, &
      msg="B-spline export after create(Xt) did not create all files", &
      group="io"); ti=ti+1

   ! 10) Getter tests
   call ut%test(ti)%check( &
      name="get_Xc(): matches input (NURBS)", &
      res=nurbs%get_Xc(), &
      expected=Xc, &
      tol=TOL, &
      msg="NURBS Xc mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Xc(): matches input (B-spline)", &
      res=bsp%get_Xc(), &
      expected=Xc, &
      tol=TOL, &
      msg="B-spline Xc mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Xc(1): matches input (NURBS)", &
      res=nurbs%get_Xc(1), &
      expected=Xc(1,:), &
      tol=TOL, &
      msg="NURBS Xc(1) mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Xc(1): matches input (B-spline)", &
      res=bsp%get_Xc(1), &
      expected=Xc(1,:), &
      tol=TOL, &
      msg="B-spline Xc(1) mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Xc(1,1): matches input (NURBS)", &
      res=nurbs%get_Xc(1,1), &
      expected=Xc(1,1), &
      tol=TOL, &
      msg="NURBS Xc(1,1) mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Xc(1,1): matches input (B-spline)", &
      res=bsp%get_Xc(1,1), &
      expected=Xc(1,1), &
      tol=TOL, &
      msg="B-spline Xc(1,1) mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Xg(): matches initial (NURBS)", &
      res=nurbs%get_Xg(), &
      expected=Xg, &
      tol=TOL, &
      msg="NURBS Xg mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Xg(): matches initial (B-spline)", &
      res=bsp%get_Xg(), &
      expected=Xgb, &
      tol=TOL, &
      msg="B-spline Xg mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Xg(1): matches initial (NURBS)", &
      res=nurbs%get_Xg(1), &
      expected=Xg(1,:), &
      tol=TOL, &
      msg="NURBS Xg(1) mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Xg(1): matches initial (B-spline)", &
      res=bsp%get_Xg(1), &
      expected=Xgb(1,:), &
      tol=TOL, &
      msg="B-spline Xg(1) mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Xg(1,1): matches initial (NURBS)", &
      res=nurbs%get_Xg(1,1), &
      expected=Xg(1,1), &
      tol=TOL, &
      msg="NURBS Xg(1,1) mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Xg(1,1): matches initial (B-spline)", &
      res=bsp%get_Xg(1,1), &
      expected=Xgb(1,1), &
      tol=TOL, &
      msg="B-spline Xg(1,1) mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Wc(): matches input (NURBS)", &
      res=nurbs%get_Wc(), &
      expected=Wc, &
      tol=TOL, &
      msg="NURBS Wc mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_Wc(1): matches input (NURBS)", &
      res=nurbs%get_Wc(1), &
      expected=Wc(1), &
      tol=TOL, &
      msg="NURBS Wc(1) mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_knot(): matches input (NURBS)", &
      res=nurbs%get_knot(), &
      expected=knot, &
      tol=TOL, &
      msg="NURBS knot mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_knot(): matches input (B-spline)", &
      res=bsp%get_knot(), &
      expected=knot, &
      tol=TOL, &
      msg="B-spline knot mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_knot(1): matches input (NURBS)", &
      res=nurbs%get_knot(1), &
      expected=knot(1), &
      tol=TOL, &
      msg="NURBS knot(1) mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_knot(1): matches input (B-spline)", &
      res=bsp%get_knot(1), &
      expected=knot(1), &
      tol=TOL, &
      msg="B-spline knot(1) mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_ng(): matches Xg size (NURBS)", &
      res=nurbs%get_ng(), &
      expected=size(Xg,1), &
      msg="NURBS ng mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_ng(): matches Xgb size (B-spline)", &
      res=bsp%get_ng(), &
      expected=size(Xgb,1), &
      msg="B-spline ng mismatch", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_degree(): matches 2 (NURBS)", &
      res=nurbs%get_degree(), &
      expected=2, &
      msg="NURBS degree not 2", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_degree(): matches 2 (B-spline)", &
      res=bsp%get_degree(), &
      expected=2, &
      msg="B-spline degree not 2", &
      group="getters"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_multiplicity(): matches [3,3] (NURBS)", &
      res=nurbs%get_multiplicity(), &
      expected=[3,3], &
      msg="NURBS multiplicity not [3,3]", &
      group="knot-ops"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_multiplicity(): matches [3,3] (B-spline)", &
      res=bsp%get_multiplicity(), &
      expected=[3,3], &
      msg="B-spline multiplicity not [3,3]", &
      group="knot-ops"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_continuity(): matches [-1,-1] (NURBS)", &
      res=nurbs%get_continuity(), &
      expected=[-1,-1], &
      msg="NURBS continuity not [-1,-1]", &
      group="knot-ops"); ti=ti+1
   call ut%test(ti)%check( &
      name="get_continuity(): matches [-1,-1] (B-spline)", &
      res=bsp%get_continuity(), &
      expected=[-1,-1], &
      msg="B-spline continuity not [-1,-1]", &
      group="knot-ops"); ti=ti+1
   call nurbs%cmp_nc()
   call bsp%cmp_nc()
   call ut%test(ti)%check( &
      name="cmp_nc(): matches nc (NURBS)", &
      res=nurbs%get_nc(), &
      expected=size(Xc,1), &
      msg="NURBS cmp_nc mismatch", &
      group="knot-ops"); ti=ti+1
   call ut%test(ti)%check( &
      name="cmp_nc(): matches nc (B-spline)", &
      res=bsp%get_nc(), &
      expected=size(Xc,1), &
      msg="B-spline cmp_nc mismatch", &
      group="knot-ops"); ti=ti+1

   ! 11) Element connectivity tests
   elemConn = nurbs%cmp_elem_Xc_vis(2)
   call nurbs%set_elem_Xc_vis(elemConn)
   call ut%test(ti)%check( &
      name="set/get elem_Xc_vis(p=2): equality (NURBS)", &
      res=nurbs%get_elem_Xc_vis(), &
      expected=elemConn, &
      msg="NURBS elem_Xc_vis(p=2) mismatch", &
      group="elements"); ti=ti+1
   deallocate(elemConn)
   elemConn = nurbs%cmp_elem_Xc_vis()
   call nurbs%set_elem_Xc_vis(elemConn)
   call ut%test(ti)%check( &
      name="set/get elem_Xc_vis(): equality (NURBS)", &
      res=nurbs%get_elem_Xc_vis(), &
      expected=elemConn, &
      msg="NURBS elem_Xc_vis mismatch", &
      group="elements"); ti=ti+1
   deallocate(elemConn)
   elemConn = bsp%cmp_elem_Xc_vis(2)
   call bsp%set_elem_Xc_vis(elemConn)
   call ut%test(ti)%check( &
      name="set/get elem_Xc_vis(p=2): equality (B-spline)", &
      res=bsp%get_elem_Xc_vis(), &
      expected=elemConn, &
      msg="B-spline elem_Xc_vis(p=2) mismatch", &
      group="elements"); ti=ti+1
   deallocate(elemConn)
   elemConn = bsp%cmp_elem_Xc_vis()
   call bsp%set_elem_Xc_vis(elemConn)
   call ut%test(ti)%check( &
      name="set/get elem_Xc_vis(): equality (B-spline)", &
      res=bsp%get_elem_Xc_vis(), &
      expected=elemConn, &
      msg="B-spline elem_Xc_vis mismatch", &
      group="elements"); ti=ti+1
   deallocate(elemConn)
   elemConn = nurbs%cmp_elem_Xg_vis(2)
   call nurbs%set_elem_Xg_vis(elemConn)
   call ut%test(ti)%check( &
      name="set/get elem_Xg_vis(p=2): equality (NURBS)", &
      res=nurbs%get_elem_Xg_vis(), &
      expected=elemConn, &
      msg="NURBS elem_Xg_vis(p=2) mismatch", &
      group="elements"); ti=ti+1
   deallocate(elemConn)
   elemConn = nurbs%cmp_elem_Xg_vis()
   call nurbs%set_elem_Xg_vis(elemConn)
   call ut%test(ti)%check( &
      name="set/get elem_Xg_vis(): equality (NURBS)", &
      res=nurbs%get_elem_Xg_vis(), &
      expected=elemConn, &
      msg="NURBS elem_Xg_vis mismatch", &
      group="elements"); ti=ti+1
   deallocate(elemConn)
   elemConn = bsp%cmp_elem_Xg_vis(2)
   call bsp%set_elem_Xg_vis(elemConn)
   call ut%test(ti)%check( &
      name="set/get elem_Xg_vis(p=2): equality (B-spline)", &
      res=bsp%get_elem_Xg_vis(), &
      expected=elemConn, &
      msg="B-spline elem_Xg_vis(p=2) mismatch", &
      group="elements"); ti=ti+1
   deallocate(elemConn)
   elemConn = bsp%cmp_elem_Xg_vis()
   call bsp%set_elem_Xg_vis(elemConn)
   call ut%test(ti)%check( &
      name="set/get elem_Xg_vis(): equality (B-spline)", &
      res=bsp%get_elem_Xg_vis(), &
      expected=elemConn, &
      msg="B-spline elem_Xg_vis mismatch", &
      group="elements"); ti=ti+1
   deallocate(elemConn)
   elemConn = nurbs%cmp_elem()
   call nurbs%set_elem(elemConn)
   call ut%test(ti)%check( &
      name="set/get elem: equality (NURBS)", &
      res=nurbs%get_elem(), &
      expected=elemConn, &
      msg="NURBS elem mismatch", &
      group="elements"); ti=ti+1
   deallocate(elemConn)
   elemConn = bsp%cmp_elem()
   call bsp%set_elem(elemConn)
   call ut%test(ti)%check( &
      name="set/get elem: equality (B-spline)", &
      res=bsp%get_elem(), &
      expected=elemConn, &
      msg="B-spline elem mismatch", &
      group="elements"); ti=ti+1
   deallocate(elemConn)

   ! 12) Modify and export
   call nurbs%modify_Xc(Xc(1,1), 1, 1)
   call bsp%modify_Xc(Xc(1,1), 1, 1)
   call nurbs%modify_Wc(Wc(1), 1)
   call nurbs%create()
   call bsp%create()
   call nurbs%export_Xc(fXc_nurbs)
   call bsp%export_Xc(fXc_bsp)
   call nurbs%export_Xg(fXg_nurbs)
   call bsp%export_Xg(fXg_bsp)
   call ut%test(ti)%check( &
      name="modify_Xc/Wc: geometry preserved (NURBS)", &
      res=nurbs%get_Xg(), &
      expected=Xg, &
      tol=TOL, &
      msg="NURBS geometry changed after modify_Xc/Wc", &
      group="modify"); ti=ti+1
   call ut%test(ti)%check( &
      name="modify_Xc: geometry preserved (B-spline)", &
      res=bsp%get_Xg(), &
      expected=Xgb, &
      tol=TOL, &
      msg="B-spline geometry changed after modify_Xc", &
      group="modify"); ti=ti+1

   ! 13) Basis and derivative tests
   call nurbs%basis(res=23, Tgc=Tgc)
   call bsp%basis(res=23, Tgc=Tgcb)
   call ut%test(ti)%check( &
      name="basis(res=23): sum(N)=1 rows (NURBS)", &
      res=all(abs(sum(Tgc,dim=2)-1.0_rk) <= TOL), &
      expected=.true., &
      msg="NURBS basis partition of unity violated", &
      group="basis"); ti=ti+1
   call ut%test(ti)%check( &
      name="basis(res=23): sum(N)=1 rows (B-spline)", &
      res=all(abs(sum(Tgcb,dim=2)-1.0_rk) <= TOL), &
      expected=.true., &
      msg="B-spline basis partition of unity violated", &
      group="basis"); ti=ti+1

   call nurbs%basis(Xt=0.0_rk, Tgc=Tgc1)
   call bsp%basis(Xt=0.0_rk, Tgc=Tgc1b)
   call ut%test(ti)%check( &
      name="basis(Xt=0): sum(N)=1 (NURBS)", &
      res=sum(Tgc1), &
      expected=1.0_rk, &
      tol=TOL, &
      msg="NURBS basis sum not 1 at t=0", &
      group="basis"); ti=ti+1
   call ut%test(ti)%check( &
      name="basis(Xt=0): sum(N)=1 (B-spline)", &
      res=sum(Tgc1b), &
      expected=1.0_rk, &
      tol=TOL, &
      msg="B-spline basis sum not 1 at t=0", &
      group="basis"); ti=ti+1

   allocate(Xt(23))
   do i = 1, 23
      Xt(i) = real(i-1,rk) / real(23-1,rk)
   end do
   call nurbs%basis(Xt=Xt, Tgc=Tgc)
   call bsp%basis(Xt=Xt, Tgc=Tgcb)
   call ut%test(ti)%check( &
      name="basis(Xt vector): sum(N)=1 rows (NURBS)", &
      res=all(abs(sum(Tgc,dim=2)-1.0_rk) <= TOL), &
      expected=.true., &
      msg="NURBS basis partition of unity violated (vector)", &
      group="basis"); ti=ti+1
   call ut%test(ti)%check( &
      name="basis(Xt vector): sum(N)=1 rows (B-spline)", &
      res=all(abs(sum(Tgcb,dim=2)-1.0_rk) <= TOL), &
      expected=.true., &
      msg="B-spline basis partition of unity violated (vector)", &
      group="basis"); ti=ti+1

   call nurbs%derivative(res=23, dTgc=dTgc, Tgc=Tgc)
   call bsp%derivative(res=23, dTgc=dTgcb, Tgc=Tgcb)
   call ut%test(ti)%check( &
      name="derivative(res=23): sum(dN)=0 rows (NURBS)", &
      res=all(abs(sum(dTgc,dim=2)) <= TOL), &
      expected=.true., &
      msg="NURBS derivative sum not zero", &
      group="derivatives"); ti=ti+1
   call ut%test(ti)%check( &
      name="derivative(res=23): sum(dN)=0 rows (B-spline)", &
      res=all(abs(sum(dTgcb,dim=2)) <= TOL), &
      expected=.true., &
      msg="B-spline derivative sum not zero", &
      group="derivatives"); ti=ti+1

   call nurbs%derivative(Xt=Xt, dTgc=dTgc, Tgc=Tgc)
   call bsp%derivative(Xt=Xt, dTgc=dTgcb, Tgc=Tgcb)
   call ut%test(ti)%check( &
      name="derivative(Xt vector): sum(dN)=0 rows (NURBS)", &
      res=all(abs(sum(dTgc,dim=2)) <= TOL), &
      expected=.true., &
      msg="NURBS derivative sum not zero (vector)", &
      group="derivatives"); ti=ti+1
   call ut%test(ti)%check( &
      name="derivative(Xt vector): sum(dN)=0 rows (B-spline)", &
      res=all(abs(sum(dTgcb,dim=2)) <= TOL), &
      expected=.true., &
      msg="B-spline derivative sum not zero (vector)", &
      group="derivatives"); ti=ti+1

   call nurbs%derivative(Xt=0.0_rk, dTgc=dTgc1, Tgc=Tgc1)
   call bsp%derivative(Xt=0.0_rk, dTgc=dTgc1b, Tgc=Tgc1b)
   call ut%test(ti)%check( &
      name="derivative(Xt=0): sum(dN)=0 (NURBS)", &
      res=sum(dTgc1), &
      expected=0.0_rk, &
      tol=TOL, &
      msg="NURBS derivative sum not zero at t=0", &
      group="derivatives"); ti=ti+1
   call ut%test(ti)%check( &
      name="derivative(Xt=0): sum(dN)=0 (B-spline)", &
      res=sum(dTgc1b), &
      expected=0.0_rk, &
      tol=TOL, &
      msg="B-spline derivative sum not zero at t=0", &
      group="derivatives"); ti=ti+1

   call nurbs%derivative(Xt=0.0_rk, dTgc=dTgc1, Tgc=Tgc1, elem=[1,2,3])
   call bsp%derivative(Xt=0.0_rk, dTgc=dTgc1b, Tgc=Tgc1b, elem=[1,2,3])
   call ut%test(ti)%check( &
      name="derivative(Xt=0, elem): sum(dN)=0 (NURBS)", &
      res=sum(dTgc1), &
      expected=0.0_rk, &
      tol=TOL, &
      msg="NURBS derivative sum not zero (elem)", &
      group="derivatives"); ti=ti+1
   call ut%test(ti)%check( &
      name="derivative(Xt=0, elem): sum(dN)=0 (B-spline)", &
      res=sum(dTgc1b), &
      expected=0.0_rk, &
      tol=TOL, &
      msg="B-spline derivative sum not zero (elem)", &
      group="derivatives"); ti=ti+1

   call nurbs%derivative2(res=23, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
   call bsp%derivative2(res=23, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)
   call ut%test(ti)%check( &
      name="derivative2(res=23): sum(d2N)=0 rows (NURBS)", &
      res=all(abs(sum(d2Tgc,dim=2)) <= TOL), &
      expected=.true., &
      msg="NURBS second derivative sum not zero", &
      group="derivatives"); ti=ti+1
   call ut%test(ti)%check( &
      name="derivative2(res=23): sum(d2N)=0 rows (B-spline)", &
      res=all(abs(sum(d2Tgcb,dim=2)) <= TOL), &
      expected=.true., &
      msg="B-spline second derivative sum not zero", &
      group="derivatives"); ti=ti+1

   call nurbs%derivative2(Xt=Xt, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)
   call bsp%derivative2(Xt=Xt, d2Tgc=d2Tgcb, dTgc=dTgcb, Tgc=Tgcb)
   call ut%test(ti)%check( &
      name="derivative2(Xt vector): sum(d2N)=0 rows (NURBS)", &
      res=all(abs(sum(d2Tgc,dim=2)) <= TOL), &
      expected=.true., &
      msg="NURBS second derivative sum not zero (vector)", &
      group="derivatives"); ti=ti+1
   call ut%test(ti)%check( &
      name="derivative2(Xt vector): sum(d2N)=0 rows (B-spline)", &
      res=all(abs(sum(d2Tgcb,dim=2)) <= TOL), &
      expected=.true., &
      msg="B-spline second derivative sum not zero (vector)", &
      group="derivatives"); ti=ti+1

   call nurbs%derivative2(Xt=0.0_rk, d2Tgc=d2Tgc1, dTgc=dTgc1, Tgc=Tgc1)
   call bsp%derivative2(Xt=0.0_rk, d2Tgc=d2Tgc1b, dTgc=dTgc1b, Tgc=Tgc1b)
   call ut%test(ti)%check( &
      name="derivative2(Xt=0): sum(d2N)=0 (NURBS)", &
      res=sum(d2Tgc1), &
      expected=0.0_rk, &
      tol=TOL, &
      msg="NURBS second derivative sum not zero at t=0", &
      group="derivatives"); ti=ti+1
   call ut%test(ti)%check( &
      name="derivative2(Xt=0): sum(d2N)=0 (B-spline)", &
      res=sum(d2Tgc1b), &
      expected=0.0_rk, &
      tol=TOL, &
      msg="B-spline second derivative sum not zero at t=0", &
      group="derivatives"); ti=ti+1

   ! 14) Rotation tests (Xc)
   call nurbs%rotate_Xc(45.0_rk, 0.0_rk, 0.0_rk)
   call nurbs%rotate_Xc(-45.0_rk, 0.0_rk, 0.0_rk)
   call bsp%rotate_Xc(45.0_rk, 0.0_rk, 0.0_rk)
   call bsp%rotate_Xc(-45.0_rk, 0.0_rk, 0.0_rk)
   call ut%test(ti)%check( &
      name="rotate_Xc(X-axis): geometry preserved (NURBS)", &
      res=nurbs%get_Xc(), &
      expected=Xc, &
      tol=TOL, &
      msg="NURBS Xc changed after X-axis rotation", &
      group="transform"); ti=ti+1
   call ut%test(ti)%check( &
      name="rotate_Xc(X-axis): geometry preserved (B-spline)", &
      res=bsp%get_Xc(), &
      expected=Xc, &
      tol=TOL, &
      msg="B-spline Xc changed after X-axis rotation", &
      group="transform"); ti=ti+1

   call nurbs%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
   call nurbs%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)
   call bsp%rotate_Xc(0.0_rk, 45.0_rk, 0.0_rk)
   call bsp%rotate_Xc(0.0_rk, -45.0_rk, 0.0_rk)
   call ut%test(ti)%check( &
      name="rotate_Xc(Y-axis): geometry preserved (NURBS)", &
      res=nurbs%get_Xc(), &
      expected=Xc, &
      tol=TOL, &
      msg="NURBS Xc changed after Y-axis rotation", &
      group="transform"); ti=ti+1
   call ut%test(ti)%check( &
      name="rotate_Xc(Y-axis): geometry preserved (B-spline)", &
      res=bsp%get_Xc(), &
      expected=Xc, &
      tol=TOL, &
      msg="B-spline Xc changed after Y-axis rotation", &
      group="transform"); ti=ti+1

   call nurbs%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
   call nurbs%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)
   call bsp%rotate_Xc(0.0_rk, 0.0_rk, 45.0_rk)
   call bsp%rotate_Xc(0.0_rk, 0.0_rk, -45.0_rk)
   call ut%test(ti)%check( &
      name="rotate_Xc(Z-axis): geometry preserved (NURBS)", &
      res=nurbs%get_Xc(), &
      expected=Xc, &
      tol=TOL, &
      msg="NURBS Xc changed after Z-axis rotation", &
      group="transform"); ti=ti+1
   call ut%test(ti)%check( &
      name="rotate_Xc(Z-axis): geometry preserved (B-spline)", &
      res=bsp%get_Xc(), &
      expected=Xc, &
      tol=TOL, &
      msg="B-spline Xc changed after Z-axis rotation", &
      group="transform"); ti=ti+1

   ! 15) Rotation tests (Xg)
   call nurbs%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
   call nurbs%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)
   call bsp%rotate_Xg(45.0_rk, 0.0_rk, 0.0_rk)
   call bsp%rotate_Xg(-45.0_rk, 0.0_rk, 0.0_rk)
   call ut%test(ti)%check( &
      name="rotate_Xg(X-axis): geometry preserved (NURBS)", &
      res=nurbs%get_Xg(), &
      expected=Xg, &
      tol=TOL, &
      msg="NURBS Xg changed after X-axis rotation", &
      group="transform"); ti=ti+1
   call ut%test(ti)%check( &
      name="rotate_Xg(X-axis): geometry preserved (B-spline)", &
      res=bsp%get_Xg(), &
      expected=Xgb, &
      tol=TOL, &
      msg="B-spline Xg changed after X-axis rotation", &
      group="transform"); ti=ti+1

   call nurbs%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
   call nurbs%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)
   call bsp%rotate_Xg(0.0_rk, 45.0_rk, 0.0_rk)
   call bsp%rotate_Xg(0.0_rk, -45.0_rk, 0.0_rk)
   call ut%test(ti)%check( &
      name="rotate_Xg(Y-axis): geometry preserved (NURBS)", &
      res=nurbs%get_Xg(), &
      expected=Xg, &
      tol=TOL, &
      msg="NURBS Xg changed after Y-axis rotation", &
      group="transform"); ti=ti+1
   call ut%test(ti)%check( &
      name="rotate_Xg(Y-axis): geometry preserved (B-spline)", &
      res=bsp%get_Xg(), &
      expected=Xgb, &
      tol=TOL, &
      msg="B-spline Xg changed after Y-axis rotation", &
      group="transform"); ti=ti+1

   call nurbs%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
   call nurbs%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)
   call bsp%rotate_Xg(0.0_rk, 0.0_rk, 45.0_rk)
   call bsp%rotate_Xg(0.0_rk, 0.0_rk, -45.0_rk)
   call ut%test(ti)%check( &
      name="rotate_Xg(Z-axis): geometry preserved (NURBS)", &
      res=nurbs%get_Xg(), &
      expected=Xg, &
      tol=TOL, &
      msg="NURBS Xg changed after Z-axis rotation", &
      group="transform"); ti=ti+1
   call ut%test(ti)%check( &
      name="rotate_Xg(Z-axis): geometry preserved (B-spline)", &
      res=bsp%get_Xg(), &
      expected=Xgb, &
      tol=TOL, &
      msg="B-spline Xg changed after Z-axis rotation", &
      group="transform"); ti=ti+1

   ! 16) Translation tests (Xc)
   call nurbs%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
   call nurbs%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])
   call bsp%translate_Xc([5.0_rk, 5.0_rk, 5.0_rk])
   call bsp%translate_Xc([-5.0_rk, -5.0_rk, -5.0_rk])
   call ut%test(ti)%check( &
      name="translate_Xc(): geometry preserved (NURBS)", &
      res=nurbs%get_Xc(), &
      expected=Xc, &
      tol=TOL, &
      msg="NURBS Xc changed after translation", &
      group="transform"); ti=ti+1
   call ut%test(ti)%check( &
      name="translate_Xc(): geometry preserved (B-spline)", &
      res=bsp%get_Xc(), &
      expected=Xc, &
      tol=TOL, &
      msg="B-spline Xc changed after translation", &
      group="transform"); ti=ti+1

   ! 17) Translation tests (Xg)
   call nurbs%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
   call nurbs%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])
   call bsp%translate_Xg([5.0_rk, 5.0_rk, 5.0_rk])
   call bsp%translate_Xg([-5.0_rk, -5.0_rk, -5.0_rk])
   call ut%test(ti)%check( &
      name="translate_Xg(): geometry preserved (NURBS)", &
      res=nurbs%get_Xg(), &
      expected=Xg, &
      tol=TOL, &
      msg="NURBS Xg changed after translation", &
      group="transform"); ti=ti+1
   call ut%test(ti)%check( &
      name="translate_Xg(): geometry preserved (B-spline)", &
      res=bsp%get_Xg(), &
      expected=Xgb, &
      tol=TOL, &
      msg="B-spline Xg changed after translation", &
      group="transform"); ti=ti+1

   ! 18) Knot insertion
   call nurbs%insert_knots([0.25_rk, 0.75_rk], [2,1])
   call bsp%insert_knots([0.25_rk, 0.75_rk], [2,1])
   call nurbs%create()
   call bsp%create()
   call ut%test(ti)%check( &
      name="insert_knots(): geometry preserved (NURBS)", &
      res=nurbs%get_Xg(), &
      expected=Xg, &
      tol=TOL, &
      msg="NURBS geometry changed after knot insertion", &
      group="knot-ops"); ti=ti+1
   call ut%test(ti)%check( &
      name="insert_knots(): geometry preserved (B-spline)", &
      res=bsp%get_Xg(), &
      expected=Xgb, &
      tol=TOL, &
      msg="B-spline geometry changed after knot insertion", &
      group="knot-ops"); ti=ti+1

   ! 19) Degree elevation
   call nurbs%elevate_degree(2)
   call bsp%elevate_degree(2)
   call nurbs%create()
   call bsp%create()
   call ut%test(ti)%check( &
      name="elevate_degree(): geometry preserved (NURBS)", &
      res=nurbs%get_Xg(), &
      expected=Xg, &
      tol=TOL, &
      msg="NURBS geometry changed after degree elevation", &
      group="knot-ops"); ti=ti+1
   call ut%test(ti)%check( &
      name="elevate_degree(): geometry preserved (B-spline)", &
      res=bsp%get_Xg(), &
      expected=Xgb, &
      tol=TOL, &
      msg="B-spline geometry changed after degree elevation", &
      group="knot-ops"); ti=ti+1

   ! 20) Knot removal
   call nurbs%remove_knots([0.25_rk, 0.75_rk], [2,1])
   call bsp%remove_knots([0.25_rk, 0.75_rk], [2,1])
   call nurbs%create()
   call bsp%create()
   call ut%test(ti)%check( &
      name="remove_knots(): geometry preserved (NURBS)", &
      res=nurbs%get_Xg(), &
      expected=Xg, &
      tol=TOL, &
      msg="NURBS geometry changed after knot removal", &
      group="knot-ops"); ti=ti+1
   call ut%test(ti)%check( &
      name="remove_knots(): geometry preserved (B-spline)", &
      res=bsp%get_Xg(), &
      expected=Xgb, &
      tol=TOL, &
      msg="B-spline geometry changed after knot removal", &
      group="knot-ops"); ti=ti+1

   ! 21) Shape functions (circle, half-circle)
   call nurbs%set_circle([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk)
   call nurbs%create(res=23)
   call nurbs%nearest_point([1.0_rk, 0.0_rk, 0.0_rk], nearest_Xg, nearest_Xt, id)
   call ut%test(ti)%check( &
      name="set_circle(): nearest point param (NURBS)", &
      res=nearest_Xt, &
      expected=0.0_rk, &
      tol=TOL, &
      msg="NURBS circle nearest point param incorrect", &
      group="shapes"); ti=ti+1

   call nurbs%set_half_circle([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk)
   call nurbs%create(res=23)
   call nurbs%nearest_point([0.0_rk, 1.0_rk, 0.0_rk], nearest_Xg, nearest_Xt, id)
   call ut%test(ti)%check( &
      name="set_half_circle(): nearest point param (NURBS)", &
      res=nearest_Xt, &
      expected=0.5_rk, &  ! Matches modern program's expectation
      tol=TOL, &
      msg="NURBS half-circle nearest point param incorrect", &
      group="shapes"); ti=ti+1

   ! 22) Least-squares B-spline fitting
   block
      type(nurbs_curve) :: bsp_fit
      integer :: n
      real(rk), allocatable :: Xt_fit(:), Xdata(:,:), Xg_eval(:,:)
      real(rk) :: err1, err2, err3, rms

      n = 42
      allocate(Xt_fit(n), Xdata(n, 3))
      do i = 1, n
         Xt_fit(i) = real(i-1, rk) / real(n-1, rk)
         Xdata(i,1) = Xt_fit(i)
         Xdata(i,2) = 0.3_rk * sin(4.0_rk * PI * Xt_fit(i))
         Xdata(i,3) = 0.3_rk * cos(4.0_rk * PI * Xt_fit(i))
      end do

      call bsp_fit%set( &
         degree=5, &
         Xth_dir=[(real(i-1, rk)/10.0_rk, i=1,11)], &
         continuity=[-1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1])
      call bsp_fit%lsq_fit_bspline(Xt_fit, Xdata, n)
      call bsp_fit%create(res=n)
      Xg_eval = bsp_fit%get_Xg()

      err1 = norm2(Xg_eval(:,1) - Xdata(:,1)) / norm2(Xdata(:,1))
      err2 = norm2(Xg_eval(:,2) - Xdata(:,2)) / norm2(Xdata(:,2))
      err3 = norm2(Xg_eval(:,3) - Xdata(:,3)) / norm2(Xdata(:,3))
      rms = sqrt((err1**2 + err2**2 + err3**2) / 3.0_rk)

      call ut%test(ti)%check( &
         name="lsq_fit_bspline(): RMS error", &
         res=rms, &
         expected=0.0_rk, &
         tol=1.0e-6_rk, &  ! Matches old program's tolerance
         msg="B-spline least-squares fit RMS error too high", &
         group="basis"); ti=ti+1

      call bsp_fit%finalize()
      deallocate(Xt_fit, Xdata, Xg_eval)
   end block

   ! Finalize
   call nurbs%finalize()
   call bsp%finalize()
   deallocate(Xc, Wc, Xg, Xgb, Xt)
   if (allocated(Tgc)) deallocate(Tgc)
   if (allocated(Tgcb)) deallocate(Tgcb)
   if (allocated(dTgc)) deallocate(dTgc)
   if (allocated(dTgcb)) deallocate(dTgcb)
   if (allocated(d2Tgc)) deallocate(d2Tgc)
   if (allocated(d2Tgcb)) deallocate(d2Tgcb)
   if (allocated(Tgc1)) deallocate(Tgc1)
   if (allocated(Tgc1b)) deallocate(Tgc1b)
   if (allocated(dTgc1)) deallocate(dTgc1)
   if (allocated(dTgc1b)) deallocate(dTgc1b)
   if (allocated(d2Tgc1)) deallocate(d2Tgc1)
   if (allocated(d2Tgc1b)) deallocate(d2Tgc1b)

   call ut%summary(verbose=3, required_score=100.0)

end program test_nurbs_curve
