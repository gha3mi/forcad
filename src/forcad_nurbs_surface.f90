!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> This module defines the 'nurbs_surface' type for representing a Non-Uniform Rational B-Spline (NURBS) surface.
module forcad_nurbs_surface

    use forcad_utils, only: rk, basis_bspline, elemConn_C0, kron, ndgrid, compute_multiplicity, compute_knot_vector, &
        basis_bspline_der, insert_knot_A_5_1, findspan, elevate_degree_A_5_9, remove_knots_A_5_8, tetragon_Xc, &
        elemConn_Cn, unique, rotation

    implicit none

    private
    public nurbs_surface

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    type nurbs_surface
        real(rk), allocatable, private :: Xc(:,:)  !! Control points (2D array: [nc(1)*nc(2), dim])
        real(rk), allocatable, private :: Xg(:,:)  !! Geometry points (2D array: [ng(1)*ng(2), dim])
        real(rk), allocatable, private :: Wc(:)    !! Weights for control points (1D array: [nc(1)*nc(2)])
        real(rk), allocatable, private :: Xt1(:)   !! Evaluation parameter values in the first direction (1D array: [ng(1)])
        real(rk), allocatable, private :: Xt2(:)   !! Evaluation parameter values in the second direction (1D array: [ng(2)])
        real(rk), allocatable, private :: Xt(:,:)  !! Evaluation parameter values (2D array: [ng(1)*ng(2), 2])
        real(rk), allocatable, private :: knot1(:) !! Knot vector in the first direction (1D array)
        real(rk), allocatable, private :: knot2(:) !! Knot vector in the second direction (1D array)
        integer, private :: degree(2)              !! Degree (order) of the surface
        integer, private :: nc(2)                  !! Number of control points in each direction
        integer, private :: ng(2)                  !! Number of geometry points in each direction
        integer, allocatable, private :: elemConn_Xc_vis(:,:) !! Connectivity for visualization of control points
        integer, allocatable, private :: elemConn_Xg_vis(:,:) !! Connectivity for visualization of geometry points
        integer, allocatable, private :: elemConn(:,:)        !! IGA element connectivity
    contains
        procedure :: set1                   !!> Set knot vectors, control points and weights for the NURBS surface object
        procedure :: set2                   !!> Set NURBS surface using nodes of parameter space, degree, continuity, control points and weights
        procedure :: set3                   !!> Set Bezier or Rational Bezier surface using control points and weights
        generic :: set => set1, set2, set3  !!> Set NURBS surface
        procedure :: create                 !!> Generate geometry points
        procedure :: cmp_Xg                !!> Compute geometry points
        procedure, private :: get_Xc_all   !!> Get all control points
        procedure, private :: get_Xci      !!> Get i-th control point
        procedure, private :: get_Xcid     !!> Get i-th control point in a specific direction
        generic :: get_Xc => get_Xc_all, get_Xci, get_Xcid !!> Get control points
        procedure, private :: get_Xg_all   !!> Get all geometry points
        procedure, private :: get_Xgi      !!> Get i-th geometry point
        procedure, private :: get_Xgid     !!> Get i-th geometry point in a specific direction
        generic :: get_Xg => get_Xg_all, get_Xgi, get_Xgid !!> Get geometry points
        procedure, private :: get_Wc_all   !!> Get all weights
        procedure, private :: get_Wci      !!> Get i-th weight
        generic :: get_Wc => get_Wc_all, get_Wci !!> Get weights
        procedure :: get_Xt                 !!> Get parameter values
        procedure, private :: get_knot_all  !!> Get all knot vectors
        procedure, private :: get_knoti     !!> Get i-th knot value
        generic :: get_knot => get_knoti, get_knot_all !!> Get knot vector
        procedure :: get_ng                 !!> Get number of geometry points
        procedure :: cmp_degree             !!> Compute degree of the NURBS surface
        procedure, private :: get_degree_all!!> Get degree of the NURBS surface in both directions
        procedure, private :: get_degree_dir!!> Get degree of the NURBS surface in a specific direction
        generic :: get_degree => get_degree_all, get_degree_dir !!> Get degree of the NURBS surface
        procedure :: finalize               !!> Finalize the NURBS surface object
        procedure :: cmp_elem_Xc_vis        !!> Generate connectivity for control points
        procedure :: cmp_elem_Xg_vis        !!> Generate connectivity for geometry points
        procedure :: cmp_elem               !!> Generate IGA element connectivity
        procedure :: get_elem_Xc_vis        !!> Get connectivity for control points
        procedure :: get_elem_Xg_vis        !!> Get connectivity for geometry points
        procedure :: get_elem               !!> Get IGA element connectivity
        procedure :: set_elem_Xc_vis        !!> Set connectivity for control points
        procedure :: set_elem_Xg_vis        !!> Set connectivity for geometry points
        procedure :: set_elem               !!> Set IGA element connectivity
        procedure :: export_Xc              !!> Export control points to VTK file
        procedure :: export_Xg              !!> Export geometry points to VTK file
        procedure :: modify_Xc              !!> Modify control points
        procedure :: modify_Wc              !!> Modify weights
        procedure :: get_multiplicity       !!> Compute and return the multiplicity of the knot vector
        procedure :: get_continuity         !!> Compute and return the continuity of the NURBS surface
        procedure :: cmp_nc                 !!> Compute number of required control points
        procedure :: get_nc                 !!> Get number of control points
        procedure, private :: basis_vector  !!> Compute the basis functions of the NURBS surface
        procedure, private :: basis_scalar  !!> Compute the basis functions of the NURBS surface
        generic :: basis => basis_vector, basis_scalar    !!> Compute the basis functions of the NURBS surface
        procedure, private :: derivative_vector      !!> Compute the derivative of the NURBS surface
        procedure, private :: derivative_scalar      !!> Compute the derivative of the NURBS surface
        generic :: derivative => derivative_vector, derivative_scalar   !!> Compute the derivative of the NURBS surface
        procedure, private :: derivative2_vector     !!> Compute the second derivative of the NURBS surface
        procedure, private :: derivative2_scalar     !!> Compute the second derivative of the NURBS surface
        generic :: derivative2 => derivative2_vector, derivative2_scalar !!> Compute the second derivative of the NURBS surface
        procedure :: insert_knots           !!> Insert knots into the knot vector
        procedure :: elevate_degree         !!> Elevate degree
        procedure :: is_rational            !!> Check if the NURBS surface is rational
        procedure :: remove_knots           !!> Remove knots from the knot vector
        procedure :: rotate_Xc              !!> Rotate control points
        procedure :: rotate_Xg              !!> Rotate geometry points
        procedure :: translate_Xc           !!> Translate control points
        procedure :: translate_Xg           !!> Translate geometry points
        procedure :: show                   !!> Show the NURBS object using PyVista
        procedure :: nearest_point          !!> Find the nearest point on the NURBS surface (Approximation)
        procedure :: nearest_point2         !!> Find the nearest point on the NURBS surface (Minimization - Newton's method)

        ! Shapes
        procedure :: set_tetragon           !!> Set a tetragon
        procedure :: set_ring               !!> Set a ring
        procedure :: set_half_ring          !!> Set a half ring
        procedure :: set_C                  !!> Set a C-shape
    end type
    !===============================================================================

    interface compute_Xg
        pure function compute_Xg_nurbs_2d(f_Xt, f_knot1, f_knot2, f_degree, f_nc, f_ng, f_Xc, f_Wc) result(f_Xg)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:,:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: f_nc(2)
            integer, intent(in) :: f_ng(2)
            real(rk), intent(in), contiguous :: f_Xc(:,:)
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable :: f_Xg(:,:)
        end function

        pure function compute_Xg_bspline_2d(f_Xt, f_knot1, f_knot2, f_degree, f_nc, f_ng, f_Xc) result(f_Xg)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:,:)
            real(rk), intent(in), contiguous :: f_knot1(:)
            real(rk), intent(in), contiguous :: f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: f_nc(2)
            integer, intent(in) :: f_ng(2)
            real(rk), intent(in), contiguous :: f_Xc(:,:)
            real(rk), allocatable :: f_Xg(:,:)
        end function

        pure function compute_Xg_nurbs_2d_1point(f_Xt, f_knot1, f_knot2, f_degree, f_nc, f_Xc, f_Wc) result(f_Xg)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: f_nc(2)
            real(rk), intent(in), contiguous :: f_Xc(:,:)
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable :: f_Xg(:)
        end function

        pure function compute_Xg_bspline_2d_1point(f_Xt, f_knot1, f_knot2, f_degree, f_nc, f_Xc) result(f_Xg)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot1(:)
            real(rk), intent(in), contiguous :: f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: f_nc(2)
            real(rk), intent(in), contiguous :: f_Xc(:,:)
            real(rk), allocatable :: f_Xg(:)
        end function
    end interface

    interface compute_dTgc
        pure subroutine compute_dTgc_nurbs_2d_vector(f_Xt, f_knot1, f_knot2, f_degree, f_nc, f_ng, f_Wc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:,:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: f_nc(2)
            integer, intent(in) :: f_ng(2)
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable, intent(out) :: f_dTgc(:,:,:)
            real(rk), allocatable, intent(out) :: f_Tgc(:,:)
        end subroutine

        pure subroutine compute_dTgc_bspline_2d_vector(f_Xt, f_knot1, f_knot2, f_degree, nc, f_ng, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:,:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: nc(2)
            integer, intent(in) :: f_ng(2)
            real(rk), allocatable, intent(out) :: f_dTgc(:,:,:)
            real(rk), allocatable, intent(out) :: f_Tgc(:,:)
        end subroutine

        pure subroutine compute_dTgc_nurbs_2d_scalar(f_Xt, f_knot1, f_knot2, f_degree, f_nc, f_Wc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: f_nc(2)
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable, intent(out) :: f_dTgc(:,:)
            real(rk), allocatable, intent(out) :: f_Tgc(:)
        end subroutine

        pure subroutine compute_dTgc_bspline_2d_scalar(f_Xt, f_knot1, f_knot2, f_degree, nc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: nc(2)
            real(rk), allocatable, intent(out) :: f_dTgc(:,:)
            real(rk), allocatable, intent(out) :: f_Tgc(:)
        end subroutine
    end interface

    interface compute_d2Tgc
        pure subroutine compute_d2Tgc_nurbs_2d_vector(f_Xt, f_knot1, f_knot2, f_degree, f_nc, f_ng, f_Wc, f_d2Tgc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:,:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: f_nc(2)
            integer, intent(in) :: f_ng(2)
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable, intent(out) :: f_d2Tgc(:,:,:)
            real(rk), allocatable, intent(out) :: f_dTgc(:,:,:)
            real(rk), allocatable, intent(out) :: f_Tgc(:,:)
        end subroutine

        pure subroutine compute_d2Tgc_bspline_2d_vector(f_Xt, f_knot1, f_knot2, f_degree, nc, f_ng, f_d2Tgc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:,:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: nc(2)
            integer, intent(in) :: f_ng(2)
            real(rk), allocatable, intent(out) :: f_d2Tgc(:,:,:)
            real(rk), allocatable, intent(out) :: f_dTgc(:,:,:)
            real(rk), allocatable, intent(out) :: f_Tgc(:,:)
        end subroutine

        pure subroutine compute_d2Tgc_nurbs_2d_scalar(f_Xt, f_knot1, f_knot2, f_degree, f_nc, f_Wc, f_d2Tgc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: f_nc(2)
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable, intent(out) :: f_d2Tgc(:,:)
            real(rk), allocatable, intent(out) :: f_dTgc(:,:)
            real(rk), allocatable, intent(out) :: f_Tgc(:)
        end subroutine

        pure subroutine compute_d2Tgc_bspline_2d_scalar(f_Xt, f_knot1, f_knot2, f_degree, nc, f_d2Tgc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: nc(2)
            real(rk), allocatable, intent(out) :: f_d2Tgc(:,:)
            real(rk), allocatable, intent(out) :: f_dTgc(:,:)
            real(rk), allocatable, intent(out) :: f_Tgc(:)
        end subroutine
    end interface

    interface compute_Tgc
        pure function compute_Tgc_nurbs_2d_vector(f_Xt, f_knot1, f_knot2, f_degree, f_nc, f_ng, f_Wc) result(f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:,:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: f_nc(2)
            integer, intent(in) :: f_ng(2)
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable :: f_Tgc(:,:)
        end function

        pure function compute_Tgc_bspline_2d_vector(f_Xt, f_knot1, f_knot2, f_degree, f_nc, f_ng) result(f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:,:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: f_nc(2)
            integer, intent(in) :: f_ng(2)
            real(rk), allocatable :: f_Tgc(:,:)
        end function

        pure function compute_Tgc_nurbs_2d_scalar(f_Xt, f_knot1, f_knot2, f_degree, f_nc, f_Wc) result(f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: f_nc(2)
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable :: f_Tgc(:)
        end function

        pure function compute_Tgc_bspline_2d_scalar(f_Xt, f_knot1, f_knot2, f_degree, f_nc) result(f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot1(:), f_knot2(:)
            integer, intent(in) :: f_degree(2)
            integer, intent(in) :: f_nc(2)
            real(rk), allocatable :: f_Tgc(:)
        end function
    end interface

    interface
        pure function nearest_point_help_2d(f_ng, f_Xg, f_point_Xg) result(f_distances)
            import :: rk
            integer, intent(in) :: f_ng(2)
            real(rk), intent(in), contiguous :: f_Xg(:,:)
            real(rk), intent(in), contiguous :: f_point_Xg(:)
            real(rk), allocatable :: f_distances(:)
        end function
    end interface

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set knot vectors, control points and weights for the NURBS surface object.
    pure subroutine set1(this, knot1, knot2, Xc, Wc)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: knot1(:)
        real(rk), intent(in), contiguous :: knot2(:)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)

        if (allocated(this%knot1)) deallocate(this%knot1)
        if (allocated(this%knot2)) deallocate(this%knot2)
        if (allocated(this%Xc)) deallocate(this%Xc)

        this%knot1 = knot1
        this%knot2 = knot2
        call this%cmp_degree()
        call this%cmp_nc()
        this%Xc = Xc
        if (present(Wc)) then
            if (size(Wc) /= this%nc(1)*this%nc(2)) then
                error stop 'Number of weights does not match the number of control points.'
            else
                if (allocated(this%Wc)) deallocate(this%Wc)
                this%Wc = Wc
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set NURBS surface using nodes of parameter space, degree, continuity, control points and weights
    pure subroutine set2(this, Xth_dir1, Xth_dir2, degree, continuity1, continuity2, Xc, Wc)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xth_dir1(:), Xth_dir2(:)
        integer, intent(in), contiguous :: degree(:)
        integer, intent(in), contiguous :: continuity1(:), continuity2(:)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)

        this%knot1 = compute_knot_vector(Xth_dir1, degree(1), continuity1)
        this%knot2 = compute_knot_vector(Xth_dir2, degree(2), continuity2)
        this%degree(1) = degree(1)
        this%degree(2) = degree(2)
        call this%cmp_nc()
        this%Xc = Xc
        if (present(Wc)) this%Wc = Wc
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set Bezier or Rational Bezier surface using control points and weights.
    pure subroutine set3(this, nc, Xc, Wc)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), contiguous :: nc(:)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)

        if (allocated(this%Xc)) deallocate(this%Xc)

        this%Xc = Xc
        this%nc = nc

        if (allocated(this%knot1)) deallocate(this%knot1)
        allocate(this%knot1(2*this%nc(1)))
        this%knot1(1:this%nc(1)) = 0.0_rk
        this%knot1(this%nc(1)+1:2*this%nc(1)) = 1.0_rk

        if (allocated(this%knot2)) deallocate(this%knot2)
        allocate(this%knot2(2*this%nc(2)))
        this%knot2(1:this%nc(2)) = 0.0_rk
        this%knot2(this%nc(2)+1:2*this%nc(2)) = 1.0_rk

        call this%cmp_degree()
        if (present(Wc)) then
            if (size(Wc) /= this%nc(1)*this%nc(2)) then
                error stop 'Number of weights does not match the number of control points.'
            else
                if (allocated(this%Wc)) deallocate(this%Wc)
                this%Wc = Wc
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine create(this, res1, res2, Xt1, Xt2, Xt)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: res1, res2
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:)
        real(rk), contiguous, intent(in), optional :: Xt(:,:)
        integer :: i

        ! check
        if (.not.allocated(this%Xc)) then
            error stop 'Control points are not set.'
        end if

        if (.not.allocated(this%knot1) .or. .not.allocated(this%knot2)) then
            error stop 'Knot vector(s) is/are not set.'
        end if


        ! Set parameter values
        if (present(Xt1)) then
            if (allocated(this%Xt1)) deallocate(this%Xt1)
            this%Xt1 = Xt1
        elseif (present(res1)) then
            if (allocated(this%Xt1)) deallocate(this%Xt1)
            allocate(this%Xt1(res1))
            this%Xt1 = [(real(i-1, rk) / real(res1-1, rk), i=1, res1)]
            ! else
            ! this%Xt1 = this%Xt1
        end if

        ! Set parameter values
        if (present(Xt2)) then
            if (allocated(this%Xt2)) deallocate(this%Xt2)
            this%Xt2 = Xt2
        elseif (present(res2)) then
            if (allocated(this%Xt2)) deallocate(this%Xt2)
            allocate(this%Xt2(res2))
            this%Xt2 = [(real(i-1, rk) / real(res2-1, rk), i=1, res2)]
            ! else
            ! this%Xt2 = this%Xt2
        end if

        if (present(Xt)) then
            this%Xt = Xt
        else

            ! Set number of geometry points
            this%ng(1) = size(this%Xt1,1)
            this%ng(2) = size(this%Xt2,1)

            call ndgrid(this%Xt1, this%Xt2, this%Xt)
        end if

        if (allocated(this%Xg)) deallocate(this%Xg)

        if (this%is_rational()) then ! NURBS
            this%Xg = compute_Xg(&
                this%Xt, this%knot1, this%knot2, this%degree, this%nc, this%ng, this%Xc, this%Wc)
        else ! B-Spline
            this%Xg = compute_Xg(&
                this%Xt, this%knot1, this%knot2, this%degree, this%nc, this%ng, this%Xc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_Xg(this, Xt) result(Xg)
        class(nurbs_surface), intent(in) :: this
        real(rk), contiguous, intent(in) :: Xt(:)
        real(rk), allocatable :: Xg(:)

        ! check
        if (.not.allocated(this%Xc)) then
            error stop 'Control points are not set.'
        end if

        if (.not.allocated(this%knot1) .or. .not.allocated(this%knot2)) then
            error stop 'Knot vector(s) is/are not set.'
        end if

        if (this%is_rational()) then ! NURBS
            Xg = compute_Xg(Xt, this%knot1, this%knot2, this%degree, this%nc, this%Xc, this%Wc)
        else ! B-Spline
            Xg = compute_Xg(Xt, this%knot1, this%knot2, this%degree, this%nc, this%Xc)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xc_all(this) result(Xc)
        class(nurbs_surface), intent(in) :: this
        real(rk), allocatable :: Xc(:,:)

        if (allocated(this%Xc)) then
            Xc = this%Xc
        else
            error stop 'Control points are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xci(this, n) result(Xc)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: n
        real(rk), allocatable :: Xc(:)

        if (allocated(this%Xc)) then
            if (n<lbound(this%Xc,1) .or. n>ubound(this%Xc,1)) then
                error stop 'Invalid index for control points.'
            end if
            Xc = this%Xc(n,:)
        else
            error stop 'Control points are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xcid(this, n, dir) result(Xc)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: n
        integer, intent(in) :: dir
        real(rk) :: Xc

        if (allocated(this%Xc)) then
            if (n<lbound(this%Xc,1) .or. n>ubound(this%Xc,1)) then
                error stop 'Invalid index for control points.'
            end if
            if (dir<lbound(this%Xc,2) .or. dir>ubound(this%Xc,2)) then
                error stop 'Invalid direction for control points.'
            end if
            Xc = this%Xc(n, dir)
        else
            error stop 'Control points are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xg_all(this) result(Xg)
        class(nurbs_surface), intent(in) :: this
        real(rk), allocatable :: Xg(:,:)

        if (allocated(this%Xg)) then
            Xg = this%Xg
        else
            error stop 'Geometry points are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xgi(this, n) result(Xg)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: n
        real(rk), allocatable :: Xg(:)

        if (allocated(this%Xg)) then
            if (n<lbound(this%Xg,1) .or. n>ubound(this%Xg,1)) then
                error stop 'Invalid index for geometry points.'
            end if
            Xg = this%Xg(n,:)
        else
            error stop 'Control points are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xgid(this, n, dir) result(Xg)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: n
        integer, intent(in) :: dir
        real(rk) :: Xg

        if (allocated(this%Xg)) then
            if (n<lbound(this%Xg,1) .or. n>ubound(this%Xg,1)) then
                error stop 'Invalid index for geometry points.'
            end if
            if (dir<lbound(this%Xg,2) .or. dir>ubound(this%Xg,2)) then
                error stop 'Invalid direction for geometry points.'
            end if
            Xg = this%Xg(n, dir)
        else
            error stop 'Control points are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Wc_all(this) result(Wc)
        class(nurbs_surface), intent(in) :: this
        real(rk), allocatable :: Wc(:)

        if (allocated(this%Wc)) then
            Wc = this%Wc
        else
            error stop 'The NURBS surface is not rational or weights are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Wci(this, n) result(Wc)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: n
        real(rk) :: Wc

        if (allocated(this%Wc)) then
            if (n<lbound(this%Wc,1) .or. n>ubound(this%Wc,1)) then
                error stop 'Invalid index for weights.'
            end if
            Wc = this%Wc(n)
        else
            error stop 'The NURBS surface is not rational or weights are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xt(this, dir) result(Xt)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), allocatable :: Xt(:)

        if (dir == 1) then
            if (allocated(this%Xt1)) then
                Xt = this%Xt1
            else
                error stop 'Parameter values are not set.'
            end if
        elseif (dir == 2) then
            if (allocated(this%Xt2)) then
                Xt = this%Xt2
            else
                error stop 'Parameter values are not set.'
            end if
        else
            error stop 'Invalid direction for parameter values.'
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_ng(this) result(ng)
        class(nurbs_surface), intent(in) :: this
        integer :: ng(2)

        ng = this%ng
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine cmp_degree(this,dir)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: dir
        integer, allocatable :: m1(:), m2(:)

        if (present(dir)) then
            if (dir == 1) then
                m1 = this%get_multiplicity(1)
                this%degree(1) = m1(1) - 1
            else if (dir == 2) then
                m2 = this%get_multiplicity(2)
                this%degree(2) = m2(1) - 1
            else
                error stop 'Invalid direction for degree.'
            end if
        else
            m1 = this%get_multiplicity(1)
            this%degree(1) = m1(1) - 1

            m2 = this%get_multiplicity(2)
            this%degree(2) = m2(1) - 1
        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_degree_all(this) result(degree)
        class(nurbs_surface), intent(in) :: this
        integer :: degree(2)

        degree = this%degree
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_degree_dir(this,dir) result(degree)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        integer :: degree

        if (dir == 1) then
            degree = this%degree(1)
        else if (dir == 2) then
            degree = this%degree(2)
        else
            error stop 'Invalid direction for degree.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_knot_all(this, dir) result(knot)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), allocatable :: knot(:)

        if (dir == 1) then
            if (allocated(this%knot1)) then
                knot = this%knot1
            else
                error stop 'Knot vector is not set.'
            end if
        elseif (dir == 2) then
            if (allocated(this%knot2)) then
                knot = this%knot2
            else
                error stop 'Knot vector is not set.'
            end if
        else
            error stop 'Invalid direction for knot vector.'
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_knoti(this, dir, i) result(knot)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        integer, intent(in) :: i
        real(rk) :: knot

        if (dir == 1) then
            if (allocated(this%knot1)) then
                if (i < 1 .or. i > size(this%knot1)) then
                    error stop 'Invalid index for knot vector.'
                else
                    knot = this%knot1(i)
                end if
            else
                error stop 'Knot vector is not set.'
            end if
        elseif (dir == 2) then
            if (allocated(this%knot2)) then
                if (i < 1 .or. i > size(this%knot2)) then
                    error stop 'Invalid index for knot vector.'
                else
                    knot = this%knot2(i)
                end if
            else
                error stop 'Knot vector is not set.'
            end if
        else
            error stop 'Invalid direction for knot vector.'
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine finalize(this)
        class(nurbs_surface), intent(inout) :: this
        if (allocated(this%Xc)) deallocate(this%Xc)
        if (allocated(this%Xg)) deallocate(this%Xg)
        if (allocated(this%Wc)) deallocate(this%Wc)
        if (allocated(this%Xt1)) deallocate(this%Xt1)
        if (allocated(this%Xt2)) deallocate(this%Xt2)
        if (allocated(this%Xt)) deallocate(this%Xt)
        if (allocated(this%knot1)) deallocate(this%knot1)
        if (allocated(this%knot2)) deallocate(this%knot2)
        if (allocated(this%elemConn_Xc_vis)) deallocate(this%elemConn_Xc_vis)
        if (allocated(this%elemConn_Xg_vis)) deallocate(this%elemConn_Xg_vis)
        if (allocated(this%elemConn)) deallocate(this%elemConn)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elem_Xc_vis(this, p) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), contiguous, optional :: p(:)

        if (present(p)) then
            elemConn = elemConn_C0(this%nc(1), this%nc(2), p(1), p(2))
        else
            elemConn = elemConn_C0(this%nc(1), this%nc(2), 1, 1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elem_Xg_vis(this, p) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), contiguous, optional :: p(:)

        if (present(p)) then
            elemConn = elemConn_C0(this%ng(1), this%ng(2), p(1), p(2))
        else
            elemConn = elemConn_C0(this%ng(1), this%ng(2), 1, 1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xc(this, filename)
        class(nurbs_surface), intent(in) :: this
        character(len=*), intent(in) :: filename
        integer :: i, nc, nunit
        integer, allocatable :: elemConn(:,:)

        ! check
        if (.not.allocated(this%Xc)) then
            error stop 'Control points are not set.'
        end if

        if (.not.allocated(this%elemConn_Xc_vis)) then
            elemConn = this%cmp_elem_Xc_vis()
        else
            elemConn = this%elemConn_Xc_vis
        end if

        nc = size(this%Xc, 1)

        open(newunit=nunit, file=filename, action='write')
        write(nunit,'(a)') '# vtk DataFile Version 2.0'
        write(nunit,'(a)') 'Generated by ForCAD'
        write(nunit,'(a)') 'ASCII'
        write(nunit,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(nunit,'(a," ",g0," ",a)') 'POINTS', nc, 'double'

        if (size(this%Xc,2) == 2) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xc(i,1), this%Xc(i,2), 0.0_rk , i = 1, nc)
        elseif (size(this%Xc,2) == 3) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xc(i,1), this%Xc(i,2), this%Xc(i,3) , i = 1, nc)
        else
            error stop 'Invalid dimension for control points.'
        end if

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*(size(elemConn,2)+1)
        write(nunit,'(g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0)')&
            (4, elemConn(i,1)-1,elemConn(i,2)-1,elemConn(i,4)-1,elemConn(i,3)-1, i = 1, size(elemConn,1))

        write(nunit,'(a," ",g0)') 'CELL_TYPES', size(elemConn,1)
        write(nunit,'(g0)') (9 , i = 1, size(elemConn,1))
        close(nunit)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xg(this, filename)
        class(nurbs_surface), intent(in) :: this
        character(len=*), intent(in) :: filename
        integer :: i, ng, nunit
        integer, allocatable :: elemConn(:,:)

        ! check
        if (.not.allocated(this%Xg)) then
            error stop 'Geometry points are not set.'
        end if

        if (.not.allocated(this%elemConn_Xg_vis)) then
            elemConn = this%cmp_elem_Xg_vis()
        else
            elemConn = this%elemConn_Xg_vis
        end if

        ng = size(this%Xg, 1)

        open(newunit=nunit, file=filename, action='write')
        write(nunit,'(a)') '# vtk DataFile Version 2.0'
        write(nunit,'(a)') 'Generated by ForCAD'
        write(nunit,'(a)') 'ASCII'
        write(nunit,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(nunit,'(a," ",g0," ",a)') 'POINTS', ng, 'double'

        if (size(this%Xg,2) == 2) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xg(i,1), this%Xg(i,2), 0.0_rk , i = 1, ng)
        elseif (size(this%Xg,2) == 3) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xg(i,1), this%Xg(i,2), this%Xg(i,3) , i = 1, ng)
        else
            error stop 'Invalid dimension for geometry points.'
        end if

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*(size(elemConn,2)+1)
        write(nunit,'(g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0)')&
            (4, elemConn(i,1)-1,elemConn(i,2)-1,elemConn(i,4)-1,elemConn(i,3)-1, i = 1, size(elemConn,1))

        write(nunit,'(a," ",g0)') 'CELL_TYPES', size(elemConn,1)
        write(nunit,'(g0)') (9 , i = 1, size(elemConn,1))
        close(nunit)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine modify_Xc(this,X,num,dir)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: X
        integer, intent(in) :: num
        integer, intent(in) :: dir

        if (allocated(this%Xc)) then
            this%Xc(num,dir) = X
            if (allocated(this%Wc)) then
                call this%set(knot1=this%get_knot(1), knot2=this%get_knot(2), Xc=this%get_Xc(), Wc=this%get_Wc())
            else
                call this%set(knot1=this%get_knot(1), knot2=this%get_knot(2), Xc=this%get_Xc())
            end if
        else
            error stop 'Control points are not set.'
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine modify_Wc(this,W,num)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: W
        integer, intent(in) :: num

        if (allocated(this%Wc)) then
            this%Wc(num) = W
            if (allocated(this%knot1) .and. allocated(this%knot2)) then
                call this%set(knot1=this%get_knot(1), knot2=this%get_knot(2), Xc=this%get_Xc(), Wc=this%get_Wc())
            else
                call this%set(nc=this%nc, Xc=this%get_Xc(), Wc=this%get_Wc())
            end if
        else
            error stop 'The NURBS surface is not rational.'
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_multiplicity(this, dir) result(m)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        integer, allocatable :: m(:)

        if (dir == 1) then

            ! check
            if (.not.allocated(this%knot1)) then
                error stop 'Knot vector is not set.'
            else
                m = compute_multiplicity(this%knot1)
            end if

        elseif (dir == 2) then

            ! check
            if (.not.allocated(this%knot2)) then
                error stop 'Knot vector is not set.'
            else
                m = compute_multiplicity(this%knot2)
            end if

        else
            error stop 'Invalid direction.'
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_continuity(this, dir) result(c)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        integer, allocatable :: c(:)

        if (dir == 1) then

            ! check
            if (.not.allocated(this%knot1)) then
                error stop 'Knot vector is not set.'
            else
                c = this%degree(1) - compute_multiplicity(this%knot1)
            end if

        elseif (dir == 2) then

            ! check
            if (.not.allocated(this%knot2)) then
                error stop 'Knot vector is not set.'
            else
                c = this%degree(2) - compute_multiplicity(this%knot2)
            end if

        else
            error stop 'Invalid direction.'
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine cmp_nc(this, dir)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: dir

        if (present(dir)) then

            if (dir == 1) then

                ! check
                if (.not.allocated(this%knot1)) then
                    error stop 'Knot vector is not set.'
                else
                    this%nc(1) = sum(compute_multiplicity(this%knot1)) - this%degree(1) - 1
                end if

            elseif (dir == 2) then

                ! check
                if (.not.allocated(this%knot2)) then
                    error stop 'Knot vector is not set.'
                else
                    this%nc(2) = sum(compute_multiplicity(this%knot2)) - this%degree(2) - 1
                end if

            else
                error stop 'Invalid direction.'
            end if

        else
            ! check
            if (.not.allocated(this%knot1)) then
                error stop 'Knot vector is not set.'
            else
                this%nc(1) = sum(compute_multiplicity(this%knot1)) - this%degree(1) - 1
            end if

            ! check
            if (.not.allocated(this%knot2)) then
                error stop 'Knot vector is not set.'
            else
                this%nc(2) = sum(compute_multiplicity(this%knot2)) - this%degree(2) - 1
            end if

        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_nc(this, dir) result(nc)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        integer :: nc

        if (dir == 1) then

            ! check
            if (.not.allocated(this%knot1)) then
                error stop 'Knot vector is not set.'
            else
                nc = sum(compute_multiplicity(this%knot1)) - this%degree(1) - 1
            end if

        elseif (dir == 2) then

            ! check
            if (.not.allocated(this%knot2)) then
                error stop 'Knot vector is not set.'
            else
                nc = sum(compute_multiplicity(this%knot2)) - this%degree(2) - 1
            end if

        else
            error stop 'Invalid direction.'
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine derivative_vector(this, res1, res2, Xt1, Xt2, dTgc, Tgc)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: res1, res2
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
        integer :: i
        real(rk), allocatable :: Xt(:,:)

        ! Set parameter values
        if (present(Xt1)) then
            if (allocated(this%Xt1)) deallocate(this%Xt1)
            this%Xt1 = Xt1
        elseif (present(res1)) then
            if (allocated(this%Xt1)) deallocate(this%Xt1)
            allocate(this%Xt1(res1))
            this%Xt1 = [(real(i-1, rk) / real(res1-1, rk), i=1, res1)]
            ! else
            ! this%Xt1 = this%Xt1
        end if

        ! Set parameter values
        if (present(Xt2)) then
            if (allocated(this%Xt2)) deallocate(this%Xt2)
            this%Xt2 = Xt2
        elseif (present(res2)) then
            if (allocated(this%Xt2)) deallocate(this%Xt2)
            allocate(this%Xt2(res2))
            this%Xt2 = [(real(i-1, rk) / real(res2-1, rk), i=1, res2)]
            ! else
            ! this%Xt2 = this%Xt2
        end if

        ! Set number of geometry points
        this%ng(1) = size(this%Xt1,1)
        this%ng(2) = size(this%Xt2,1)

        call ndgrid(this%Xt1, this%Xt2, Xt)

        if (this%is_rational()) then ! NURBS
            call compute_dTgc(Xt, this%knot1, this%knot2, this%degree, this%nc, this%ng, this%Wc, dTgc, Tgc)
        else ! B-Spline
            call compute_dTgc(Xt, this%knot1, this%knot2, this%degree, this%nc, this%ng, dTgc, Tgc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine derivative_scalar(this, Xt, dTgc, Tgc)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out), optional :: Tgc(:)

        if (this%is_rational()) then ! NURBS
            call compute_dTgc(Xt, this%knot1, this%knot2, this%degree, this%nc, this%Wc, dTgc, Tgc)
        else ! B-Spline
            call compute_dTgc(Xt, this%knot1, this%knot2, this%degree, this%nc, dTgc, Tgc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine derivative2_vector(this, res1, res2, Xt1, Xt2, d2Tgc, dTgc, Tgc)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: res1, res2
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:,:)
        real(rk), allocatable, intent(out), optional :: dTgc(:,:,:)
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
        integer :: i
        real(rk), allocatable :: Xt(:,:)

        ! Set parameter values
        if (present(Xt1)) then
            if (allocated(this%Xt1)) deallocate(this%Xt1)
            this%Xt1 = Xt1
        elseif (present(res1)) then
            if (allocated(this%Xt1)) deallocate(this%Xt1)
            allocate(this%Xt1(res1))
            this%Xt1 = [(real(i-1, rk) / real(res1-1, rk), i=1, res1)]
            ! else
            ! this%Xt1 = this%Xt1
        end if

        ! Set parameter values
        if (present(Xt2)) then
            if (allocated(this%Xt2)) deallocate(this%Xt2)
            this%Xt2 = Xt2
        elseif (present(res2)) then
            if (allocated(this%Xt2)) deallocate(this%Xt2)
            allocate(this%Xt2(res2))
            this%Xt2 = [(real(i-1, rk) / real(res2-1, rk), i=1, res2)]
            ! else
            ! this%Xt2 = this%Xt2
        end if

        ! Set number of geometry points
        this%ng(1) = size(this%Xt1,1)
        this%ng(2) = size(this%Xt2,1)

        call ndgrid(this%Xt1, this%Xt2, Xt)

        if (this%is_rational()) then ! NURBS
            call compute_d2Tgc(Xt, this%knot1, this%knot2, this%degree, this%nc, this%ng, this%Wc, d2Tgc, dTgc, Tgc)
        else ! B-Spline
            call compute_d2Tgc(Xt, this%knot1, this%knot2, this%degree, this%nc, this%ng, d2Tgc, dTgc, Tgc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine derivative2_scalar(this, Xt, d2Tgc, dTgc, Tgc)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
        real(rk), allocatable, intent(out), optional :: dTgc(:,:)
        real(rk), allocatable, intent(out), optional :: Tgc(:)

        if (this%is_rational()) then ! NURBS
            call compute_d2Tgc(Xt, this%knot1, this%knot2, this%degree, this%nc, this%Wc, d2Tgc, dTgc, Tgc)
        else ! B-Spline
            call compute_d2Tgc(Xt, this%knot1, this%knot2, this%degree, this%nc, d2Tgc, dTgc, Tgc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine basis_vector(this, res1, res2, Xt1, Xt2, Tgc)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: res1, res2
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        integer :: i
        real(rk), allocatable :: Xt(:,:)

        ! Set parameter values
        if (present(Xt1)) then
            if (allocated(this%Xt1)) deallocate(this%Xt1)
            this%Xt1 = Xt1
        elseif (present(res1)) then
            if (allocated(this%Xt1)) deallocate(this%Xt1)
            allocate(this%Xt1(res1))
            this%Xt1 = [(real(i-1, rk) / real(res1-1, rk), i=1, res1)]
            ! else
            ! this%Xt1 = this%Xt1
        end if

        ! Set parameter values
        if (present(Xt2)) then
            if (allocated(this%Xt2)) deallocate(this%Xt2)
            this%Xt2 = Xt2
        elseif (present(res2)) then
            if (allocated(this%Xt2)) deallocate(this%Xt2)
            allocate(this%Xt2(res2))
            this%Xt2 = [(real(i-1, rk) / real(res2-1, rk), i=1, res2)]
            ! else
            ! this%Xt2 = this%Xt2
        end if

        ! Set number of geometry points
        this%ng(1) = size(this%Xt1,1)
        this%ng(2) = size(this%Xt2,1)

        call ndgrid(this%Xt1, this%Xt2, Xt)

        if (this%is_rational()) then ! NURBS
            Tgc = compute_Tgc(Xt, this%knot1, this%knot2, this%degree, this%nc, this%ng, this%Wc)
        else ! B-Spline
            Tgc = compute_Tgc(Xt, this%knot1, this%knot2, this%degree, this%nc, this%ng)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine basis_scalar(this, Xt, Tgc)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: Xt(:)
        real(rk), allocatable, intent(out) :: Tgc(:)

        if (this%is_rational()) then ! NURBS
            Tgc = compute_Tgc(Xt, this%knot1, this%knot2, this%degree, this%nc, this%Wc)
        else ! B-Spline
            Tgc = compute_Tgc(Xt, this%knot1, this%knot2, this%degree, this%nc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine insert_knots(this, dir ,Xth,r)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in) :: dir
        real(rk), intent(in), contiguous :: Xth(:)
        integer, intent(in), contiguous :: r(:)
        integer :: k, i, s, dim, j, n_new
        real(rk), allocatable :: Xc(:,:), Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:)
        real(rk), allocatable:: Xc3(:,:,:)


        if (dir == 1) then ! direction 1

            if(this%is_rational()) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(1)-1,this%degree(1),Xth(i),this%knot1)
                    if (this%knot1(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot1, Xth(i))
                    else
                        s = 0
                    end if

                    dim = size(this%Xc,2)
                    allocate(Xcw(size(this%Xc,1),dim+1))
                    do j = 1, size(this%Xc,1)
                        Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                    end do
                    Xcw(:,dim+1) = this%Wc(:)

                    Xcw = reshape(Xcw,[this%nc(1),this%nc(2)*(dim+1)])

                    call insert_knot_A_5_1(&
                        this%degree(1),&
                        this%knot1,&
                        Xcw,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        n_new,&
                        knot_new,&
                        Xcw_new)

                    Xcw_new = reshape(Xcw_new,[this%nc(2)*(n_new+1),dim+1])

                    allocate(Xc_new(1:this%nc(2)*(n_new+1),1:dim))
                    allocate(Wc_new(1:this%nc(2)*(n_new+1)))
                    do j = 1, this%nc(2)*(n_new+1)
                        Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                    end do
                    Wc_new(:) = Xcw_new(:,dim+1)

                    call this%set(knot1=knot_new, knot2=this%get_knot(2), Xc=Xc_new, Wc=Wc_new)
                    deallocate(Xcw, Xcw_new, Xc_new, Wc_new)
                end do


            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(1)-1,this%degree(1),Xth(i),this%knot1)
                    if (this%knot1(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot1, Xth(i))
                    else
                        s = 0
                    end if

                    dim = size(this%Xc,2)

                    Xc = reshape(this%Xc,[this%nc(1),this%nc(2)*dim])

                    call insert_knot_A_5_1(&
                        this%degree(1),&
                        this%knot1,&
                        Xc,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        n_new,&
                        knot_new,&
                        Xc_new)

                    Xc_new = reshape(Xc_new,[(this%nc(2))*(n_new+1),dim])

                    call this%set(knot1=knot_new, knot2=this%get_knot(2), Xc=Xc_new)
                end do

            end if


        elseif (dir == 2) then! direction 2

            if(this%is_rational()) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(2)-1,this%degree(2),Xth(i),this%knot2)
                    if (this%knot2(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot2, Xth(i))
                    else
                        s = 0
                    end if

                    dim = size(this%Xc,2)
                    allocate(Xcw(size(this%Xc,1),dim+1))
                    do j = 1, size(this%Xc,1)
                        Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                    end do
                    Xcw(:,dim+1) = this%Wc(:)

                    Xc3 = reshape(Xcw, [this%nc(1),this%nc(2),dim+1])
                    Xc3 = reshape(Xc3, [this%nc(2),this%nc(1),dim+1], order=[2,1,3])
                    Xcw = reshape(Xc3,[this%nc(2),this%nc(1)*(dim+1)])

                    call insert_knot_A_5_1(&
                        this%degree(2),&
                        this%knot2,&
                        Xcw,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        n_new,&
                        knot_new,&
                        Xcw_new)

                    Xc3 = reshape(Xcw_new, [n_new+1,this%nc(1),dim+1])
                    Xc3 = reshape(Xc3, [this%nc(1),n_new+1,dim+1], order=[2,1,3])
                    Xcw_new = reshape(Xc3,[(this%nc(1))*(n_new+1),dim+1])

                    allocate(Xc_new(1:(n_new+1)*this%nc(1),1:dim))
                    allocate(Wc_new(1:(n_new+1)*this%nc(1)))
                    do j = 1, (n_new+1)*this%nc(1)
                        Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                    end do
                    Wc_new(:) = Xcw_new(:,dim+1)

                    call this%set(knot1=this%get_knot(1), knot2=knot_new, Xc=Xc_new, Wc=Wc_new)
                    deallocate(Xcw, Xcw_new, Xc_new, Wc_new)
                end do

            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(2)-1,this%degree(2),Xth(i),this%knot2)
                    if (this%knot2(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot2, Xth(i))
                    else
                        s = 0
                    end if

                    dim = size(this%Xc,2)

                    Xc3 = reshape(this%Xc, [this%nc(1),this%nc(2),dim])
                    Xc3 = reshape(Xc3, [this%nc(2),this%nc(1),dim], order=[2,1,3])
                    Xc = reshape(Xc3,[this%nc(2),this%nc(1)*dim])

                    call insert_knot_A_5_1(&
                        this%degree(2),&
                        this%knot2,&
                        Xc,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        n_new,&
                        knot_new,&
                        Xc_new)

                    Xc3 = reshape(Xc_new, [n_new+1,this%nc(1),dim])
                    Xc3 = reshape(Xc3, [this%nc(1),n_new+1,dim], order=[2,1,3])
                    Xc_new = reshape(Xc3,[(this%nc(1))*(n_new+1),dim])

                    call this%set(knot1=this%get_knot(1), knot2=knot_new, Xc=Xc_new)
                end do


            end if

        else
            error stop 'Invalid direction.'
        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine elevate_degree(this, dir, t)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in) :: dir
        integer, intent(in) :: t
        real(rk), allocatable :: Xc(:,:), Xcw(:,:), Xcw_new(:,:), knot_new(:), Xc_new(:,:), Wc_new(:)
        integer :: dim, j, nc_new
        real(rk), allocatable:: Xc3(:,:,:)


        if (dir == 1) then ! direction 1

            if(this%is_rational()) then ! NURBS

                dim = size(this%Xc,2)
                allocate(Xcw(size(this%Xc,1),dim+1))
                do j = 1, size(this%Xc,1)
                    Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                end do
                Xcw(:,dim+1) = this%Wc(:)

                Xcw = reshape(Xcw,[this%nc(1),this%nc(2)*(dim+1)],order=[1,2])

                call elevate_degree_A_5_9(t, this%knot1, this%degree(1), Xcw, nc_new, knot_new, Xcw_new)

                Xcw_new = reshape(Xcw_new,[this%nc(2)*nc_new,dim+1],order=[1,2])

                allocate(Xc_new(1:this%nc(2)*nc_new,1:dim))
                allocate(Wc_new(1:this%nc(2)*nc_new))
                do j = 1, this%nc(2)*nc_new
                    Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                end do

                Wc_new(:) = Xcw_new(:,dim+1)

                call this%set(knot1=knot_new, knot2=this%get_knot(2), Xc=Xc_new, Wc=Wc_new)
                deallocate(Xcw, Xcw_new, Xc_new, Wc_new)

            else ! B-Spline

                dim = size(this%Xc,2)
                Xc = reshape(this%Xc,[this%nc(1),this%nc(2)*(dim)],order=[1,2])

                call elevate_degree_A_5_9(t, this%knot1, this%degree(1), Xc, nc_new, knot_new, Xc_new)

                Xc_new = reshape(Xc_new,[this%nc(2)*nc_new,dim],order=[1,2])

                call this%set(knot1=knot_new, knot2=this%get_knot(2), Xc=Xc_new)
                deallocate(Xc, Xc_new)

            end if

        elseif (dir == 2) then ! direction 2

            if(this%is_rational()) then ! NURBS

                dim = size(this%Xc,2)
                allocate(Xcw(size(this%Xc,1),dim+1))
                do j = 1, size(this%Xc,1)
                    Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                end do

                Xcw(:,dim+1) = this%Wc(:)

                Xc3 = reshape(Xcw, [this%nc(1),this%nc(2),dim+1])
                Xc3 = reshape(Xc3, [this%nc(2),this%nc(1),dim+1], order=[2,1,3])
                Xcw = reshape(Xc3,[this%nc(2),this%nc(1)*(dim+1)])

                call elevate_degree_A_5_9(t, this%knot2, this%degree(2), Xcw, nc_new, knot_new, Xcw_new)

                Xc3 = reshape(Xcw_new, [nc_new,this%nc(1),dim+1])
                Xc3 = reshape(Xc3, [this%nc(1),nc_new,dim+1], order=[2,1,3])
                Xcw_new = reshape(Xc3,[(this%nc(1))*nc_new,dim+1])

                allocate(Xc_new(1:nc_new*this%nc(1),1:dim))
                allocate(Wc_new(1:nc_new*this%nc(1)))
                do j = 1, nc_new*this%nc(1)
                    Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                end do

                Wc_new(:) = Xcw_new(:,dim+1)

                call this%set(knot1=this%get_knot(1), knot2=knot_new, Xc=Xc_new, Wc=Wc_new)
                deallocate(Xcw, Xcw_new, Xc_new, Wc_new)

            else ! B-Spline

                dim = size(this%Xc,2)

                Xc3 = reshape(this%Xc, [this%nc(1),this%nc(2),dim])
                Xc3 = reshape(Xc3, [this%nc(2),this%nc(1),dim], order=[2,1,3])
                Xc = reshape(Xc3,[this%nc(2),this%nc(1)*dim])

                call elevate_degree_A_5_9(t, this%knot2, this%degree(2), Xc, nc_new, knot_new, Xc_new)

                Xc3 = reshape(Xc_new, [nc_new,this%nc(1),dim])
                Xc3 = reshape(Xc3, [this%nc(1),nc_new,dim], order=[2,1,3])
                Xc_new = reshape(Xc3,[(this%nc(1))*nc_new,dim])

                call this%set(knot1=this%get_knot(1), knot2=knot_new, Xc=Xc_new)

            end if

        else
            error stop 'Invalid direction.'
        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function is_rational(this) result(r)
        class(nurbs_surface), intent(in) :: this
        logical :: r

        r = .false.
        if(allocated(this%Wc)) then
            if (any(this%Wc /= this%Wc(1))) then
                r = .true.
            end if
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_elem_Xc_vis(this, elemConn)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (allocated(this%elemConn_Xc_vis)) deallocate(this%elemConn_Xc_vis)
        this%elemConn_Xc_vis = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_elem_Xg_vis(this, elemConn)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (allocated(this%elemConn_Xg_vis)) deallocate(this%elemConn_Xg_vis)
        this%elemConn_Xg_vis = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_elem(this, elemConn)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (allocated(this%elemConn)) deallocate(this%elemConn)
        this%elemConn = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem_Xc_vis(this) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        elemConn = this%elemConn_Xc_vis
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem_Xg_vis(this) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        elemConn = this%elemConn_Xg_vis
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem(this) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        elemConn = this%elemConn
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine remove_knots(this, dir ,Xth,r)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in) :: dir
        real(rk), intent(in), contiguous :: Xth(:)
        integer, intent(in), contiguous :: r(:)
        integer :: k, i, s, dim, j, nc_new, t
        real(rk), allocatable :: Xc(:,:), Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:)
        real(rk), allocatable:: Xc3(:,:,:)


        if (dir == 1) then ! direction 1

            if(this%is_rational()) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(1)-1,this%degree(1),Xth(i),this%knot1)
                    if (this%knot1(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot1, Xth(i))
                    else
                        s = 0
                    end if
                    k = k + 1

                    dim = size(this%Xc,2)
                    allocate(Xcw(size(this%Xc,1),dim+1))
                    do j = 1, size(this%Xc,1)
                        Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                    end do
                    Xcw(:,dim+1) = this%Wc(:)

                    Xcw = reshape(Xcw,[this%nc(1),this%nc(2)*(dim+1)],order=[1,2])

                    call remove_knots_A_5_8(&
                        this%degree(1),&
                        this%knot1,&
                        Xcw,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        t,&
                        knot_new,&
                        Xcw_new)

                    if (allocated(Xcw)) deallocate(Xcw)

                    if (t == 0) then
                        ! no change
                    else
                        nc_new = size(Xcw_new,1)
                        Xcw_new = reshape(Xcw_new,[this%nc(2)*(nc_new),dim+1],order=[1,2])

                        allocate(Xc_new(1:this%nc(2)*(nc_new),1:dim))
                        allocate(Wc_new(1:this%nc(2)*(nc_new)))
                        do j = 1, this%nc(2)*(nc_new)
                            Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                        end do

                        Wc_new(:) = Xcw_new(:,dim+1)

                        call this%set(knot1=knot_new, knot2=this%get_knot(2), Xc=Xc_new, Wc=Wc_new)
                        deallocate(Xcw_new, Xc_new, Wc_new)
                    end if
                end do


            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(1)-1,this%degree(1),Xth(i),this%knot1)
                    if (this%knot1(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot1, Xth(i))
                    else
                        s = 0
                    end if
                    k = k + 1

                    dim = size(this%Xc,2)

                    Xc = reshape(this%Xc,[this%nc(1),this%nc(2)*dim],order=[1,2])

                    call remove_knots_A_5_8(&
                        this%degree(1),&
                        this%knot1,&
                        Xc,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        t,&
                        knot_new,&
                        Xc_new)

                    if (allocated(Xc)) deallocate(Xc)

                    if (t == 0) then
                        ! no change
                    else
                        nc_new = size(Xc_new,1)
                        Xc_new = reshape(Xc_new,[(this%nc(2))*(nc_new),dim],order=[1,2])

                        call this%set(knot1=knot_new, knot2=this%get_knot(2), Xc=Xc_new)
                    end if
                end do

            end if


        elseif (dir == 2) then! direction 2

            if(this%is_rational()) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(2)-1,this%degree(2),Xth(i),this%knot2)
                    if (this%knot2(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot2, Xth(i))
                    else
                        s = 0
                    end if
                    k = k + 1

                    dim = size(this%Xc,2)
                    allocate(Xcw(size(this%Xc,1),dim+1))
                    do j = 1, size(this%Xc,1)
                        Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                    end do

                    Xcw(:,dim+1) = this%Wc(:)

                    Xc3 = reshape(Xcw, [this%nc(1),this%nc(2),dim+1])
                    Xc3 = reshape(Xc3, [this%nc(2),this%nc(1),dim+1], order=[2,1,3])
                    Xcw = reshape(Xc3, [this%nc(2),this%nc(1)*(dim+1)])

                    call remove_knots_A_5_8(&
                        this%degree(2),&
                        this%knot2,&
                        Xcw,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        t,&
                        knot_new,&
                        Xcw_new)

                    if (allocated(Xcw)) deallocate(Xcw)

                    if (t == 0) then
                        ! no change
                    else

                        nc_new = size(Xcw_new,1)
                        Xc3 = reshape(Xcw_new, [nc_new,this%nc(1),dim+1])
                        Xc3 = reshape(Xc3, [this%nc(1),nc_new,dim+1], order=[2,1,3])
                        Xcw_new = reshape(Xc3,[(this%nc(1))*(nc_new),dim+1])

                        allocate(Xc_new(1:(nc_new)*this%nc(1),1:dim))
                        allocate(Wc_new(1:(nc_new)*this%nc(1)))
                        do j = 1, (nc_new)*this%nc(1)
                            Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                        end do

                        Wc_new(:) = Xcw_new(:,dim+1)

                        call this%set(knot1=this%get_knot(1), knot2=knot_new, Xc=Xc_new, Wc=Wc_new)
                        deallocate(Xcw_new, Xc_new, Wc_new)
                    end if

                end do

            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(2)-1,this%degree(2),Xth(i),this%knot2)
                    if (this%knot2(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot2, Xth(i))
                    else
                        s = 0
                    end if
                    k = k + 1

                    dim = size(this%Xc,2)

                    Xc3 = reshape(this%Xc, [this%nc(1),this%nc(2),dim])
                    Xc3 = reshape(Xc3, [this%nc(2),this%nc(1),dim], order=[2,1,3])
                    Xc = reshape(Xc3,[this%nc(2),this%nc(1)*dim])

                    call remove_knots_A_5_8(&
                        this%degree(2),&
                        this%knot2,&
                        Xc,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        t,&
                        knot_new,&
                        Xc_new)

                    if (allocated(Xc)) deallocate(Xc)

                    if (t == 0) then
                        ! no change
                    else
                        nc_new = size(Xc_new,1)

                        Xc3 = reshape(Xc_new, [nc_new,this%nc(1),dim])
                        Xc3 = reshape(Xc3, [this%nc(1),nc_new,dim], order=[2,1,3])
                        Xc_new = reshape(Xc3,[(this%nc(1))*(nc_new),dim])

                        call this%set(knot1=this%get_knot(1), knot2=knot_new, Xc=Xc_new)
                    end if

                end do


            end if

        else
            error stop 'Invalid direction.'
        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_tetragon(this, L, nc, Wc)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: L(2)
        integer, intent(in) :: nc(2)
        real(rk), intent(in), contiguous, optional :: Wc(:)

        call this%set(nc = nc, Xc = tetragon_Xc(L, nc), Wc = Wc)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elem(this) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        call elemConn_Cn(this%nc(1), this%nc(2),&
            this%degree(1),this%degree(2),&
            unique(this%knot1),unique(this%knot2),&
            this%get_multiplicity(1),this%get_multiplicity(2),&
            elemConn)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine rotate_Xc(this, alpha, beta, theta)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: alpha, beta, theta
        integer :: i

        do i = 1, this%nc(1)*this%nc(2)
            this%Xc(i, :) = matmul(rotation(alpha,beta,theta), this%Xc(i, :))
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine rotate_Xg(this, alpha, beta, theta)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: alpha, beta, theta
        integer :: i

        do i = 1, this%ng(1)*this%ng(2)
            this%Xg(i, :) = matmul(rotation(alpha,beta,theta), this%Xg(i, :))
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine translate_Xc(this, vec)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: vec(:)
        integer :: i

        do i = 1, this%nc(1)*this%nc(2)
            this%Xc(i, :) = this%Xc(i, :) + vec
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine translate_Xg(this, vec)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: vec(:)
        integer :: i

        do i = 1, this%ng(1)*this%ng(2)
            this%Xg(i, :) = this%Xg(i, :) + vec
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine show(this, vtkfile_Xc, vtkfile_Xg)
        class(nurbs_surface), intent(inout) :: this
        character(len=*), intent(in) :: vtkfile_Xc, vtkfile_Xg
        character(len=3000) :: pyvista_script

        pyvista_script = &
            "import pyvista as pv"//achar(10)//&
            "pv.global_theme.color = 'white'"//achar(10)//&
            "Xc = pv.read('"//trim(vtkfile_Xc)//"')"//achar(10)//&
            "Xg = pv.read('"//trim(vtkfile_Xg)//"')"//achar(10)//&
            "p = pv.Plotter(lighting='light kit')"//achar(10)//&
            "actor_Xcp = p.add_mesh("//achar(10)//&
            "    Xc,"//achar(10)//&
            "    style='points',"//achar(10)//&
            "    point_size=10,"//achar(10)//&
            "    color='red',"//achar(10)//&
            "    render_points_as_spheres=True,"//achar(10)//&
            "    opacity=0.5,"//achar(10)//&
            ")"//achar(10)//&
            "actor_Xcw = p.add_mesh("//achar(10)//&
            "    Xc,"//achar(10)//&
            "    show_edges=True,"//achar(10)//&
            "    color='yellow',"//achar(10)//&
            "    line_width=3,"//achar(10)//&
            "    style='wireframe',"//achar(10)//&
            "    opacity=0.2"//achar(10)//&
            ")"//achar(10)//&
            "actor_Xg = p.add_mesh("//achar(10)//&
            "    Xg,"//achar(10)//&
            "    show_edges=False,"//achar(10)//&
            "    color='cyan',"//achar(10)//&
            "    line_width=7,"//achar(10)//&
            "    metallic=0.6,"//achar(10)//&
            "    pbr=True,"//achar(10)//&
            "    split_sharp_edges=True,"//achar(10)//&
            ")"//achar(10)//&
            "p.add_axes(interactive=False)"//achar(10)//&
            "def point_picker_callback(point):"//achar(10)//&
            "    mesh = Xc"//achar(10)//&
            "    point_id = mesh.find_closest_point(point)"//achar(10)//&
            "    point_coords = mesh.points[point_id]"//achar(10)//&
            "    label = f'ID: {point_id + 1}\n({point_coords[0]:.3f}, {point_coords[1]:.3f}, {point_coords[2]:.3f})'"//achar(10)//&
            "    p.add_point_labels("//achar(10)//&
            "        [point_coords],"//achar(10)//&
            "        [label],"//achar(10)//&
            "        font_size=14,"//achar(10)//&
            "        text_color='black',"//achar(10)//&
            "        show_points=False,"//achar(10)//&
            "        fill_shape=False,"//achar(10)//&
            "        shape=None,"//achar(10)//&
            "    )"//achar(10)//&
            "picker = p.enable_point_picking(callback=point_picker_callback, show_message=False)"//achar(10)//&
            "window_size = p.window_size"//achar(10)//&
            "y_pos = window_size[1]"//achar(10)//&
            "def Xcp_toggle_vis(flag):"//achar(10)//&
            "    actor_Xcp.SetVisibility(flag)"//achar(10)//&
            "def Xcw_toggle_vis(flag):"//achar(10)//&
            "    actor_Xcw.SetVisibility(flag)"//achar(10)//&
            "def Xg_toggle_vis(flag):"//achar(10)//&
            "    actor_Xg.SetVisibility(flag)"//achar(10)//&
            "p.add_checkbox_button_widget("//achar(10)//&
            "    Xcp_toggle_vis,"//achar(10)//&
            "    value=True,"//achar(10)//&
            "    color_on='red',"//achar(10)//&
            "    size=25,"//achar(10)//&
            "    position=(0, y_pos - 1 * 25),"//achar(10)//&
            ")"//achar(10)//&
            "p.add_checkbox_button_widget("//achar(10)//&
            "    Xcw_toggle_vis,"//achar(10)//&
            "    value=True,"//achar(10)//&
            "    color_on='yellow',"//achar(10)//&
            "    size=25,"//achar(10)//&
            "    position=(0, y_pos - 2 * 25),"//achar(10)//&
            ")"//achar(10)//&
            "p.add_checkbox_button_widget("//achar(10)//&
            "    Xg_toggle_vis,"//achar(10)//&
            "    value=True,"//achar(10)//&
            "    color_on='cyan',"//achar(10)//&
            "    size=25,"//achar(10)//&
            "    position=(0, y_pos - 3 * 25),"//achar(10)//&
            ")"//achar(10)//&
            "p.add_text("//achar(10)//&
            "    'Xc (Points)',"//achar(10)//&
            "    position=(25 + 3, y_pos - 1 * 25),"//achar(10)//&
            "    font_size=8,"//achar(10)//&
            "    color='black',"//achar(10)//&
            "    font='times',"//achar(10)//&
            ")"//achar(10)//&
            "p.add_text("//achar(10)//&
            "    'Xc (Control geometry)',"//achar(10)//&
            "    position=(25 + 3, y_pos - 2 * 25),"//achar(10)//&
            "    font_size=8,"//achar(10)//&
            "    color='black',"//achar(10)//&
            "    font='times',"//achar(10)//&
            ")"//achar(10)//&
            "p.add_text("//achar(10)//&
            "    'Xg (Geometry)',"//achar(10)//&
            "    position=(25 + 3, y_pos - 3 * 25),"//achar(10)//&
            "    font_size=8,"//achar(10)//&
            "    color='black',"//achar(10)//&
            "    font='times',"//achar(10)//&
            ")"//achar(10)//&
            "p.add_text('ForCAD', position=(0.0, 10.0), font_size=14, color='black', font='times')"//achar(10)//&
            "p.add_text("//achar(10)//&
            "    'https://github.com/gha3mi/forcad',"//achar(10)//&
            "    position=(0.0, 0.0),"//achar(10)//&
            "    font_size=7,"//achar(10)//&
            "    color='blue',"//achar(10)//&
            "    font='times',"//achar(10)//&
            ")"//achar(10)//&
            "p.show(title='ForCAD', interactive=True)"

        call execute_command_line('python -c "'//trim(adjustl(pyvista_script))//'"')
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_ring(this, center, radius1, radius2)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: center(:)
        real(rk), intent(in) :: radius1, radius2
        real(rk), allocatable :: Xc(:,:), Wc(:), knot1(:), knot2(:)
        integer :: i

        ! Define control points for ring
        allocate(Xc(14, 3))
        Xc(1,:) = [ 1.0_rk,  0.0_rk,              0.0_rk]
        Xc(2,:) = [ 1.0_rk,  sqrt(3.0_rk),        0.0_rk]
        Xc(3,:) = [-0.5_rk,  sqrt(3.0_rk)/2.0_rk, 0.0_rk]
        Xc(4,:) = [-2.0_rk,  0.0_rk,              0.0_rk]
        Xc(5,:) = [-0.5_rk, -sqrt(3.0_rk)/2.0_rk, 0.0_rk]
        Xc(6,:) = [ 1.0_rk, -sqrt(3.0_rk),        0.0_rk]
        Xc(7,:) = [ 1.0_rk,  0.0_rk,              0.0_rk]

        Xc(1:7,1:2) = Xc(1:7,1:2) * radius1

        Xc(8,:) = [ 1.0_rk,  0.0_rk,              0.0_rk]
        Xc(9,:) = [ 1.0_rk,  sqrt(3.0_rk),        0.0_rk]
        Xc(10,:)= [-0.5_rk,  sqrt(3.0_rk)/2.0_rk, 0.0_rk]
        Xc(11,:)= [-2.0_rk,  0.0_rk,              0.0_rk]
        Xc(12,:)= [-0.5_rk, -sqrt(3.0_rk)/2.0_rk, 0.0_rk]
        Xc(13,:)= [ 1.0_rk, -sqrt(3.0_rk),        0.0_rk]
        Xc(14,:)= [ 1.0_rk,  0.0_rk,              0.0_rk]

        Xc(8:14,1:2) = Xc(8:14,1:2) * radius2

        ! Translate the control points
        do i = 1, size(Xc, 1)
            Xc(i,:) = center + Xc(i,:)
        end do

        ! Define weights for the control points
        Wc = [1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk,&
            1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk]

        ! Define knot vector
        knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/3.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot2 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

        ! Set knot vector, control points, and weights
        call this%set(knot1, knot2, Xc, Wc)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_C(this, center, radius1, radius2)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: center(:)
        real(rk), intent(in) :: radius1, radius2
        real(rk), allocatable :: Xc(:,:), Wc(:), knot1(:), knot2(:)
        integer :: i

        ! Define control points for C-shape
        allocate(Xc(10, 3))
        Xc(1,:)= [ 1.0_rk,  0.0_rk,              0.0_rk]
        Xc(2,:)= [ 1.0_rk,  sqrt(3.0_rk),        0.0_rk]
        Xc(3,:)= [-0.5_rk,  sqrt(3.0_rk)/2.0_rk, 0.0_rk]
        Xc(4,:)= [-2.0_rk,  0.0_rk,              0.0_rk]
        Xc(5,:)= [-0.5_rk, -sqrt(3.0_rk)/2.0_rk, 0.0_rk]

        Xc(1:5,1:2) = Xc(1:5,1:2) * radius1

        Xc(6,:)= [ 1.0_rk,  0.0_rk,              0.0_rk]
        Xc(7,:)= [ 1.0_rk,  sqrt(3.0_rk),        0.0_rk]
        Xc(8,:)= [-0.5_rk,  sqrt(3.0_rk)/2.0_rk, 0.0_rk]
        Xc(9,:)= [-2.0_rk,  0.0_rk,              0.0_rk]
        Xc(10,:)=[-0.5_rk, -sqrt(3.0_rk)/2.0_rk, 0.0_rk]

        Xc(6:10,1:2) = Xc(6:10,1:2) * radius2

        ! Translate the control points
        do i = 1, size(Xc, 1)
            Xc(i,:) = center + Xc(i,:)
        end do

        ! Define weights for the control points
        Wc = [1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk,&
            1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk]

        ! Define knot vector
        knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/2.0_rk, 1.0_rk/2.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot2 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

        ! Set knot vector, control points, and weights
        call this%set(knot1, knot2, Xc, Wc)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_half_ring(this, center, radius1, radius2)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: center(:)
        real(rk), intent(in) :: radius1, radius2
        real(rk), allocatable :: Xc(:,:), Wc(:), knot1(:), knot2(:)
        integer :: i

        ! Define control points for half ring
        allocate(Xc(10, 3))
        Xc(1,:)  = [ 0.5_rk, 0.0_rk, 0.0_rk]
        Xc(2,:)  = [ 0.5_rk, 0.5_rk, 0.0_rk]
        Xc(3,:)  = [ 0.0_rk, 0.5_rk, 0.0_rk]
        Xc(4,:)  = [-0.5_rk, 0.5_rk, 0.0_rk]
        Xc(5,:)  = [-0.5_rk, 0.0_rk, 0.0_rk]

        Xc(1:5,1:2) = Xc(1:5,1:2) * radius1

        Xc(6,:)  = [ 0.5_rk, 0.0_rk, 0.0_rk]
        Xc(7,:)  = [ 0.5_rk, 0.5_rk, 0.0_rk]
        Xc(8,:)  = [ 0.0_rk, 0.5_rk, 0.0_rk]
        Xc(9,:)  = [-0.5_rk, 0.5_rk, 0.0_rk]
        Xc(10,:) = [-0.5_rk, 0.0_rk, 0.0_rk]

        Xc(6:10,1:2) = Xc(6:10,1:2) * radius2

        ! Translate the control points
        do i = 1, size(Xc, 1)
            Xc(i,:) = center + Xc(i,:)
        end do

        ! Define weights for the control points
        Wc = [1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk,&
            1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk]

        ! Define knot vector
        knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/2.0_rk, &
            1.0_rk/2.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot2 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

        ! Set knot vector, control points, and weights
        call this%set(knot1, knot2, Xc, Wc)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine nearest_point(this, point_Xg, nearest_Xg, nearest_Xt, id)
        class(nurbs_surface), intent(in) :: this
        real(rk), intent(in) :: point_Xg(:)
        real(rk), intent(out), allocatable, optional :: nearest_Xg(:)
        real(rk), intent(out), allocatable, optional :: nearest_Xt(:)
        integer, intent(out), optional :: id
        integer :: id_
        real(rk), allocatable :: distances(:)

        allocate(distances(this%ng(1)*this%ng(2)))
        distances = nearest_point_help_2d(this%ng, this%Xg, point_Xg)

        id_ = minloc(distances, dim=1)
        if (present(id)) id = id_
        if (present(nearest_Xg)) nearest_Xg = this%Xg(id_,:)
        if (present(nearest_Xt)) nearest_Xt = this%Xt(id_,:)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine nearest_point2(this, point_Xg, tol, maxit, nearest_Xt, nearest_Xg)

        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: point_Xg(:)
        real(rk), intent(in) :: tol
        integer, intent(in) :: maxit
        real(rk), intent(out) :: nearest_Xt(2)
        real(rk), allocatable, intent(out), optional :: nearest_Xg(:)
        real(rk):: xk(2), obj, grad(2), hess(2,2), dk(2), alphak, tau, beta, det_inv, Ainv(2,2), lower_bounds(2), upper_bounds(2)
        real(rk), allocatable :: Xg(:), Tgc(:), dTgc(:,:), d2Tgc(:,:), distances(:)
        integer :: k, l
        logical :: convergenz
        type(nurbs_surface) :: copy_this

        k = 0

        ! lower and upper bounds
        lower_bounds = [minval(this%knot1), minval(this%knot2)]
        upper_bounds = [maxval(this%knot1), maxval(this%knot2)]

        ! guess initial point
        copy_this = this
        call this%create(10,10)
        allocate(distances(copy_this%ng(1)*copy_this%ng(2)))
        distances = nearest_point_help_2d(copy_this%ng, copy_this%Xg, point_Xg)
        xk = copy_this%Xt(minloc(distances, dim=1),:)
        call copy_this%finalize()

        ! Check if xk is within the knot vector range
        if (xk(1) < minval(this%knot1)) then
            xk(1) = minval(this%knot1)
        else if (xk(1) > maxval(this%knot1)) then
            xk(1) = maxval(this%knot1)
        end if

        if (xk(2) < minval(this%knot2)) then
            xk(2) = minval(this%knot2)
        else if (xk(2) > maxval(this%knot2)) then
            xk(2) = maxval(this%knot2)
        end if

        convergenz = .false.

        allocate(Xg(size(this%Xc,2)))
        ! allocate(dTgc(size(this%Xc,1), 2))
        ! allocate(d2Tgc(size(this%Xc,1), 2))

        do while (.not. convergenz .and. k < maxit)

            ! objection, gradient and hessian
            Xg = this%cmp_Xg(xk)
            call this%derivative2(Xt=xk, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc) ! Tgc is not needed

            obj = norm2(Xg - point_Xg) + 0.001_rk ! add a small number to avoid division by zero

            grad(1) = dot_product((Xg-point_Xg)/obj, matmul(dTgc(:,1),this%Xc))
            grad(2) = dot_product((Xg-point_Xg)/obj, matmul(dTgc(:,2),this%Xc))

            hess(1,1) = ( dot_product(matmul(dTgc(:,1),this%Xc), matmul(dTgc(:,1),this%Xc)) + &
                dot_product((Xg-point_Xg), matmul(d2Tgc(1:this%nc(1)*this%nc(2)                        ,1),this%Xc)) ) /obj &
                - ( dot_product(Xg-point_Xg, matmul(dTgc(:,1), this%Xc))*grad(1) ) / obj**2
            hess(2,1) = ( dot_product(matmul(dTgc(:,1),this%Xc), matmul(dTgc(:,2),this%Xc)) + &
                dot_product((Xg-point_Xg), matmul(d2Tgc(this%nc(1)*this%nc(2)+1:2*this%nc(1)*this%nc(2),1),this%Xc)) ) /obj &
                - ( dot_product(Xg-point_Xg, matmul(dTgc(:,2), this%Xc))*grad(1) ) / obj**2
            hess(1,2) = ( dot_product(matmul(dTgc(:,2),this%Xc), matmul(dTgc(:,1),this%Xc)) + &
                dot_product((Xg-point_Xg), matmul(d2Tgc(1:this%nc(1)*this%nc(2)                        ,2),this%Xc)) ) /obj &
                - ( dot_product(Xg-point_Xg, matmul(dTgc(:,1), this%Xc))*grad(2) ) / obj**2
            hess(2,2) = ( dot_product(matmul(dTgc(:,2),this%Xc), matmul(dTgc(:,2),this%Xc)) + &
                dot_product((Xg-point_Xg), matmul(d2Tgc(this%nc(1)*this%nc(2)+1:2*this%nc(1)*this%nc(2),2),this%Xc)) ) /obj &
                - ( dot_product(Xg-point_Xg, matmul(dTgc(:,2), this%Xc))*grad(2) ) / obj**2

            ! debug
            print '(i3,1x,2e20.10,1x,e20.10)', k, xk, norm2(grad)

            if (norm2(grad) <= tol) then
                convergenz = .true.
                nearest_Xt = xk
                if (present(nearest_Xg)) nearest_Xg = this%cmp_Xg(nearest_Xt)
            else

                ! Inverse of Hessian
                det_inv = 1.0_rk/(hess(1,1)*hess(2,2) - hess(1,2)*hess(2,1))
                Ainv(1,1) =  hess(2,2)
                Ainv(2,1) = -hess(2,1)
                Ainv(1,2) = -hess(1,2)
                Ainv(2,2) =  hess(1,1)
                Ainv = Ainv * det_inv

                dk = - matmul(Ainv, grad)

                ! Backtracking-Armijo Line Search
                alphak = 1.0_rk
                tau = 0.5_rk     ! 0 < tau  < 1
                beta = 1.0e-4_rk ! 0 < beta < 1
                l = 0
                do while (.not. norm2(this%cmp_Xg(xk + alphak*dk) - point_Xg) <= obj + alphak*beta*dot_product(grad,dk) .and. l<50)
                    alphak = tau * alphak
                    l = l + 1
                end do

                xk = xk + alphak*dk
                ! Check if xk is within the knot vector range
                xk = max(min(xk, upper_bounds), lower_bounds)
                k = k + 1
            end if
        end do

    end subroutine
    !===============================================================================

end module forcad_nurbs_surface

!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Xg_nurbs_2d(Xt, knot1, knot2, degree, nc, ng, Xc, Wc) result(Xg)
    use forcad_utils, only: rk, basis_bspline, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:,:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    integer, intent(in) :: ng(2)
    real(rk), intent(in), contiguous :: Xc(:,:)
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable :: Xg(:,:)
    real(rk), allocatable :: Tgc(:)
    integer :: i

    allocate(Xg(ng(1)*ng(2), size(Xc,2)))
    allocate(Tgc(nc(1)*nc(2)))

    !$OMP PARALLEL DO PRIVATE(Tgc)
    do i = 1, ng(1)*ng(2)
        Tgc = kron(&
            basis_bspline(Xt(i,2), knot2, nc(2), degree(2)),&
            basis_bspline(Xt(i,1), knot1, nc(1), degree(1)))
        Tgc = Tgc*(Wc/(dot_product(Tgc,Wc)))
        Xg(i,:) = matmul(Tgc, Xc)
    end do
    !$OMP END PARALLEL DO
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Xg_nurbs_2d_1point(Xt, knot1, knot2, degree, nc, Xc, Wc) result(Xg)
    use forcad_utils, only: rk, basis_bspline, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    real(rk), intent(in), contiguous :: Xc(:,:)
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable :: Xg(:)
    real(rk), allocatable :: Tgc(:)

    allocate(Xg(size(Xc,2)))
    allocate(Tgc(nc(1)*nc(2)))

    Tgc = kron(&
        basis_bspline(Xt(2), knot2, nc(2), degree(2)),&
        basis_bspline(Xt(1), knot1, nc(1), degree(1)))
    Tgc = Tgc*(Wc/(dot_product(Tgc,Wc)))
    Xg= matmul(Tgc, Xc)
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Xg_bspline_2d(Xt, knot1, knot2, degree, nc, ng, Xc) result(Xg)
    use forcad_utils, only: rk, basis_bspline, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:,:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    integer, intent(in) :: ng(2)
    real(rk), intent(in), contiguous :: Xc(:,:)
    real(rk), allocatable :: Xg(:,:)
    integer :: i

    allocate(Xg(ng(1)*ng(2), size(Xc,2)))
    !$OMP PARALLEL DO
    do i = 1, ng(1)*ng(2)
        Xg(i,:) = matmul(kron(&
            basis_bspline(Xt(i,2), knot2, nc(2), degree(2)),&
            basis_bspline(Xt(i,1), knot1, nc(1), degree(1))),&
            Xc)
    end do
    !$OMP END PARALLEL DO
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Xg_bspline_2d_1point(Xt, knot1, knot2, degree, nc, Xc) result(Xg)
    use forcad_utils, only: rk, basis_bspline, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    real(rk), intent(in), contiguous :: Xc(:,:)
    real(rk), allocatable :: Xg(:)

    allocate(Xg(size(Xc,2)))
    Xg= matmul(kron(&
        basis_bspline(Xt(2), knot2, nc(2), degree(2)),&
        basis_bspline(Xt(1), knot1, nc(1), degree(1))),&
        Xc)
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_dTgc_nurbs_2d_vector(Xt, knot1, knot2, degree, nc, ng, Wc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_der, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:,:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    integer, intent(in) :: ng(2)
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable, intent(out) :: dTgc(:,:,:)
    real(rk), allocatable, intent(out) :: Tgc(:,:)
    real(rk), allocatable :: dBi(:,:), dB1(:), dB2(:)
    real(rk), allocatable :: Bi(:), B1(:), B2(:)
    integer :: i

    allocate(dTgc(ng(1)*ng(2), nc(1)*nc(2), 2), Tgc(ng(1)*ng(2), nc(1)*nc(2)))
    allocate(Bi(nc(1)*nc(2)), dBi(nc(1)*nc(2), 2))
    do i = 1, size(Xt, 1)
        call basis_bspline_der(Xt(i,1), knot1, nc(1), degree(1), dB1, B1)
        call basis_bspline_der(Xt(i,2), knot2, nc(2), degree(2), dB2, B2)

        Bi = kron(B2, B1)

        Tgc(i,:) = Bi*(Wc/(dot_product(Bi,Wc)))

        dBi(:,1) = kron(B2, dB1)
        dBi(:,2) = kron(dB2, B1)

        dTgc(i,:,1) = ( dBi(:,1)*Wc - Tgc(i,:)*dot_product(dBi(:,1),Wc) ) / dot_product(Bi,Wc)
        dTgc(i,:,2) = ( dBi(:,2)*Wc - Tgc(i,:)*dot_product(dBi(:,2),Wc) ) / dot_product(Bi,Wc)
    end do
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_dTgc_nurbs_2d_scalar(Xt, knot1, knot2, degree, nc, Wc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_der, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable, intent(out) :: dTgc(:,:)
    real(rk), allocatable, intent(out) :: Tgc(:)
    real(rk), allocatable :: dB1(:), dB2(:), dBi(:,:)
    real(rk), allocatable :: B1(:), B2(:), Bi(:)

    allocate(dTgc(nc(1)*nc(2), 2), Tgc(nc(1)*nc(2)))
    allocate(dBi(nc(1)*nc(2), 2), Bi(nc(1)*nc(2)))

    call basis_bspline_der(Xt(1), knot1, nc(1), degree(1), dB1, B1)
    call basis_bspline_der(Xt(2), knot2, nc(2), degree(2), dB2, B2)

    Bi = kron(B2, B1)
    Tgc = Bi*(Wc/(dot_product(Bi,Wc)))

    dBi(:,1) = kron(B2, dB1)
    dBi(:,2) = kron(dB2, B1)

    dTgc(:,1) = ( dBi(:,1)*Wc - Tgc*dot_product(dBi(:,1),Wc) ) / dot_product(Bi,Wc)
    dTgc(:,2) = ( dBi(:,2)*Wc - Tgc*dot_product(dBi(:,2),Wc) ) / dot_product(Bi,Wc)
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_dTgc_bspline_2d_vector(Xt, knot1, knot2, degree, nc, ng, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_der, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:,:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    integer, intent(in) :: ng(2)
    real(rk), allocatable, intent(out) :: dTgc(:,:,:)
    real(rk), allocatable, intent(out) :: Tgc(:,:)
    real(rk), allocatable :: dB1(:), dB2(:)
    real(rk), allocatable :: B1(:), B2(:)
    integer :: i

    allocate(dTgc(ng(1)*ng(2), nc(1)*nc(2), 2))
    !$OMP PARALLEL DO PRIVATE(dB1, dB2)
    do i = 1, size(Xt, 1)
        call basis_bspline_der(Xt(i,1), knot1, nc(1), degree(1), dB1, B1)
        call basis_bspline_der(Xt(i,2), knot2, nc(2), degree(2), dB2, B2)
        Tgc(i,:) = kron(B2, B1)

        dTgc(i,:,1) = kron(B2, dB1)
        dTgc(i,:,2) = kron(dB2, B1)
    end do
    !$OMP END PARALLEL DO
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_dTgc_bspline_2d_scalar(Xt, knot1, knot2, degree, nc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_der, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    real(rk), allocatable, intent(out) :: dTgc(:,:)
    real(rk), allocatable, intent(out) :: Tgc(:)
    real(rk), allocatable :: dTgc1(:), dTgc2(:)
    real(rk), allocatable :: Tgc1(:), Tgc2(:)

    allocate(dTgc(nc(1)*nc(2), 2))
    call basis_bspline_der(Xt(1), knot1, nc(1), degree(1), dTgc1, Tgc1)
    call basis_bspline_der(Xt(2), knot2, nc(2), degree(2), dTgc2, Tgc2)
    Tgc = kron(Tgc2, Tgc1)

    dTgc(:,1) = kron(Tgc2, dTgc1)
    dTgc(:,2) = kron(dTgc2, Tgc1)
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_d2Tgc_nurbs_2d_vector(Xt, knot1, knot2, degree, nc, ng, Wc, d2Tgc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_2der, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:,:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    integer, intent(in) :: ng(2)
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable, intent(out) :: d2Tgc(:,:,:)
    real(rk), allocatable, intent(out) :: dTgc(:,:,:)
    real(rk), allocatable, intent(out) :: Tgc(:,:)
    real(rk), allocatable :: d2Bi(:,:), d2B1(:), d2B2(:)
    real(rk), allocatable :: dBi(:,:), dB1(:), dB2(:)
    real(rk), allocatable :: Bi(:), B1(:), B2(:)
    real(rk), allocatable :: Tgci(:), dTgci(:)

    integer :: i

    allocate(d2Tgc(ng(1)*ng(2), 2*nc(1)*nc(2), 2))

    allocate(B1(nc(1)), B2(nc(2)))
    allocate(dB1(nc(1)), dB2(nc(2)))
    allocate(d2B1(nc(1)), d2B2(nc(2)))
    allocate(Bi(nc(1)*nc(2)), dBi(nc(1)*nc(2), 2), d2Bi(2*nc(1)*nc(2), 2))

    allocate(Tgci(nc(1)*nc(2)), dTgci(nc(1)*nc(2)))
    allocate(Tgc(ng(1)*ng(2), nc(1)*nc(2)), dTgc(ng(1)*ng(2), nc(1)*nc(2), 2))
    do i = 1, size(Xt, 1)
        call basis_bspline_2der(Xt(i,1), knot1, nc(1), degree(1), d2B1, dB1, B1)
        call basis_bspline_2der(Xt(i,2), knot2, nc(2), degree(2), d2B2, dB2, B2)

        Bi = kron(B2, B1)

        Tgc(i,:) = Bi*(Wc/(dot_product(Bi,Wc)))

        dBi(:,1) = kron(B2, dB1)
        dBi(:,2) = kron(dB2, B1)

        dTgc(i,:,1) = ( dBi(:,1)*Wc - Tgc(i,:)*dot_product(dBi(:,1),Wc) ) / dot_product(Bi,Wc)
        dTgc(i,:,2) = ( dBi(:,2)*Wc - Tgc(i,:)*dot_product(dBi(:,2),Wc) ) / dot_product(Bi,Wc)

        d2Bi(1:nc(1)*nc(2)              ,1) = kron(B2, d2B1)
        d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),1) = kron(dB2, dB1)
        d2Bi(1:nc(1)*nc(2)              ,2) = kron(dB2, dB1)
        d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),2) = kron(d2B2, B1)

        d2Tgc(i,1:nc(1)*nc(2)              ,1) = &
            (d2Bi(1:nc(1)*nc(2)              ,1)*Wc - 2.0_rk*dTgc(i,:,1)*dot_product(dBi(:,1),Wc)                                 &
            - Tgc(i,:)*dot_product(d2Bi(1:nc(1)*nc(2)              ,1),Wc)) / dot_product(Bi,Wc)
        d2Tgc(i,nc(1)*nc(2)+1:2*nc(1)*nc(2),1) = &
            (d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),1)*Wc - dTgc(i,:,1)*dot_product(dBi(:,2),Wc) - dTgc(i,:,2)*dot_product(dBi(:,1),Wc) &
            - Tgc(i,:)*dot_product(d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),1),Wc)) / dot_product(Bi,Wc)
        d2Tgc(i,1:nc(1)*nc(2)              ,2) = &
            (d2Bi(1:nc(1)*nc(2)              ,2)*Wc - dTgc(i,:,1)*dot_product(dBi(:,2),Wc) - dTgc(i,:,2)*dot_product(dBi(:,1),Wc) &
            - Tgc(i,:)*dot_product(d2Bi(1:nc(1)*nc(2)              ,2),Wc)) / dot_product(Bi,Wc)
        d2Tgc(i,nc(1)*nc(2)+1:2*nc(1)*nc(2),2) = &
            (d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),2)*Wc - 2.0_rk*dTgc(i,:,2)*dot_product(dBi(:,2),Wc)                                 &
            - Tgc(i,:)*dot_product(d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),2),Wc)) / dot_product(Bi,Wc)
    end do
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_d2Tgc_nurbs_2d_scalar(Xt, knot1, knot2, degree, nc, Wc, d2Tgc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_2der, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable, intent(out) :: d2Tgc(:,:)
    real(rk), allocatable, intent(out) :: dTgc(:,:)
    real(rk), allocatable, intent(out) :: Tgc(:)
    real(rk), allocatable :: d2Bi(:,:), d2B1(:), d2B2(:)
    real(rk), allocatable :: dBi(:,:), dB1(:), dB2(:)
    real(rk), allocatable :: Bi(:), B1(:), B2(:)

    allocate(B1(nc(1)), B2(nc(2)))
    allocate(dB1(nc(1)), dB2(nc(2)))
    allocate(d2B1(nc(1)), d2B2(nc(2)))
    allocate(Bi(nc(1)*nc(2)), dBi(nc(1)*nc(2), 2), d2Bi(2*nc(1)*nc(2), 2))

    allocate(Tgc(nc(1)*nc(2)), dTgc(nc(1)*nc(2), 2), d2Tgc(2*nc(1)*nc(2), 2))

    call basis_bspline_2der(Xt(1), knot1, nc(1), degree(1), d2B1, dB1, B1)
    call basis_bspline_2der(Xt(2), knot2, nc(2), degree(2), d2B2, dB2, B2)

    Bi = kron(B2, B1)

    Tgc = Bi*(Wc/(dot_product(Bi,Wc)))

    dBi(:,1) = kron(B2, dB1)
    dBi(:,2) = kron(dB2, B1)

    dTgc(:,1) = ( dBi(:,1)*Wc - Tgc*dot_product(dBi(:,1),Wc) ) / dot_product(Bi,Wc)
    dTgc(:,2) = ( dBi(:,2)*Wc - Tgc*dot_product(dBi(:,2),Wc) ) / dot_product(Bi,Wc)

    d2Bi(1:nc(1)*nc(2)              ,1) = kron(B2, d2B1)
    d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),1) = kron(dB2, dB1)
    d2Bi(1:nc(1)*nc(2)              ,2) = kron(dB2, dB1)
    d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),2) = kron(d2B2, B1)

    d2Tgc(1:nc(1)*nc(2)              ,1) = &
        (d2Bi(1:nc(1)*nc(2)              ,1)*Wc - 2.0_rk*dTgc(:,1)*dot_product(dBi(:,1),Wc)                               &
        - Tgc*dot_product(d2Bi(1:nc(1)*nc(2)              ,1),Wc)) / dot_product(Bi,Wc)
    d2Tgc(nc(1)*nc(2)+1:2*nc(1)*nc(2),1) = &
        (d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),1)*Wc - dTgc(:,1)*dot_product(dBi(:,2),Wc) - dTgc(:,2)*dot_product(dBi(:,1),Wc) &
        - Tgc*dot_product(d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),1),Wc)) / dot_product(Bi,Wc)
    d2Tgc(1:nc(1)*nc(2)              ,2) = &
        (d2Bi(1:nc(1)*nc(2)              ,2)*Wc - dTgc(:,1)*dot_product(dBi(:,2),Wc) - dTgc(:,2)*dot_product(dBi(:,1),Wc) &
        - Tgc*dot_product(d2Bi(1:nc(1)*nc(2)              ,2),Wc)) / dot_product(Bi,Wc)
    d2Tgc(nc(1)*nc(2)+1:2*nc(1)*nc(2),2) = &
        (d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),2)*Wc - 2.0_rk*dTgc(:,2)*dot_product(dBi(:,2),Wc)                               &
        - Tgc*dot_product(d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),2),Wc)) / dot_product(Bi,Wc)
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_d2Tgc_bspline_2d_vector(Xt, knot1, knot2, degree, nc, ng, d2Tgc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_2der, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:,:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    integer, intent(in) :: ng(2)
    real(rk), allocatable, intent(out) :: d2Tgc(:,:,:)
    real(rk), allocatable, intent(out) :: dTgc(:,:,:)
    real(rk), allocatable, intent(out) :: Tgc(:,:)
    real(rk), allocatable :: d2B1(:), d2B2(:)
    real(rk), allocatable :: dB1(:), dB2(:)
    real(rk), allocatable :: B1(:), B2(:)
    integer :: i

    allocate(d2Tgc(ng(1)*ng(2), 2*nc(1)*nc(2), 2))
    allocate(dTgc(ng(1)*ng(2), nc(1)*nc(2), 2))
    allocate(Tgc(ng(1)*ng(2), nc(1)*nc(2)))
    do i = 1, size(Xt, 1)
        call basis_bspline_2der(Xt(i,1), knot1, nc(1), degree(1), d2B1, dB1, B1)
        call basis_bspline_2der(Xt(i,2), knot2, nc(2), degree(2), d2B2, dB2, B2)

        Tgc(i,:) = kron(B2, B1)

        dTgc(i,:,1) = kron(B2, dB1)
        dTgc(i,:,2) = kron(dB2, B1)

        d2Tgc(i,1:nc(1)*nc(2)              ,1) = kron(B2, d2B1)
        d2Tgc(i,nc(1)*nc(2)+1:2*nc(1)*nc(2),1) = kron(dB2, dB1)
        d2Tgc(i,1:nc(1)*nc(2)              ,2) = kron(dB2, dB1)
        d2Tgc(i,nc(1)*nc(2)+1:2*nc(1)*nc(2),2) = kron(d2B2, B1)
    end do
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_d2Tgc_bspline_2d_scalar(Xt, knot1, knot2, degree, nc, d2Tgc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_2der, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    real(rk), allocatable, intent(out) :: d2Tgc(:,:)
    real(rk), allocatable, intent(out) :: dTgc(:,:)
    real(rk), allocatable, intent(out) :: Tgc(:)
    real(rk), allocatable :: d2B1(:), d2B2(:)
    real(rk), allocatable :: dB1(:), dB2(:)
    real(rk), allocatable :: B1(:), B2(:)

    allocate(d2Tgc(2*nc(1)*nc(2), 2))
    allocate(dTgc(nc(1)*nc(2), 2))
    allocate(Tgc(nc(1)*nc(2)))
    call basis_bspline_2der(Xt(1), knot1, nc(1), degree(1), d2B1, dB1, B1)
    call basis_bspline_2der(Xt(2), knot2, nc(2), degree(2), d2B2, dB2, B2)
    Tgc = kron(B2, B1)

    dTgc(:,1) = kron(B2, dB1)
    dTgc(:,2) = kron(dB2, B1)

    d2Tgc(1:nc(1)*nc(2)              ,1) = kron(B2, d2B1)
    d2Tgc(nc(1)*nc(2)+1:2*nc(1)*nc(2),1) = kron(dB2, dB1)
    d2Tgc(1:nc(1)*nc(2)              ,2) = kron(dB2, dB1)
    d2Tgc(nc(1)*nc(2)+1:2*nc(1)*nc(2),2) = kron(d2B2, B1)
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Tgc_nurbs_2d_vector(Xt, knot1, knot2, degree, nc, ng, Wc) result(Tgc)
    use forcad_utils, only: rk, basis_bspline, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:,:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    integer, intent(in) :: ng(2)
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable :: Tgc(:,:)
    real(rk), allocatable :: Tgci(:)
    integer :: i

    allocate(Tgc(ng(1)*ng(2), nc(1)*nc(2)))
    allocate(Tgci(nc(1)*nc(2)))
    !$OMP PARALLEL DO PRIVATE(Tgci)
    do i = 1, size(Xt, 1)
        Tgci = kron(&
            basis_bspline(Xt(i,2), knot2, nc(2), degree(2)),&
            basis_bspline(Xt(i,1), knot1, nc(1), degree(1)))
        Tgc(i,:) = Tgci*(Wc/(dot_product(Tgci,Wc)))
    end do
    !$OMP END PARALLEL DO
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Tgc_nurbs_2d_scalar(Xt, knot1, knot2, degree, nc, Wc) result(Tgc)
    use forcad_utils, only: rk, basis_bspline, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable :: Tgc(:)

    allocate(Tgc(nc(1)*nc(2)))
    Tgc = kron(&
        basis_bspline(Xt(2), knot2, nc(2), degree(2)),&
        basis_bspline(Xt(1), knot1, nc(1), degree(1)))
    Tgc = Tgc*(Wc/(dot_product(Tgc,Wc)))
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Tgc_bspline_2d_vector(Xt, knot1, knot2, degree, nc, ng) result(Tgc)
    use forcad_utils, only: rk, basis_bspline, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:,:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    integer, intent(in) :: ng(2)
    real(rk), allocatable :: Tgc(:,:)
    integer :: i

    allocate(Tgc(ng(1)*ng(2), nc(1)*nc(2)))
    !$OMP PARALLEL DO
    do i = 1, size(Xt, 1)
        Tgc(i,:) = kron(&
            basis_bspline(Xt(i,2), knot2, nc(2), degree(2)),&
            basis_bspline(Xt(i,1), knot1, nc(1), degree(1)))
    end do
    !$OMP END PARALLEL DO
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Tgc_bspline_2d_scalar(Xt, knot1, knot2, degree, nc) result(Tgc)
    use forcad_utils, only: rk, basis_bspline, kron

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot1(:), knot2(:)
    integer, intent(in) :: degree(2)
    integer, intent(in) :: nc(2)
    real(rk), allocatable :: Tgc(:)

    allocate(Tgc(nc(1)*nc(2)))
    Tgc = kron(&
        basis_bspline(Xt(2), knot2, nc(2), degree(2)),&
        basis_bspline(Xt(1), knot1, nc(1), degree(1)))
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function nearest_point_help_2d(ng, Xg, point_Xg) result(distances)
    use forcad_utils, only: rk

    implicit none
    integer, intent(in) :: ng(2)
    real(rk), intent(in), contiguous :: Xg(:,:)
    real(rk), intent(in), contiguous :: point_Xg(:)
    real(rk), allocatable :: distances(:)
    integer :: i

    allocate(distances(ng(1)*ng(2)))
    !$OMP PARALLEL DO
    do i = 1, ng(1)*ng(2)
        distances(i) = norm2(Xg(i,:) - point_Xg)
    end do
    !$OMP END PARALLEL DO

end function
!===============================================================================
