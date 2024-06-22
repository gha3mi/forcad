!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> This module defines the 'nurbs_curve' type for representing a Non-Uniform Rational B-Spline (NURBS) curve.
module forcad_nurbs_curve

    use forcad_utils, only: rk, basis_bspline, elemConn_C0, compute_multiplicity, compute_knot_vector, basis_bspline_der,&
        insert_knot_A_5_1, findspan, elevate_degree_A_5_9, remove_knots_A_5_8, &
        elemConn_Cn, unique, rotation

    implicit none

    private
    public nurbs_curve

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    type nurbs_curve
        real(rk), allocatable, private :: Xc(:,:) !! Control points (2D array: [nc, dim])
        real(rk), allocatable, private :: Xg(:,:) !! Geometry points (2D array: [ng, dim])
        real(rk), allocatable, private :: Wc(:)   !! Weights for control points (1D array: [nc])
        real(rk), allocatable, private :: Xt(:)   !! Evaluation points (1D array: [ng])
        real(rk), allocatable, private :: knot(:) !! Knot vector (1D array)
        integer, private :: degree                !! Degree (order) of the curve
        integer, private :: nc                    !! Number of control points
        integer, private :: ng                    !! Number of geometry points
        integer, allocatable, private :: elemConn_Xc_vis(:,:) !! Connectivity for visualization of control points
        integer, allocatable, private :: elemConn_Xg_vis(:,:) !! Connectivity for visualization of geometry points
        integer, allocatable, private :: elemConn(:,:)        !! IGA element connectivity
    contains
        procedure :: set1                  !!> Set knot vector, control points and weights for the NURBS curve object
        procedure :: set2                  !!> Set NURBS curve using nodes of parameter space, degree, continuity, control points and weights
        procedure :: set3                  !!> Set Bezier or Rational Bezier curve using control points and weights
        generic :: set => set1, set2, set3 !!> Set NURBS curve
        procedure :: create                !!> Generate geometry points
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
        procedure :: get_Xt                !!> Get parameter values
        procedure, private :: get_knot_all !!> Get all knot vectors
        procedure, private :: get_knoti    !!> Get i-th knot value
        generic :: get_knot => get_knoti, get_knot_all !!> Get knot vector
        procedure :: get_ng                !!> Get number of geometry points
        procedure :: cmp_degree            !!> Compute degree of the NURBS curve
        procedure :: get_degree            !!> Get degree of the NURBS curve
        procedure :: finalize              !!> Finalize the NURBS curve object
        procedure :: cmp_elem_Xc_vis       !!> Generate connectivity for control points
        procedure :: cmp_elem_Xg_vis       !!> Generate connectivity for geometry points
        procedure :: cmp_elem              !!> Generate IGA element connectivity
        procedure :: get_elem_Xc_vis       !!> Get connectivity for control points
        procedure :: get_elem_Xg_vis       !!> Get connectivity for geometry points
        procedure :: get_elem              !!> Get IGA element connectivity
        procedure :: set_elem_Xc_vis       !!> Set connectivity for control points
        procedure :: set_elem_Xg_vis       !!> Set connectivity for geometry points
        procedure :: set_elem              !!> Set IGA element connectivity
        procedure :: export_Xc             !!> Export control points to VTK file
        procedure :: export_Xg             !!> Export geometry points to VTK file
        procedure :: modify_Xc             !!> Modify control points
        procedure :: modify_Wc             !!> Modify weights
        procedure :: get_multiplicity      !!> Compute and return the multiplicity of the knots
        procedure :: get_continuity        !!> Compute and return the continuity of the curve
        procedure :: cmp_nc                !!> Compute number of required control points
        procedure :: get_nc                !!> Get number of control points
        procedure :: insert_knots          !!> Insert knots into the knot vector
        procedure :: elevate_degree        !!> Elevate the degree of the curve
        procedure, private :: basis_vector !!> Compute the basis functions of the NURBS curve
        procedure, private :: basis_scalar !!> Compute the basis functions of the NURBS curve
        generic :: basis => basis_vector, basis_scalar   !!> Compute the basis functions of the NURBS curve
        procedure, private :: derivative_vector      !!> Compute the derivative of the NURBS curve
        procedure, private :: derivative_scalar      !!> Compute the derivative of the NURBS curve
        generic :: derivative => derivative_vector, derivative_scalar   !!> Compute the derivative of the NURBS curve
        procedure, private :: derivative2_vector     !!> Compute the second derivative of the NURBS curve
        procedure, private :: derivative2_scalar     !!> Compute the second derivative of the NURBS curve
        generic :: derivative2 => derivative2_vector, derivative2_scalar !!> Compute the second derivative of the NURBS curve
        procedure :: is_rational           !!> Check if the NURBS curve is rational
        procedure :: remove_knots          !!> Remove knots from the knot vector
        procedure :: rotate_Xc             !!> Rotate control points
        procedure :: rotate_Xg             !!> Rotate geometry points
        procedure :: translate_Xc          !!> Translate control points
        procedure :: translate_Xg          !!> Translate geometry points
        procedure :: show                  !!> Show the NURBS object using PyVista
        procedure :: nearest_point         !!> Find the nearest point on the NURBS curve (Approximation)
        procedure :: nearest_point2        !!> Find the nearest point on the NURBS curve (Minimization - Newton's method)

        ! Shapes
        procedure :: set_circle            !!> Set a circle
        procedure :: set_half_circle       !!> Set a half circle
        procedure :: set_C                 !!> Set a C-shape
    end type
    !===============================================================================


    interface compute_Xg
        pure function compute_Xg_nurbs_1d(f_Xt, f_knot, f_degree, f_nc, f_ng, f_Xc, f_Wc) result(f_Xg)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            integer, intent(in) :: f_ng
            real(rk), intent(in), contiguous :: f_Xc(:,:)
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable :: f_Xg(:,:)
        end function

        pure function compute_Xg_bspline_1d(f_Xt, f_knot, f_degree, f_nc, f_ng, f_Xc) result(f_Xg)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            integer, intent(in) :: f_ng
            real(rk), intent(in), contiguous :: f_Xc(:,:)
            real(rk), allocatable :: f_Xg(:,:)
        end function

        pure function compute_Xg_nurbs_1d_1point(f_Xt, f_knot, f_degree, f_nc, f_Xc, f_Wc) result(f_Xg)
            import :: rk
            real(rk), intent(in) :: f_Xt
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            real(rk), intent(in), contiguous :: f_Xc(:,:)
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable :: f_Xg(:)
        end function

        pure function compute_Xg_bspline_1d_1point(f_Xt, f_knot, f_degree, f_nc, f_Xc) result(f_Xg)
            import :: rk
            real(rk), intent(in) :: f_Xt
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            real(rk), intent(in), contiguous :: f_Xc(:,:)
            real(rk), allocatable :: f_Xg(:)
        end function
    end interface

    interface compute_d2Tgc
        pure subroutine compute_d2Tgc_nurbs_1d_vector(f_Xt, f_knot, f_degree, f_nc, f_ng, f_Wc, f_d2Tgc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            integer, intent(in) :: f_ng
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable, intent(out) :: f_d2Tgc(:,:)
            real(rk), allocatable, intent(out) :: f_dTgc(:,:)
            real(rk), allocatable, intent(out) :: f_Tgc(:,:)
        end subroutine

        pure subroutine compute_d2Tgc_bspline_1d_vector(f_Xt, f_knot, f_degree, f_nc, f_ng, f_d2Tgc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            integer, intent(in) :: f_ng
            real(rk), allocatable, intent(out) :: f_d2Tgc(:,:)
            real(rk), allocatable, intent(out) :: f_dTgc(:,:)
            real(rk), allocatable, intent(out) :: f_Tgc(:,:)
        end subroutine

        pure subroutine compute_d2Tgc_nurbs_1d_scalar(f_Xt, f_knot, f_degree, f_nc, f_Wc, f_d2Tgc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in) :: f_Xt
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable, intent(out) :: f_d2Tgc(:)
            real(rk), allocatable, intent(out) :: f_dTgc(:)
            real(rk), allocatable, intent(out) :: f_Tgc(:)
        end subroutine

        pure subroutine compute_d2Tgc_bspline_1d_scalar(f_Xt, f_knot, f_degree, f_nc, f_d2Tgc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in) :: f_Xt
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            real(rk), allocatable, intent(out) :: f_d2Tgc(:)
            real(rk), allocatable, intent(out) :: f_dTgc(:)
            real(rk), allocatable, intent(out) :: f_Tgc(:)
        end subroutine
    end interface

    interface compute_dTgc
        pure subroutine compute_dTgc_nurbs_1d_vector(f_Xt, f_knot, f_degree, f_nc, f_ng, f_Wc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            integer, intent(in) :: f_ng
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable, intent(out) :: f_dTgc(:,:)
            real(rk), allocatable, intent(out) :: f_Tgc(:,:)
        end subroutine

        pure subroutine compute_dTgc_bspline_1d_vector(f_Xt, f_knot, f_degree, f_nc, f_ng, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            integer, intent(in) :: f_ng
            real(rk), allocatable, intent(out) :: f_dTgc(:,:)
            real(rk), allocatable, intent(out) :: f_Tgc(:,:)
        end subroutine

        pure subroutine compute_dTgc_nurbs_1d_scalar(f_Xt, f_knot, f_degree, f_nc, f_Wc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in) :: f_Xt
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable, intent(out) :: f_dTgc(:)
            real(rk), allocatable, intent(out) :: f_Tgc(:)
        end subroutine

        pure subroutine compute_dTgc_bspline_1d_scalar(f_Xt, f_knot, f_degree, f_nc, f_dTgc, f_Tgc)
            import :: rk
            real(rk), intent(in) :: f_Xt
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            real(rk), allocatable, intent(out) :: f_dTgc(:)
            real(rk), allocatable, intent(out) :: f_Tgc(:)
        end subroutine
    end interface

    interface compute_Tgc
        pure function compute_Tgc_nurbs_1d_vector(f_Xt, f_knot, f_degree, f_nc, f_ng, f_Wc) result(f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            integer, intent(in) :: f_ng
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable :: f_Tgc(:,:)
        end function

        pure function compute_Tgc_bspline_1d_vector(f_Xt, f_knot, f_degree, f_nc, f_ng) result(f_Tgc)
            import :: rk
            real(rk), intent(in), contiguous :: f_Xt(:)
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            integer, intent(in) :: f_ng
            real(rk), allocatable :: f_Tgc(:,:)
        end function

        pure function compute_Tgc_nurbs_1d_scalar(f_Xt, f_knot, f_degree, f_nc, f_Wc) result(f_Tgc)
            import :: rk
            real(rk), intent(in) :: f_Xt
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            real(rk), intent(in), contiguous :: f_Wc(:)
            real(rk), allocatable :: f_Tgc(:)
        end function

        pure function compute_Tgc_bspline_1d_scalar(f_Xt, f_knot, f_degree, f_nc) result(f_Tgc)
            import :: rk
            real(rk), intent(in) :: f_Xt
            real(rk), intent(in), contiguous :: f_knot(:)
            integer, intent(in) :: f_degree
            integer, intent(in) :: f_nc
            real(rk), allocatable :: f_Tgc(:)
        end function
    end interface

    interface
        pure function nearest_point_help_1d(f_ng, f_Xg, f_point_Xg) result(f_distances)
            import :: rk
            integer, intent(in) :: f_ng
            real(rk), intent(in), contiguous :: f_Xg(:,:)
            real(rk), intent(in), contiguous :: f_point_Xg(:)
            real(rk), allocatable :: f_distances(:)
        end function
    end interface

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set knot vector, control points and weights for the NURBS curve object.
    pure subroutine set1(this, knot, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: knot(:)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)

        if (allocated(this%knot)) deallocate(this%knot)
        if (allocated(this%Xc)) deallocate(this%Xc)

        this%knot = knot
        call this%cmp_degree()
        this%Xc = Xc
        this%nc = size(this%Xc, 1)
        if (present(Wc)) then
            if (size(Wc) /= this%nc) then
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
    !> Set NURBS curve using nodes of parameter space (Xth), degree, continuity, control points and weights.
    pure subroutine set2(this, Xth_dir, degree, continuity, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xth_dir(:)
        integer, intent(in) :: degree
        integer, intent(in), contiguous :: continuity(:)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)

        if (allocated(this%knot)) deallocate(this%knot)
        if (allocated(this%Xc)) deallocate(this%Xc)

        this%knot = compute_knot_vector(Xth_dir, degree, continuity)
        this%degree = degree
        this%Xc = Xc
        this%nc = size(this%Xc, 1)
        if (present(Wc)) then
            if (size(Wc) /= this%nc) then
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
    !> Set Bezier or Rational Bezier curve using control points and weights.
    pure subroutine set3(this, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)

        if (allocated(this%knot)) deallocate(this%knot)
        if (allocated(this%Xc)) deallocate(this%Xc)

        this%Xc = Xc
        this%nc = size(this%Xc, 1)

        if (allocated(this%knot)) deallocate(this%knot)
        allocate(this%knot(2*this%nc))
        this%knot(1:this%nc) = 0.0_rk
        this%knot(this%nc+1:2*this%nc) = 1.0_rk

        call this%cmp_degree()
        if (present(Wc)) then
            if (size(Wc) /= this%nc) then
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
    pure subroutine create(this, res, Xt)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: res
        real(rk), intent(in), contiguous, optional :: Xt(:)
        integer :: i

        ! check
        if (.not.allocated(this%Xc)) then
            error stop 'Control points are not set.'
        end if

        if (.not.allocated(this%knot)) then
            error stop 'Knot vector is not set.'
        end if

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            this%Xt = Xt
        elseif (present(res)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            allocate(this%Xt(res))
            this%Xt = [(real(i-1, rk) / real(res-1, rk), i=1, res)]
            ! else
            ! this%Xt = this%Xt
        end if

        ! Set number of geometry points
        this%ng = size(this%Xt)

        ! Allocate memory for geometry points
        if (allocated(this%Xg)) deallocate(this%Xg)

        if (this%is_rational()) then ! NURBS
            this%Xg = compute_Xg(&
                this%Xt, this%knot, this%degree, this%nc, this%ng, this%Xc, this%Wc)
        else ! B-Spline
            this%Xg = compute_Xg(&
                this%Xt, this%knot, this%degree, this%nc, this%ng, this%Xc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_Xg(this, Xt) result(Xg)
        class(nurbs_curve), intent(in) :: this
        real(rk), intent(in) :: Xt
        real(rk), allocatable :: Xg(:)

        ! check
        if (.not.allocated(this%Xc)) then
            error stop 'Control points are not set.'
        end if

        if (.not.allocated(this%knot)) then
            error stop 'Knot vector is not set.'
        end if

        if (this%is_rational()) then ! NURBS
            Xg = compute_Xg(Xt, this%knot, this%degree, this%nc, this%Xc, this%Wc)
        else ! B-Spline
            Xg = compute_Xg(Xt, this%knot, this%degree, this%nc, this%Xc)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xc_all(this) result(Xc)
        class(nurbs_curve), intent(in) :: this
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
        class(nurbs_curve), intent(in) :: this
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
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: n
        integer, intent(in) :: dir
        real(rk) :: Xc

        if (allocated(this%Xc)) then
            if (n<lbound(this%Xc,1) .or. n>ubound(this%Xc,1)) then
                error stop 'Invalid index for control points.'
            end if
            if (dir<lbound(this%Xc,2) .or. dir>ubound(this%Xc,2)) then
                error stop 'Invalid index for control points.'
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
        class(nurbs_curve), intent(in) :: this
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
        class(nurbs_curve), intent(in) :: this
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
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: n
        integer, intent(in) :: dir
        real(rk) :: Xg

        if (allocated(this%Xg)) then
            if (n<lbound(this%Xg,1) .or. n>ubound(this%Xg,1)) then
                error stop 'Invalid index for geometry points.'
            end if
            if (dir<lbound(this%Xg,2) .or. dir>ubound(this%Xg,2)) then
                error stop 'Invalid index for geometry points.'
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
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: Wc(:)

        if (allocated(this%Wc)) then
            Wc = this%Wc
        else
            error stop 'The NURBS curve is not rational or weights are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Wci(this, n) result(Wc)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: n
        real(rk) :: Wc

        if (allocated(this%Wc)) then
            if (n<lbound(this%Wc,1) .or. n>ubound(this%Wc,1)) then
                error stop 'Invalid index for weights.'
            end if
            Wc = this%Wc(n)
        else
            error stop 'The NURBS curve is not rational or weights are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xt(this) result(Xt)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: Xt(:)

        if (allocated(this%Xt)) then
            Xt = this%Xt
        else
            error stop 'Parameter values are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_ng(this) result(ng)
        class(nurbs_curve), intent(in) :: this
        integer :: ng

        ng = this%ng
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine cmp_degree(this)
        class(nurbs_curve), intent(inout) :: this
        integer, allocatable :: m(:)

        m = this%get_multiplicity()

        this%degree = m(1) - 1
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_degree(this) result(degree)
        class(nurbs_curve), intent(in) :: this
        integer :: degree

        degree = this%degree
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_knot_all(this) result(knot)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: knot(:)

        if (allocated(this%knot)) then
            knot = this%knot
        else
            error stop 'Knot vector is not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_knoti(this,i) result(knot)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: i
        real(rk) :: knot

        if (allocated(this%knot)) then
            if (i < 1 .or. i > size(this%knot)) then
                error stop 'Invalid index for knot vector.'
            else
                knot = this%knot(i)
            end if
        else
            error stop 'Knot vector is not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine finalize(this)
        class(nurbs_curve), intent(inout) :: this
        if (allocated(this%Xc)) deallocate(this%Xc)
        if (allocated(this%Xg)) deallocate(this%Xg)
        if (allocated(this%Wc)) deallocate(this%Wc)
        if (allocated(this%Xt)) deallocate(this%Xt)
        if (allocated(this%knot)) deallocate(this%knot)
        if (allocated(this%elemConn_Xc_vis)) deallocate(this%elemConn_Xc_vis)
        if (allocated(this%elemConn_Xg_vis)) deallocate(this%elemConn_Xg_vis)
        if (allocated(this%elemConn)) deallocate(this%elemConn)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elem_Xc_vis(this, p) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), optional :: p

        if (present(p)) then
            elemConn = elemConn_C0(this%nc,p)
        else
            elemConn = elemConn_C0(this%nc,1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elem_Xg_vis(this, p) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), optional :: p

        if (present(p)) then
            elemConn = elemConn_C0(this%ng,p)
        else
            elemConn = elemConn_C0(this%ng,1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xc(this, filename)
        class(nurbs_curve), intent(in) :: this
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

        if (size(this%Xc, 2) == 2) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xc(i,1), this%Xc(i,2), 0.0_rk , i = 1, nc)
        elseif (size(this%Xc, 2) == 3) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xc(i,1), this%Xc(i,2), this%Xc(i,3) , i = 1, nc)
        else
            error stop 'Invalid dimension of the control points.'
        end if

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*(size(elemConn,2)+1)
        write(nunit,'(g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0)')&
            (2, elemConn(i,1)-1,elemConn(i,2)-1, i = 1, size(elemConn,1))

        write(nunit,'(a," ",g0)') 'CELL_TYPES', size(elemConn,1)
        write(nunit,'(g0)') (3 , i = 1, size(elemConn,1))
        close(nunit)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xg(this, filename)
        class(nurbs_curve), intent(in) :: this
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

        if (size(this%Xg, 2) == 2) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xg(i,1), this%Xg(i,2), 0.0_rk , i = 1, ng)
        elseif (size(this%Xg, 2) == 3) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xg(i,1), this%Xg(i,2), this%Xg(i,3) , i = 1, ng)
        else
            error stop 'Invalid dimension of the geometry points.'
        end if

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*(size(elemConn,2)+1)
        write(nunit,'(g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0)')&
            (2, elemConn(i,1)-1,elemConn(i,2)-1, i = 1, size(elemConn,1))

        write(nunit,'(a," ",g0)') 'CELL_TYPES', size(elemConn,1)
        write(nunit,'(g0)') (3 , i = 1, size(elemConn,1))
        close(nunit)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine modify_Xc(this,X,num,dir)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: X
        integer, intent(in) :: num
        integer, intent(in) :: dir

        if (allocated(this%Xc)) then
            this%Xc(num,dir) = X
            if (allocated(this%Wc)) then
                call this%set(knot=this%get_knot(), Xc=this%get_Xc(), Wc=this%get_Wc())
            else
                call this%set(knot=this%get_knot(), Xc=this%get_Xc())
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
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: W
        integer, intent(in) :: num

        if (allocated(this%Wc)) then
            this%Wc(num) = W
            if (allocated(this%knot)) then
                call this%set(knot=this%get_knot(), Xc=this%get_Xc(), Wc=this%get_Wc())
            else
                call this%set(Xc=this%get_Xc(), Wc=this%get_Wc())
            end if
        else
            error stop 'The NURBS curve is not rational.'
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_multiplicity(this) result(m)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: m(:)

        ! check
        if (.not.allocated(this%knot)) then
            error stop 'Knot vector is not set.'
        else
            m = compute_multiplicity(this%knot)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_continuity(this) result(c)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: c(:)

        ! check
        if (.not.allocated(this%knot)) then
            error stop 'Knot vector is not set.'
        else
            c = this%degree - compute_multiplicity(this%knot)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine cmp_nc(this)
        class(nurbs_curve), intent(inout) :: this

        this%nc = sum(compute_multiplicity(this%knot)) - this%degree - 1
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_nc(this) result(nc)
        class(nurbs_curve), intent(in) :: this
        integer :: nc
        nc = this%nc
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine insert_knots(this,Xth,r)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xth(:)
        integer, intent(in), contiguous :: r(:)
        integer :: k, i, s, dim, j, n_new
        real(rk), allocatable :: Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:)

        if (this%is_rational()) then ! NURBS

            do i = 1, size(Xth)
                k = findspan(this%nc-1,this%degree,Xth(i),this%knot)
                if (this%knot(k+1) == Xth(i)) then
                    s = compute_multiplicity(this%knot,Xth(i))
                else
                    s = 0
                end if

                dim = size(this%Xc,2)
                allocate(Xcw(size(this%Xc,1),dim+1))

                do j = 1, size(this%Xc,1)
                    Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                end do
                Xcw(:,dim+1) = this%Wc(:)

                call insert_knot_A_5_1(&
                    this%degree,&
                    this%knot,&
                    Xcw,&
                    Xth(i),&
                    k,&
                    s,&
                    r(i),&
                    n_new,&
                    knot_new,&
                    Xcw_new)

                allocate(Xc_new(1:n_new+1,1:dim))
                allocate(Wc_new(1:n_new+1))
                do j = 1, n_new+1
                    Xc_new(j,1:dim) = Xcw_new(j-1,1:dim)/Xcw_new(j-1,dim+1)
                    Wc_new(j) = Xcw_new(j-1,dim+1)
                end do

                call this%set(knot=knot_new, Xc=Xc_new, Wc=Wc_new)
                deallocate(Xcw, Xcw_new, Xc_new, Wc_new)
            end do

        else ! B-Spline

            do i = 1, size(Xth)
                k = findspan(this%nc-1,this%degree,Xth(i),this%knot)
                if (this%knot(k+1) == Xth(i)) then
                    s = compute_multiplicity(this%knot,Xth(i))
                else
                    s = 0
                end if

                call insert_knot_A_5_1(&
                    this%degree,&
                    this%knot,&
                    this%Xc,&
                    Xth(i),&
                    k,&
                    s,&
                    r(i),&
                    n_new,&
                    knot_new,&
                    Xc_new)

                call this%set(knot=knot_new, Xc=Xc_new)
            end do

        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine elevate_degree(this, t)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in) :: t
        real(rk), allocatable :: Xcw(:,:), Xcw_new(:,:), knot_new(:), Xc_new(:,:), Wc_new(:)
        integer :: dim, j, nc_new

        if (this%is_rational()) then ! NURBS

            dim = size(this%Xc,2)
            allocate(Xcw(size(this%Xc,1),dim+1))
            do j = 1, size(this%Xc,1)
                Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                Xcw(j,dim+1) = this%Wc(j)
            end do

            call elevate_degree_A_5_9(t, this%knot, this%degree, Xcw, nc_new, knot_new, Xcw_new)

            allocate(Xc_new(1:nc_new,1:dim))
            allocate(Wc_new(1:nc_new))
            do j = 1, nc_new
                Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
            end do
            Wc_new(:) = Xcw_new(:,dim+1)

            call this%set(knot=knot_new, Xc=Xc_new, Wc=Wc_new)
            deallocate(Xcw, Xcw_new, Xc_new, Wc_new)

        else ! B-Spline

            dim = size(this%Xc,2)

            call elevate_degree_A_5_9(t, this%knot, this%degree, this%Xc, nc_new, knot_new, Xc_new)

            call this%set(knot=knot_new, Xc=Xc_new)
            deallocate(Xc_new)

        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine derivative_vector(this, res, Xt, dTgc, Tgc)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: res
        real(rk), intent(in), contiguous, optional :: Xt(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
        integer :: i

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            this%Xt = Xt
        elseif (present(res)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            allocate(this%Xt(res))
            this%Xt = [(real(i-1, rk) / real(res-1, rk), i=1, res)]
            ! else
            ! this%Xt = this%Xt
        end if

        ! Set number of geometry points
        this%ng = size(this%Xt,1)

        if (this%is_rational()) then ! NURBS
            call compute_dTgc(this%Xt, this%knot, this%degree, this%nc, this%ng, this%Wc, dTgc, Tgc)
        else ! B-Spline
            call compute_dTgc(this%Xt, this%knot, this%degree, this%nc, this%ng, dTgc, Tgc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine derivative_scalar(this, Xt, dTgc, Tgc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xt
        real(rk), allocatable, intent(out) :: dTgc(:)
        real(rk), allocatable, intent(out), optional :: Tgc(:)

        if (this%is_rational()) then ! NURBS
            call compute_dTgc(Xt, this%knot, this%degree, this%nc, this%Wc, dTgc, Tgc)
        else ! B-Spline
            call compute_dTgc(Xt, this%knot, this%degree, this%nc, dTgc, Tgc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine derivative2_vector(this, res, Xt, d2Tgc, dTgc, Tgc)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: res
        real(rk), intent(in), contiguous, optional :: Xt(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
        real(rk), allocatable, intent(out), optional :: dTgc(:,:)
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
        integer :: i

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            this%Xt = Xt
        elseif (present(res)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            allocate(this%Xt(res))
            this%Xt = [(real(i-1, rk) / real(res-1, rk), i=1, res)]
            ! else
            ! this%Xt = this%Xt
        end if

        ! Set number of geometry points
        this%ng = size(this%Xt,1)

        if (this%is_rational()) then ! NURBS
            call compute_d2Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, this%Wc, d2Tgc, dTgc, Tgc)
        else ! B-Spline
            call compute_d2Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, d2Tgc, dTgc, Tgc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine derivative2_scalar(this, Xt, d2Tgc, dTgc, Tgc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xt
        real(rk), allocatable, intent(out) :: d2Tgc(:)
        real(rk), allocatable, intent(out), optional :: dTgc(:)
        real(rk), allocatable, intent(out), optional :: Tgc(:)

        if (this%is_rational()) then ! NURBS
            call compute_d2Tgc(Xt, this%knot, this%degree, this%nc, this%Wc, d2Tgc, dTgc, Tgc)
        else ! B-Spline
            call compute_d2Tgc(Xt, this%knot, this%degree, this%nc, d2Tgc, dTgc, Tgc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine basis_vector(this, res, Xt, Tgc)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: res
        real(rk), intent(in), contiguous, optional :: Xt(:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        integer :: i

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            this%Xt = Xt
        elseif (present(res)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            allocate(this%Xt(res))
            this%Xt = [(real(i-1, rk) / real(res-1, rk), i=1, res)]
            ! else
            ! this%Xt = this%Xt
        end if

        ! Set number of geometry points
        this%ng = size(this%Xt,1)

        if (this%is_rational()) then ! NURBS
            Tgc = compute_Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, this%Wc)
        else ! B-Spline
            Tgc = compute_Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine basis_scalar(this, Xt, Tgc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xt
        real(rk), allocatable, intent(out) :: Tgc(:)

        if (this%is_rational()) then ! NURBS
            Tgc = compute_Tgc(Xt, this%knot, this%degree, this%nc, this%Wc)
        else ! B-Spline
            Tgc = compute_Tgc(Xt, this%knot, this%degree, this%nc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function is_rational(this) result(r)
        class(nurbs_curve), intent(in) :: this
        logical :: r

        r = .false.
        if (allocated(this%Wc)) then
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
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (allocated(this%elemConn_Xc_vis)) deallocate(this%elemConn_Xc_vis)
        this%elemConn_Xc_vis = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_elem_Xg_vis(this, elemConn)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (allocated(this%elemConn_Xg_vis)) deallocate(this%elemConn_Xg_vis)
        this%elemConn_Xg_vis = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_elem(this, elemConn)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (allocated(this%elemConn)) deallocate(this%elemConn)
        this%elemConn = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem_Xc_vis(this) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        elemConn = this%elemConn_Xc_vis
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem_Xg_vis(this) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        elemConn = this%elemConn_Xg_vis
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem(this) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        elemConn = this%elemConn
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine remove_knots(this,Xth,r)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xth(:)
        integer, intent(in), contiguous :: r(:)
        integer :: k, i, s, dim, j, nc_new, t
        real(rk), allocatable :: Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:)

        if (this%is_rational()) then ! NURBS

            do i = 1, size(Xth)
                k = findspan(this%nc-1,this%degree,Xth(i),this%knot)
                if (this%knot(k+1) == Xth(i)) then
                    s = compute_multiplicity(this%knot,Xth(i))
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

                call remove_knots_A_5_8(&
                    this%degree,&
                    this%knot,&
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
                    allocate(Xc_new(nc_new,dim))
                    allocate(Wc_new(nc_new))
                    do j = 1, nc_new
                        Xc_new(j,:) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                    end do
                    Wc_new(:) = Xcw_new(:,dim+1)

                    call this%set(knot=knot_new, Xc=Xc_new, Wc=Wc_new)
                    if (allocated(Xcw_new)) deallocate(Xcw_new)
                    if (allocated(Xc_new)) deallocate(Xc_new)
                    if (allocated(Wc_new)) deallocate(Wc_new)
                end if
            end do

        else ! B-Spline

            do i = 1, size(Xth)
                k = findspan(this%nc-1,this%degree,Xth(i),this%knot)
                if (this%knot(k+1) == Xth(i)) then
                    s = compute_multiplicity(this%knot,Xth(i))
                else
                    s = 0
                end if
                k = k + 1

                call remove_knots_A_5_8(&
                    this%degree,&
                    this%knot,&
                    this%Xc,&
                    Xth(i),&
                    k,&
                    s,&
                    r(i),&
                    t,&
                    knot_new,&
                    Xc_new)

                if (t == 0) then
                    ! no change
                else
                    call this%set(knot=knot_new, Xc=Xc_new)
                end if
            end do

        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_circle(this, center, radius)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: center(:)
        real(rk), intent(in) :: radius
        real(rk), allocatable :: Xc(:,:), Wc(:), knot(:)
        integer :: i

        ! Define control points for circle
        allocate(Xc(7, 3))
        Xc(1,:)= [ 1.0_rk,  0.0_rk,              0.0_rk]
        Xc(2,:)= [ 1.0_rk,  sqrt(3.0_rk),        0.0_rk]
        Xc(3,:)= [-0.5_rk,  sqrt(3.0_rk)/2.0_rk, 0.0_rk]
        Xc(4,:)= [-2.0_rk,  0.0_rk,              0.0_rk]
        Xc(5,:)= [-0.5_rk, -sqrt(3.0_rk)/2.0_rk, 0.0_rk]
        Xc(6,:)= [ 1.0_rk, -sqrt(3.0_rk),        0.0_rk]
        Xc(7,:)= [ 1.0_rk,  0.0_rk,              0.0_rk]

        ! Scale and translate the control points
        do i = 1, size(Xc, 1)
            Xc(i,:) = center + Xc(i,:) * radius
        end do

        ! Define weights for the control points
        Wc = [1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk]

        ! Define knot vector
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/3.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

        ! Set knot vector, control points, and weights
        call this%set(knot, Xc, Wc)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_C(this, center, radius)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: center(:)
        real(rk), intent(in) :: radius
        real(rk), allocatable :: Xc(:,:), Wc(:), knot(:)
        integer :: i

        ! Define control points for C-shape
        allocate(Xc(5, 3))
        Xc(1,:)= [ 1.0_rk,  0.0_rk,              0.0_rk]
        Xc(2,:)= [ 1.0_rk,  sqrt(3.0_rk),        0.0_rk]
        Xc(3,:)= [-0.5_rk,  sqrt(3.0_rk)/2.0_rk, 0.0_rk]
        Xc(4,:)= [-2.0_rk,  0.0_rk,              0.0_rk]
        Xc(5,:)= [-0.5_rk, -sqrt(3.0_rk)/2.0_rk, 0.0_rk]

        ! Scale and translate the control points
        do i = 1, size(Xc, 1)
            Xc(i,:) = center + Xc(i,:) * radius
        end do

        ! Define weights for the control points
        Wc = [1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk]

        ! Define knot vector
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/2.0_rk, 1.0_rk/2.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

        ! Set knot vector, control points, and weights
        call this%set(knot, Xc, Wc)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elem(this) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        call elemConn_Cn(this%nc, this%degree, unique(this%knot), this%get_multiplicity(),&
            elemConn)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine rotate_Xc(this, alpha, beta, theta)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: alpha, beta, theta
        integer :: i

        do i = 1, this%nc
            this%Xc(i, :) = matmul(rotation(alpha,beta,theta), this%Xc(i, :))
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine rotate_Xg(this, alpha, beta, theta)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: alpha, beta, theta
        integer :: i

        do i = 1, this%ng
            this%Xg(i, :) = matmul(rotation(alpha,beta,theta), this%Xg(i, :))
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine translate_Xc(this, vec)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: vec(:)
        integer :: i

        do i = 1, this%nc
            this%Xc(i, :) = this%Xc(i, :) + vec
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine translate_Xg(this, vec)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: vec(:)
        integer :: i

        do i = 1, this%nc
            this%Xg(i, :) = this%Xg(i, :) + vec
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine show(this, vtkfile_Xc, vtkfile_Xg)
        class(nurbs_curve), intent(inout) :: this
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
    pure subroutine set_half_circle(this, center, radius)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: center(:)
        real(rk), intent(in) :: radius
        real(rk), allocatable :: Xc(:,:), Wc(:), knot(:)
        integer :: i

        ! Define control points for half circle
        allocate(Xc(5, 3))
        Xc(1,:) = [ 0.5_rk, 0.0_rk, 0.0_rk]
        Xc(2,:) = [ 0.5_rk, 0.5_rk, 0.0_rk]
        Xc(3,:) = [ 0.0_rk, 0.5_rk, 0.0_rk]
        Xc(4,:) = [-0.5_rk, 0.5_rk, 0.0_rk]
        Xc(5,:) = [-0.5_rk, 0.0_rk, 0.0_rk]

        ! Scale and translate the control points
        do i = 1, size(Xc, 1)
            Xc(i,:) = center + Xc(i,:) * radius
        end do

        ! Define weights for the control points
        Wc = [1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk]

        ! Define knot vector
        knot = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/2.0_rk, &
            1.0_rk/2.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]

        ! Set knot vector, control points, and weights
        call this%set(knot, Xc, Wc)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine nearest_point(this, point_Xg, nearest_Xg, nearest_Xt, id)
        class(nurbs_curve), intent(in) :: this
        real(rk), intent(in) :: point_Xg(:)
        real(rk), intent(out), allocatable, optional :: nearest_Xg(:)
        real(rk), intent(out), optional :: nearest_Xt
        integer, intent(out), optional :: id
        integer :: id_
        real(rk), allocatable :: distances(:)

        allocate(distances(this%ng))
        distances = nearest_point_help_1d(this%ng, this%Xg, point_Xg)

        id_ = minloc(distances, dim=1)
        if (present(id)) id = id_
        if (present(nearest_Xg)) nearest_Xg = this%Xg(id_,:)
        if (present(nearest_Xt)) nearest_Xt = this%Xt(id_)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine nearest_point2(this, point_Xg, tol, maxit, nearest_Xt, nearest_Xg)

        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: point_Xg(:)
        real(rk), intent(in) :: tol
        integer, intent(in) :: maxit
        real(rk), intent(out) :: nearest_Xt
        real(rk), allocatable, intent(out), optional :: nearest_Xg(:)
        real(rk):: xk, obj, grad, hess, dk, alphak, tau, beta, lower_bounds, upper_bounds
        real(rk), allocatable :: Xg(:), Tgc(:), dTgc(:), d2Tgc(:), distances(:)
        integer :: k, l
        logical :: convergenz
        type(nurbs_curve) :: copy_this

        k = 0

        ! lower and upper bounds
        lower_bounds = minval(this%knot)
        upper_bounds = maxval(this%knot)

        ! guess initial point
        copy_this = this
        call this%create(10)
        allocate(distances(copy_this%ng))
        distances = nearest_point_help_1d(copy_this%ng, copy_this%Xg, point_Xg)
        xk = copy_this%Xt(minloc(distances, dim=1))
        call copy_this%finalize()

        ! Check if xk is within the knot vector range
        if (xk < minval(this%knot)) then
            xk = minval(this%knot)
        else if (xk > maxval(this%knot)) then
            xk = maxval(this%knot)
        end if

        convergenz = .false.

        allocate(Xg(size(this%Xc,2)))
        ! allocate(dTgc(size(this%Xc,1)))
        ! allocate(d2Tgc(size(this%Xc,1)))

        do while (.not. convergenz .and. k < maxit)

            ! objection, gradient and hessian
            Xg = this%cmp_Xg(xk)
            call this%derivative2(Xt=xk, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc) ! Tgc is not needed

            obj  = norm2(Xg - point_Xg) + 0.001_rk ! add a small number to avoid division by zero
            grad = dot_product((Xg-point_Xg)/obj, matmul(dTgc,this%Xc))
            hess = dot_product(matmul(dTgc,this%Xc) - (Xg-point_Xg)/obj*grad, matmul(dTgc,this%Xc))/obj&
                + dot_product((Xg-point_Xg)/obj, matmul(d2Tgc,this%Xc))

            ! debug
            print '(i3,1x,e20.10,1x,e20.10)', k, xk, abs(grad)

            if (abs(grad) <= tol) then
                convergenz = .true.
                nearest_Xt = xk
                if (present(nearest_Xg)) nearest_Xg = this%cmp_Xg(nearest_Xt)
            else
                dk = - grad / hess

                ! Backtracking-Armijo Line Search
                alphak = 1.0_rk
                tau = 0.5_rk     ! 0 < tau  < 1
                beta = 1.0e-4_rk ! 0 < beta < 1
                l = 0
                do while (.not. norm2(this%cmp_Xg(xk + alphak*dk) - point_Xg)  <= obj + alphak*beta*grad*dk .and. l<50)
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

end module forcad_nurbs_curve

!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Xg_nurbs_1d(Xt, knot, degree, nc, ng, Xc, Wc) result(Xg)
    use forcad_utils, only: rk, basis_bspline

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    integer, intent(in) :: ng
    real(rk), intent(in), contiguous :: Xc(:,:)
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable :: Xg(:,:)
    real(rk), allocatable :: Tgc(:)
    integer :: i

    allocate(Xg(ng, size(Xc,2)))
    allocate(Tgc(nc))
    !$OMP PARALLEL DO PRIVATE(Tgc)
    do i = 1, ng
        Tgc = basis_bspline(Xt(i), knot, nc, degree)
        Tgc = Tgc*(Wc/(dot_product(Tgc,Wc)))
        Xg(i,:) = matmul(Tgc, Xc)
    end do
    !$OMP END PARALLEL DO
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Xg_nurbs_1d_1point(Xt, knot, degree, nc, Xc, Wc) result(Xg)
    use forcad_utils, only: rk, basis_bspline

    implicit none
    real(rk), intent(in) :: Xt
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    real(rk), intent(in), contiguous :: Xc(:,:)
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable :: Xg(:)
    real(rk), allocatable :: Tgc(:)

    allocate(Xg(size(Xc)))
    allocate(Tgc(nc))
    Tgc = basis_bspline(Xt, knot, nc, degree)
    Tgc = Tgc*(Wc/(dot_product(Tgc,Wc)))
    Xg = matmul(Tgc, Xc)
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Xg_bspline_1d(Xt, knot, degree, nc, ng, Xc) result(Xg)
    use forcad_utils, only: rk, basis_bspline

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    integer, intent(in) :: ng
    real(rk), intent(in), contiguous :: Xc(:,:)
    real(rk), allocatable :: Xg(:,:)
    integer :: i

    allocate(Xg(ng, size(Xc,2)))
    !$OMP PARALLEL DO
    do i = 1, ng
        Xg(i,:) = matmul(basis_bspline(Xt(i), knot, nc, degree), Xc)
    end do
    !$OMP END PARALLEL DO
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Xg_bspline_1d_1point(Xt, knot, degree, nc, Xc) result(Xg)
    use forcad_utils, only: rk, basis_bspline

    implicit none
    real(rk), intent(in) :: Xt
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    real(rk), intent(in), contiguous :: Xc(:,:)
    real(rk), allocatable :: Xg(:)

    allocate(Xg(size(Xc)))
    Xg = matmul(basis_bspline(Xt, knot, nc, degree), Xc)
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_dTgc_nurbs_1d_vector(Xt, knot, degree, nc, ng, Wc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_der

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    integer, intent(in) :: ng
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable, intent(out) :: dTgc(:,:)
    real(rk), allocatable, intent(out) :: Tgc(:,:)
    real(rk), allocatable :: dBi(:), Bi(:)
    integer :: i

    allocate(dTgc(ng, nc), Tgc(ng, nc), dBi(nc), Bi(nc))
    do i = 1, size(Xt)
        call basis_bspline_der(Xt(i), knot, nc, degree, dBi, Bi)
        Tgc(i,:) = Bi*(Wc/(dot_product(Bi,Wc)))
        dTgc(i,:) = ( dBi*Wc - Tgc(i,:)*dot_product(dBi,Wc) ) / dot_product(Bi,Wc)
    end do
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_dTgc_nurbs_1d_scalar(Xt, knot, degree, nc, Wc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_der

    implicit none
    real(rk), intent(in) :: Xt
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable, intent(out) :: dTgc(:)
    real(rk), allocatable, intent(out) :: Tgc(:)
    real(rk), allocatable :: dBi(:), Bi(:)

    allocate(dTgc(nc), Tgc(nc))
    call basis_bspline_der(Xt, knot, nc, degree, dBi, Bi)
    Tgc = Bi*(Wc/(dot_product(Bi,Wc)))
    dTgc = ( dBi*Wc - Tgc*dot_product(dBi,Wc) ) / dot_product(Bi,Wc)
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_dTgc_bspline_1d_vector(Xt, knot, degree, nc, ng, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_der

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    integer, intent(in) :: ng
    real(rk), allocatable, intent(out) :: dTgc(:,:)
    real(rk), allocatable, intent(out) :: Tgc(:,:)
    real(rk), allocatable :: dBi(:), Bi(:)
    integer :: i

    allocate(dTgc(ng, nc), Tgc(ng, nc), dBi(nc), Bi(nc))
    do i = 1, size(Xt)
        call basis_bspline_der(Xt(i), knot, nc, degree, dBi, Bi)
        Tgc(i,:) = Bi
        dTgc(i,:) = dBi
    end do
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_dTgc_bspline_1d_scalar(Xt, knot, degree, nc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_der

    implicit none
    real(rk), intent(in) :: Xt
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    real(rk), allocatable, intent(out) :: dTgc(:)
    real(rk), allocatable, intent(out) :: Tgc(:)

    allocate(dTgc(nc), Tgc(nc))
    call basis_bspline_der(Xt, knot, nc, degree, dTgc, Tgc)
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_d2Tgc_nurbs_1d_vector(Xt, knot, degree, nc, ng, Wc, d2Tgc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_2der

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    integer, intent(in) :: ng
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable, intent(out) :: d2Tgc(:,:)
    real(rk), allocatable, intent(out) :: dTgc(:,:)
    real(rk), allocatable, intent(out) :: Tgc(:,:)
    real(rk), allocatable :: d2Bi(:), dBi(:), Tgci(:), dTgci(:), Bi(:)
    integer :: i

    allocate(d2Tgc(ng, nc), dTgc(ng, nc), Tgc(ng, nc), d2Bi(nc), dTgci(nc), dBi(nc), Tgci(nc), Bi(nc))

    do i = 1, size(Xt)
        call basis_bspline_2der(Xt(i), knot, nc, degree, d2Bi, dBi, Bi)
        Tgci = Bi*(Wc/(dot_product(Bi,Wc)))
        Tgc(i,:) = Tgci

        dTgci = ( dBi*Wc - Tgci*dot_product(dBi,Wc) ) / dot_product(Bi,Wc)
        dTgc(i,:) = dTgci

        d2Tgc(i,:) = (d2Bi*Wc - 2.0_rk*dTgci*dot_product(dBi,Wc) - Tgci*dot_product(d2Bi,Wc)) / dot_product(Bi,Wc)
    end do
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_d2Tgc_nurbs_1d_scalar(Xt, knot, degree, nc, Wc, d2Tgc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_2der

    implicit none
    real(rk), intent(in) :: Xt
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable, intent(out) :: d2Tgc(:)
    real(rk), allocatable, intent(out) :: dTgc(:)
    real(rk), allocatable, intent(out) :: Tgc(:)
    real(rk), allocatable :: d2Bi(:), dBi(:), Bi(:)

    allocate(d2Tgc(nc), dTgc(nc), Tgc(nc), d2Bi(nc), dBi(nc), Bi(nc))

    call basis_bspline_2der(Xt, knot, nc, degree, d2Bi, dBi, Bi)
    Tgc = Bi*(Wc/(dot_product(Bi,Wc)))
    dTgc = ( dBi*Wc - Tgc*dot_product(dBi,Wc) ) / dot_product(Bi,Wc)
    d2Tgc = (d2Bi*Wc - 2.0_rk*dTgc*dot_product(dBi,Wc) - Tgc*dot_product(d2Bi,Wc)) / dot_product(Bi,Wc)
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_d2Tgc_bspline_1d_vector(Xt, knot, degree, nc, ng, d2Tgc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_2der

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    integer, intent(in) :: ng
    real(rk), allocatable, intent(out) :: d2Tgc(:,:)
    real(rk), allocatable, intent(out) :: dTgc(:,:)
    real(rk), allocatable, intent(out) :: Tgc(:,:)
    real(rk), allocatable :: d2Tgci(:), dTgci(:), Tgci(:)
    integer :: i

    allocate(d2Tgc(ng, nc), dTgc(ng, nc), Tgc(ng, nc), dTgci(nc), Tgci(nc), d2Tgci(nc))
    do i = 1, size(Xt)
        call basis_bspline_2der(Xt(i), knot, nc, degree, d2Tgci, dTgci, Tgci)
        Tgc(i,:) = Tgci
        dTgc(i,:) = dTgci
        d2Tgc(i,:) = d2Tgci
    end do
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure subroutine compute_d2Tgc_bspline_1d_scalar(Xt, knot, degree, nc, d2Tgc, dTgc, Tgc)
    use forcad_utils, only: rk, basis_bspline_2der

    implicit none
    real(rk), intent(in) :: Xt
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    real(rk), allocatable, intent(out) :: d2Tgc(:)
    real(rk), allocatable, intent(out) :: dTgc(:)
    real(rk), allocatable, intent(out) :: Tgc(:)

    allocate(d2Tgc(nc), dTgc(nc), Tgc(nc))
    call basis_bspline_2der(Xt, knot, nc, degree, d2Tgc, dTgc, Tgc)
end subroutine
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Tgc_nurbs_1d_vector(Xt, knot, degree, nc, ng, Wc) result(Tgc)
    use forcad_utils, only: rk, basis_bspline

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    integer, intent(in) :: ng
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable :: Tgc(:,:)
    real(rk), allocatable :: Tgci(:)
    integer :: i

    allocate(Tgc(ng, nc), Tgci(nc))
    !$OMP PARALLEL DO PRIVATE(Tgci)
    do i = 1, size(Xt,1)
        Tgci = basis_bspline(Xt(i), knot, nc, degree)
        Tgc(i,:) = Tgci*(Wc/(dot_product(Tgci,Wc)))
    end do
    !$OMP END PARALLEL DO
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Tgc_nurbs_1d_scalar(Xt, knot, degree, nc, Wc) result(Tgc)
    use forcad_utils, only: rk, basis_bspline

    implicit none
    real(rk), intent(in) :: Xt
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    real(rk), intent(in), contiguous :: Wc(:)
    real(rk), allocatable :: Tgc(:)

    allocate(Tgc(nc))
    Tgc = basis_bspline(Xt, knot, nc, degree)
    Tgc = Tgc*(Wc/(dot_product(Tgc,Wc)))
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Tgc_bspline_1d_vector(Xt, knot, degree, nc, ng) result(Tgc)
    use forcad_utils, only: rk, basis_bspline

    implicit none
    real(rk), intent(in), contiguous :: Xt(:)
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    integer, intent(in) :: ng
    real(rk), allocatable :: Tgc(:,:)
    integer :: i

    allocate(Tgc(ng, nc))
    !$OMP PARALLEL DO
    do i = 1, size(Xt,1)
        Tgc(i,:) = basis_bspline(Xt(i), knot, nc, degree)
    end do
    !$OMP END PARALLEL DO
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function compute_Tgc_bspline_1d_scalar(Xt, knot, degree, nc) result(Tgc)
    use forcad_utils, only: rk, basis_bspline

    implicit none
    real(rk), intent(in) :: Xt
    real(rk), intent(in), contiguous :: knot(:)
    integer, intent(in) :: degree
    integer, intent(in) :: nc
    real(rk), allocatable :: Tgc(:)

    allocate(Tgc(nc))
    Tgc = basis_bspline(Xt, knot, nc, degree)
end function
!===============================================================================


!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
impure function nearest_point_help_1d(ng, Xg, point_Xg) result(distances)
    use forcad_utils, only: rk

    implicit none
    integer, intent(in) :: ng
    real(rk), intent(in), contiguous :: Xg(:,:)
    real(rk), intent(in), contiguous :: point_Xg(:)
    real(rk), allocatable :: distances(:)
    integer :: i

    allocate(distances(ng))
    !$OMP PARALLEL DO
    do i = 1, ng
        distances(i) = norm2(Xg(i,:) - point_Xg)
    end do
    !$OMP END PARALLEL DO

end function
!===============================================================================
