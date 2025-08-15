!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> This module defines the 'nurbs_curve' type for representing a Non-Uniform Rational B-Spline (NURBS) curve.
module forcad_nurbs_curve

    use forcad_kinds, only: rk
    use forcad_utils, only: basis_bspline, elemConn_C0, ndgrid, compute_multiplicity, compute_knot_vector, basis_bspline_der,&
        insert_knot_A_5_1, findspan, elevate_degree_A_5_9, remove_knots_A_5_8, &
        elemConn_Cn, unique, rotation, dyad, gauss_leg, export_vtk_legacy, basis_bspline_2der
    use fordebug, only: debug

    implicit none

    private
    public nurbs_curve, compute_Tgc, compute_dTgc

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

        type(debug) :: err !! 101: size mismatch (weights vs control points), 102: missing control points, 103: missing knot vector, 104: missing geometry points, 105: missing weights, 106: lsq fit underdetermined
    contains
        procedure, private :: set1                  !!> Set knot vector, control points and weights for the NURBS curve object
        procedure, private :: set1a
        procedure, private :: set2                  !!> Set NURBS curve using nodes of parameter space, degree, continuity, control points and weights
        procedure, private :: set3                  !!> Set Bezier or Rational Bezier curve using control points and weights
        procedure, private :: set4                  !!> Set NURBS curve using degree, number of control points, control points and weights
        generic :: set => set1, set1a, set2, set3, set4 !!> Set NURBS curve
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
        procedure :: cmp_elem_Xth          !!> Generate connectivity for parameter points
        procedure :: cmp_elem              !!> Generate IGA element connectivity
        procedure :: get_elem_Xc_vis       !!> Get connectivity for control points
        procedure :: get_elem_Xg_vis       !!> Get connectivity for geometry points
        procedure :: get_elem              !!> Get IGA element connectivity
        procedure :: set_elem_Xc_vis       !!> Set connectivity for control points
        procedure :: set_elem_Xg_vis       !!> Set connectivity for geometry points
        procedure :: set_elem              !!> Set IGA element connectivity
        procedure :: export_Xc             !!> Export control points to VTK file
        procedure :: export_Xg             !!> Export geometry points to VTK file
        procedure :: export_Xth            !!> Export parameter space to VTK file
        procedure :: export_iges           !!> Export the NURBS curve to an IGES file
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
        procedure :: nearest_point2        !!> Find the nearest point on the NURBS curve (Minimization - Newtons method)
        procedure :: ansatz                !!> Compute the shape functions, derivative of shape functions and dL
        procedure :: cmp_length            !!> Compute the length of the NURBS curve
        procedure :: lsq_fit_bspline       !!> Fit B-spline curve to structured data points using least squares

        ! Shapes
        procedure :: set_circle            !!> Set a circle
        procedure :: set_half_circle       !!> Set a half circle
        procedure :: set_C                 !!> Set a C-shape
    end type
    !===============================================================================

    interface compute_Xg
        module procedure compute_Xg_nurbs_1d
        module procedure compute_Xg_bspline_1d
        module procedure compute_Xg_nurbs_1d_1point
        module procedure compute_Xg_bspline_1d_1point
    end interface

    interface compute_Tgc
        module procedure compute_Tgc_nurbs_1d_vector
        module procedure compute_Tgc_bspline_1d_vector
        module procedure compute_Tgc_nurbs_1d_scalar
        module procedure compute_Tgc_bspline_1d_scalar
    end interface

    interface compute_dTgc
        module procedure compute_dTgc_nurbs_1d_vector
        module procedure compute_dTgc_bspline_1d_vector
        module procedure compute_dTgc_nurbs_1d_scalar
        module procedure compute_dTgc_bspline_1d_scalar
    end interface

    interface compute_d2Tgc
        module procedure compute_d2Tgc_nurbs_1d_vector
        module procedure compute_d2Tgc_bspline_1d_vector
        module procedure compute_d2Tgc_nurbs_1d_scalar
        module procedure compute_d2Tgc_bspline_1d_scalar
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

        if (.not. this%err%ok) return

        if (allocated(this%knot)) then
            if (size(this%knot) /= size(knot)) deallocate(this%knot)
        end if

        if (allocated(this%Xc)) then
            if (size(this%Xc, 1) /= size(Xc, 1) .or. size(this%Xc, 2) /= size(Xc, 2)) deallocate(this%Xc)
        end if

        this%knot = knot
        call this%cmp_degree()
        this%Xc = Xc
        this%nc = size(this%Xc, 1)
        if (present(Wc)) then
            if (size(Wc) /= this%nc) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_curve',&
                    message    = 'Weights length mismatch: size(Wc) must equal number of control points.',&
                    location   = 'set1',&
                    suggestion = 'Provide Wc with size(Wc) == size(Xc,1).')
                return
            else
                if (allocated(this%Wc)) then
                    if (size(this%Wc) /= size(Wc)) deallocate(this%Wc)
                end if
                this%Wc = Wc
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set knot vector, control points and weights for the NURBS curve object.
    pure subroutine set1a(this, knot, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: knot(:)
        real(rk), intent(in), contiguous :: Xc(:)
        real(rk), intent(in), contiguous, optional :: Wc(:)

        if (.not. this%err%ok) return

        if (allocated(this%knot)) then
            if (size(this%knot) /= size(knot)) deallocate(this%knot)
        end if

        if (allocated(this%Xc)) then
            if (size(this%Xc, 1) /= size(Xc) .or. size(this%Xc, 2) /= 3) then
                deallocate(this%Xc)
                allocate(this%Xc(size(Xc), 3), source = 0.0_rk)
            end if
        else
            allocate(this%Xc(size(Xc), 3), source = 0.0_rk)
        end if

        this%knot = knot
        call this%cmp_degree()
        this%Xc(:,1) = Xc
        this%nc = size(this%Xc, 1)
        if (present(Wc)) then
            if (size(Wc) /= this%nc) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_curve',&
                    message    = 'Weights length mismatch: size(Wc) must equal number of control points.',&
                    location   = 'set1a',&
                    suggestion = 'Provide Wc with size(Wc) == size(Xc,1).')
                return
            else
                if (allocated(this%Wc)) then
                    if (size(this%Wc) /= size(Wc)) deallocate(this%Wc)
                end if
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
        real(rk), intent(in), contiguous, optional :: Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)

        if (.not. this%err%ok) return

        if (allocated(this%knot)) deallocate(this%knot)
        this%knot   = compute_knot_vector(Xth_dir, degree, continuity)
        this%degree = degree
        call this%cmp_nc()

        if (present(Xc)) then
            if (size(Xc,1) /= this%nc) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_curve', &
                    message    = 'Control points size mismatch in set2',&
                    location   = 'set2', &
                    suggestion = 'size(Xc,1) must equal computed nc.')
                return
            end if
        end if

        if (present(Wc)) then
            if (size(Wc) /= this%nc) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_curve', &
                    message    = 'Weights size mismatch in set2',&
                    location   = 'set2', &
                    suggestion = 'size(Wc) must equal computed nc.')
                return
            end if
        end if

        if (present(Xc)) then
            if (allocated(this%Xc)) then
                if (size(this%Xc, 1) /= size(Xc, 1) .or. size(this%Xc, 2) /= size(Xc, 2)) deallocate(this%Xc)
            end if
            this%Xc = Xc
        end if

        if (present(Wc)) then
            if (allocated(this%Wc)) then
                if (size(this%Wc) /= size(Wc)) deallocate(this%Wc)
            end if
            this%Wc = Wc
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

        if (.not. this%err%ok) return

        if (allocated(this%Xc)) then
            if (size(this%Xc, 1) /= size(Xc, 1) .or. size(this%Xc, 2) /= size(Xc, 2)) deallocate(this%Xc)
        end if

        this%Xc = Xc
        this%nc = size(this%Xc, 1)

        if (allocated(this%knot)) then
            if (size(this%knot) /= 2*this%nc) then
                deallocate(this%knot)
                allocate(this%knot(2*this%nc))
            end if
        else
            allocate(this%knot(2*this%nc))
        end if
        this%knot(1:this%nc) = 0.0_rk
        this%knot(this%nc+1:2*this%nc) = 1.0_rk

        call this%cmp_degree()
        if (present(Wc)) then
            if (size(Wc) /= this%nc) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_curve',&
                    message    = 'Weights length mismatch: size(Wc) must equal number of control points.',&
                    location   = 'set3',&
                    suggestion = 'Provide Wc with size(Wc) == size(Xc,1).')
                return
            else
                if (allocated(this%Wc)) then
                    if (size(this%Wc) /= size(Wc)) deallocate(this%Wc)
                end if
                this%Wc = Wc
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set4(this, degree, nc, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        integer :: m, i

        if (.not. this%err%ok) return

        if (allocated(this%Xc)) then
            if (size(this%Xc, 1) /= size(Xc, 1) .or. size(this%Xc, 2) /= size(Xc, 2)) deallocate(this%Xc)
        end if

        this%Xc = Xc
        this%nc = nc
        this%degree = degree

        ! Size of knot vectors
        m = nc + degree + 1

        if (allocated(this%knot)) then
            if (size(this%knot) /= m) then
                deallocate(this%knot)
                allocate(this%knot(m))
            end if
        else
            allocate(this%knot(m))
        end if
        this%knot(1:degree+1) = 0.0_rk
        this%knot(degree+2:m-degree-1) = [(real(i, rk)/(m-2*degree-1), i=1, m-2*degree-2)]
        this%knot(m-degree:m) = 1.0_rk

        if (present(Wc)) then
            if (size(Wc) /= nc) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_curve',&
                    message    = 'Weights length mismatch: size(Wc) must equal number of control points.',&
                    location   = 'set4',&
                    suggestion = 'Provide Wc with size(Wc) == nc.')
                return
            else
                if (allocated(this%Wc)) then
                    if (size(this%Wc) /= size(Wc)) deallocate(this%Wc)
                end if
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

        if (.not. this%err%ok) return

        ! check
        if (.not.allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = 'forcad_nurbs_curve',&
                message    = 'Control points are not set.',&
                location   = 'create',&
                suggestion = 'Call set(...) first before create().')
            return
        end if

        if (.not.allocated(this%knot)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = 'forcad_nurbs_curve',&
                message    = 'Knot vector is not set.',&
                location   = 'create',&
                suggestion = 'Call set(...) first before create().')
            return
        end if

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) then
                if (size(this%Xt) /= size(Xt)) deallocate(this%Xt)
            end if
            this%Xt = Xt
        elseif (present(res)) then
            if (allocated(this%Xt)) then
                if (size(this%Xt) /= res) then
                    deallocate(this%Xt)
                    allocate(this%Xt(res))
                end if
            else
                allocate(this%Xt(res))
            end if
            this%Xt = [(this%knot(1)+(this%knot(size(this%knot))-this%knot(1))*real(i-1,rk)/real(res-1,rk), i=1, res)]
            ! else
            ! this%Xt = this%Xt
        end if

        ! Set number of geometry points
        this%ng = size(this%Xt)

        ! Allocate memory for geometry points
        if (allocated(this%Xg)) then
            if (size(this%Xg, 1) /= this%ng .or. size(this%Xg, 2) /= size(this%Xc, 2)) deallocate(this%Xg)
        end if

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

        ng = this%ng
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine cmp_degree(this)
        class(nurbs_curve), intent(inout) :: this
        integer, allocatable :: m(:)

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

        degree = this%degree
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_knot_all(this) result(knot)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: knot(:)

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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
    pure function cmp_elem_Xth(this, p) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), optional :: p

        if (.not. this%err%ok) return

        if (present(p)) then
            elemConn = elemConn_C0(size(unique(this%knot)), p)
        else
            elemConn = elemConn_C0(size(unique(this%knot)), 1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xc(this, filename, point_data, field_names, encoding)
        class(nurbs_curve), intent(inout) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), optional :: point_data(:,:)
        character(len=*), intent(in), optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        ! check
        if (.not.allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = 'forcad_nurbs_curve',&
                message    = 'Control points are not set.',&
                location   = 'export_Xc',&
                suggestion = 'Call set(...) first before exporting.')
            return
        end if

        if (.not.allocated(this%elemConn_Xc_vis)) then
            elemConn = this%cmp_elem_Xc_vis()
        else
            elemConn = this%elemConn_Xc_vis
        end if

        call export_vtk_legacy(filename=filename, points=this%Xc, elemConn=elemConn, vtkCellType=3, &
                               point_data=point_data, field_names=field_names, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xg(this, filename, point_data, field_names, encoding)
        class(nurbs_curve), intent(inout) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), optional :: point_data(:,:)
        character(len=*), intent(in), optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        ! check
        if (.not.allocated(this%Xg)) then
            call this%err%set(&
                code       = 104,&
                severity   = 1,&
                category   = 'forcad_nurbs_curve',&
                message    = 'Geometry points are not set.',&
                location   = 'export_Xg',&
                suggestion = 'Generate Xg by calling create(...) before exporting.')
            return
        end if

        if (.not.allocated(this%elemConn_Xg_vis)) then
            elemConn = this%cmp_elem_Xg_vis()
        else
            elemConn = this%elemConn_Xg_vis
        end if

        call export_vtk_legacy(filename=filename, points=this%Xg, elemConn=elemConn, vtkCellType=3, &
                               point_data=point_data, field_names=field_names, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xth(this, filename, point_data, field_names, encoding)
        class(nurbs_curve), intent(in) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), optional :: point_data(:,:)
        character(len=*), intent(in), optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Xth(:,:), Xth1(:), Xth2(:), Xth3(:)

        if (.not. this%err%ok) return

        elemConn = this%cmp_elem_Xth()
        Xth1 = unique(this%knot)
        Xth2 = [0.0_rk]
        Xth3 = [0.0_rk]
        call ndgrid(Xth1, Xth2, Xth3, Xth)

        call export_vtk_legacy(filename=filename, points=Xth, elemConn=elemConn, vtkCellType=3, &
                               point_data=point_data, field_names=field_names, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_iges(this, filename)
        use forIGES, only: Gsection_t, Dentry_t, entity126_t, DElist_t, PElist_t,&
                           makeSsection, makeGsection, makeDPsections, writeIGESfile, wp

        class(nurbs_curve), intent(inout) :: this
        character(len=*), intent(in)      :: filename

        type(Gsection_t)  :: G
        type(Dentry_t)    :: D
        type(entity126_t) :: curve126
        type(DElist_t)    :: Dlist
        type(PElist_t)    :: Plist
        character(80), allocatable :: Ssection(:), Gsection(:), Dsection(:), Psection(:), Ssec_out(:)
        integer :: i, K, M, N, prop3
        real(wp) :: T(-this%degree:1+this%degree), X(0:this%degree), Y(0:this%degree), Z(0:this%degree), W(0:this%degree), V(0:1)
        real(wp) :: XNORM, YNORM, ZNORM

        if (.not. this%err%ok) return

        ! Parameters for IGES knot vector
        K = this%degree
        M = this%degree
        N = 1 + K - M

        ! Copy your knot vector to IGES indexing
        do i = -M, N + K
            T(i) = real(this%knot(i + M + 1), kind=wp)
        end do

        ! Copy control points
        if (this%is_rational()) then
            do i = 0, K
                X(i) = real(this%Xc(i+1, 1), kind=wp)
                Y(i) = real(this%Xc(i+1, 2), kind=wp)
                Z(i) = real(this%Xc(i+1, 3), kind=wp)
                W(i) = real(this%Wc(i+1), kind=wp)
            end do
            prop3 = 1
        else
            do i = 0, K
                X(i) = real(this%Xc(i+1, 1), kind=wp)
                Y(i) = real(this%Xc(i+1, 2), kind=wp)
                Z(i) = real(this%Xc(i+1, 3), kind=wp)
                W(i) = real(1.0_rk, kind=wp)
            end do
            prop3 = 0
        end if

        XNORM = real(0.0_rk, kind=wp)
        YNORM = real(0.0_rk, kind=wp)
        ZNORM = real(0.0_rk, kind=wp)

        V = real([minval(this%knot), maxval(this%knot)], kind=wp)

        ! Initialize IGES entity126 (Rational B-spline Curve)
        call curve126%init(&
            DEP   = 1,&
            form  = 0,&
            K     = K,&
            M     = M,&
            PROP1 = 0,&
            PROP2 = 0,&
            PROP3 = prop3,&
            PROP4 = 0,&
            T     = T,&
            W     = W,&
            X     = X,&
            Y     = Y,&
            Z     = Z,&
            V     = V,&
            XNORM = XNORM,&
            YNORM = YNORM,&
            ZNORM = ZNORM)

        ! Directory entry
        call D%init(entity_type=126, param_data=1, transformation_matrix=0, form_number=0)

        ! Entity and directory lists
        call Dlist%init()
        call Plist%init()
        call Dlist%append(D)
        call Plist%append(curve126)

        ! Global section
        call G%init(filename=filename)

        ! S-section description
        allocate(Ssection(1))
        Ssection(1) = 'ForCAD'

        ! Create IGES sections
        call makeSsection(Ssection, Ssec_out)
        call makeGsection(G, Gsection)
        call makeDPsections(Dlist, Plist, Dsection, Psection)

        ! Write IGES file
        call writeIGESfile(filename, Ssec_out, Gsection, Dsection, Psection)

        ! Cleanup
        call Dlist%delete()
        call Plist%delete()
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

        if (.not. this%err%ok) return

        if (allocated(this%Xc)) then
            this%Xc(num,dir) = X
            if (allocated(this%Wc)) then
                call this%set(knot=this%get_knot(), Xc=this%get_Xc(), Wc=this%get_Wc())
            else
                call this%set(knot=this%get_knot(), Xc=this%get_Xc())
            end if
        else
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = 'forcad_nurbs_curve',&
                message    = 'Control points are not set.',&
                location   = 'modify_Xc',&
                suggestion = 'Call set(...) before modifying it.')
            return
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

        if (.not. this%err%ok) return

        if (allocated(this%Wc)) then
            this%Wc(num) = W
            if (allocated(this%knot)) then
                call this%set(knot=this%get_knot(), Xc=this%get_Xc(), Wc=this%get_Wc())
            else
                call this%set(Xc=this%get_Xc(), Wc=this%get_Wc())
            end if
        else
            call this%err%set(&
                code       = 105,&
                severity   = 1,&
                category   = 'forcad_nurbs_curve',&
                message    = 'Weights are not set.',&
                location   = 'modify_Wc',&
                suggestion = 'Pass Wc when calling set(...), before modifying weights.')
            return
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_multiplicity(this) result(m)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: m(:)

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

        this%nc = sum(compute_multiplicity(this%knot)) - this%degree - 1
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_nc(this) result(nc)
        class(nurbs_curve), intent(in) :: this
        integer :: nc

        if (.not. this%err%ok) return

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
        integer :: k, i, s, d, j, n_new
        real(rk), allocatable :: Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:)

        if (.not. this%err%ok) return

        if (this%is_rational()) then ! NURBS

            do i = 1, size(Xth)
                k = findspan(this%nc-1,this%degree,Xth(i),this%knot)
                ! if (this%knot(k+1) == Xth(i)) then
                if (abs(this%knot(k+1) - Xth(i)) < 2.0_rk*epsilon(0.0_rk)) then
                    s = compute_multiplicity(this%knot,Xth(i))
                else
                    s = 0
                end if

                d = size(this%Xc,2)
                allocate(Xcw(size(this%Xc,1),d+1))

                do j = 1, size(this%Xc,1)
                    Xcw(j,1:d) = this%Xc(j,1:d)*this%Wc(j)
                end do
                Xcw(:,d+1) = this%Wc(:)

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

                allocate(Xc_new(1:n_new+1,1:d))
                allocate(Wc_new(1:n_new+1))
                do j = 1, n_new+1
                    Xc_new(j,1:d) = Xcw_new(j-1,1:d)/Xcw_new(j-1,d+1)
                    Wc_new(j) = Xcw_new(j-1,d+1)
                end do

                call this%set(knot=knot_new, Xc=Xc_new, Wc=Wc_new)
                deallocate(Xcw, Xcw_new, Xc_new, Wc_new)
            end do

        else ! B-Spline

            do i = 1, size(Xth)
                k = findspan(this%nc-1,this%degree,Xth(i),this%knot)
                ! if (this%knot(k+1) == Xth(i)) then
                if (abs(this%knot(k+1) - Xth(i)) < 2.0_rk*epsilon(0.0_rk)) then
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
        integer :: d, j, nc_new

        if (.not. this%err%ok) return

        if (this%is_rational()) then ! NURBS

            d = size(this%Xc,2)
            allocate(Xcw(size(this%Xc,1),d+1))
            do j = 1, size(this%Xc,1)
                Xcw(j,1:d) = this%Xc(j,1:d)*this%Wc(j)
                Xcw(j,d+1) = this%Wc(j)
            end do

            call elevate_degree_A_5_9(t, this%knot, this%degree, Xcw, nc_new, knot_new, Xcw_new)

            allocate(Xc_new(1:nc_new,1:d))
            allocate(Wc_new(1:nc_new))
            do j = 1, nc_new
                Xc_new(j,1:d) = Xcw_new(j,1:d)/Xcw_new(j,d+1)
            end do
            Wc_new(:) = Xcw_new(:,d+1)

            call this%set(knot=knot_new, Xc=Xc_new, Wc=Wc_new)
            deallocate(Xcw, Xcw_new, Xc_new, Wc_new)

        else ! B-Spline

            d = size(this%Xc,2)

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

        if (.not. this%err%ok) return

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) then
                if (size(Xt,1) /= size(this%Xt,1)) deallocate(this%Xt)
            end if
            this%Xt = Xt
        elseif (present(res)) then
            if (allocated(this%Xt)) then
                if (size(this%Xt,1) /= res) then
                    deallocate(this%Xt)
                    allocate(this%Xt(res))
                end if
            else
                allocate(this%Xt(res))
            end if
            this%Xt = [(this%knot(1)+(this%knot(size(this%knot))-this%knot(1))*real(i-1,rk)/real(res-1,rk), i=1, res)]
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
    pure subroutine derivative_scalar(this, Xt, dTgc, Tgc, elem)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xt
        integer, intent(in), optional :: elem(:)
        real(rk), allocatable, intent(out) :: dTgc(:)
        real(rk), allocatable, intent(out), optional :: Tgc(:)

        if (.not. this%err%ok) return

        if (this%is_rational()) then ! NURBS
            if (present(elem)) then
                associate(Wce => this%Wc(elem))
                    call compute_dTgc(Xt, this%knot, this%degree, this%nc, Wce, dTgc, Tgc, elem)
                end associate
            else
                call compute_dTgc(Xt, this%knot, this%degree, this%nc, this%Wc, dTgc, Tgc, elem)
            end if
        else ! B-Spline
            call compute_dTgc(Xt, this%knot, this%degree, this%nc, dTgc, Tgc, elem)
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

        if (.not. this%err%ok) return

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) then
                if (size(Xt,1) /= size(this%Xt,1)) deallocate(this%Xt)
            end if
            this%Xt = Xt
        elseif (present(res)) then
            if (allocated(this%Xt)) then
                if (size(this%Xt,1) /= res) then
                    deallocate(this%Xt)
                    allocate(this%Xt(res))
                end if
            else
                allocate(this%Xt(res))
            end if
            this%Xt = [(this%knot(1)+(this%knot(size(this%knot))-this%knot(1))*real(i-1,rk)/real(res-1,rk), i=1, res)]
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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) then
                if (size(Xt,1) /= size(this%Xt,1)) deallocate(this%Xt)
            end if
            this%Xt = Xt
        elseif (present(res)) then
            if (allocated(this%Xt)) then
                if (size(this%Xt,1) /= res) then
                    deallocate(this%Xt)
                    allocate(this%Xt(res))
                end if
            else
                allocate(this%Xt(res))
            end if
            this%Xt = [(this%knot(1)+(this%knot(size(this%knot))-this%knot(1))*real(i-1,rk)/real(res-1,rk), i=1, res)]
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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

        r = .false.
        if (allocated(this%Wc)) then
            ! if (any(this%Wc /= this%Wc(1))) then
            if (any(abs(this%Wc - this%Wc(1)) > 2.0_rk*epsilon(0.0_rk))) then
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

        if (.not. this%err%ok) return

        if (allocated(this%elemConn_Xc_vis)) then
            if (size(this%elemConn_Xc_vis,1) /= size(elemConn,1) .or. size(this%elemConn_Xc_vis,2) /= size(elemConn,2)) then
                deallocate(this%elemConn_Xc_vis)
            end if
        end if
        this%elemConn_Xc_vis = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_elem_Xg_vis(this, elemConn)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (.not. this%err%ok) return

        if (allocated(this%elemConn_Xg_vis)) then
            if (size(this%elemConn_Xg_vis,1) /= size(elemConn,1) .or. size(this%elemConn_Xg_vis,2) /= size(elemConn,2)) then
                deallocate(this%elemConn_Xg_vis)
            end if
        end if
        this%elemConn_Xg_vis = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_elem(this, elemConn)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (.not. this%err%ok) return

        if (allocated(this%elemConn)) then
            if (size(this%elemConn,1) /= size(elemConn,1) .or. size(this%elemConn,2) /= size(elemConn,2)) then
                deallocate(this%elemConn)
            end if
        end if
        this%elemConn = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem_Xc_vis(this) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        elemConn = this%elemConn_Xc_vis
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem_Xg_vis(this) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        elemConn = this%elemConn_Xg_vis
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem(this) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

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
        integer :: k, i, s, d, j, nc_new, t
        real(rk), allocatable :: Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:)

        if (.not. this%err%ok) return

        if (this%is_rational()) then ! NURBS

            do i = 1, size(Xth)
                k = findspan(this%nc-1,this%degree,Xth(i),this%knot)
                ! if (this%knot(k+1) == Xth(i)) then
                if (abs(this%knot(k+1) - Xth(i)) < 2.0_rk*epsilon(0.0_rk)) then
                    s = compute_multiplicity(this%knot,Xth(i))
                else
                    s = 0
                end if
                k = k + 1

                d = size(this%Xc,2)
                allocate(Xcw(size(this%Xc,1),d+1))
                do j = 1, size(this%Xc,1)
                    Xcw(j,1:d) = this%Xc(j,1:d)*this%Wc(j)
                end do
                Xcw(:,d+1) = this%Wc(:)

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
                    allocate(Xc_new(nc_new,d))
                    allocate(Wc_new(nc_new))
                    do j = 1, nc_new
                        Xc_new(j,:) = Xcw_new(j,1:d)/Xcw_new(j,d+1)
                    end do
                    Wc_new(:) = Xcw_new(:,d+1)

                    call this%set(knot=knot_new, Xc=Xc_new, Wc=Wc_new)
                    if (allocated(Xcw_new)) deallocate(Xcw_new)
                    if (allocated(Xc_new)) deallocate(Xc_new)
                    if (allocated(Wc_new)) deallocate(Wc_new)
                end if
            end do

        else ! B-Spline

            do i = 1, size(Xth)
                k = findspan(this%nc-1,this%degree,Xth(i),this%knot)
                ! if (this%knot(k+1) == Xth(i)) then
                if (abs(this%knot(k+1) - Xth(i)) < 2.0_rk*epsilon(0.0_rk)) then
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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

        do i = 1, this%ng
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

        if (.not. this%err%ok) return

#ifndef NOSHOW_PYVISTA
        block
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
            "p.show(title='ForCAD', interactive=True)"//achar(10)//&
            "p.deep_clean()"//achar(10)//&
            "del p"
        call execute_command_line('python -c "'//trim(adjustl(pyvista_script))//'"')
        end block
#endif
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

        if (.not. this%err%ok) return

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
        real(rk), intent(in), contiguous :: point_Xg(:)
        real(rk), intent(out), optional :: nearest_Xg(size(point_Xg))
        real(rk), intent(out), optional :: nearest_Xt
        integer, intent(out), optional :: id
        integer :: id_, i
        real(rk), allocatable :: distances(:)

        if (.not. this%err%ok) return

        allocate(distances(this%ng))

#if defined(__NVCOMPILER)
        do i = 1, this%ng
#else
        do concurrent (i = 1: this%ng)
#endif
            distances(i) = norm2(this%Xg(i,:) - point_Xg)
        end do

        ! replaced minloc due to NVFortran bug
#if defined(__NVCOMPILER)
        id_ = 1
        do i = 2, size(distances)
            if (distances(i) < distances(id_)) id_ = i
        end do
#else
        id_ = minloc(distances, dim=1)
#endif
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
        real(rk), intent(in), contiguous :: point_Xg(:)
        real(rk), intent(in) :: tol
        integer, intent(in) :: maxit
        real(rk), intent(out) :: nearest_Xt
        real(rk), intent(out), optional :: nearest_Xg(size(this%Xc,2))

        real(rk) :: xk, xkn, obj, obj_trial, grad, hess, dk, alphak
        real(rk) :: tau, beta, eps, lower_bounds, upper_bounds, alpha_max, alpha_i, xt
        real(rk) :: Xg(size(this%Xc,2))
        real(rk), allocatable :: Tgc(:), dTgc(:), d2Tgc(:)
        integer :: k, l
        logical :: convergenz
        type(nurbs_curve) :: copy_this

        if (.not. this%err%ok) return

        dk  = 0.0_rk
        k   = 0
        eps = 10.0_rk*tiny(1.0_rk)

        ! bounds
        lower_bounds = minval(this%knot)
        upper_bounds = maxval(this%knot)

        ! initial guess (coarse search)
        copy_this = this
        call copy_this%create(10)
        call copy_this%nearest_point(point_Xg=point_Xg, nearest_Xt=xk)
        call copy_this%finalize()

        ! clamp initial guess to bounds
        xk = max(min(xk, upper_bounds), lower_bounds)

        xkn = xk
        convergenz = .false.

        do while (.not. convergenz .and. k < maxit)

            ! objective, gradient, hessian
            Xg = this%cmp_Xg(xk)
            call this%derivative2(Xt=xk, d2Tgc=d2Tgc, dTgc=dTgc, Tgc=Tgc)  ! Tgc unused

            obj  = norm2(Xg - point_Xg) + 0.001_rk  ! small epsilon to avoid divide-by-zero
            grad = dot_product((Xg-point_Xg)/obj, matmul(dTgc, this%Xc))
            hess = dot_product(matmul(dTgc,this%Xc) - (Xg-point_Xg)/obj*grad, matmul(dTgc,this%Xc))/obj &
                + dot_product((Xg-point_Xg)/obj, matmul(d2Tgc,this%Xc))

            ! debug
            print '(i3,1x,e20.10,1x,e20.10)', k, xk, abs(grad)

            if (abs(grad) <= tol .or. (k>0 .and. abs(xk-xkn) <= tol)) then
                convergenz = .true.
                nearest_Xt = xk
                if (present(nearest_Xg)) nearest_Xg = this%cmp_Xg(nearest_Xt)
            else
                ! Newton step
                dk = - grad / hess

                ! Backtracking-Armijo with feasibility
                tau  = 0.5_rk
                beta = 1.0e-4_rk

                ! maximum feasible step so xk + alpha*dk stays within [lower_bounds, upper_bounds]
                if (dk > 0.0_rk) then
                    if (upper_bounds > xk) then
                        alpha_i  = (upper_bounds - xk) / dk
                        alpha_max = max(0.0_rk, alpha_i)
                    else
                        alpha_max = 0.0_rk
                    end if
                else if (dk < 0.0_rk) then
                    if (lower_bounds < xk) then
                        alpha_i  = (lower_bounds - xk) / dk
                        alpha_max = max(0.0_rk, alpha_i)
                    else
                        alpha_max = 0.0_rk
                    end if
                else
                    alpha_max = 0.0_rk
                end if

                if (alpha_max <= eps) then
                    convergenz = .true.
                    nearest_Xt = xk
                    if (present(nearest_Xg)) nearest_Xg = this%cmp_Xg(nearest_Xt)
                    exit
                end if

                alphak = min(1.0_rk, alpha_max)
                l = 0
                do
                    if (alphak <= eps .or. l >= 50) exit
                    xt = xk + alphak*dk        ! feasible since alphak ≤ alpha_max
                    obj_trial = norm2(this%cmp_Xg(xt) - point_Xg) + 0.001_rk
                    if (obj_trial <= obj + alphak*beta*grad*dk) exit
                    alphak = min(tau*alphak, alpha_max)  ! shrink but stay feasible
                    l = l + 1
                end do

                xkn = xk
                if (alphak > eps) then
                    xk = xk + alphak*dk
                end if

                ! clamp updated iterate
                xk = max(min(xk, upper_bounds), lower_bounds)

                k = k + 1
            end if
        end do

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine ansatz(this, ie, ig, Tgc, dTgc_dXg, dL, ngauss)
        class(nurbs_curve), intent(inout) :: this

        integer, intent(in) :: ie, ig
        real(rk), intent(out) :: dL
        real(rk), allocatable, intent(out) :: Tgc(:), dTgc_dXg(:,:)
        integer, intent(in), optional :: ngauss
        real(rk), allocatable :: Xth(:), Xth_e(:), Xc_eT(:,:), Xksi(:), Wksi(:)
        integer, allocatable :: elem_th(:,:), elem_c(:,:), elem_ce(:)
        type(nurbs_curve) :: th, th_e
        real(rk), allocatable :: dTtth_dXksi(:), Ttth(:), dTgc_dXt(:), dXg_dXt(:)
        real(rk) :: Xt, dXt_dXksi
        real(rk), allocatable :: dXg_dXksi(:) !! Jacobian matrix
        real(rk) :: det_dXg_dXksi !! Determinant of the Jacobian matrix

        if (.not. this%err%ok) return

        if (present(ngauss)) then
            call gauss_leg([0.0_rk, 1.0_rk], ngauss-1, Xksi, Wksi)
        else
            call gauss_leg([0.0_rk, 1.0_rk], this%degree, Xksi, Wksi)
        end if

        Xth = unique(this%knot)

        call th%set([0.0_rk,Xth,1.0_rk], Xth)
        elem_th = th%cmp_elem()

        elem_c = this%cmp_elem()

        Xth_e = Xth(elem_th(ie,:))
        call th_e%set([0.0_rk,0.0_rk,1.0_rk,1.0_rk], Xth_e)

        elem_ce = elem_c(ie,:)
        Xc_eT = transpose(this%Xc(elem_ce,:))

        call th_e%derivative(Xksi(ig), dTtth_dXksi, Ttth)
        Xt = dot_product(Xth_e, Ttth)
        dXt_dXksi = dot_product(Xth_e, dTtth_dXksi)

        call this%derivative(Xt, dTgc_dXt, Tgc, elem_ce)
        dXg_dXt = matmul(Xc_eT, dTgc_dXt)

        dTgc_dXg = dyad(dTgc_dXt, dXg_dXt)/norm2(dXg_dXt)

        dXg_dXksi = dXg_dXt*dXt_dXksi
        det_dXg_dXksi = norm2(dXg_dXksi)

        dL = det_dXg_dXksi*Wksi(ig)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine cmp_length(this, length, ngauss)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(out) :: length
        integer, intent(in), optional :: ngauss
        real(rk), allocatable :: Tgc(:), dTgc_dXg(:,:)
        integer :: ie, ig, ngauss_
        real(rk) :: dL, dL_ig

        if (.not. this%err%ok) return

        if (present(ngauss)) then
            ngauss_ = ngauss
        else
            ngauss_ = this%degree + 1
        end if

        length = 0.0_rk
#if defined(__NVCOMPILER) || (defined(__GFORTRAN__) && (__GNUC__ < 15 || (__GNUC__ == 15 && __GNUC_MINOR__ < 1)))
        do ie = 1, size(this%cmp_elem(),1)
#else
        do concurrent (ie = 1: size(this%cmp_elem(),1)) reduce(+:length)
#endif
            dL = 0.0_rk
            do ig = 1, ngauss_
                call this%ansatz(ie, ig, Tgc, dTgc_dXg, dL_ig, ngauss_)
                dL = dL + dL_ig
            end do
            length = length + dL
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_Tgc_1d(Xti, knot, nc, degree, Wc) result(Tgc)
        real(rk), intent(in) :: Xti
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree, nc
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk) :: Tgc(nc)
        real(rk) :: tmp
        integer :: i

        Tgc = basis_bspline(Xti, knot, nc, degree)
        tmp = dot_product(Tgc, Wc)
        do concurrent (i = 1: nc)
           Tgc(i) = (Tgc(i) * Wc(i)) / tmp
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_Xg_nurbs_1d(Xt, knot, degree, nc, ng, Xc, Wc) result(Xg)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in) :: ng
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable :: Xg(:,:)
        integer :: i

        allocate(Xg(ng, size(Xc,2)), source = 0.0_rk)
#if defined(__NVCOMPILER)
        do i = 1, ng
#else
        do concurrent (i = 1: ng)
#endif
            Xg(i,:) = matmul(cmp_Tgc_1d(Xt(i), knot, nc, degree, Wc), Xc)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_Xg_nurbs_1d_1point(Xt, knot, degree, nc, Xc, Wc) result(Xg)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk) :: Xg(size(Xc,2))
        real(rk), allocatable :: Tgc(:)

        allocate(Tgc(nc))
        Tgc = basis_bspline(Xt, knot, nc, degree)
        Tgc = Tgc*(Wc/(dot_product(Tgc,Wc)))
        Xg = matmul(Tgc, Xc)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_Xg_bspline_1d(Xt, knot, degree, nc, ng, Xc) result(Xg)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in) :: ng
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), allocatable :: Xg(:,:)
        integer :: i

        allocate(Xg(ng, size(Xc,2)))
#if defined(__NVCOMPILER)
        do i = 1, ng
#else
        do concurrent (i = 1: ng)
#endif
            Xg(i,:) = matmul(basis_bspline(Xt(i), knot, nc, degree), Xc)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_Xg_bspline_1d_1point(Xt, knot, degree, nc, Xc) result(Xg)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk) :: Xg(size(Xc,2))

        Xg = matmul(basis_bspline(Xt, knot, nc, degree), Xc)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine compute_d2Tgc_nurbs_1d_vector(Xt, knot, degree, nc, ng, Wc, d2Tgc, dTgc, Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in) :: ng
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        real(rk) :: d2Bi(nc), dBi(nc), Bi(nc)
        integer :: i

        allocate(d2Tgc(ng, nc), dTgc(ng, nc), Tgc(ng, nc))

#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, size(Xt)
#else
        do concurrent (i = 1: size(Xt)) local(d2Bi, dBi, Bi)
#endif
            call basis_bspline_2der(Xt(i), knot, nc, degree, d2Bi, dBi, Bi)
            Tgc(i,:) = Bi*(Wc/(dot_product(Bi,Wc)))
            dTgc(i,:) = ( dBi*Wc - Tgc(i,:)*dot_product(dBi,Wc) ) / dot_product(Bi,Wc)
            d2Tgc(i,:) = (d2Bi*Wc - 2.0_rk*dTgc(i,:)*dot_product(dBi,Wc) - Tgc(i,:)*dot_product(d2Bi,Wc)) / dot_product(Bi,Wc)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine compute_d2Tgc_nurbs_1d_scalar(Xt, knot, degree, nc, Wc, d2Tgc, dTgc, Tgc)
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
    pure subroutine compute_d2Tgc_bspline_1d_vector(Xt, knot, degree, nc, ng, d2Tgc, dTgc, Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in) :: ng
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        integer :: i
        real(rk) :: Xti, d2Tgci(nc), dTgci(nc), Tgci(nc)

        allocate(d2Tgc(ng, nc), dTgc(ng, nc), Tgc(ng, nc))
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, size(Xt)
#else
        do concurrent (i = 1: size(Xt)) local(Xti, d2Tgci, dTgci, Tgci)
#endif
            Xti = Xt(i)
            call basis_bspline_2der(Xti, knot, nc, degree, d2Tgci , dTgci, Tgci)
            d2Tgc(i,:) = d2Tgci
            dTgc(i,:) = dTgci
            Tgc(i,:) = Tgci
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine compute_d2Tgc_bspline_1d_scalar(Xt, knot, degree, nc, d2Tgc, dTgc, Tgc)
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
    pure subroutine compute_dTgc_nurbs_1d_vector(Xt, knot, degree, nc, ng, Wc, dTgc, Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in) :: ng
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        real(rk) :: dBi(nc), Bi(nc)
        integer :: i

        allocate(dTgc(ng, nc), Tgc(ng, nc))
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, size(Xt)
#else
        do concurrent (i = 1: size(Xt)) local(dBi, Bi)
#endif
            call basis_bspline_der(Xt(i), knot, nc, degree, dBi, Bi)
            Tgc(i,:) = Bi*(Wc/(dot_product(Bi,Wc)))
            dTgc(i,:) = ( dBi*Wc - Tgc(i,:)*dot_product(dBi,Wc) ) / dot_product(Bi,Wc)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine compute_dTgc_nurbs_1d_scalar(Xt, knot, degree, nc, Wc, dTgc, Tgc, elem)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        real(rk), intent(in), contiguous :: Wc(:)
        integer, intent(in), optional :: elem(:)
        real(rk), allocatable, intent(out) :: dTgc(:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        real(rk) :: dBi(nc), Bi(nc)

        call basis_bspline_der(Xt, knot, nc, degree, dBi, Bi)

        if (.not. present(elem)) then
            allocate(dTgc(nc), Tgc(nc))
            Tgc = Bi*(Wc/(dot_product(Bi,Wc)))
            dTgc = ( dBi*Wc - Tgc*dot_product(dBi,Wc) ) / dot_product(Bi,Wc)
        else
            allocate(dTgc(size(elem)), Tgc(size(elem)))
            Tgc = Bi(elem)*(Wc/(dot_product(Bi(elem),Wc)))
            dTgc = ( dBi(elem)*Wc - Tgc*dot_product(dBi(elem),Wc) ) / dot_product(Bi(elem),Wc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine compute_dTgc_bspline_1d_vector(Xt, knot, degree, nc, ng, dTgc, Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in) :: ng
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        integer :: i
        real(rk) :: Xti, dTgci(nc), Tgci(nc)

        allocate(dTgc(ng, nc), Tgc(ng, nc))
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, size(Xt)
#else
        do concurrent (i = 1: size(Xt)) local(Xti, dTgci, Tgci)
#endif
            Xti = Xt(i)
            call basis_bspline_der(Xti, knot, nc, degree, dTgci, Tgci)
            dTgc(i,:) = dTgci
            Tgc(i,:) = Tgci
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine compute_dTgc_bspline_1d_scalar(Xt, knot, degree, nc, dTgc, Tgc, elem)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in), optional :: elem(:)
        real(rk), allocatable, intent(out) :: dTgc(:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        real(rk), allocatable :: dB(:), B(:)

        if (.not. present(elem)) then
            allocate(dTgc(nc), Tgc(nc))
            call basis_bspline_der(Xt, knot, nc, degree, dTgc, Tgc)
        else
            allocate(dB(size(elem)), B(size(elem)))
            call basis_bspline_der(Xt, knot, nc, degree, dB, B)
            Tgc = B(elem)
            dTgc = dB(elem)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_Tgc_nurbs_1d_vector(Xt, knot, degree, nc, ng, Wc) result(Tgc)
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
#if defined(__NVCOMPILER)
        do i = 1, size(Xt,1)
#else
        do concurrent (i = 1: size(Xt,1))
#endif
            Tgci = basis_bspline(Xt(i), knot, nc, degree)
            Tgc(i,:) = Tgci*(Wc/(dot_product(Tgci,Wc)))
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_Tgc_nurbs_1d_scalar(Xt, knot, degree, nc, Wc) result(Tgc)
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
    pure function compute_Tgc_bspline_1d_vector(Xt, knot, degree, nc, ng) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in) :: ng
        real(rk), allocatable :: Tgc(:,:)
        integer :: i

        allocate(Tgc(ng, nc))
#if defined(__NVCOMPILER)
        do i = 1, size(Xt,1)
#else
        do concurrent (i = 1: size(Xt,1))
#endif
            Tgc(i,:) = basis_bspline(Xt(i), knot, nc, degree)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_Tgc_bspline_1d_scalar(Xt, knot, degree, nc) result(Tgc)
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
    pure subroutine lsq_fit_bspline(this, Xt, Xdata, ndata)
        use forcad_interface, only: solve
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:), Xdata(:,:)
        integer, intent(in) :: ndata
        real(rk), allocatable :: T(:,:), Tt(:,:), TtT(:,:), TtX(:,:)
        integer :: i

        if (.not. this%err%ok) return

        if (this%nc > ndata) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = 'forcad_nurbs_curve',&
                message    = 'Too few data points for the requested number of control points.',&
                location   = 'lsq_fit_bspline',&
                suggestion = 'Use nc <= ndata: reduce nc or increase the number of data points.')
            return
        end if

        allocate(T(ndata, this%nc))
#if defined(__NVCOMPILER)
        do i = 1, ndata
#else
        do concurrent (i = 1: ndata)
#endif
            T(i,:) = basis_bspline(Xt(i), this%knot, this%nc, this%degree)
        end do
        Tt = transpose(T)
        TtT = matmul(Tt, T)
        TtX = matmul(Tt, Xdata)
        this%Xc = solve(TtT, TtX)
    end subroutine
   !===============================================================================

end module forcad_nurbs_curve