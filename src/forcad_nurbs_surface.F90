!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> This module defines the 'nurbs_surface' type for representing a Non-Uniform Rational B-Spline (NURBS) surface.
module forcad_nurbs_surface

    use forcad_kinds, only: rk
    use forcad_utils, only: basis_bspline, elemConn_C0, kron, ndgrid, compute_multiplicity, compute_knot_vector, &
        basis_bspline_der, insert_knot_A_5_1, findspan, elevate_degree_A_5_9, remove_knots_A_5_8, tetragon_Xc, &
        elemConn_Cn, unique, rotation, det, inv, gauss_leg, export_vtk_legacy, basis_bspline_2der
    use fordebug, only: debug

    implicit none

    private
    public nurbs_surface, compute_Tgc, compute_dTgc

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

        type(debug) :: err !! 101: size mismatch (weights vs control points), 102: missing control points, 103: missing knot vector, 104: missing geometry points, 105: missing weights, 106: lsq fit underdetermined
    contains
        procedure, private :: set1                   !!> Set knot vectors, control points and weights for the NURBS surface object
        procedure, private :: set2                   !!> Set NURBS surface using nodes of parameter space, degree, continuity, control points and weights
        procedure, private :: set3                   !!> Set Bezier or Rational Bezier surface using control points and weights
        procedure, private :: set4                   !!> Set NURBS surface using degree, number of control points, control points and weights
        generic :: set => set1, set2, set3, set4  !!> Set NURBS surface
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
        procedure, private :: get_nc_dir             !!> Get number of control points in a specific direction
        procedure, private :: get_nc_all             !!> Get number of control points in all directions
        generic :: get_nc => get_nc_all, get_nc_dir !!> Get number of control points
        procedure :: cmp_degree             !!> Compute degree of the NURBS surface
        procedure, private :: get_degree_all!!> Get degree of the NURBS surface in both directions
        procedure, private :: get_degree_dir!!> Get degree of the NURBS surface in a specific direction
        generic :: get_degree => get_degree_all, get_degree_dir !!> Get degree of the NURBS surface
        procedure :: finalize               !!> Finalize the NURBS surface object
        procedure :: cmp_elem_Xc_vis        !!> Generate connectivity for control points
        procedure :: cmp_elem_Xg_vis        !!> Generate connectivity for geometry points
        procedure :: cmp_elem_Xth           !!> Generate connectivity for parameter points
        procedure :: cmp_elem               !!> Generate IGA element connectivity
        procedure :: get_elem_Xc_vis        !!> Get connectivity for control points
        procedure :: get_elem_Xg_vis        !!> Get connectivity for geometry points
        procedure :: get_elem               !!> Get IGA element connectivity
        procedure :: set_elem_Xc_vis        !!> Set connectivity for control points
        procedure :: set_elem_Xg_vis        !!> Set connectivity for geometry points
        procedure :: set_elem               !!> Set IGA element connectivity
        procedure :: export_Xc              !!> Export control points to VTK file
        procedure :: export_Xg              !!> Export geometry points to VTK file
        procedure :: export_Xth             !!> Export parameter space to VTK file
        procedure :: export_Xth_in_Xg       !!> Export parameter space in geometry points to VTK file
        procedure :: export_iges            !!> Export the NURBS surface to IGES format
        procedure :: modify_Xc              !!> Modify control points
        procedure :: modify_Wc              !!> Modify weights
        procedure :: get_multiplicity       !!> Compute and return the multiplicity of the knot vector
        procedure :: get_continuity         !!> Compute and return the continuity of the NURBS surface
        procedure :: cmp_nc                 !!> Compute number of required control points
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
        procedure :: nearest_point2         !!> Find the nearest point on the NURBS surface (Minimization - Newtons method)
        procedure :: ansatz                 !!> Compute the shape functions, derivative of shape functions and dA
        procedure :: cmp_area               !!> Compute the area of the NURBS surface
        procedure :: lsq_fit_bspline        !!> Fit B-spline volume to structured data points using least squares

        ! Shapes
        procedure :: set_tetragon           !!> Set a tetragon
        procedure :: set_ring               !!> Set a ring
        procedure :: set_half_ring          !!> Set a half ring
        procedure :: set_C                  !!> Set a C-shape
    end type
    !===============================================================================

    interface compute_Xg
        module procedure compute_Xg_nurbs_2d
        module procedure compute_Xg_bspline_2d
        module procedure compute_Xg_nurbs_2d_1point
        module procedure compute_Xg_bspline_2d_1point
    end interface

    interface compute_Tgc
        module procedure compute_Tgc_nurbs_2d_vector
        module procedure compute_Tgc_bspline_2d_vector
        module procedure compute_Tgc_nurbs_2d_scalar
        module procedure compute_Tgc_bspline_2d_scalar
    end interface

    interface compute_dTgc
        module procedure compute_dTgc_nurbs_2d_vector
        module procedure compute_dTgc_bspline_2d_vector
        module procedure compute_dTgc_nurbs_2d_scalar
        module procedure compute_dTgc_bspline_2d_scalar
    end interface

    interface compute_d2Tgc
        module procedure compute_d2Tgc_nurbs_2d_vector
        module procedure compute_d2Tgc_bspline_2d_vector
        module procedure compute_d2Tgc_nurbs_2d_scalar
        module procedure compute_d2Tgc_bspline_2d_scalar
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

        if (.not. this%err%ok) return

        if (allocated(this%knot1)) then
            if (size(this%knot1) /= size(knot1)) deallocate(this%knot1)
        end if
        if (allocated(this%knot2)) then
            if (size(this%knot2) /= size(knot2)) deallocate(this%knot2)
        end if
        if (allocated(this%Xc)) then
            if (size(this%Xc,1) /= size(Xc,1) .or. size(this%Xc,2) /= size(Xc,2)) deallocate(this%Xc)
        end if

        this%knot1 = knot1
        this%knot2 = knot2
        call this%cmp_degree()
        call this%cmp_nc()
        this%Xc = Xc

        if (present(Wc)) then
            if (size(Wc) /= this%nc(1)*this%nc(2)) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_surface',&
                    message    = 'Weights length mismatch: size(Wc) must equal number of control points.',&
                    location   = 'set1',&
                    suggestion = 'Provide Wc with size(Wc) == nc(1)*nc(2).')
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
    !> Set NURBS surface using nodes of parameter space, degree, continuity, control points and weights
    pure subroutine set2(this, Xth_dir1, Xth_dir2, degree, continuity1, continuity2, Xc, Wc)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xth_dir1(:), Xth_dir2(:)
        integer, intent(in), contiguous :: degree(:)
        integer, intent(in), contiguous :: continuity1(:), continuity2(:)
        real(rk), intent(in), contiguous, optional :: Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)

        if (.not. this%err%ok) return

        if (allocated(this%knot1)) deallocate(this%knot1)
        if (allocated(this%knot2)) deallocate(this%knot2)
        this%knot1    = compute_knot_vector(Xth_dir1, degree(1), continuity1)
        this%knot2    = compute_knot_vector(Xth_dir2, degree(2), continuity2)
        this%degree   = degree
        call this%cmp_nc()

        if (present(Xc)) then
            if (size(Xc,1) /= this%nc(1)*this%nc(2)) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_surface', &
                    message    = 'Control points size mismatch in set2',&
                    location   = 'set2', &
                    suggestion = 'size(Xc,1) must equal nc(1)*nc(2).' )
                return
            end if
        end if
        if (present(Wc)) then
            if (size(Wc) /= this%nc(1)*this%nc(2)) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_surface', &
                    message    = 'Weights size mismatch in set2',&
                    location   = 'set2', &
                    suggestion = 'size(Wc) must equal nc(1)*nc(2).' )
                return
            end if
        end if

        if (present(Xc)) then
            if (allocated(this%Xc)) then
                if (size(this%Xc,1) /= size(Xc,1) .or. size(this%Xc,2) /= size(Xc,2)) deallocate(this%Xc)
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
    !> Set Bezier or Rational Bezier surface using control points and weights.
    pure subroutine set3(this, nc, Xc, Wc)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), contiguous :: nc(:)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)

        if (.not. this%err%ok) return

        if (allocated(this%Xc)) then
            if (size(this%Xc,1) /= size(Xc,1) .or. size(this%Xc,2) /= size(Xc,2)) deallocate(this%Xc)
        end if

        this%Xc = Xc
        this%nc = nc

        if (allocated(this%knot1)) then
            if (size(this%knot1) /= 2*this%nc(1)) then
                deallocate(this%knot1)
                allocate(this%knot1(2*this%nc(1)))
            end if
        else
            allocate(this%knot1(2*this%nc(1)))
        end if
        this%knot1(1:this%nc(1)) = 0.0_rk
        this%knot1(this%nc(1)+1:2*this%nc(1)) = 1.0_rk

        if (allocated(this%knot2)) then
            if (size(this%knot2) /= 2*this%nc(2)) then
                deallocate(this%knot2); allocate(this%knot2(2*this%nc(2)))
            end if
        else
            allocate(this%knot2(2*this%nc(2)))
        end if
        this%knot2(1:this%nc(2)) = 0.0_rk
        this%knot2(this%nc(2)+1:2*this%nc(2)) = 1.0_rk

        call this%cmp_degree()
        if (present(Wc)) then
            if (size(Wc) /= this%nc(1)*this%nc(2)) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_surface',&
                    message    = 'Weights length mismatch: size(Wc) must equal number of control points.',&
                    location   = 'set3',&
                    suggestion = 'Provide Wc with size(Wc) == nc(1)*nc(2).')
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
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), contiguous :: degree(:)
        integer, intent(in), contiguous :: nc(:)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        integer :: m(2), i

        if (.not. this%err%ok) return

        if (allocated(this%Xc)) then
            if (size(this%Xc,1) /= size(Xc,1) .or. size(this%Xc,2) /= size(Xc,2)) deallocate(this%Xc)
        end if

        this%Xc = Xc
        this%nc = nc
        this%degree = degree

        m = nc + degree + 1

        if (allocated(this%knot1)) then
            if (size(this%knot1) /= m(1)) then
                deallocate(this%knot1)
                allocate(this%knot1(m(1)))
            end if
        else
            allocate(this%knot1(m(1)))
        end if
        this%knot1(1:degree(1)+1) = 0.0_rk
        this%knot1(degree(1)+2:m(1)-degree(1)-1) = [(real(i, rk)/(m(1)-2*degree(1)-1), i=1, m(1)-2*degree(1)-2)]
        this%knot1(m(1)-degree(1):m(1)) = 1.0_rk

        if (allocated(this%knot2)) then
            if (size(this%knot2) /= m(2)) then
                deallocate(this%knot2); allocate(this%knot2(m(2)))
            end if
        else
            allocate(this%knot2(m(2)))
        end if
        this%knot2(1:degree(2)+1) = 0.0_rk
        this%knot2(degree(2)+2:m(2)-degree(2)-1) = [(real(i, rk)/(m(2)-2*degree(2)-1), i=1, m(2)-2*degree(2)-2)]
        this%knot2(m(2)-degree(2):m(2)) = 1.0_rk

        if (present(Wc)) then
            if (size(Wc) /= nc(1)*nc(2)) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_surface',&
                    message    = 'Weights length mismatch: size(Wc) must equal number of control points.',&
                    location   = 'set4',&
                    suggestion = 'Provide Wc with size(Wc) == nc(1)*nc(2).')
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
    pure subroutine create(this, res1, res2, Xt1, Xt2, Xt)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: res1, res2
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:)
        real(rk), contiguous, intent(in), optional :: Xt(:,:)
        integer :: i

        if (.not. this%err%ok) return

        ! check
        if (.not.allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'Control points are not set.',&
                location   = 'create',&
                suggestion = 'Call set(...) first before create().')
            return
        end if

        if (.not.allocated(this%knot1) .or. .not.allocated(this%knot2)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'Knot vector is not set.',&
                location   = 'create',&
                suggestion = 'Call set(...) first before create().')
            return
        end if

        if (present(Xt1)) then
            if (allocated(this%Xt1)) then
                if (size(this%Xt1) /= size(Xt1)) deallocate(this%Xt1)
            end if
            this%Xt1 = Xt1
        elseif (present(res1)) then
            if (allocated(this%Xt1)) then
                if (size(this%Xt1) /= res1) then
                    deallocate(this%Xt1)
                    allocate(this%Xt1(res1))
                end if
            else
                allocate(this%Xt1(res1))
            end if
            this%Xt1 = [(this%knot1(1)+(this%knot1(size(this%knot1))-this%knot1(1))*real(i-1,rk)/real(res1-1,rk), i=1, res1)]
        end if

        if (present(Xt2)) then
            if (allocated(this%Xt2)) then
                if (size(this%Xt2) /= size(Xt2)) deallocate(this%Xt2)
            end if
            this%Xt2 = Xt2
        elseif (present(res2)) then
            if (allocated(this%Xt2)) then
                if (size(this%Xt2) /= res2) then
                    deallocate(this%Xt2)
                    allocate(this%Xt2(res2))
                end if
            else
                allocate(this%Xt2(res2))
            end if
            this%Xt2 = [(this%knot2(1)+(this%knot2(size(this%knot2))-this%knot2(1))*real(i-1,rk)/real(res2-1,rk), i=1, res2)]
        end if

        if (present(Xt)) then
            this%Xt = Xt
        else
            this%ng(1) = size(this%Xt1,1)
            this%ng(2) = size(this%Xt2,1)
            call ndgrid(this%Xt1, this%Xt2, this%Xt)
        end if

        if (allocated(this%Xg)) then
            if (size(this%Xg,1) /= this%ng(1)*this%ng(2) .or. size(this%Xg,2) /= size(this%Xc,2)) then
                deallocate(this%Xg)
            end if
        end if

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: n
        integer, intent(in) :: dir
        real(rk) :: Xc

        if (.not. this%err%ok) return

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
        class(nurbs_surface), intent(in) :: this
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
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: n
        integer, intent(in) :: dir
        real(rk) :: Xg

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

        if (present(dir)) then
            if (dir == 1) then
                m1 = this%get_multiplicity(1)
                this%degree(1) = m1(1) - 1
            else if (dir == 2) then
                m2 = this%get_multiplicity(2)
                this%degree(2) = m2(1) - 1
            else
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_surface',&
                    message    = 'Invalid direction for degree.',&
                    location   = 'cmp_degree',&
                    suggestion = 'Check the direction argument.')
                return
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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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
    pure function cmp_elem_Xth(this, p) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), contiguous, optional :: p(:)

        if (.not. this%err%ok) return

        if (present(p)) then
            elemConn = elemConn_C0(size(unique(this%knot1)), size(unique(this%knot2)), p(1), p(2))
        else
            elemConn = elemConn_C0(size(unique(this%knot1)), size(unique(this%knot2)), 1, 1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xc(this, filename, point_data, field_names, encoding)
        class(nurbs_surface), intent(inout) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), optional :: point_data(:,:)
        character(len=*), intent(in), optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        if (.not.allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
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

        call export_vtk_legacy(filename=filename, points=this%Xc, elemConn=elemConn, vtkCellType=9, &
                               point_data=point_data, field_names=field_names, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xg(this, filename, point_data, field_names, encoding)
        class(nurbs_surface), intent(inout) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), optional :: point_data(:,:)
        character(len=*), intent(in), optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        if (.not.allocated(this%Xg)) then
            call this%err%set(&
                code       = 104,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
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

        call export_vtk_legacy(filename=filename, points=this%Xg, elemConn=elemConn, vtkCellType=9, &
                               point_data=point_data, field_names=field_names, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xth(this, filename, point_data, field_names, encoding)
        class(nurbs_surface), intent(in) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), optional :: point_data(:,:)
        character(len=*), intent(in), optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Xth(:,:), Xth1(:), Xth2(:)

        if (.not. this%err%ok) return

        elemConn = this%cmp_elem_Xth()
        Xth1 = unique(this%knot1)
        Xth2 = unique(this%knot2)
        call ndgrid(Xth1, Xth2, Xth)

        call export_vtk_legacy(filename=filename, points=Xth, elemConn=elemConn, vtkCellType=9, &
                               point_data=point_data, field_names=field_names, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xth_in_Xg(this, filename, res, encoding)
        class(nurbs_surface), intent(inout) :: this
        character(len=*),     intent(in) :: filename
        integer, intent(in),  optional   :: res
        character(len=*), intent(in), optional :: encoding

        integer :: ne_u, ne_v, ne_total, np, j, i, m, s, r, o, a, b, t, g, offsetP, line_nodes
        integer :: res_min, dim, N1sp, N2sp, L, N, res1, res2

        real(rk), allocatable :: U1(:), U2(:)
        real(rk), allocatable :: U1r(:), U2r(:)
        real(rk), allocatable :: Xt_all(:,:), Xg_all(:,:)
        integer,  allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        if (.not. allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'Control points are not set.',&
                location   = 'export_Xth_in_Xg',&
                suggestion = 'Call set(...) first before exporting.')
            return
        end if

        if (.not. allocated(this%knot1) .or. .not. allocated(this%knot2)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'Knot vector is not set.',&
                location   = 'export_Xth_in_Xg',&
                suggestion = 'Call set(...) first before exporting.')
            return
        end if

        res_min = 10
        if (present(res)) res_min = max(2, res)

        U1 = unique(this%knot1)
        if (size(U1) < 2) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'knot1 needs >= 2 unique values.',&
                location   = 'export_Xth_in_Xg',&
                suggestion = 'Check the knot vector for sufficient unique values.')
            return
        end if
        U2 = unique(this%knot2)
        if (size(U2) < 2) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'knot2 needs >= 2 unique values.',&
                location   = 'export_Xth_in_Xg',&
                suggestion = 'Check the knot vector for sufficient unique values.')
            return
        end if

        N1sp = size(U1) - 1
        N2sp = size(U2) - 1

        L = N1sp
        if (N2sp > 0) then
            a = L; b = N2sp
            do; t = mod(a,b); if (t==0) exit; a=b; b=t; end do
            g = b
            L = (L/g) * N2sp
        end if

        L = L * max(1, res_min - 1)

        N    = L + 1
        res1 = L / N1sp + 1
        res2 = L / N2sp + 1

        dim = size(this%Xc,2)
        if (dim < 2 .or. dim > 3) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'Invalid geometry dimension.',&
                location   = 'export_Xth_in_Xg',&
                suggestion = 'Check the geometry dimension before exporting the NURBS surface.')
            return
        end if

        ! Allocate refined knot vectors
        allocate(U1r( (size(U1)-1)*(res1-1) + 1 ))
        allocate(U2r( (size(U2)-1)*(res2-1) + 1 ))

        do s = 1, size(U1)-1
            o = (s-1)*(res1-1)
            do r = 1, res1
                U1r(o+r) = U1(s) + (U1(s+1)-U1(s)) * real(r-1,rk)/real(res1-1,rk)
            end do
        end do
        do s = 1, size(U2)-1
            o = (s-1)*(res2-1)
            do r = 1, res2
                U2r(o+r) = U2(s) + (U2(s+1)-U2(s)) * real(r-1,rk)/real(res2-1,rk)
            end do
        end do
        if (size(U1r)/=N .or. size(U2r)/=N) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'Refinement size mismatch.',&
                location   = 'export_Xth_in_Xg',&
                suggestion = 'Check the refinement process for consistency.')
            return
        end if

        ! total element count and node count
        ne_u = size(U2)
        ne_v = size(U1)
        ne_total = size(U2) + size(U1)
        np = ne_total * N

        ! Allocate global arrays
        allocate(Xt_all(np,2), Xg_all(np,dim), elemConn(ne_total, N))

        ! build all parametric points
        offsetP    = 0
        line_nodes = N

        ! dir-1: u varies (v=U2(j), w=U3(k))
        do concurrent (j = 1:size(U2))
            Xt_all(offsetP + (j-1)*line_nodes + 1 : offsetP + j*line_nodes, 1) = U1r
            Xt_all(offsetP + (j-1)*line_nodes + 1 : offsetP + j*line_nodes, 2) = U2(j)
        end do
        offsetP = offsetP + ne_u*line_nodes

        ! dir-2: v varies (u=U1(i), w=U3(k))
        do concurrent (i = 1:size(U1))
            Xt_all(offsetP + (i-1)*line_nodes + 1 : offsetP + i*line_nodes, 1) = U1(i)
            Xt_all(offsetP + (i-1)*line_nodes + 1 : offsetP + i*line_nodes, 2) = U2r
        end do

        ! compute global points
        if (this%is_rational()) then
            Xg_all = compute_Xg(Xt_all, this%knot1, this%knot2, this%degree, this%nc, [np,1], this%Xc, this%Wc)
        else
            Xg_all = compute_Xg(Xt_all, this%knot1, this%knot2, this%degree, this%nc, [np,1], this%Xc)
        end if

        ! connectivity
        do concurrent (l = 1:ne_total, m = 1:N)
            elemConn(l, m) = (l-1)*N + m
        end do

        ! write VTK file
        call export_vtk_legacy(filename=filename, points=Xg_all, elemConn=elemConn, vtkCellType=4, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_iges(this, filename)
        use forIGES, only: Gsection_t, Dentry_t, entity128_t, DElist_t, PElist_t,&
                           makeSsection, makeGsection, makeDPsections, writeIGESfile, wp

        class(nurbs_surface), intent(inout) :: this
        character(len=*),     intent(in)    :: filename

        type(Gsection_t)  :: G
        type(Dentry_t)    :: D
        type(entity128_t) :: surf128
        type(DElist_t)    :: Dlist
        type(PElist_t)    :: Plist
        character(80), allocatable :: Ssection(:), Gsection(:), Dsection(:), Psection(:), Ssec_out(:)
        real(wp) :: X(0:this%degree(1), 0:this%degree(2)), Y(0:this%degree(1), 0:this%degree(2)), Z(0:this%degree(1), 0:this%degree(2)), W(0:this%degree(1), 0:this%degree(2)), S(-this%degree(1):1+this%degree(1)), T(-this%degree(2):1+this%degree(2))
        integer :: i, j, idx
        integer :: K1, K2, M1, M2, N1, N2, prop3
        real(wp) :: U(0:1), V(0:1)

        if (.not. this%err%ok) return

        ! Parameters consistent with the IGES definition
        K1 = this%degree(1)
        K2 = this%degree(2)
        M1 = this%degree(1)
        M2 = this%degree(2)

        ! Compute required N1 and N2 based on IGES standard
        N1 = 1 + K1 - M1
        N2 = 1 + K2 - M2

        ! Copy knots explicitly, matching IGES indexing exactly
        do i = -M1, N1 + K1
            S(i) = real(this%knot1(i + M1 + 1), kind=wp)
        end do

        do i = -M2, N2 + K2
            T(i) = real(this%knot2(i + M2 + 1), kind=wp)
        end do

        ! Correctly map control points and weights
        if (this%is_rational()) then
            do j = 0, K2
                do i = 0, K1
                    idx = j * this%nc(1) + i + 1
                    X(i,j) = real(this%Xc(idx,1), kind=wp)
                    Y(i,j) = real(this%Xc(idx,2), kind=wp)
                    Z(i,j) = real(this%Xc(idx,3), kind=wp)
                    W(i,j) = real(this%Wc(idx), kind=wp)
                end do
            end do
            prop3 = 1  ! Rational surface
        else
            do j = 0, K2
                do i = 0, K1
                    idx = j * this%nc(1) + i + 1
                    X(i,j) = real(this%Xc(idx,1), kind=wp)
                    Y(i,j) = real(this%Xc(idx,2), kind=wp)
                    Z(i,j) = real(this%Xc(idx,3), kind=wp)
                    W(i,j) = real(1.0_rk, kind=wp)
                end do
            end do
            prop3 = 0  ! b-Spline surface
        end if

        U = real([minval(this%knot1), maxval(this%knot1)], kind=wp)
        V = real([minval(this%knot2), maxval(this%knot2)], kind=wp)

        ! Initialize IGES entity 128 (Rational B-spline Surface)
        call surf128%init(&
            DEP   = 1,&
            form  = 0,&
            K1    = K1,&
            K2    = K2,&
            M1    = M1,&
            M2    = M2,&
            PROP1 = 0,&
            PROP2 = 0,&
            PROP3 = prop3,&
            PROP4 = 0,&
            PROP5 = 0,&
            S     = S,&
            T     = T,&
            W     = W,&
            X     = X,&
            Y     = Y,&
            Z     = Z,&
            U     = U,&
            V     = V)

        ! Directory entry
        call D%init(entity_type=128, param_data=1, transformation_matrix=0, form_number=0)

        ! Create entity and directory lists
        call Dlist%init()
        call Plist%init()
        call Dlist%append(D)
        call Plist%append(surf128)

        ! Global section initialization
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
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: X
        integer, intent(in) :: num
        integer, intent(in) :: dir

        if (.not. this%err%ok) return

        if (allocated(this%Xc)) then
            this%Xc(num,dir) = X
            if (allocated(this%Wc)) then
                call this%set(knot1=this%get_knot(1), knot2=this%get_knot(2), Xc=this%get_Xc(), Wc=this%get_Wc())
            else
                call this%set(knot1=this%get_knot(1), knot2=this%get_knot(2), Xc=this%get_Xc())
            end if
        else
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
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
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: W
        integer, intent(in) :: num

        if (.not. this%err%ok) return

        if (allocated(this%Wc)) then
            this%Wc(num) = W
            if (allocated(this%knot1) .and. allocated(this%knot2)) then
                call this%set(knot1=this%get_knot(1), knot2=this%get_knot(2), Xc=this%get_Xc(), Wc=this%get_Wc())
            else
                call this%set(nc=this%nc, Xc=this%get_Xc(), Wc=this%get_Wc())
            end if
        else
            call this%err%set(&
                code       = 105,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
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
    pure function get_multiplicity(this, dir) result(m)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        integer, allocatable :: m(:)

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

        if (present(dir)) then
            if (dir == 1) then
                if (.not.allocated(this%knot1)) then
                    call this%err%set(&
                        code       = 103,&
                        severity   = 1,&
                        category   = 'forcad_nurbs_surface',&
                        message    = 'Knot vector is not set.',&
                        location   = 'cmp_nc',&
                        suggestion = 'Call set(...) first before computing nc.')
                    return
                else
                    this%nc(1) = sum(compute_multiplicity(this%knot1)) - this%degree(1) - 1
                end if
            elseif (dir == 2) then
                if (.not.allocated(this%knot2)) then
                    call this%err%set(&
                        code       = 103,&
                        severity   = 1,&
                        category   = 'forcad_nurbs_surface',&
                        message    = 'Knot vector is not set.',&
                        location   = 'cmp_nc',&
                        suggestion = 'Call set(...) first before computing nc.')
                    return
                else
                    this%nc(2) = sum(compute_multiplicity(this%knot2)) - this%degree(2) - 1
                end if
            else
                call this%err%set(&
                    code       = 103,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_surface',&
                    message    = 'Invalid direction for computing number of control points.',&
                    location   = 'cmp_nc',&
                    suggestion = 'Use dir=1 or dir=2 to specify the direction.')
                return
            end if
        else
            if (.not.allocated(this%knot1)) then
                call this%err%set(&
                    code       = 103,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_surface',&
                    message    = 'Knot vector is not set.',&
                    location   = 'cmp_nc',&
                    suggestion = 'Call set(...) first before computing nc.')
                return
            else
                this%nc(1) = sum(compute_multiplicity(this%knot1)) - this%degree(1) - 1
            end if

            if (.not.allocated(this%knot2)) then
                call this%err%set(&
                    code       = 103,&
                    severity   = 1,&
                    category   = 'forcad_nurbs_surface',&
                    message    = 'Knot vector is not set.',&
                    location   = 'cmp_nc',&
                    suggestion = 'Call set(...) first before computing nc.')
                return
            else
                this%nc(2) = sum(compute_multiplicity(this%knot2)) - this%degree(2) - 1
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_nc_all(this) result(nc)
        class(nurbs_surface), intent(in) :: this
        integer :: nc(2)

        if (.not. this%err%ok) return

        nc = this%nc
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_nc_dir(this, dir) result(nc)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        integer :: nc

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

        ! Set parameter values
        if (present(Xt1)) then
            if (allocated(this%Xt1)) then
                if (size(this%Xt1) /= size(Xt1)) then
                    deallocate(this%Xt1)
                end if
            end if
            this%Xt1 = Xt1
        elseif (present(res1)) then
            if (allocated(this%Xt1)) then
                if (size(this%Xt1) /= res1) then
                    deallocate(this%Xt1)
                    allocate(this%Xt1(res1))
                end if
            else
                allocate(this%Xt1(res1))
            end if
            this%Xt1 = [(this%knot1(1)+(this%knot1(size(this%knot1))-this%knot1(1))*real(i-1,rk)/real(res1-1,rk), i=1, res1)]
            ! else
            ! this%Xt1 = this%Xt1
        end if

        ! Set parameter values
        if (present(Xt2)) then
            if (allocated(this%Xt2)) then
                if (size(this%Xt2) /= size(Xt2)) then
                    deallocate(this%Xt2)
                end if
            end if
            this%Xt2 = Xt2
        elseif (present(res2)) then
            if (allocated(this%Xt2)) then
                if (size(this%Xt2) /= res2) then
                    deallocate(this%Xt2)
                    allocate(this%Xt2(res2))
                end if
            else
                allocate(this%Xt2(res2))
            end if
            this%Xt2 = [(this%knot2(1)+(this%knot2(size(this%knot2))-this%knot2(1))*real(i-1,rk)/real(res2-1,rk), i=1, res2)]
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
    pure subroutine derivative_scalar(this, Xt, dTgc, Tgc, elem)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:)
        integer, intent(in), optional :: elem(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out), optional :: Tgc(:)

        if (.not. this%err%ok) return

        if (this%is_rational()) then ! NURBS
            if (present(elem)) then
                associate(Wce => this%Wc(elem))
                    call compute_dTgc(Xt, this%knot1, this%knot2, this%degree, this%nc, Wce, dTgc, Tgc, elem)
                end associate
            else
                call compute_dTgc(Xt, this%knot1, this%knot2, this%degree, this%nc, this%Wc, dTgc, Tgc)
            end if
        else ! B-Spline
            call compute_dTgc(Xt, this%knot1, this%knot2, this%degree, this%nc, dTgc, Tgc, elem)
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

        if (.not. this%err%ok) return

        ! Set parameter values
        if (present(Xt1)) then
            if (allocated(this%Xt1)) then
                if (size(this%Xt1) /= size(Xt1)) then
                    deallocate(this%Xt1)
                end if
            end if
            this%Xt1 = Xt1
        elseif (present(res1)) then
            if (allocated(this%Xt1)) then
                if (size(this%Xt1) /= res1) then
                    deallocate(this%Xt1)
                    allocate(this%Xt1(res1))
                end if
            else
                allocate(this%Xt1(res1))
            end if
            this%Xt1 = [(this%knot1(1)+(this%knot1(size(this%knot1))-this%knot1(1))*real(i-1,rk)/real(res1-1,rk), i=1, res1)]
            ! else
            ! this%Xt1 = this%Xt1
        end if

        ! Set parameter values
        if (present(Xt2)) then
            if (allocated(this%Xt2)) then
                if (size(this%Xt2) /= size(Xt2)) then
                    deallocate(this%Xt2)
                end if
            end if
            this%Xt2 = Xt2
        elseif (present(res2)) then
            if (allocated(this%Xt2)) then
                if (size(this%Xt2) /= res2) then
                    deallocate(this%Xt2)
                    allocate(this%Xt2(res2))
                end if
            else
                allocate(this%Xt2(res2))
            end if
            this%Xt2 = [(this%knot2(1)+(this%knot2(size(this%knot2))-this%knot2(1))*real(i-1,rk)/real(res2-1,rk), i=1, res2)]
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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

        ! Set parameter values
        if (present(Xt1)) then
            if (allocated(this%Xt1)) then
                if (size(this%Xt1) /= size(Xt1)) then
                    deallocate(this%Xt1)
                end if
            end if
            this%Xt1 = Xt1
        elseif (present(res1)) then
            if (allocated(this%Xt1)) then
                if (size(this%Xt1) /= res1) then
                    deallocate(this%Xt1)
                    allocate(this%Xt1(res1))
                end if
            else
                allocate(this%Xt1(res1))
            end if
            this%Xt1 = [(this%knot1(1)+(this%knot1(size(this%knot1))-this%knot1(1))*real(i-1,rk)/real(res1-1,rk), i=1, res1)]
            ! else
            ! this%Xt1 = this%Xt1
        end if

        ! Set parameter values
        if (present(Xt2)) then
            if (allocated(this%Xt2)) then
                if (size(this%Xt2) /= size(Xt2)) then
                    deallocate(this%Xt2)
                end if
            end if
            this%Xt2 = Xt2
        elseif (present(res2)) then
            if (allocated(this%Xt2)) then
                if (size(this%Xt2) /= res2) then
                    deallocate(this%Xt2)
                    allocate(this%Xt2(res2))
                end if
            else
                allocate(this%Xt2(res2))
            end if
            this%Xt2 = [(this%knot2(1)+(this%knot2(size(this%knot2))-this%knot2(1))*real(i-1,rk)/real(res2-1,rk), i=1, res2)]
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
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), allocatable, intent(out) :: Tgc(:)

        if (.not. this%err%ok) return

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
        integer :: k, i, s, d, j, n_new
        real(rk), allocatable :: Xc(:,:), Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:)
        real(rk), allocatable:: Xc3(:,:,:)

        if (.not. this%err%ok) return

        if (dir == 1) then ! direction 1

            if(this%is_rational()) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(1)-1,this%degree(1),Xth(i),this%knot1)
                    ! if (this%knot1(k+1) == Xth(i)) then
                    if (abs(this%knot1(k+1) - Xth(i)) < 2.0_rk*epsilon(0.0_rk)) then
                        s = compute_multiplicity(this%knot1, Xth(i))
                    else
                        s = 0
                    end if

                    d = size(this%Xc,2)
                    allocate(Xcw(size(this%Xc,1),d+1))
                    do j = 1, size(this%Xc,1)
                        Xcw(j,1:d) = this%Xc(j,1:d)*this%Wc(j)
                    end do
                    Xcw(:,d+1) = this%Wc(:)

                    Xcw = reshape(Xcw,[this%nc(1),this%nc(2)*(d+1)])

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

                    Xcw_new = reshape(Xcw_new,[this%nc(2)*(n_new+1),d+1])

                    allocate(Xc_new(1:this%nc(2)*(n_new+1),1:d))
                    allocate(Wc_new(1:this%nc(2)*(n_new+1)))
                    do j = 1, this%nc(2)*(n_new+1)
                        Xc_new(j,1:d) = Xcw_new(j,1:d)/Xcw_new(j,d+1)
                    end do
                    Wc_new(:) = Xcw_new(:,d+1)

                    call this%set(knot1=knot_new, knot2=this%get_knot(2), Xc=Xc_new, Wc=Wc_new)
                    deallocate(Xcw, Xcw_new, Xc_new, Wc_new)
                end do


            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(1)-1,this%degree(1),Xth(i),this%knot1)
                    ! if (this%knot1(k+1) == Xth(i)) then
                    if (abs(this%knot1(k+1) - Xth(i)) < 2.0_rk*epsilon(0.0_rk)) then
                        s = compute_multiplicity(this%knot1, Xth(i))
                    else
                        s = 0
                    end if

                    d = size(this%Xc,2)

                    Xc = reshape(this%Xc,[this%nc(1),this%nc(2)*d])

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

                    Xc_new = reshape(Xc_new,[(this%nc(2))*(n_new+1),d])

                    call this%set(knot1=knot_new, knot2=this%get_knot(2), Xc=Xc_new)
                end do

            end if


        elseif (dir == 2) then! direction 2

            if(this%is_rational()) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(2)-1,this%degree(2),Xth(i),this%knot2)
                    ! if (this%knot2(k+1) == Xth(i)) then
                    if (abs(this%knot2(k+1) - Xth(i)) < 2.0_rk*epsilon(0.0_rk)) then
                        s = compute_multiplicity(this%knot2, Xth(i))
                    else
                        s = 0
                    end if

                    d = size(this%Xc,2)
                    allocate(Xcw(size(this%Xc,1),d+1))
                    do j = 1, size(this%Xc,1)
                        Xcw(j,1:d) = this%Xc(j,1:d)*this%Wc(j)
                    end do
                    Xcw(:,d+1) = this%Wc(:)

                    Xc3 = reshape(Xcw, [this%nc(1),this%nc(2),d+1])
                    Xc3 = reshape(Xc3, [this%nc(2),this%nc(1),d+1], order=[2,1,3])
                    Xcw = reshape(Xc3,[this%nc(2),this%nc(1)*(d+1)])

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

                    Xc3 = reshape(Xcw_new, [n_new+1,this%nc(1),d+1])
                    Xc3 = reshape(Xc3, [this%nc(1),n_new+1,d+1], order=[2,1,3])
                    Xcw_new = reshape(Xc3,[(this%nc(1))*(n_new+1),d+1])

                    allocate(Xc_new(1:(n_new+1)*this%nc(1),1:d))
                    allocate(Wc_new(1:(n_new+1)*this%nc(1)))
                    do j = 1, (n_new+1)*this%nc(1)
                        Xc_new(j,1:d) = Xcw_new(j,1:d)/Xcw_new(j,d+1)
                    end do
                    Wc_new(:) = Xcw_new(:,d+1)

                    call this%set(knot1=this%get_knot(1), knot2=knot_new, Xc=Xc_new, Wc=Wc_new)
                    deallocate(Xcw, Xcw_new, Xc_new, Wc_new)
                end do

            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(2)-1,this%degree(2),Xth(i),this%knot2)
                    ! if (this%knot2(k+1) == Xth(i)) then
                    if (abs(this%knot2(k+1) - Xth(i)) < 2.0_rk*epsilon(0.0_rk)) then
                        s = compute_multiplicity(this%knot2, Xth(i))
                    else
                        s = 0
                    end if

                    d = size(this%Xc,2)

                    Xc3 = reshape(this%Xc, [this%nc(1),this%nc(2),d])
                    Xc3 = reshape(Xc3, [this%nc(2),this%nc(1),d], order=[2,1,3])
                    Xc = reshape(Xc3,[this%nc(2),this%nc(1)*d])

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

                    Xc3 = reshape(Xc_new, [n_new+1,this%nc(1),d])
                    Xc3 = reshape(Xc3, [this%nc(1),n_new+1,d], order=[2,1,3])
                    Xc_new = reshape(Xc3,[(this%nc(1))*(n_new+1),d])

                    call this%set(knot1=this%get_knot(1), knot2=knot_new, Xc=Xc_new)
                end do


            end if

        else
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'Invalid direction for inserting knots.',&
                location   = 'insert_knots',&
                suggestion = 'Use dir=1 or dir=2 to specify the direction.')
            return
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
        integer :: d, j, nc_new
        real(rk), allocatable:: Xc3(:,:,:)

        if (.not. this%err%ok) return

        if (dir == 1) then ! direction 1

            if(this%is_rational()) then ! NURBS

                d = size(this%Xc,2)
                allocate(Xcw(size(this%Xc,1),d+1))
                do j = 1, size(this%Xc,1)
                    Xcw(j,1:d) = this%Xc(j,1:d)*this%Wc(j)
                end do
                Xcw(:,d+1) = this%Wc(:)

                Xcw = reshape(Xcw,[this%nc(1),this%nc(2)*(d+1)],order=[1,2])

                call elevate_degree_A_5_9(t, this%knot1, this%degree(1), Xcw, nc_new, knot_new, Xcw_new)

                Xcw_new = reshape(Xcw_new,[this%nc(2)*nc_new,d+1],order=[1,2])

                allocate(Xc_new(1:this%nc(2)*nc_new,1:d))
                allocate(Wc_new(1:this%nc(2)*nc_new))
                do j = 1, this%nc(2)*nc_new
                    Xc_new(j,1:d) = Xcw_new(j,1:d)/Xcw_new(j,d+1)
                end do

                Wc_new(:) = Xcw_new(:,d+1)

                call this%set(knot1=knot_new, knot2=this%get_knot(2), Xc=Xc_new, Wc=Wc_new)
                deallocate(Xcw, Xcw_new, Xc_new, Wc_new)

            else ! B-Spline

                d = size(this%Xc,2)
                Xc = reshape(this%Xc,[this%nc(1),this%nc(2)*(d)],order=[1,2])

                call elevate_degree_A_5_9(t, this%knot1, this%degree(1), Xc, nc_new, knot_new, Xc_new)

                Xc_new = reshape(Xc_new,[this%nc(2)*nc_new,d],order=[1,2])

                call this%set(knot1=knot_new, knot2=this%get_knot(2), Xc=Xc_new)
                deallocate(Xc, Xc_new)

            end if

        elseif (dir == 2) then ! direction 2

            if(this%is_rational()) then ! NURBS

                d = size(this%Xc,2)
                allocate(Xcw(size(this%Xc,1),d+1))
                do j = 1, size(this%Xc,1)
                    Xcw(j,1:d) = this%Xc(j,1:d)*this%Wc(j)
                end do

                Xcw(:,d+1) = this%Wc(:)

                Xc3 = reshape(Xcw, [this%nc(1),this%nc(2),d+1])
                Xc3 = reshape(Xc3, [this%nc(2),this%nc(1),d+1], order=[2,1,3])
                Xcw = reshape(Xc3,[this%nc(2),this%nc(1)*(d+1)])

                call elevate_degree_A_5_9(t, this%knot2, this%degree(2), Xcw, nc_new, knot_new, Xcw_new)

                Xc3 = reshape(Xcw_new, [nc_new,this%nc(1),d+1])
                Xc3 = reshape(Xc3, [this%nc(1),nc_new,d+1], order=[2,1,3])
                Xcw_new = reshape(Xc3,[(this%nc(1))*nc_new,d+1])

                allocate(Xc_new(1:nc_new*this%nc(1),1:d))
                allocate(Wc_new(1:nc_new*this%nc(1)))
                do j = 1, nc_new*this%nc(1)
                    Xc_new(j,1:d) = Xcw_new(j,1:d)/Xcw_new(j,d+1)
                end do

                Wc_new(:) = Xcw_new(:,d+1)

                call this%set(knot1=this%get_knot(1), knot2=knot_new, Xc=Xc_new, Wc=Wc_new)
                deallocate(Xcw, Xcw_new, Xc_new, Wc_new)

            else ! B-Spline

                d = size(this%Xc,2)

                Xc3 = reshape(this%Xc, [this%nc(1),this%nc(2),d])
                Xc3 = reshape(Xc3, [this%nc(2),this%nc(1),d], order=[2,1,3])
                Xc = reshape(Xc3,[this%nc(2),this%nc(1)*d])

                call elevate_degree_A_5_9(t, this%knot2, this%degree(2), Xc, nc_new, knot_new, Xc_new)

                Xc3 = reshape(Xc_new, [nc_new,this%nc(1),d])
                Xc3 = reshape(Xc3, [this%nc(1),nc_new,d], order=[2,1,3])
                Xc_new = reshape(Xc3,[(this%nc(1))*nc_new,d])

                call this%set(knot1=this%get_knot(1), knot2=knot_new, Xc=Xc_new)

            end if

        else
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'Invalid direction for elevating degree.',&
                location   = 'elevate_degree',&
                suggestion = 'Use dir=1 or dir=2 to specify the direction.')
            return
        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function is_rational(this) result(r)
        class(nurbs_surface), intent(in) :: this
        logical :: r

        if (.not. this%err%ok) return

        r = .false.
        if(allocated(this%Wc)) then
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
        class(nurbs_surface), intent(inout) :: this
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
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (.not. this%err%ok) return

        if (allocated(this%elemConn_Xc_vis)) then
            if (size(this%elemConn_Xc_vis,1) /= size(elemConn,1) .or. size(this%elemConn_Xc_vis,2) /= size(elemConn,2)) then
                deallocate(this%elemConn_Xc_vis)
            end if
        end if
        this%elemConn_Xg_vis = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_elem(this, elemConn)
        class(nurbs_surface), intent(inout) :: this
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
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        elemConn = this%elemConn_Xc_vis
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem_Xg_vis(this) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        elemConn = this%elemConn_Xg_vis
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem(this) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

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
        integer :: k, i, s, d, j, nc_new, t
        real(rk), allocatable :: Xc(:,:), Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:)
        real(rk), allocatable:: Xc3(:,:,:)

        if (.not. this%err%ok) return

        if (dir == 1) then ! direction 1

            if(this%is_rational()) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(1)-1,this%degree(1),Xth(i),this%knot1)
                    ! if (this%knot1(k+1) == Xth(i)) then
                    if (abs(this%knot1(k+1) - Xth(i)) < 2.0_rk*epsilon(0.0_rk)) then
                        s = compute_multiplicity(this%knot1, Xth(i))
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

                    Xcw = reshape(Xcw,[this%nc(1),this%nc(2)*(d+1)],order=[1,2])

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
                        Xcw_new = reshape(Xcw_new,[this%nc(2)*(nc_new),d+1],order=[1,2])

                        allocate(Xc_new(1:this%nc(2)*(nc_new),1:d))
                        allocate(Wc_new(1:this%nc(2)*(nc_new)))
                        do j = 1, this%nc(2)*(nc_new)
                            Xc_new(j,1:d) = Xcw_new(j,1:d)/Xcw_new(j,d+1)
                        end do

                        Wc_new(:) = Xcw_new(:,d+1)

                        call this%set(knot1=knot_new, knot2=this%get_knot(2), Xc=Xc_new, Wc=Wc_new)
                        deallocate(Xcw_new, Xc_new, Wc_new)
                    end if
                end do


            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(1)-1,this%degree(1),Xth(i),this%knot1)
                    ! if (this%knot1(k+1) == Xth(i)) then
                    if (abs(this%knot1(k+1) - Xth(i)) < 2.0_rk*epsilon(0.0_rk)) then
                        s = compute_multiplicity(this%knot1, Xth(i))
                    else
                        s = 0
                    end if
                    k = k + 1

                    d = size(this%Xc,2)

                    Xc = reshape(this%Xc,[this%nc(1),this%nc(2)*d],order=[1,2])

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
                        Xc_new = reshape(Xc_new,[(this%nc(2))*(nc_new),d],order=[1,2])

                        call this%set(knot1=knot_new, knot2=this%get_knot(2), Xc=Xc_new)
                    end if
                end do

            end if


        elseif (dir == 2) then! direction 2

            if(this%is_rational()) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(2)-1,this%degree(2),Xth(i),this%knot2)
                    ! if (this%knot2(k+1) == Xth(i)) then
                    if (abs(this%knot2(k+1) - Xth(i)) < 2.0_rk*epsilon(0.0_rk)) then
                        s = compute_multiplicity(this%knot2, Xth(i))
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

                    Xc3 = reshape(Xcw, [this%nc(1),this%nc(2),d+1])
                    Xc3 = reshape(Xc3, [this%nc(2),this%nc(1),d+1], order=[2,1,3])
                    Xcw = reshape(Xc3, [this%nc(2),this%nc(1)*(d+1)])

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
                        Xc3 = reshape(Xcw_new, [nc_new,this%nc(1),d+1])
                        Xc3 = reshape(Xc3, [this%nc(1),nc_new,d+1], order=[2,1,3])
                        Xcw_new = reshape(Xc3,[(this%nc(1))*(nc_new),d+1])

                        allocate(Xc_new(1:(nc_new)*this%nc(1),1:d))
                        allocate(Wc_new(1:(nc_new)*this%nc(1)))
                        do j = 1, (nc_new)*this%nc(1)
                            Xc_new(j,1:d) = Xcw_new(j,1:d)/Xcw_new(j,d+1)
                        end do

                        Wc_new(:) = Xcw_new(:,d+1)

                        call this%set(knot1=this%get_knot(1), knot2=knot_new, Xc=Xc_new, Wc=Wc_new)
                        deallocate(Xcw_new, Xc_new, Wc_new)
                    end if

                end do

            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(2)-1,this%degree(2),Xth(i),this%knot2)
                    ! if (this%knot2(k+1) == Xth(i)) then
                    if (abs(this%knot2(k+1) - Xth(i)) < 2.0_rk*epsilon(0.0_rk)) then
                        s = compute_multiplicity(this%knot2, Xth(i))
                    else
                        s = 0
                    end if
                    k = k + 1

                    d = size(this%Xc,2)

                    Xc3 = reshape(this%Xc, [this%nc(1),this%nc(2),d])
                    Xc3 = reshape(Xc3, [this%nc(2),this%nc(1),d], order=[2,1,3])
                    Xc = reshape(Xc3,[this%nc(2),this%nc(1)*d])

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

                        Xc3 = reshape(Xc_new, [nc_new,this%nc(1),d])
                        Xc3 = reshape(Xc3, [this%nc(1),nc_new,d], order=[2,1,3])
                        Xc_new = reshape(Xc3,[(this%nc(1))*(nc_new),d])

                        call this%set(knot1=this%get_knot(1), knot2=knot_new, Xc=Xc_new)
                    end if

                end do


            end if

        else
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'Invalid direction for removing knots.',&
                location   = 'remove_knots',&
                suggestion = 'Use dir=1 or dir=2 to specify the direction.')
            return
        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_tetragon(this, L, nc, Wc)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: L(:)
        integer, intent(in), contiguous :: nc(:)
        real(rk), intent(in), contiguous, optional :: Wc(:)

        if (.not. this%err%ok) return

        call this%set(nc = nc, Xc = tetragon_Xc(L, nc), Wc = Wc)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elem(this) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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
#ifndef NOSHOW_PYVISTA
        block
        character(len=3000) :: pyvista_script

        if (.not. this%err%ok) return

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
    pure subroutine set_ring(this, center, radius1, radius2)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: center(:)
        real(rk), intent(in) :: radius1, radius2
        real(rk), allocatable :: Xc(:,:), Wc(:), knot1(:), knot2(:)
        integer :: i

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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

        if (.not. this%err%ok) return

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
        real(rk), intent(in), contiguous :: point_Xg(:)
        real(rk), intent(out), optional :: nearest_Xg(size(point_Xg))
        real(rk), intent(out), optional :: nearest_Xt(2)
        integer, intent(out), optional :: id
        integer :: id_, i
        real(rk), allocatable :: distances(:)

        if (.not. this%err%ok) return

        allocate(distances(this%ng(1)*this%ng(2)))

#if defined(__NVCOMPILER)
        do i = 1, this%ng(1)*this%ng(2)
#else
        do concurrent (i = 1: this%ng(1)*this%ng(2))
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
        if (present(nearest_Xt)) nearest_Xt = this%Xt(id_,:)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine nearest_point2(this, point_Xg, tol, maxit, nearest_Xt, nearest_Xg)

        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: point_Xg(:)
        real(rk), intent(in) :: tol
        integer, intent(in) :: maxit
        real(rk), intent(out) :: nearest_Xt(2)
        real(rk), intent(out), optional :: nearest_Xg(size(this%Xc,2))

        real(rk) :: obj, obj_trial, grad(2), hess(2,2), dk(2)
        real(rk) :: alphak, alpha_max, alpha_i, tau, beta, eps
        real(rk) :: lower_bounds(2), upper_bounds(2), xt(2)
        real(rk), allocatable :: Tgc(:), dTgc(:,:), d2Tgc(:,:)
        real(rk) :: Xg(size(this%Xc,2)), xk(2), xkn(2)
        integer :: k, l, i
        logical :: convergenz
        type(nurbs_surface) :: copy_this

        if (.not. this%err%ok) return

        alphak = 0.0_rk
        dk     = 0.0_rk
        k      = 0
        eps    = 10.0_rk*tiny(1.0_rk)

        ! bounds
        lower_bounds = [minval(this%knot1), minval(this%knot2)]
        upper_bounds = [maxval(this%knot1), maxval(this%knot2)]

        ! initial guess (coarse search)
        copy_this = this
        call copy_this%create(10,10)
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

            obj = norm2(Xg - point_Xg) + 0.001_rk  ! small epsilon to avoid divide-by-zero

            grad(1) = dot_product((Xg-point_Xg)/obj, matmul(dTgc(:,1), this%Xc))
            grad(2) = dot_product((Xg-point_Xg)/obj, matmul(dTgc(:,2), this%Xc))

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

            if (norm2(grad) <= tol .or. (k>0 .and. norm2(xk-xkn) <= tol)) then
                convergenz  = .true.
                nearest_Xt  = xk
                if (present(nearest_Xg)) nearest_Xg = this%cmp_Xg(nearest_Xt)
            else
                ! Newton step
                dk = - matmul(inv(hess), grad)

                ! Backtracking-Armijo with feasibility (box constraints)
                tau  = 0.5_rk
                beta = 1.0e-4_rk

                ! compute maximum feasible step so xk + alpha*dk stays in [lower_bounds, upper_bounds]
                alpha_max = 1.0_rk
                do i = 1, 2
                    if (dk(i) > 0.0_rk) then
                        if (upper_bounds(i) > xk(i)) then
                            alpha_i = (upper_bounds(i) - xk(i)) / dk(i)
                            alpha_max = min(alpha_max, max(0.0_rk, alpha_i))
                        else
                            alpha_max = 0.0_rk
                        end if
                    else if (dk(i) < 0.0_rk) then
                        if (lower_bounds(i) < xk(i)) then
                            alpha_i = (lower_bounds(i) - xk(i)) / dk(i)
                            alpha_max = min(alpha_max, max(0.0_rk, alpha_i))
                        else
                            alpha_max = 0.0_rk
                        end if
                    end if
                end do

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
                    if (obj_trial <= obj + alphak*beta*dot_product(grad, dk)) exit
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
    pure subroutine ansatz(this, ie, ig, Tgc, dTgc_dXg, dA, ngauss)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in) :: ie, ig
        real(rk), intent(out) :: dA
        real(rk), allocatable, intent(out) :: Tgc(:), dTgc_dXg(:,:)
        integer, intent(in), optional :: ngauss(2)
        real(rk), allocatable :: Xth(:,:), Xth_e(:,:), Xth_eT(:,:), Xc_eT(:,:), Xth1(:), Xth2(:), Xksi(:,:), Wksi(:)
        integer, allocatable :: elem_th(:,:), elem_c(:,:), elem_ce(:)
        type(nurbs_surface) :: th, th_e
        real(rk), allocatable :: dTtth_dXksi(:,:), Ttth(:), dTgc_dXt(:,:), Xt(:), dXt_dXksi(:,:), dXg_dXt(:,:)
        real(rk), allocatable :: dXg_dXksi(:,:) !! Jacobian matrix
        real(rk) :: det_dXg_dXksi !! Determinant of the Jacobian matrix
        real(rk) :: Xksii(2)

        if (.not. this%err%ok) return

        if (present(ngauss)) then
            call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], ngauss-1, Xksi, Wksi)
        else
            call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], this%degree, Xksi, Wksi)
        end if

        Xth1 = unique(this%knot1)
        Xth2 = unique(this%knot2)
        call ndgrid(Xth1, Xth2, Xth)

        call th%set([0.0_rk,Xth1,1.0_rk], [0.0_rk,Xth2,1.0_rk], Xth)
        elem_th = th%cmp_elem()

        elem_c = this%cmp_elem()

        Xth_e = Xth(elem_th(ie,:),:)
        call th_e%set([0.0_rk,0.0_rk,1.0_rk,1.0_rk], [0.0_rk,0.0_rk,1.0_rk,1.0_rk], Xth_e)

        Xth_eT = transpose(Xth_e)
        elem_ce = elem_c(ie,:)
        Xc_eT = transpose(this%Xc(elem_ce,:))

        Xksii = Xksi(ig,:)
        call th_e%derivative(Xksii, dTtth_dXksi, Ttth)
        Xt = matmul(Xth_eT, Ttth)
        dXt_dXksi = matmul(Xth_eT, dTtth_dXksi)

        call this%derivative(Xt, dTgc_dXt, Tgc, elem_ce)
        dXg_dXt = matmul(Xc_eT, dTgc_dXt)

        dTgc_dXg = matmul(dTgc_dXt, inv(dXg_dXt))

        dXg_dXksi = matmul(dXg_dXt, dXt_dXksi)
        det_dXg_dXksi = det(dXg_dXksi)

        dA = det_dXg_dXksi*Wksi(ig)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine cmp_area(this, area, ngauss)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(out) :: area
        integer, intent(in), optional :: ngauss(2)
        real(rk), allocatable :: Tgc(:), dTgc_dXg(:,:)
        integer :: ie, ig
        integer :: ngauss_(2)
        real(rk) :: dA, dA_ig

        if (.not. this%err%ok) return

        if (present(ngauss)) then
            ngauss_ = ngauss
        else
            ngauss_ = this%degree + 1
        end if

        area = 0.0_rk
#if defined(__NVCOMPILER)
        do ie = 1, size(this%cmp_elem(),1)
#else
        do concurrent (ie = 1: size(this%cmp_elem(),1)) reduce(+:area)
#endif
            dA = 0.0_rk
            do ig = 1, product(ngauss_)
                call this%ansatz(ie, ig, Tgc, dTgc_dXg, dA_ig, ngauss_)
                dA = dA + dA_ig
            end do
            area = area + dA
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_Tgc_2d(Xti, knot1, knot2, nc, degree, Wc) result(Tgc)
        real(rk), intent(in), contiguous :: Xti(:)
        real(rk), intent(in), contiguous :: knot1(:)
        real(rk), intent(in), contiguous :: knot2(:)
        integer, intent(in), contiguous :: degree(:), nc(:)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk) :: Tgc(nc(1)*nc(2))
        real(rk) :: tmp
        integer :: i

        Tgc = kron(&
            basis_bspline(Xti(2), knot2, nc(2), degree(2)),&
            basis_bspline(Xti(1), knot1, nc(1), degree(1)))
        tmp = dot_product(Tgc, Wc)
        do concurrent (i = 1: nc(1)*nc(2))
           Tgc(i) = (Tgc(i) * Wc(i)) / tmp
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_Xg_nurbs_2d(Xt, knot1, knot2, degree, nc, ng, Xc, Wc) result(Xg)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in), contiguous :: degree(:)
        integer, intent(in), contiguous :: nc(:)
        integer, intent(in), contiguous :: ng(:)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable :: Xg(:,:)
        integer :: i

        allocate(Xg(ng(1)*ng(2), size(Xc,2)))
#if defined(__NVCOMPILER)
        do i = 1, ng(1)*ng(2)
#else
        do concurrent (i = 1: ng(1)*ng(2))
#endif
            Xg(i,:) = matmul(cmp_Tgc_2d(Xt(i,:), knot1, knot2, nc, degree, Wc), Xc)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_Xg_nurbs_2d_1point(Xt, knot1, knot2, degree, nc, Xc, Wc) result(Xg)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk) :: Xg(size(Xc,2))
        real(rk), allocatable :: Tgc(:)

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
    pure function compute_Xg_bspline_2d(Xt, knot1, knot2, degree, nc, ng, Xc) result(Xg)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in) :: ng(2)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), allocatable :: Xg(:,:)
        integer :: i

        allocate(Xg(ng(1)*ng(2), size(Xc,2)))
#if defined(__NVCOMPILER)
        do i = 1, ng(1)*ng(2)
#else
        do concurrent (i = 1: ng(1)*ng(2))
#endif
            Xg(i,:) = matmul(kron(&
                basis_bspline(Xt(i,2), knot2, nc(2), degree(2)),&
                basis_bspline(Xt(i,1), knot1, nc(1), degree(1))),&
                Xc)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_Xg_bspline_2d_1point(Xt, knot1, knot2, degree, nc, Xc) result(Xg)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk) :: Xg(size(Xc,2))

        Xg = matmul(kron(&
            basis_bspline(Xt(2), knot2, nc(2), degree(2)),&
            basis_bspline(Xt(1), knot1, nc(1), degree(1))),&
            Xc)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine compute_dTgc_nurbs_2d_vector(Xt, knot1, knot2, degree, nc, ng, Wc, dTgc, Tgc)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in) :: ng(2)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        real(rk) :: dB1(nc(1)), dB2(nc(2))
        real(rk) :: B1(nc(1)), B2(nc(2))
        real(rk) :: dBi(nc(1)*nc(2), 2), Bi(nc(1)*nc(2))
        integer :: i

        allocate(dTgc(ng(1)*ng(2), nc(1)*nc(2), 2), Tgc(ng(1)*ng(2), nc(1)*nc(2)))
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, size(Xt, 1)
#else
        do concurrent (i = 1: size(Xt, 1)) local(B1, B2, dB1, dB2, Bi, dBi)
#endif
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
    !> If `elem` is not present: `Wc` refers to the full weight vector.
    !> If `elem` is present:     `Wc` refers to the element-local weight vector (`Wce`).
    pure subroutine compute_dTgc_nurbs_2d_scalar(Xt, knot1, knot2, degree, nc, Wc, dTgc, Tgc, elem)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        real(rk), intent(in), contiguous :: Wc(:)
        integer, intent(in), optional :: elem(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        real(rk) :: dB1(nc(1)), dB2(nc(2))
        real(rk) :: B1(nc(1)), B2(nc(2))
        real(rk), allocatable :: dBi(:,:), Bi(:)


        call basis_bspline_der(Xt(1), knot1, nc(1), degree(1), dB1, B1)
        call basis_bspline_der(Xt(2), knot2, nc(2), degree(2), dB2, B2)

        if (.not. present(elem)) then
            allocate(dTgc(nc(1)*nc(2), 2), Tgc(nc(1)*nc(2)))
            allocate(dBi(nc(1)*nc(2), 2), Bi(nc(1)*nc(2)))

            Bi = kron(B2, B1)
            Tgc = Bi*(Wc/(dot_product(Bi,Wc)))

            dBi(:,1) = kron(B2, dB1)
            dBi(:,2) = kron(dB2, B1)

            dTgc(:,1) = ( dBi(:,1)*Wc - Tgc*dot_product(dBi(:,1),Wc) ) / dot_product(Bi,Wc)
            dTgc(:,2) = ( dBi(:,2)*Wc - Tgc*dot_product(dBi(:,2),Wc) ) / dot_product(Bi,Wc)
        else
            allocate(dTgc(size(elem), 2), Tgc(size(elem)))
            allocate(dBi(size(elem), 2), Bi(size(elem)))

            associate(Biall => kron(B2, B1))
                Bi = Biall(elem)
                Tgc = Bi*(Wc/(dot_product(Bi,Wc)))
            end associate

            associate(dB1all => kron(B2, dB1), dB2all => kron(dB2, B1))
                dBi(:,1) = dB1all(elem)
                dBi(:,2) = dB2all(elem)
            end associate

            dTgc(:,1) = ( dBi(:,1)*Wc - Tgc*dot_product(dBi(:,1),Wc) ) / dot_product(Bi,Wc)
            dTgc(:,2) = ( dBi(:,2)*Wc - Tgc*dot_product(dBi(:,2),Wc) ) / dot_product(Bi,Wc)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine compute_dTgc_bspline_2d_vector(Xt, knot1, knot2, degree, nc, ng, dTgc, Tgc)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in) :: ng(2)
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        real(rk) :: dB1(nc(1)), dB2(nc(2))
        real(rk) :: B1(nc(1)), B2(nc(2))
        integer :: i

        allocate(dTgc(ng(1)*ng(2), nc(1)*nc(2), 2), Tgc(ng(1)*ng(2), nc(1)*nc(2)))

#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, size(Xt, 1)
#else
        do concurrent (i = 1: size(Xt, 1)) local(B1, B2, dB1, dB2)
#endif
            call basis_bspline_der(Xt(i,1), knot1, nc(1), degree(1), dB1, B1)
            call basis_bspline_der(Xt(i,2), knot2, nc(2), degree(2), dB2, B2)
            Tgc(i,:) = kron(B2, B1)

            dTgc(i,:,1) = kron(B2, dB1)
            dTgc(i,:,2) = kron(dB2, B1)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine compute_dTgc_bspline_2d_scalar(Xt, knot1, knot2, degree, nc, dTgc, Tgc, elem)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in), optional :: elem(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        real(rk) :: dTgc1(nc(1)), dTgc2(nc(2))
        real(rk) :: Tgc1(nc(1)), Tgc2(nc(2))

        call basis_bspline_der(Xt(1), knot1, nc(1), degree(1), dTgc1, Tgc1)
        call basis_bspline_der(Xt(2), knot2, nc(2), degree(2), dTgc2, Tgc2)

        if (.not. present(elem)) then
            allocate(dTgc(nc(1)*nc(2), 2))
            Tgc = kron(Tgc2, Tgc1)

            dTgc(:,1) = kron(Tgc2, dTgc1)
            dTgc(:,2) = kron(dTgc2, Tgc1)
        else
            allocate(dTgc(size(elem), 2))

            associate(B => kron(Tgc2, Tgc1))
                Tgc = B(elem)
            end associate

            associate(dB1 => kron(Tgc2, dTgc1), dB2 => kron(dTgc2, Tgc1))
                dTgc(:,1) = dB1(elem)
                dTgc(:,2) = dB2(elem)
            end associate
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine compute_d2Tgc_nurbs_2d_vector(Xt, knot1, knot2, degree, nc, ng, Wc, d2Tgc, dTgc, Tgc)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in) :: ng(2)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:,:)
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        real(rk) :: d2B1(nc(1)), d2B2(nc(2))
        real(rk) :: dB1(nc(1)), dB2(nc(2))
        real(rk) :: B1(nc(1)), B2(nc(2))
        real(rk), allocatable :: Tgci(:), dTgci(:)
        real(rk) :: d2Bi(2*nc(1)*nc(2),2), dBi(nc(1)*nc(2),2), Bi(nc(1)*nc(2))

        integer :: i

        allocate(d2Tgc(ng(1)*ng(2), 2*nc(1)*nc(2), 2))

        ! allocate(Bi(nc(1)*nc(2)), dBi(nc(1)*nc(2), 2), d2Bi(2*nc(1)*nc(2), 2))

        allocate(Tgci(nc(1)*nc(2)), dTgci(nc(1)*nc(2)))
        allocate(Tgc(ng(1)*ng(2), nc(1)*nc(2)), dTgc(ng(1)*ng(2), nc(1)*nc(2), 2))
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, size(Xt, 1)
#else
        do concurrent (i = 1: size(Xt, 1)) local(B1, B2, dB1, dB2, Bi, dBi, d2Bi)
#endif
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

            d2Tgc(i,1:nc(1)*nc(2) ,1) = &
                (d2Bi(1:nc(1)*nc(2) ,1)*Wc - 2.0_rk*dTgc(i,:,1)*dot_product(dBi(:,1),Wc) &
                - Tgc(i,:)*dot_product(d2Bi(1:nc(1)*nc(2) ,1),Wc)) / dot_product(Bi,Wc)
            d2Tgc(i,nc(1)*nc(2)+1:2*nc(1)*nc(2),1) = &
                (d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),1)*Wc &
                - dTgc(i,:,1)*dot_product(dBi(:,2),Wc) - dTgc(i,:,2)*dot_product(dBi(:,1),Wc) &
                - Tgc(i,:)*dot_product(d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),1),Wc)) / dot_product(Bi,Wc)
            d2Tgc(i,1:nc(1)*nc(2) ,2) = &
                (d2Bi(1:nc(1)*nc(2) ,2)*Wc - dTgc(i,:,1)*dot_product(dBi(:,2),Wc) - dTgc(i,:,2)*dot_product(dBi(:,1),Wc) &
                - Tgc(i,:)*dot_product(d2Bi(1:nc(1)*nc(2) ,2),Wc)) / dot_product(Bi,Wc)
            d2Tgc(i,nc(1)*nc(2)+1:2*nc(1)*nc(2),2) = &
                (d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),2)*Wc - 2.0_rk*dTgc(i,:,2)*dot_product(dBi(:,2),Wc) &
                - Tgc(i,:)*dot_product(d2Bi(nc(1)*nc(2)+1:2*nc(1)*nc(2),2),Wc)) / dot_product(Bi,Wc)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine compute_d2Tgc_nurbs_2d_scalar(Xt, knot1, knot2, degree, nc, Wc, d2Tgc, dTgc, Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        real(rk) :: d2B1(nc(1)), d2B2(nc(2))
        real(rk) :: dB1(nc(1)), dB2(nc(2))
        real(rk) :: B1(nc(1)), B2(nc(2))
        real(rk), allocatable :: d2Bi(:,:), dBi(:,:), Bi(:)

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
    pure subroutine compute_d2Tgc_bspline_2d_vector(Xt, knot1, knot2, degree, nc, ng, d2Tgc, dTgc, Tgc)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in) :: ng(2)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:,:)
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        real(rk) :: d2B1(nc(1)), d2B2(nc(2))
        real(rk) :: dB1(nc(1)), dB2(nc(2))
        real(rk) :: B1(nc(1)), B2(nc(2))
        integer :: i

        allocate(d2Tgc(ng(1)*ng(2), 2*nc(1)*nc(2), 2))
        allocate(dTgc(ng(1)*ng(2), nc(1)*nc(2), 2))
        allocate(Tgc(ng(1)*ng(2), nc(1)*nc(2)))
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, size(Xt, 1)
#else
        do concurrent (i = 1: size(Xt, 1)) local(B1, B2, dB1, dB2, d2B1, d2B2)
#endif
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
    pure subroutine compute_d2Tgc_bspline_2d_scalar(Xt, knot1, knot2, degree, nc, d2Tgc, dTgc, Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        real(rk) :: d2B1(nc(1)), d2B2(nc(2))
        real(rk) :: dB1(nc(1)), dB2(nc(2))
        real(rk) :: B1(nc(1)), B2(nc(2))

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
    pure function compute_Tgc_nurbs_2d_vector(Xt, knot1, knot2, degree, nc, ng, Wc) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in) :: ng(2)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable :: Tgc(:,:)
        real(rk) :: Tgci(nc(1)*nc(2))
        integer :: i

        allocate(Tgc(ng(1)*ng(2), nc(1)*nc(2)))
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, size(Xt, 1)
#else
        do concurrent (i = 1: size(Xt, 1)) local(Tgci)
#endif
            Tgci = kron(&
                basis_bspline(Xt(i,2), knot2, nc(2), degree(2)),&
                basis_bspline(Xt(i,1), knot1, nc(1), degree(1)))
            Tgc(i,:) = Tgci*(Wc/(dot_product(Tgci,Wc)))
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_Tgc_nurbs_2d_scalar(Xt, knot1, knot2, degree, nc, Wc) result(Tgc)
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
    pure function compute_Tgc_bspline_2d_vector(Xt, knot1, knot2, degree, nc, ng) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in) :: ng(2)
        real(rk), allocatable :: Tgc(:,:)
        integer :: i

        allocate(Tgc(ng(1)*ng(2), nc(1)*nc(2)))
#if defined(__NVCOMPILER)
        do i = 1, size(Xt, 1)
#else
        do concurrent (i = 1: size(Xt, 1))
#endif
            Tgc(i,:) = kron(&
                basis_bspline(Xt(i,2), knot2, nc(2), degree(2)),&
                basis_bspline(Xt(i,1), knot1, nc(1), degree(1)))
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_Tgc_bspline_2d_scalar(Xt, knot1, knot2, degree, nc) result(Tgc)
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
    pure subroutine lsq_fit_bspline(this, Xt, Xdata, ndata)
        use forcad_interface, only: solve
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:,:), Xdata(:,:)
        integer, intent(in) :: ndata(2)
        real(rk), allocatable :: T(:,:), Tt(:,:), TtT(:,:), TtX(:,:)
        integer :: i, n

        if (.not. this%err%ok) return

        if (this%nc(1) > ndata(1)) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'Invalid number of control points in the first direction.',&
                location   = 'lsq_fit_bspline',&
                suggestion = 'Ensure that the number of control points does not exceed the number of data points.')
            return
        end if

        if (this%nc(1) > ndata(1)) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = 'forcad_nurbs_surface',&
                message    = 'Invalid number of control points in the second direction.',&
                location   = 'lsq_fit_bspline',&
                suggestion = 'Ensure that the number of control points does not exceed the number of data points.')
            return
        end if

        n = ndata(1)*ndata(2)

        allocate(T(n, this%nc(1)*this%nc(2)))
#if defined(__NVCOMPILER) || (defined(__GFORTRAN__) && (__GNUC__ < 15 || (__GNUC__ == 15 && __GNUC_MINOR__ < 1)))
        do i = 1, n
#else
        do concurrent (i = 1: n)
#endif
            T(i,:) = kron(&
            basis_bspline(Xt(i,2), this%knot2, this%nc(2), this%degree(2)),&
            basis_bspline(Xt(i,1), this%knot1, this%nc(1), this%degree(1)))
        end do
        Tt = transpose(T)
        TtT = matmul(Tt, T)
        TtX = matmul(Tt, Xdata)
        this%Xc = solve(TtT, TtX)
    end subroutine
    !===============================================================================

end module forcad_nurbs_surface