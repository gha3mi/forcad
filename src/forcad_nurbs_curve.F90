!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Curves represented by B-spline and NURBS bases.
!!
!! A curve of degree \(p\), control points \(\mathbf{P}_i\), positive weights
!! \(w_i\), and B-spline basis functions \(N_{i,p}\) is evaluated as
!!
!! \[\mathbf{C}(u)=
!!   \frac{\sum_{i=1}^{n_c}N_{i,p}(u)w_i\mathbf{P}_i}
!!        {\sum_{i=1}^{n_c}N_{i,p}(u)w_i}.\]
!!
!! Control subscripts in this module's formulas are one-based to match Fortran
!! storage; the conventional zero-based \(i=0,\ldots,n\) notation is equivalent
!! under \(n_c=n+1\).
!!
!! With \(W(u)=\sum_jN_{j,p}(u)w_j\), the rational basis is
!! \(R_i(u)=N_{i,p}(u)w_i/W(u)\). Positive weights imply \(W(u)>0\) on
!! the active domain. At parameter values where the requested derivatives
!! exist,
!!
!! \[
!! R_i(u)\ge 0,\qquad \sum_iR_i(u)=1,\qquad
!! \sum_iR_i^{(k)}(u)=0\quad(k>0).
!! \]
!!
!! Polynomial basis derivatives vanish for \(k>p\), but rational derivatives
!! generally do not: differentiating the quotient introduces derivatives of
!! \(1/W\) at all orders supported by the spline's local smoothness.
!! Rational evaluation is invariant under the projective gauge
!! \(w_i\mapsto\alpha w_i\), \(\alpha>0\). Before each local rational sum, the
!! active weights are multiplied by one exact radix power chosen from their
!! largest model exponent. Homogeneous refinement applies the same normalization
!! to the complete weight vector and restores the original gauge afterward.
!! Rational geometric derivatives use the nested quotient form
!! \((S'-\mathbf C W')/W\), avoiding an explicit \(W^2\).
!!
!! Omitting the weights selects the polynomial B-spline form. Supplied
!! weights are also classified as uniform when their maximum deviation from
!! the first weight is at most `32*epsilon(rk)*maxval(Wc)`; in that case the
!! evaluation path treats them as exactly equal and omits rational
!! normalization. Control points use `[control point, physical coordinate]`
!! storage. The active parameter interval is
!! `[knot(degree+1), knot(nc+1)]` in Fortran indexing.
!!
!! Any finite, nondecreasing knot vector satisfying the spline-space size,
!! multiplicity, and positive-active-interval invariants is accepted. This
!! includes uniform or nonuniform, clamped, one-sided clamped, unclamped, and
!! periodic-form vectors. A complete cyclic knot extension together with
!! repeated control points and weights is verified by
!! [[nurbs_curve:get_parameter_topology]]; malformed periodic-form data remain
!! valid bounded splines rather than being mislabeled as periodic. Parameter
!! wrapping is an independent evaluation policy and does not impose topology.
!!
!! At an interior knot of multiplicity \(s\), a degree-\(p\) spline space has
!! guaranteed continuity \(C^{p-s}\); special control data may yield smoother
!! geometry. A full multiplicity \(p+1\) separates adjacent pieces. Endpoint
!! interpolation is guaranteed only when the corresponding endpoint is clamped
!! with multiplicity \(p+1\). At a knot where a derivative is discontinuous,
!! evaluation returns the derivative associated with the span selected by the
!! right-continuous basis convention (with explicit final-endpoint handling).
!!
!! Refinement operates on homogeneous controls
!! \(\mathbf H_i=(w_i\mathbf P_i,w_i)\) and preserves
!! \(\mathbf C_{\mathrm{old}}(u)=\mathbf C_{\mathrm{new}}(u)\) to working
!! precision on the active domain. Knot insertion and positive degree elevation
!! preserve a verified periodic topology, including its exterior extension and
!! repeated control rows. Other non-clamped representations are converted by
!! degree elevation to an equivalent active-domain open-clamped form. Periodic
!! knot removal is rejected because the bounded removal kernel cannot preserve
!! the cyclic identities. For curve integration,
!!
!! \[
!! \mathrm dL=\left\lVert\frac{\partial\mathbf C}{\partial u}\right\rVert
!! \mathrm du .
!! \]
!!
!! The [[nurbs_curve]] API includes evaluation, arbitrary-order basis
!! derivatives, IGA element data, geometry-preserving knot insertion and degree elevation,
!! tolerance-controlled knot removal, quadrature, fitting, discrete and local
!! continuous nearest-point searches, transformations, and VTK/IGES export.
!! Recoverable failures are reported in [[nurbs_curve:err]]. Unless a procedure
!! explicitly documents an output on failure, callers must inspect that
!! diagnostic before using returned data.
!!
!! **Primary reference:** L. Piegl and W. Tiller, *The NURBS Book*, 2nd ed.,
!! Springer, 1997.
!! [doi:10.1007/978-3-642-59223-2](https://doi.org/10.1007/978-3-642-59223-2).
!! **Finite-range reference:** IEEE, *IEEE Standard for Floating-Point
!! Arithmetic*, IEEE Std 754-2019, 2019.
!! [IEEE 754-2019](https://standards.ieee.org/ieee/754/6210/).
module forcad_nurbs_curve

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    use forcad_kinds, only: rk
    use forcad_utils, only: basis_bspline, elemConn_C0, ndgrid, compute_multiplicity, compute_knot_vector, basis_bspline_der,&
        insert_knot_A_5_1, insert_knot_periodic_A_5_1, findspan, elevate_degree_A_5_9, remove_knots_A_5_8, &
        elemConn_Cn, unique, rotate_points, gauss_leg, export_vtk_legacy, basis_bspline_2der, fill_uniform, &
        valid_knot_vector, knot_start, knot_end, knot_tolerance, active_knots, active_span_count, active_knot_multiplicity, &
        sparse_left_matmul, show_pyvista_singlepatch, basis_bspline_der_all_active, tensor_basis_derivative_local, &
        tensor_basis_derivatives2_local, map_parameter, periodic_topology
    use fordebug, only: debug

    implicit none

    private
    public nurbs_curve, compute_Tgc, compute_dTgc

    !===============================================================================
    !> Mutable B-spline or NURBS curve.
    !!
    !! Construct a valid curve with [[nurbs_curve:set]]. Sampling with
    !! [[nurbs_curve:create]] caches geometry points and their parameters;
    !! [[nurbs_curve:cmp_Xg]] evaluates requested parameters without requiring
    !! that cache. Refinement preserves the spline mapping where the underlying
    !! operation is exact, but changes control and element topology.
    !!
    !! `degree` means polynomial degree; the corresponding order is
    !! `degree+1`. Knot and control-point indices exposed by this API are
    !! one-based. Weights must be finite and strictly positive. Control
    !! coordinates are required to be finite and have at least one component;
    !! coordinate finiteness is a caller precondition in setter overloads that
    !! do not explicitly diagnose it.
    !!
    !! Diagnostic codes use the following categories: 100 invalid argument or
    !! index; 101 shape mismatch; 102--105 missing geometry state; 106
    !! underdetermined fitting problem; 107 invalid parameters; 108 singular
    !! system or mapping; and 109 iterative nonconvergence.
    !!
    !! @note Query methods return allocatable copies of owned arrays. Private
    !! storage prevents direct external mutation; mutating methods validate the
    !! core state they replace, while cache refresh remains an explicit caller
    !! operation as described below.
    !! @endnote
    !! @warning Core-data setters, fitting, control-point transformations, and
    !! weight changes do not automatically clear or recompute every sampled or
    !! connectivity cache. Call `create` after a geometry change before using
    !! cached `Xg`, and rebuild connectivity after a control-count or spline-
    !! space change.
    !! @endwarning
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    type nurbs_curve
        real(rk), allocatable, private :: Xc(:,:) !! Control points `[nc,ncoord]`.
        real(rk), allocatable, private :: Xg(:,:) !! Cached sampled points `[ng,ncoord]`.
        real(rk), allocatable, private :: Wc(:)   !! Positive control weights `[nc]`.
        real(rk), allocatable, private :: Xt(:)   !! Parameters of cached samples `[ng]`.
        real(rk), allocatable, private :: knot(:) !! Nondecreasing knot vector `[nc+degree+1]`.
        integer, private :: degree = 0            !! Polynomial degree \(p\).
        integer, private :: nc = 0                !! Number of control points.
        integer, private :: ng = 0                !! Number of cached geometry points.
        logical, private :: rational = .false.    !! True when nonuniform weights affect the basis.
        logical, private :: wrap_parameters = .false. !! Apply modulo mapping on evaluation.
        integer, allocatable, private :: elemConn_Xc_vis(:,:) !! Control-polygon VTK connectivity.
        integer, allocatable, private :: elemConn_Xg_vis(:,:) !! Sampled-curve VTK connectivity.
        integer, allocatable, private :: elemConn(:,:)        !! IGA element-to-control-point map.

        type(debug) :: err !! Recoverable diagnostic state for the most recent failed operation.
    contains
        procedure, private :: set1                  !!> Set an explicit knot vector, control points, and optional weights.
        procedure, private :: set1a                 !!> Set an explicit knot vector with inferred degree.
        procedure, private :: set2                  !!> Build a knot vector from breakpoints, degree, and continuity.
        procedure, private :: set3                  !!> Construct a Bezier or rational Bezier curve.
        procedure, private :: set4                  !!> Construct an open uniform curve from degree and size.
        generic :: set => set1, set1a, set2, set3, set4 !!> Replace curve geometry or initialize a fitting space.
        procedure :: create                !!> Sample the curve and cache parameters, points, and visualization connectivity.
        procedure :: cmp_Xg                !!> Evaluate physical points at one or more parameters.
        procedure, private :: get_Xc_all   !!> Return all control points.
        procedure, private :: get_Xci      !!> Return one control point.
        procedure, private :: get_Xcid     !!> Return one coordinate of one control point.
        generic :: get_Xc => get_Xc_all, get_Xci, get_Xcid !!> Query control-point data.
        procedure, private :: get_Xg_all   !!> Return all cached geometry points.
        procedure, private :: get_Xgi      !!> Return one cached geometry point.
        procedure, private :: get_Xgid     !!> Return one coordinate of one cached point.
        generic :: get_Xg => get_Xg_all, get_Xgi, get_Xgid !!> Query cached geometry points.
        procedure, private :: get_Wc_all   !!> Return all control weights.
        procedure, private :: get_Wci      !!> Return one control weight.
        generic :: get_Wc => get_Wc_all, get_Wci !!> Query control weights.
        procedure :: get_Xt                !!> Return parameters associated with cached geometry points.
        procedure, private :: get_knot_all !!> Return the complete knot vector.
        procedure, private :: get_knoti    !!> Return one knot value.
        generic :: get_knot => get_knoti, get_knot_all !!> Query knot-vector data.
        procedure :: get_ng                !!> Return the number of cached geometry points.
        procedure :: cmp_degree            !!> Infer the degree from knot and control-point counts.
        procedure, private :: get_degree_all !!> Return the polynomial degree.
        procedure, private :: get_degree_dir !!> Return the degree after validating direction 1.
        generic :: get_degree => get_degree_all, get_degree_dir !!> Query the polynomial degree.
        procedure :: finalize              !!> Release all owned and derived curve data while preserving diagnostics.
        procedure :: cmp_elem_Xc_vis       !!> Build line connectivity for the control polygon.
        procedure :: cmp_elem_Xg_vis       !!> Build line connectivity for cached curve samples.
        procedure :: cmp_elem_Xth          !!> Build line connectivity for the sampled parameter curve.
        procedure :: cmp_elem              !!> Build nonzero-span IGA element connectivity.
        procedure :: cmp_dof_map           !!> Map stored periodic controls to independent cyclic degrees of freedom.
        procedure :: get_elem_Xc_vis       !!> Return control-polygon visualization connectivity.
        procedure :: get_elem_Xg_vis       !!> Return sampled-curve visualization connectivity.
        procedure :: get_elem              !!> Return IGA element-to-control-point connectivity.
        procedure :: set_elem_Xc_vis       !!> Replace control-polygon visualization connectivity.
        procedure :: set_elem_Xg_vis       !!> Replace sampled-curve visualization connectivity.
        procedure :: set_elem              !!> Replace IGA element connectivity after shape validation.
        procedure :: export_Xc             !!> Export the control polygon in legacy VTK format.
        procedure :: export_Xg             !!> Export cached curve samples in legacy VTK format.
        procedure :: export_Xth            !!> Export the sampled parameter curve in legacy VTK format.
        procedure :: export_Xth_in_Xg      !!> Export parameter-grid lines embedded in physical space.
        procedure :: export_iges           !!> Export the curve as IGES entity 126.
        procedure :: modify_Xc             !!> Replace one control coordinate without resampling cached geometry.
        procedure :: modify_Wc             !!> Replace one positive weight and update rational state.
        procedure :: get_multiplicity      !!> Return multiplicities of all distinct stored knots.
        procedure :: get_continuity        !!> Return \(p-s\) metadata for all distinct stored knots.
        procedure :: cmp_nc                !!> Compute the control-point count implied by knots and degree.
        procedure, private :: get_nc_dir   !!> Return the control-point count after validating direction 1.
        procedure, private :: get_nc_all   !!> Return the control-point count.
        generic :: get_nc => get_nc_all, get_nc_dir !!> Query the control-point count.
        procedure :: insert_knots          !!> Insert knots while preserving active-domain geometry.
        procedure :: elevate_degree        !!> Elevate degree while preserving the active-domain geometry.
        procedure :: is_valid              !!> Report whether complete, internally consistent curve geometry is stored.
        procedure, private :: basis_vector !!> Evaluate all rational basis functions at many parameters.
        procedure, private :: basis_scalar !!> Evaluate all rational basis functions at one parameter.
        generic :: basis => basis_vector, basis_scalar !!> Evaluate the dense rational or polynomial basis.
        procedure, private :: derivative_vector      !!> Evaluate first parametric basis derivatives at many points.
        procedure, private :: derivative_scalar      !!> Evaluate first parametric basis derivatives at one point.
        generic :: derivative => derivative_vector, derivative_scalar !!> Evaluate first parametric basis derivatives.
        procedure, private :: derivative2_vector     !!> Evaluate second parametric basis derivatives at many points.
        procedure, private :: derivative2_scalar     !!> Evaluate second parametric basis derivatives at one point.
        generic :: derivative2 => derivative2_vector, derivative2_scalar !!> Evaluate second parametric basis derivatives.
        procedure, private :: derivative_order_vector !!> Evaluate a requested basis derivative order at many points.
        procedure, private :: derivative_order_scalar !!> Evaluate a requested basis derivative order at one point.
        generic :: derivative_order => derivative_order_vector, derivative_order_scalar !!> Evaluate dense derivatives of any order.
        procedure, private :: derivative_order_active_vector !!> Evaluate only locally active derivatives at many points.
        procedure, private :: derivative_order_active_scalar !!> Evaluate only locally active derivatives at one point.
        !> Evaluate active derivatives of any order.
        generic :: derivative_order_active => &
            derivative_order_active_vector, derivative_order_active_scalar
        procedure :: is_rational           !!> Report whether nonuniform weights affect the current basis.
        procedure :: get_parameter_wrapping !!> Report whether modulo parameter mapping is enabled.
        procedure :: get_parameter_topology !!> Classify the stored parameter topology.
        procedure :: put_to_nurbs          !!> Map rows by their first coordinate to the parameter interval and cache the curve.
        procedure :: remove_knots          !!> Remove knots that pass an internal roundoff-scaled geometry test.
        procedure :: rotate_Xc             !!> Rotate control points in three-dimensional embedding space.
        procedure :: rotate_Xg             !!> Rotate cached geometry points in three-dimensional embedding space.
        procedure :: translate_Xc          !!> Translate control points by a physical-space vector.
        procedure :: translate_Xg          !!> Translate cached geometry points by a physical-space vector.
        procedure :: show                  !!> Display previously exported representations with PyVista.
        procedure :: nearest_point         !!> Select the nearest cached geometry sample.
        procedure :: nearest_point2        !!> Minimize squared distance with bounded Newton iterations and backtracking.
        procedure :: ansatz                !!> Evaluate element-local physical basis data and differential length.
        procedure :: cmp_length            !!> Integrate curve length with elementwise Gauss-Legendre quadrature.
        procedure :: lsq_fit_bspline       !!> Fit polynomial control points by structured least squares.
        procedure :: lsq_fit_nurbs         !!> Fit control points and positive weights by damped nonlinear least squares.

        ! Shapes
        procedure :: set_circle            !!> Construct an exact rational circle.
        procedure :: set_half_circle       !!> Construct an exact rational semicircle.
        procedure :: set_C                 !!> Construct the predefined planar C-shaped curve.
    end type
    !===============================================================================

    !> Evaluate polynomial or rational curve points from explicit spline data.
    !! This implementation generic is private; use [[nurbs_curve:cmp_Xg]] in
    !! user code so validation and parameter-mapping policy remain attached to
    !! the object.
    interface compute_Xg
        module procedure compute_Xg_nurbs_1d
        module procedure compute_Xg_bspline_1d
        module procedure compute_Xg_nurbs_1d_1point
        module procedure compute_Xg_bspline_1d_1point
    end interface

    !> Evaluate dense curve basis values from explicit spline data.
    !!
    !! Scalar overloads return one vector of length `nc`; vector overloads
    !! return one basis row per parameter. Supplying `Wc` selects rational
    !! normalization. Values outside the active domain are zero unless wrapping
    !! is enabled.
    !!
    !! For the rational overload,
    !! \(R_i=N_{i,p}w_i/W\), with \(W=\sum_jN_{j,p}w_j\). The result therefore
    !! has partition of unity wherever \(W\ne0\).
    interface compute_Tgc
        module procedure compute_Tgc_nurbs_1d_vector
        module procedure compute_Tgc_bspline_1d_vector
        module procedure compute_Tgc_nurbs_1d_scalar
        module procedure compute_Tgc_bspline_1d_scalar
    end interface

    !> Evaluate dense first parametric derivatives of the curve basis.
    !! Optional `Tgc` output returns basis values computed in the same pass.
    !! Rational overloads apply the quotient rule consistently to values and
    !! derivatives.
    !!
    !! \[
    !! R_i'=\frac{w_i(N_i'W-N_iW')}{W^2},\qquad
    !! W'=\sum_jN_j'w_j .
    !! \]
    interface compute_dTgc
        module procedure compute_dTgc_nurbs_1d_vector
        module procedure compute_dTgc_bspline_1d_vector
        module procedure compute_dTgc_nurbs_1d_scalar
        module procedure compute_dTgc_bspline_1d_scalar
    end interface

    !> Evaluate dense second parametric derivatives of the curve basis.
    !! This implementation generic is private; the public object API is
    !! [[nurbs_curve:derivative2]].
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
    !> Test whether every supplied NURBS weight is finite and strictly positive.
    !! An empty array satisfies the elemental conjunction and therefore returns
    !! true; public setters separately require the expected control count.
    pure logical function valid_weights(Wc) result(ok)
        real(rk), intent(in), contiguous :: Wc(:)
            !! Candidate control weights.

        ok = all(ieee_is_finite(Wc) .and. Wc > 0.0_rk)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report whether the defining curve geometry is complete and consistent.
    !!
    !! A valid curve has \(n_c\) finite control points with at least one
    !! physical coordinate, a degree \(0\le p<n_c\), and a finite
    !! nondecreasing knot vector satisfying
    !! \(\operatorname{size}(U)=n_c+p+1\). If explicit weights are stored,
    !! exactly \(n_c\) finite positive values are required and the cached
    !! rational classification must agree with them. Sample and connectivity
    !! caches are derived state and are not inspected.
    pure function is_valid(this) result(ok)
        class(nurbs_curve), intent(in) :: this
            !! Curve geometry to inspect.
        logical :: ok

        ok = .false.
        if (.not. this%err%ok) return
        if (.not. allocated(this%knot) .or. .not. allocated(this%Xc)) return
        if (.not. valid_knot_vector(this%knot, this%nc, this%degree)) return
        if (size(this%Xc,1) /= this%nc .or. size(this%Xc,2) < 1) return
        if (.not. all(ieee_is_finite(this%Xc))) return

        if (allocated(this%Wc)) then
            if (size(this%Wc) /= this%nc .or. .not. valid_weights(this%Wc)) return
            if (this%rational .neqv. &
                max(maxval(this%Wc)-this%Wc(1),this%Wc(1)-minval(this%Wc)) > &
                32.0_rk*epsilon(1.0_rk)*maxval(this%Wc)) return
        else if (this%rational) then
            return
        end if

        ok = .true.
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Replace the curve with explicit knot and control data.
    !!
    !! The knot vector must be finite, nondecreasing, have size
    !! `size(Xc,1)+degree+1`, have no multiplicity above `degree+1`, and define
    !! a positive active interval. If `degree` is absent it is inferred from
    !! this size relation. At least one finite coordinate per control point is
    !! a caller precondition; this overload validates the control-row count but
    !! does not diagnose a zero coordinate dimension or nonfinite coordinates.
    !! The operation validates knots and weights before replacing core state.
    !!
    pure subroutine set1(this, knot, Xc, Wc, degree, wrap_parameters)
        class(nurbs_curve), intent(inout) :: this
            !! Curve to replace; an existing diagnostic must be cleared first.
        real(rk), intent(in), contiguous :: knot(:)
            !! Complete knot vector, including knots outside the active interval.
        real(rk), intent(in), contiguous :: Xc(:,:)
            !! Control points `[nc,ncoord]`, with `ncoord>=1`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional finite positive weights `[nc]`.
        integer, intent(in), optional :: degree
            !! Optional polynomial degree; inferred when absent.
        logical, intent(in), optional :: wrap_parameters
            !! Enable half-open modulo evaluation on the active interval.
        integer :: degree_, nc_
        logical :: wrap_parameters_

        if (.not. this%err%ok) return

        nc_ = size(Xc, 1)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
        if (present(degree)) then
            degree_ = degree
        else
            degree_ = size(knot) - nc_ - 1
        end if

        if (.not. valid_knot_vector(knot, nc_, degree_)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Invalid knot vector for the supplied control points.",&
                location   = "set1",&
                suggestion = "Use finite nondecreasing knots with size(knot) == nc + degree + 1 and degree < nc.")
            return
        end if
        if (present(Wc)) then
            if (size(Wc) /= nc_) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Weights length mismatch: size(Wc) must equal number of control points.",&
                    location   = "set1",&
                    suggestion = "Provide Wc with size(Wc) == size(Xc,1).")
                return
            else if (.not. valid_weights(Wc)) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Invalid rational weights: all weights must be finite and positive.",&
                    location   = "set1",&
                    suggestion = "Use strictly positive finite Wc values for rational NURBS geometry.")
                return
            end if
        end if

        this%Xc = Xc
        this%nc = nc_
        this%knot = knot
        this%degree = degree_
        this%wrap_parameters = wrap_parameters_
        if (present(Wc)) then
            this%Wc = Wc
            this%rational = max(maxval(this%Wc)-this%Wc(1),this%Wc(1)-minval(this%Wc)) > &
                32.0_rk*epsilon(1.0_rk)*maxval(this%Wc)
        else
            if (allocated(this%Wc)) deallocate(this%Wc)
            this%rational = .false.
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export knot-span parameter lines evaluated in physical space.
    !! Every nonzero active span receives at least `res` samples, with shared
    !! endpoints emitted once. The curve embedding must have two or three
    !! coordinates. `encoding` is forwarded to the legacy VTK writer.
    impure subroutine export_Xth_in_Xg(this, filename, res, encoding)
        class(nurbs_curve), intent(inout) :: this
        character(len=*), intent(in) :: filename
            !! Destination VTK path.
        integer, intent(in), optional :: res
            !! Minimum points per active span; default 10, minimum 2.
        character(len=*), intent(in), optional :: encoding
            !! Optional legacy VTK encoding selector.
        real(rk), allocatable :: U(:), Ur(:), Ur_eval(:), Xg(:,:)
        integer, allocatable :: elemConn(:,:)
        integer :: nspan, res_min, n, s, r, o, ncoord

        if (.not. this%err%ok) return

        if (.not. allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Control points are not set.",&
                location   = "export_Xth_in_Xg",&
                suggestion = "Call set(...) first before exporting.")
            return
        end if

        if (.not. allocated(this%knot)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Knot vector is not set.",&
                location   = "export_Xth_in_Xg",&
                suggestion = "Call set(...) first before exporting.")
            return
        end if

        ncoord = size(this%Xc, 2)
        if (ncoord < 2 .or. ncoord > 3) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Invalid geometry dimension.",&
                location   = "export_Xth_in_Xg",&
                suggestion = "Use two- or three-dimensional control points before exporting.")
            return
        end if

        res_min = 10
        if (present(res)) res_min = max(2, res)

        U = active_knots(this%knot, this%nc, this%degree)
        if (size(U) < 2) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Knot vector needs at least two active unique values.",&
                location   = "export_Xth_in_Xg",&
                suggestion = "Check the knot vector for sufficient unique values.")
            return
        end if

        nspan = size(U) - 1
        n = nspan*(res_min - 1) + 1
        allocate(Ur(n), elemConn(1,n))

        do concurrent (s = 1:nspan, r = 1:res_min-1) local(o)
            o = (s - 1)*(res_min - 1)
            Ur(o+r) = U(s) + (U(s+1) - U(s))*real(r - 1, rk)/real(res_min - 1, rk)
        end do
        Ur(n) = U(size(U))
        Ur_eval = Ur

        if (this%is_rational()) then
            Xg = compute_Xg(Ur_eval, this%knot, this%degree, this%nc, n, this%Xc, this%Wc, this%wrap_parameters)
        else
            Xg = compute_Xg(Ur_eval, this%knot, this%degree, this%nc, n, this%Xc, this%wrap_parameters)
        end if

        do concurrent (r = 1:n)
            elemConn(1,r) = r
        end do

        call export_vtk_legacy(filename=filename, points=Xg, elemConn=elemConn, vtkCellType=4, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Replace the curve with scalar control ordinates embedded on the x-axis.
    !! The one-dimensional `Xc` input becomes a three-coordinate control array
    !! `[Xc,0,0]`. Knot, degree, weight, and wrapping rules are identical to the
    !! rank-two-control-array `set` overload.
    pure subroutine set1a(this, knot, Xc, Wc, degree, wrap_parameters)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: knot(:)
            !! Complete knot vector.
        real(rk), intent(in), contiguous :: Xc(:)
            !! Scalar control ordinates `[nc]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional finite positive weights `[nc]`.
        integer, intent(in), optional :: degree
            !! Optional polynomial degree; inferred when absent.
        logical, intent(in), optional :: wrap_parameters
            !! Enable modulo evaluation on the active interval.
        integer :: degree_, nc_
        logical :: wrap_parameters_

        if (.not. this%err%ok) return

        nc_ = size(Xc)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
        if (present(degree)) then
            degree_ = degree
        else
            degree_ = size(knot) - nc_ - 1
        end if

        if (.not. valid_knot_vector(knot, nc_, degree_)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Invalid knot vector for the supplied control points.",&
                location   = "set1a",&
                suggestion = "Use finite nondecreasing knots with size(knot) == nc + degree + 1 and degree < nc.")
            return
        end if
        if (present(Wc)) then
            if (size(Wc) /= nc_) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Weights length mismatch: size(Wc) must equal number of control points.",&
                    location   = "set1a",&
                    suggestion = "Provide Wc with size(Wc) == size(Xc,1).")
                return
            else if (.not. valid_weights(Wc)) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Invalid rational weights: all weights must be finite and positive.",&
                    location   = "set1a",&
                    suggestion = "Use strictly positive finite Wc values for rational NURBS geometry.")
                return
            end if
        end if

        if (allocated(this%Xc)) then
            if (size(this%Xc, 1) /= nc_ .or. size(this%Xc, 2) /= 3) deallocate(this%Xc)
        end if
        if (.not. allocated(this%Xc)) allocate(this%Xc(nc_, 3))
        this%Xc = 0.0_rk
        this%Xc(:,1) = Xc
        this%nc = nc_
        this%knot = knot
        this%degree = degree_
        this%wrap_parameters = wrap_parameters_
        if (present(Wc)) then
            this%Wc = Wc
            this%rational = max(maxval(this%Wc)-this%Wc(1),this%Wc(1)-minval(this%Wc)) > &
                32.0_rk*epsilon(1.0_rk)*maxval(this%Wc)
        else
            if (allocated(this%Wc)) deallocate(this%Wc)
            this%rational = .false.
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct a spline space from breakpoints and requested continuities.
    !!
    !! `continuity(i)` is converted to knot multiplicity `degree-continuity(i)`
    !! for each parameter node according to [[compute_knot_vector]]. Optional
    !! control data must match the resulting count. Breakpoints must be finite
    !! and strictly increasing. Values
    !! \(-1\leq\mathtt{continuity(i)}\leq p-1\) retain a breakpoint with
    !! multiplicity \(p-\mathtt{continuity(i)}\); `-1` therefore gives full
    !! multiplicity \(p+1\). The value `continuity(i)=p` omits that breakpoint
    !! from the knot vector. The generated vector is validated before object
    !! state changes. If `Xc` is absent, any previously stored control points
    !! are released and the result is a spline space for subsequent fitting,
    !! not complete curve geometry. Optional weights are retained only when
    !! explicitly supplied.
    pure subroutine set2(this, Xth_dir, degree, continuity, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xth_dir(:)
            !! Finite, strictly increasing breakpoint values.
        integer, intent(in) :: degree
            !! Polynomial degree.
        integer, intent(in), contiguous :: continuity(:)
            !! One value per breakpoint in `-1:degree`; `degree` omits the breakpoint.
        real(rk), intent(in), contiguous, optional :: Xc(:,:)
            !! Optional control points `[computed_nc,ncoord]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional finite positive weights `[computed_nc]`.
        real(rk), allocatable :: knot_(:)
        integer :: nc_

        if (.not. this%err%ok) return

        knot_ = compute_knot_vector(Xth_dir, degree, continuity)
        nc_ = size(knot_) - degree - 1
        if (.not. valid_knot_vector(knot_, nc_, degree)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Invalid knot vector generated from Xth_dir, degree, and continuity.",&
                location   = "set2",&
                suggestion = "Use nondecreasing Xth_dir and continuity values that give knot multiplicity <= degree+1.")
            return
        end if

        if (present(Xc)) then
            if (size(Xc,1) /= nc_) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve", &
                    message    = "Control points size mismatch in set2",&
                    location   = "set2", &
                    suggestion = "size(Xc,1) must equal computed nc.")
                return
            end if
        end if

        if (present(Wc)) then
            if (size(Wc) /= nc_) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve", &
                    message    = "Weights size mismatch in set2",&
                    location   = "set2", &
                    suggestion = "size(Wc) must equal computed nc.")
                return
            else if (.not. valid_weights(Wc)) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Invalid rational weights: all weights must be finite and positive.",&
                    location   = "set2",&
                    suggestion = "Use strictly positive finite Wc values for rational NURBS geometry.")
                return
            end if
        end if

        this%knot = knot_
        this%degree = degree
        this%nc = nc_
        if (present(Xc)) then
            this%Xc = Xc
        else
            if (allocated(this%Xc)) deallocate(this%Xc)
        end if

        if (present(Wc)) then
            this%Wc = Wc
            this%rational = max(maxval(this%Wc)-this%Wc(1),this%Wc(1)-minval(this%Wc)) > &
                32.0_rk*epsilon(1.0_rk)*maxval(this%Wc)
        else
            if (allocated(this%Wc)) deallocate(this%Wc)
            this%rational = .false.
        end if
        this%wrap_parameters = .false.
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct a Bezier or rational Bezier curve on \([0,1]\).
    !! The degree is `size(Xc,1)-1` and both endpoint knots have multiplicity
    !! `degree+1`.
    pure subroutine set3(this, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xc(:,:)
            !! Bezier control points `[degree+1,ncoord]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional finite positive Bezier weights `[degree+1]`.
        integer :: nc_

        if (.not. this%err%ok) return

        nc_ = size(Xc, 1)
        if (present(Wc)) then
            if (size(Wc) /= nc_) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Weights length mismatch: size(Wc) must equal number of control points.",&
                    location   = "set3",&
                    suggestion = "Provide Wc with size(Wc) == size(Xc,1).")
                return
            else if (.not. valid_weights(Wc)) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Invalid rational weights: all weights must be finite and positive.",&
                    location   = "set3",&
                    suggestion = "Use strictly positive finite Wc values for rational NURBS geometry.")
                return
            end if
        end if

        this%Xc = Xc
        this%nc = nc_

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
        this%wrap_parameters = .false.
        if (present(Wc)) then
            this%Wc = Wc
            this%rational = max(maxval(this%Wc)-this%Wc(1),this%Wc(1)-minval(this%Wc)) > &
                32.0_rk*epsilon(1.0_rk)*maxval(this%Wc)
        else
            if (allocated(this%Wc)) deallocate(this%Wc)
            this%rational = .false.
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct a clamped open-uniform B-spline or NURBS curve on \([0,1]\).
    !! Interior knots are uniformly spaced and simple.
    pure subroutine set4(this, degree, nc, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in) :: degree
            !! Polynomial degree satisfying `0<=degree<nc`.
        integer, intent(in) :: nc
            !! Number of control points.
        real(rk), intent(in), contiguous :: Xc(:,:)
            !! Control points `[nc,ncoord]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional finite positive weights `[nc]`.
        integer :: m, i

        if (.not. this%err%ok) return

        if (degree < 0 .or. nc <= 0 .or. degree >= nc .or. size(Xc, 1) /= nc .or. size(Xc, 2) < 1) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Invalid curve dimensions for generated knot-vector construction.",&
                location   = "set4",&
                suggestion = "Require degree >= 0, nc > degree, size(Xc,1) == nc, and at least one coordinate dimension.")
            return
        end if

        if (present(Wc)) then
            if (size(Wc) /= nc) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Weights length mismatch: size(Wc) must equal number of control points.",&
                    location   = "set4",&
                    suggestion = "Provide Wc with size(Wc) == nc.")
                return
            else if (.not. valid_weights(Wc)) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Invalid rational weights: all weights must be finite and positive.",&
                    location   = "set4",&
                    suggestion = "Use strictly positive finite Wc values for rational NURBS geometry.")
                return
            end if
        end if

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
        this%wrap_parameters = .false.

        if (present(Wc)) then
            if (allocated(this%Wc)) then
                if (size(this%Wc) /= size(Wc)) deallocate(this%Wc)
            end if
            this%Wc = Wc
            this%rational = max(maxval(this%Wc)-this%Wc(1),this%Wc(1)-minval(this%Wc)) > &
                32.0_rk*epsilon(1.0_rk)*maxval(this%Wc)
        else
            if (allocated(this%Wc)) deallocate(this%Wc)
            this%rational = .false.
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Sample and cache curve geometry.
    !!
    !! With explicit `Xt`, those parameters are evaluated. Otherwise, when
    !! `res` is present, uniformly spaced samples cover the active interval
    !! (both ends are included when `res>1`). If both are omitted, an existing
    !! cached parameter vector is reused; absence of all three is an error.
    !! Geometry samples and visualization connectivity are replaced.
    pure subroutine create(this, res, Xt)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: res
            !! Number of generated samples when `Xt` is absent.
        real(rk), intent(in), contiguous, optional :: Xt(:)
            !! Optional explicit parameter vector.
        real(rk), allocatable :: Xg_new(:,:)

        if (.not. this%err%ok) return

        ! check
        if (.not.allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Control points are not set.",&
                location   = "create",&
                suggestion = "Call set(...) first before create().")
            return
        end if

        if (.not.allocated(this%knot)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Knot vector is not set.",&
                location   = "create",&
                suggestion = "Call set(...) first before create().")
            return
        end if

        ! Set parameter values
        if (present(Xt)) then
            if (size(Xt) < 1) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Parameter values are empty.",&
                    location   = "create",&
                    suggestion = "Pass a nonempty Xt array or use res >= 1.")
                return
            end if
            if (allocated(this%Xt)) then
                if (size(this%Xt) /= size(Xt)) deallocate(this%Xt)
            end if
            this%Xt = Xt
        else if (present(res)) then
            if (res < 1) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Resolution must be at least one.",&
                    location   = "create",&
                    suggestion = "Use res >= 1 or pass explicit nonempty Xt values.")
                return
            end if
            if (allocated(this%Xt)) then
                if (size(this%Xt) /= res) then
                    deallocate(this%Xt)
                    allocate(this%Xt(res))
                end if
            else
                allocate(this%Xt(res))
            end if
            call fill_uniform(knot_start(this%knot, this%nc, this%degree), knot_end(this%knot, this%nc, this%degree), this%Xt)
        end if
        if (.not. allocated(this%Xt)) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Parameter values are not set.",&
                location   = "create",&
                suggestion = "Pass Xt or res when calling create(...).")
            return
        end if

        ! Set number of geometry points
        this%ng = size(this%Xt)

        if (this%is_rational()) then ! NURBS
            Xg_new = compute_Xg(this%Xt, this%knot, this%degree, this%nc, this%ng, this%Xc, this%Wc, this%wrap_parameters)
        else ! B-Spline
            Xg_new = compute_Xg(this%Xt, this%knot, this%degree, this%nc, this%ng, this%Xc, this%wrap_parameters)
        end if
        call move_alloc(Xg_new, this%Xg)
        this%elemConn_Xg_vis = this%cmp_elem_Xg_vis()
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate one physical curve point.
    !! Returns a vector of physical coordinate dimension. Outside-domain values
    !! evaluate to the zero basis unless parameter wrapping is enabled.
    pure function cmp_Xg(this, Xt) result(Xg)
        class(nurbs_curve), intent(in) :: this
        real(rk), intent(in) :: Xt
            !! Parameter to evaluate.
        real(rk), allocatable :: Xg(:)

        if (.not. this%err%ok .or. .not. allocated(this%Xc) .or. .not. allocated(this%knot)) then
            allocate(Xg(0))
            return
        end if

        if (this%is_rational()) then ! NURBS
            Xg = compute_Xg(map_parameter(Xt, this%knot, this%nc, this%degree, this%wrap_parameters), &
                this%knot, this%degree, this%nc, this%Xc, this%Wc)
        else ! B-Spline
            Xg = compute_Xg(map_parameter(Xt, this%knot, this%nc, this%degree, this%wrap_parameters), &
                this%knot, this%degree, this%nc, this%Xc)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Parameterize rows by their first coordinate and cache the mapped curve.
    !!
    !! The first column of `X` is affinely mapped from its data range
    !! \([x_{\min},x_{\max}]\) to the active knot interval \([u_a,u_b]\):
    !!
    !! \[
    !! u_i=u_a+\frac{X_{i1}-x_{\min}}{x_{\max}-x_{\min}}(u_b-u_a).
    !! \]
    !!
    !! Other columns of `X` do not affect the parameterization. The routine
    !! atomically replaces the cached parameter-point pairs
    !! \((u_i,\mathbf C(u_i))\), the sample count, and sampled visualization
    !! connectivity. When `elemConn` is absent, consecutive input rows are
    !! connected as a piecewise-linear curve.
    pure subroutine put_to_nurbs(this, X, elemConn)
        class(nurbs_curve), intent(inout) :: this
            !! Initialized curve whose sampled geometry cache is replaced.
        real(rk), intent(in), contiguous :: X(:,:)
            !! Ordering data `[npoint,ncolumn]`; column 1 must have nonzero range.
        integer, intent(in), contiguous, optional :: elemConn(:,:)
            !! Optional one-based visualization connectivity for the cached points.
        real(rk), allocatable :: Xt(:)
        real(rk), allocatable :: Xg_new(:,:)
        real(rk) :: min_X1, max_X1, start1, end1, scale1
        integer :: i

        if (.not. this%err%ok) return

        if (.not. allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Control points are not set.",&
                location   = "put_to_nurbs",&
                suggestion = "Call set(...) first before mapping points to the NURBS curve.")
            return
        end if
        if (.not. allocated(this%knot)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Knot vector is not set.",&
                location   = "put_to_nurbs",&
                suggestion = "Call set(...) first before mapping points to the NURBS curve.")
            return
        end if
        if (size(X,1) < 1 .or. size(X,2) < 1) then
            call this%err%set(&
                code       = 101,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Input points must have at least one point and one coordinate.",&
                location   = "put_to_nurbs",&
                suggestion = "Provide X(:,:) with size(X,1) >= 1 and size(X,2) >= 1.")
            return
        end if

        min_X1 = X(1,1)
        max_X1 = X(1,1)
        do i = 2, size(X,1)
            min_X1 = min(min_X1, X(i,1))
            max_X1 = max(max_X1, X(i,1))
        end do
        if (max_X1 <= min_X1) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Input points have zero extent in the curve parameter direction.",&
                location   = "put_to_nurbs",&
                suggestion = "Provide points with distinct first-coordinate values.")
            return
        end if

        start1 = knot_start(this%knot, this%nc, this%degree)
        end1 = knot_end(this%knot, this%nc, this%degree)
        allocate(Xt(size(X,1)))
        scale1 = (end1 - start1)/(max_X1 - min_X1)
        do concurrent (i = 1:size(X,1))
            Xt(i) = (X(i,1) - min_X1)*scale1 + start1
        end do

        if (this%is_rational()) then
            Xg_new = compute_Xg(Xt, this%knot, this%degree, this%nc, size(X,1), this%Xc, this%Wc, this%wrap_parameters)
        else
            Xg_new = compute_Xg(Xt, this%knot, this%degree, this%nc, size(X,1), this%Xc, this%wrap_parameters)
        end if
        call move_alloc(Xt, this%Xt)
        call move_alloc(Xg_new, this%Xg)
        this%ng = size(this%Xt)

        if (present(elemConn)) then
            call this%set_elem_Xg_vis(elemConn)
        else
            this%elemConn_Xg_vis = this%cmp_elem_Xg_vis()
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocatable copy of the complete control polygon.
    !! The result has shape `[nc,ncoord]`, or `[0,0]` when unavailable.
    pure function get_Xc_all(this) result(Xc)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: Xc(:,:)

        if (.not. this%err%ok .or. .not. allocated(this%Xc)) then
            allocate(Xc(0,0))
            return
        end if

        Xc = this%Xc
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocatable copy of control point `n`.
    !! The result has size `ncoord`, or zero for an invalid index or state.
    pure function get_Xci(this, n) result(Xc)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: n
        real(rk), allocatable :: Xc(:)

        if (.not. this%err%ok .or. .not. allocated(this%Xc)) then
            allocate(Xc(0))
            return
        end if
        if (n<lbound(this%Xc,1) .or. n>ubound(this%Xc,1)) then
            allocate(Xc(0))
            return
        end if

        Xc = this%Xc(n,:)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return coordinate `dir` of control point `n`.
    !! A quiet NaN denotes an invalid index or unavailable control data.
    pure elemental function get_Xcid(this, n, dir) result(Xc)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: n
        integer, intent(in) :: dir
        real(rk) :: Xc

        Xc = 0.0_rk
        if (.not. this%err%ok .or. .not. allocated(this%Xc)) return
        if (n<lbound(this%Xc,1) .or. n>ubound(this%Xc,1)) return
        if (dir<lbound(this%Xc,2) .or. dir>ubound(this%Xc,2)) return

        Xc = this%Xc(n, dir)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocatable copy of all cached geometry points.
    !! The result has shape `[ng,ncoord]`, or `[0,0]` when unavailable.
    pure function get_Xg_all(this) result(Xg)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: Xg(:,:)

        if (.not. this%err%ok .or. .not. allocated(this%Xg)) then
            allocate(Xg(0,0))
            return
        end if

        Xg = this%Xg
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocatable copy of cached geometry point `n`.
    !! The result has size `ncoord`, or zero for an invalid index or state.
    pure function get_Xgi(this, n) result(Xg)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: n
        real(rk), allocatable :: Xg(:)

        if (.not. this%err%ok .or. .not. allocated(this%Xg)) then
            allocate(Xg(0))
            return
        end if
        if (n<lbound(this%Xg,1) .or. n>ubound(this%Xg,1)) then
            allocate(Xg(0))
            return
        end if

        Xg = this%Xg(n,:)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return coordinate `dir` of cached geometry point `n`.
    !! A quiet NaN denotes an invalid index or unavailable sampled data.
    pure elemental function get_Xgid(this, n, dir) result(Xg)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: n
        integer, intent(in) :: dir
        real(rk) :: Xg

        Xg = 0.0_rk
        if (.not. this%err%ok .or. .not. allocated(this%Xg)) return
        if (n<lbound(this%Xg,1) .or. n>ubound(this%Xg,1)) return
        if (dir<lbound(this%Xg,2) .or. dir>ubound(this%Xg,2)) return

        Xg = this%Xg(n, dir)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocatable copy of all explicit control weights.
    !! Polynomial curves without stored weights return a zero-length array.
    pure function get_Wc_all(this) result(Wc)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: Wc(:)

        if (.not. this%err%ok .or. .not. allocated(this%Wc)) then
            allocate(Wc(0))
            return
        end if

        Wc = this%Wc
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return explicit control weight `n`.
    !! A quiet NaN denotes an invalid index or unavailable weight data.
    pure elemental function get_Wci(this, n) result(Wc)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: n
        real(rk) :: Wc

        Wc = 0.0_rk
        if (.not. this%err%ok .or. .not. allocated(this%Wc)) return
        if (n<lbound(this%Wc,1) .or. n>ubound(this%Wc,1)) return

        Wc = this%Wc(n)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocatable copy of parameters associated with cached points.
    !! The result has size `ng`, or zero when no parameter cache exists.
    pure function get_Xt(this) result(Xt)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: Xt(:)

        if (.not. this%err%ok .or. .not. allocated(this%Xt)) then
            allocate(Xt(0))
            return
        end if

        Xt = this%Xt
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the number of cached geometry points.
    !! Zero denotes an empty or invalid sampled state.
    pure elemental function get_ng(this) result(ng)
        class(nurbs_curve), intent(in) :: this
        integer :: ng

        ng = 0
        if (.not. this%err%ok) return

        ng = this%ng
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Infer the polynomial degree from current knot and control counts.
    !! For a complete knot vector, \(p=\operatorname{size}(U)-n_c-1\).
    pure subroutine cmp_degree(this)
        class(nurbs_curve), intent(inout) :: this
        integer, allocatable :: m(:)

        if (.not. this%err%ok) return

        if (allocated(this%knot) .and. this%nc > 0) then
            this%degree = size(this%knot) - this%nc - 1
        else
            m = this%get_multiplicity()
            this%degree = m(1) - 1
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the current polynomial degree.
    !! The default value zero is returned when the object diagnostic is active.
    pure elemental function get_degree_all(this) result(degree)
        class(nurbs_curve), intent(in) :: this
        integer :: degree

        degree = 0
        if (.not. this%err%ok) return

        degree = this%degree
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the polynomial degree after validating curve direction 1.
    !! Directions other than one return zero.
    pure elemental function get_degree_dir(this, dir) result(degree)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: dir
        integer :: degree

        degree = 0
        if (.not. this%err%ok) return

        if (dir == 1) degree = this%degree
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocatable copy of the complete knot vector.
    !! The result is zero length when no valid knot vector is available.
    pure function get_knot_all(this) result(knot)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: knot(:)

        if (.not. this%err%ok .or. .not. allocated(this%knot)) then
            allocate(knot(0))
            return
        end if

        knot = this%knot
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return knot `i`.
    !! A quiet NaN denotes an invalid index or object state.
    pure elemental function get_knoti(this,i) result(knot)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: i
        real(rk) :: knot

        knot = 0.0_rk
        if (.not. this%err%ok .or. .not. allocated(this%knot)) return
        if (i < 1 .or. i > size(this%knot)) return

        knot = this%knot(i)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Release all curve-owned arrays and restore default scalar state.
    !! The diagnostic object is retained and can be cleared through its own API.
    pure elemental subroutine finalize(this)
        class(nurbs_curve), intent(inout) :: this
        if (allocated(this%Xc)) deallocate(this%Xc)
        if (allocated(this%Xg)) deallocate(this%Xg)
        if (allocated(this%Wc)) deallocate(this%Wc)
        if (allocated(this%Xt)) deallocate(this%Xt)
        if (allocated(this%knot)) deallocate(this%knot)
        if (allocated(this%elemConn_Xc_vis)) deallocate(this%elemConn_Xc_vis)
        if (allocated(this%elemConn_Xg_vis)) deallocate(this%elemConn_Xg_vis)
        if (allocated(this%elemConn)) deallocate(this%elemConn)
        this%rational = .false.
        this%wrap_parameters = .false.
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build one-based polyline connectivity for the control polygon.
    !! Optional `p` is the polynomial degree per visualization segment and
    !! defaults to one.
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
    !> Build one-based polyline connectivity for cached geometry points.
    !! Optional `p` defaults to one.
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
    !> Build polyline connectivity over distinct active knot values.
    !! Repeated knots do not add zero-length visualization cells.
    pure function cmp_elem_Xth(this, p) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), optional :: p

        if (.not. this%err%ok) return

        if (present(p)) then
            elemConn = elemConn_C0(active_span_count(this%knot, this%nc, this%degree) + 1, p)
        else
            elemConn = elemConn_C0(active_span_count(this%knot, this%nc, this%degree) + 1, 1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export the control polygon and optional point fields to legacy VTK.
    !! Encoding is `BINARY` by default; `ASCII` is also accepted.
    impure subroutine export_Xc(this, filename, point_data, field_names, encoding)
        class(nurbs_curve), intent(inout) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), contiguous, optional :: point_data(:,:)
        character(len=*), intent(in), contiguous, optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        ! check
        if (.not.allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Control points are not set.",&
                location   = "export_Xc",&
                suggestion = "Call set(...) first before exporting.")
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
    !> Export cached curve samples and optional point fields to legacy VTK.
    !! Call [[nurbs_curve:create]] before this routine.
    impure subroutine export_Xg(this, filename, point_data, field_names, encoding)
        class(nurbs_curve), intent(inout) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), contiguous, optional :: point_data(:,:)
        character(len=*), intent(in), contiguous, optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        ! check
        if (.not.allocated(this%Xg)) then
            call this%err%set(&
                code       = 104,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Geometry points are not set.",&
                location   = "export_Xg",&
                suggestion = "Generate Xg by calling create(...) before exporting.")
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
    !> Export the one-dimensional active knot mesh to legacy VTK.
    !! Knot coordinates are embedded as \((u,0,0)\).
    impure subroutine export_Xth(this, filename, point_data, field_names, encoding)
        class(nurbs_curve), intent(in) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), contiguous, optional :: point_data(:,:)
        character(len=*), intent(in), contiguous, optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Xth(:,:), Xth1(:), Xth2(:), Xth3(:)

        if (.not. this%err%ok) return

        elemConn = this%cmp_elem_Xth()
        Xth1 = active_knots(this%knot, this%nc, this%degree)
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
    !> Export the current curve as IGES entity 126.
    !!
    !! Knot values, degree, control points, weights, and the active parameter
    !! range are transferred to a Rational B-Spline Curve entity. One- and
    !! two-dimensional embeddings are padded with zero coordinates; for
    !! embeddings above three dimensions, only the first three coordinates are
    !! written because entity 126 stores three Cartesian coordinates.
    !! The closed and periodic properties are classified independently of the
    !! parameter-wrapping evaluation option. The curve is closed when
    !! \(\mathbf C(u_s)=\mathbf C(u_e)\) at the exported parameter limits.
    !! It is periodic only when it is closed and the final \(p\) control data
    !! repeat the first \(p\) control data while the knot sequence repeats after
    !! \(n_c-p\) intervals. Rational periodicity additionally requires the
    !! corresponding weights to repeat. Comparisons use roundoff-scaled
    !! tolerances; wrapping is deliberately disabled for endpoint evaluation.
    !! The polynomial/rational property follows [[nurbs_curve:is_rational]],
    !! including its roundoff-scale uniform-weight classification.
    !!
    !! The routine requires initialized curve data. File-system and writer
    !! failures are delegated to `forIGES` and are not converted to
    !! [[nurbs_curve:err]].
    !!
    !! **Format reference:** U.S. National Bureau of Standards, *Initial
    !! Graphics Exchange Specification (IGES), Version 4.0*, NBSIR 88-3813
    !! (1988), Rational B-Spline Curve Entity 126.
    !! [NIST publication](https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir88-3813.pdf).
    impure subroutine export_iges(this, filename)
        use forIGES, only: Gsection_t, Dentry_t, entity126_t, DElist_t, PElist_t,&
                           makeSsection, makeGsection, makeDPsections, writeIGESfile, wp

        class(nurbs_curve), intent(inout) :: this
            !! Valid initialized curve. Existing diagnostic errors suppress export.
        character(len=*), intent(in)      :: filename
            !! Nonempty destination path for the IGES text file.

        type(Gsection_t)  :: G
        type(Dentry_t)    :: D
        type(entity126_t) :: curve126
        type(DElist_t)    :: Dlist
        type(PElist_t)    :: Plist
        character(80), allocatable :: Ssection(:), Gsection(:), Dsection(:), Psection(:), Ssec_out(:)
        integer :: i, K, M, N, prop3, dim_, dim_export
        real(wp), allocatable :: T(:), X(:), Y(:), Z(:), W(:)
        real(wp) :: V(0:1)
        real(wp) :: XNORM, YNORM, ZNORM
        real(rk), allocatable :: Xg_start(:), Xg_end(:)
        real(rk) :: geometry_tol
        logical :: closed, periodic

        if (.not. this%err%ok) return

        ! Parameters for IGES knot vector
        K = this%nc - 1
        M = this%degree
        N = 1 + K - M
        dim_ = size(this%Xc, 2)
        dim_export = min(dim_, 3)

        allocate(T(-M:N+M))
        allocate(X(0:K), Y(0:K), Z(0:K), W(0:K))

        ! Copy your knot vector to IGES indexing
        do i = -M, K + 1
            T(i) = real(this%knot(i + M + 1), kind=wp)
        end do

        ! Copy control points
        X = real(0.0_rk, kind=wp)
        Y = real(0.0_rk, kind=wp)
        Z = real(0.0_rk, kind=wp)
        if (this%is_rational()) then
            do i = 0, K
                if (dim_ >= 1) X(i) = real(this%Xc(i+1, 1), kind=wp)
                if (dim_ >= 2) Y(i) = real(this%Xc(i+1, 2), kind=wp)
                if (dim_ >= 3) Z(i) = real(this%Xc(i+1, 3), kind=wp)
                W(i) = real(this%Wc(i+1), kind=wp)
            end do
            prop3 = 0
        else
            do i = 0, K
                if (dim_ >= 1) X(i) = real(this%Xc(i+1, 1), kind=wp)
                if (dim_ >= 2) Y(i) = real(this%Xc(i+1, 2), kind=wp)
                if (dim_ >= 3) Z(i) = real(this%Xc(i+1, 3), kind=wp)
                W(i) = real(1.0_rk, kind=wp)
            end do
            prop3 = 1
        end if

        XNORM = real(0.0_rk, kind=wp)
        YNORM = real(0.0_rk, kind=wp)
        ZNORM = real(0.0_rk, kind=wp)

        V = real([knot_start(this%knot, this%nc, this%degree), knot_end(this%knot, this%nc, this%degree)], kind=wp)

        allocate(Xg_start(dim_), Xg_end(dim_))
        if (prop3 == 0) then
            Xg_start = compute_Xg(&
                knot_start(this%knot, this%nc, this%degree),&
                this%knot,&
                this%degree,&
                this%nc,&
                this%Xc,&
                this%Wc)
            Xg_end = compute_Xg(&
                knot_end(this%knot, this%nc, this%degree),&
                this%knot,&
                this%degree,&
                this%nc,&
                this%Xc,&
                this%Wc)
        else
            Xg_start = compute_Xg(&
                knot_start(this%knot, this%nc, this%degree),&
                this%knot,&
                this%degree,&
                this%nc,&
                this%Xc)
            Xg_end = compute_Xg(&
                knot_end(this%knot, this%nc, this%degree),&
                this%knot,&
                this%degree,&
                this%nc,&
                this%Xc)
        end if

        geometry_tol = 64.0_rk*real(this%degree+1, rk)*epsilon(1.0_rk)*&
            max(1.0_rk, maxval(abs(this%Xc(:,1:dim_export))))
        closed = all(abs(Xg_start(1:dim_export) - Xg_end(1:dim_export)) <= geometry_tol)

        if (prop3 == 0) then
            periodic = closed .and. periodic_topology(this%knot, this%degree, this%Xc, this%Wc)
        else
            periodic = closed .and. periodic_topology(this%knot, this%degree, this%Xc)
        end if

        ! Initialize IGES entity126 (Rational B-spline Curve).
        curve126%entity_type = 126
        curve126%DEP         = 1
        curve126%form        = 0
        curve126%K           = K
        curve126%M           = M
        curve126%PROP1       = 0
        curve126%PROP2       = merge(1, 0, closed)
        curve126%PROP3       = prop3
        curve126%PROP4       = merge(1, 0, periodic)
        allocate(curve126%T(lbound(T, 1):ubound(T, 1)))
        allocate(curve126%W(0:K), curve126%X(0:K), curve126%Y(0:K), curve126%Z(0:K))
        curve126%T     = T
        curve126%W     = W
        curve126%X     = X
        curve126%Y     = Y
        curve126%Z     = Z
        curve126%V     = V
        curve126%XNORM = XNORM
        curve126%YNORM = YNORM
        curve126%ZNORM = ZNORM

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
        Ssection(1) = "ForCAD"

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
    !> Replace one control point or one of its physical coordinates.
    !! `X` is expected to be finite but is not checked. Cached geometry remains
    !! unchanged; call [[nurbs_curve:create]] to resample after mutation.
    pure subroutine modify_Xc(this,X,num,dir)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: X
        integer, intent(in) :: num
        integer, intent(in) :: dir

        if (.not. this%err%ok) return

        if (allocated(this%Xc)) then
            if (num < lbound(this%Xc,1) .or. num > ubound(this%Xc,1) .or. &
                dir < lbound(this%Xc,2) .or. dir > ubound(this%Xc,2)) then
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Invalid control point index or direction.",&
                    location   = "modify_Xc",&
                    suggestion = "Use 1 <= num <= nc and 1 <= dir <= spatial dimension.")
                return
            end if
            this%Xc(num,dir) = X
        else
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Control points are not set.",&
                location   = "modify_Xc",&
                suggestion = "Call set(...) before modifying it.")
            return
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Replace one finite positive control weight.
    !! The rational-state flag is recomputed. Cached geometry remains unchanged;
    !! call [[nurbs_curve:create]] to resample after mutation.
    pure subroutine modify_Wc(this,W,num)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: W
        integer, intent(in) :: num

        if (.not. this%err%ok) return

        if (allocated(this%Wc)) then
            if (num < lbound(this%Wc,1) .or. num > ubound(this%Wc,1)) then
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Invalid weight index.",&
                    location   = "modify_Wc",&
                    suggestion = "Use 1 <= num <= number of control points.")
                return
            end if
            if (.not. ieee_is_finite(W) .or. W <= 0.0_rk) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Invalid rational weight: weight must be finite and positive.",&
                    location   = "modify_Wc",&
                    suggestion = "Use a strictly positive finite weight value.")
                return
            end if
            this%Wc(num) = W
            this%rational = max(maxval(this%Wc)-this%Wc(1),this%Wc(1)-minval(this%Wc)) > &
                32.0_rk*epsilon(1.0_rk)*maxval(this%Wc)
        else
            call this%err%set(&
                code       = 105,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Weights are not set.",&
                location   = "modify_Wc",&
                suggestion = "Pass Wc when calling set(...), before modifying weights.")
            return
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return multiplicities of every distinct value in the stored knot vector.
    !! Entries follow the nondecreasing order of `unique(this%get_knot())` and
    !! include knots outside the active interval of unclamped or periodic-form
    !! representations. Invalid or unavailable state returns an empty array.
    pure function get_multiplicity(this) result(m)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: m(:)

        if (.not. this%err%ok .or. .not. allocated(this%knot)) then
            allocate(m(0))
            return
        end if

        m = compute_multiplicity(this%knot)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return \(p-s_i\) for every distinct value in the stored knot vector.
    !! At an interior active knot this is the guaranteed continuity of the
    !! spline space; a particular curve can be smoother. Entries at active
    !! endpoints or outside the active interval describe representation
    !! multiplicity, not continuity between two active spans.
    pure function get_continuity(this) result(c)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: c(:)

        if (.not. this%err%ok .or. .not. allocated(this%knot)) then
            allocate(c(0))
            return
        end if

        c = this%degree - compute_multiplicity(this%knot)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Recompute the control count implied by knot-vector size and degree.
    !! The relation is \(n_c=\operatorname{size}(U)-p-1\).
    pure subroutine cmp_nc(this)
        class(nurbs_curve), intent(inout) :: this

        if (.not. this%err%ok) return

        if (.not. allocated(this%knot)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Knot vector is not set.",&
                location   = "cmp_nc",&
                suggestion = "Call set(...) first before computing nc.")
            return
        end if
        this%nc = size(this%knot) - this%degree - 1
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the current control-point count.
    pure elemental function get_nc_all(this) result(nc)
        class(nurbs_curve), intent(in) :: this
        integer :: nc

        nc = 0
        if (.not. this%err%ok) return

        nc = this%nc
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the control count after validating curve direction 1.
    pure elemental function get_nc_dir(this, dir) result(nc)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: dir
        integer :: nc

        nc = 0
        if (.not. this%err%ok) return

        if (dir == 1 .and. allocated(this%knot)) nc = this%nc
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Insert one or more knots without changing the represented curve.
    !!
    !! Requests are applied in array order. `r(i)` copies of `Xth(i)` are added,
    !! and the resulting multiplicity may not exceed `degree+1`; multiplicity
    !! `degree+1` represents a \(C^{-1}\) break. For a verified periodic
    !! topology, the maximum is `degree` and insertion also updates the cyclic
    !! exterior knots and repeated control rows. Let \([a,b]\) be the current
    !! active interval and \(\tau\) the value returned by `knot_tolerance`.
    !! Writing \(u_i\) for `Xth(i)`, validation accepts
    !! \(u_i\in[a-\tau,b+\tau]\), but inserts the supplied value exactly: it is
    !! not snapped to an endpoint and multiplicity uses exact stored-knot
    !! equality.
    !!
    !! Rational curves are refined in homogeneous coordinates. Writing
    !! \(\mathbf P\) and \(\mathbf Q\) for the old and new Euclidean control
    !! arrays, after the new weights are fixed the optional matrix `Bs` satisfies
    !! \(\mathbf Q_{:,c}=\mathrm{Bs}\,\mathbf P_{:,c}\) for every coordinate
    !! \(c\), and has shape `[new_nc,old_nc]`. `B` is the corresponding
    !! coordinate-block operator of shape
    !! `[new_nc*ncoord,old_nc*ncoord]` for control-point-major vectors. For a
    !! rational curve these operators therefore depend on both old and new
    !! weights; they are not the homogeneous refinement operator.
    !!
    !! **Algorithm:** W. Boehm, "Inserting New Knots into B-Spline Curves,"
    !! *Computer-Aided Design* 12 (1980), 199--201.
    !! [doi:10.1016/0010-4485(80)90154-2](https://doi.org/10.1016/0010-4485(80)90154-2).
    !! Implemented in the equivalent form of Piegl and Tiller, Algorithm A5.1.
    pure subroutine insert_knots(this, Xth, r, B, Bs)
        class(nurbs_curve),               intent(inout) :: this
        real(rk), contiguous,             intent(in)    :: Xth(:)
            !! Finite knot values in the tolerance-expanded closed active interval; values are not snapped.
        integer,  contiguous,             intent(in)    :: r(:)
            !! Nonnegative insertion counts corresponding to `Xth`.
        real(rk), allocatable, optional,  intent(out)   :: B(:,:)
            !! Optional coordinate-expanded refinement operator.
        real(rk), allocatable, optional,  intent(out)   :: Bs(:,:)
            !! Optional scalar refinement operator.

        integer :: i, j, k, s, n_new, ncoord, n_old, mS, nS, c, nc_work, weight_exponent
        real(rk) :: lo, hi, ktol
        real(rk), allocatable :: Xc(:,:), Xcw(:,:), Xcw_new(:,:), knot_new(:), knot_work(:)
        real(rk), allocatable :: Wc_new(:), Wc_old(:), S_loc(:,:), A1(:,:), A_re_loc(:,:), A_work(:,:)
        logical :: periodic

        if (.not. this%err%ok) return
        if (size(Xth) /= size(r)) then
            call this%err%set(&
                code       = 101,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Knot insertion input sizes do not match.",&
                location   = "insert_knots",&
                suggestion = "Pass Xth and r arrays with the same size.")
            return
        end if
        if (any(.not. ieee_is_finite(Xth)) .or. any(r < 0)) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Invalid knot insertion request.",&
                location   = "insert_knots",&
                suggestion = "Use finite knot values and nonnegative insertion counts.")
            return
        end if

        ncoord = size(this%Xc,2)

        if (present(B) .or. present(Bs)) then
            n_old = size(this%Xc,1)
            if (this%is_rational()) then
                allocate(Wc_old(n_old)); Wc_old = this%Wc
            end if
            allocate(A1(n_old,n_old), source=0.0_rk)
            do concurrent (j = 1:n_old)
                A1(j,j) = 1.0_rk
            end do
        end if

        if (this%is_rational()) then
            weight_exponent = exponent(maxval(this%Wc))
            allocate(Xcw(size(this%Xc,1), ncoord+1))
            do concurrent (j = 1: size(this%Xc,1))
                Xcw(j,1:ncoord) = this%Xc(j,1:ncoord)*scale(this%Wc(j), -weight_exponent)
            end do
            Xcw(:,ncoord+1) = scale(this%Wc(:), -weight_exponent)
            knot_work = this%knot
            nc_work = size(Xcw,1)
            periodic = periodic_topology(knot_work, this%degree, Xcw)

            do i = 1, size(Xth)
                if (r(i) == 0) cycle
                lo = knot_start(knot_work, nc_work, this%degree)
                hi = knot_end(knot_work, nc_work, this%degree)
                ktol = knot_tolerance(knot_work, this%degree+1, nc_work+1)
                if (Xth(i) < lo - ktol .or. Xth(i) > hi + ktol) then
                    call this%err%set(&
                        code       = 100,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Knot insertion value is outside the active parameter interval.",&
                        location   = "insert_knots",&
                        suggestion = "Insert knots only inside [knot_start, knot_end].")
                    return
                end if
                s = compute_multiplicity(knot_work, Xth(i))
                if (s + r(i) > this%degree + merge(0,1,periodic)) then
                    call this%err%set(&
                        code       = 100,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Knot insertion would exceed the supported knot multiplicity.",&
                        location   = "insert_knots",&
                        suggestion = "Keep the final multiplicity within the current topology's limit.")
                    return
                end if
                k = findspan(nc_work-1, this%degree, Xth(i), knot_work)

                if (periodic) then
                    if (present(B) .or. present(Bs)) then
                        call insert_knot_periodic_A_5_1(this%degree, knot_work, Xcw, Xth(i), k, s, r(i), &
                            n_new, knot_new, Xcw_new, A_re_loc)
                    else
                        call insert_knot_periodic_A_5_1(this%degree, knot_work, Xcw, Xth(i), k, s, r(i), &
                            n_new, knot_new, Xcw_new)
                    end if
                else
                    if (present(B) .or. present(Bs)) then
                        call insert_knot_A_5_1(this%degree, knot_work, Xcw, Xth(i), k, s, r(i), &
                            n_new, knot_new, Xcw_new, A_re_loc)
                    else
                        call insert_knot_A_5_1(this%degree, knot_work, Xcw, Xth(i), k, s, r(i), &
                            n_new, knot_new, Xcw_new)
                    end if
                end if
                if (n_new < 0) then
                    call this%err%set(&
                        code       = 108,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Knot insertion failed to preserve the spline topology.",&
                        location   = "insert_knots",&
                        suggestion = "Check the knot value, multiplicity, and periodic representation.")
                    return
                end if
                if (present(B) .or. present(Bs)) then
                    call sparse_left_matmul(A_re_loc, A1, A_work)
                    call move_alloc(A_work, A1)
                end if

                call move_alloc(Xcw_new, Xcw)
                call move_alloc(knot_new, knot_work)
                nc_work = n_new + 1
            end do
            allocate(Xc(nc_work, ncoord), Wc_new(nc_work))
            do concurrent (i = 1:nc_work)
                Wc_new(i) = Xcw(lbound(Xcw, 1)+i-1, ncoord+1)
            end do
            do concurrent (j = 1:ncoord, i = 1:nc_work)
                Xc(i,j) = Xcw(lbound(Xcw, 1)+i-1, j)/Wc_new(i)
            end do
            Wc_new = scale(Wc_new, weight_exponent)
            call this%set(knot=knot_work, Xc=Xc, Wc=Wc_new, wrap_parameters=this%wrap_parameters)

        else
            Xc = this%Xc
            knot_work = this%knot
            nc_work = size(Xc,1)
            periodic = periodic_topology(knot_work, this%degree, Xc)
            do i = 1, size(Xth)
                if (r(i) == 0) cycle
                lo = knot_start(knot_work, nc_work, this%degree)
                hi = knot_end(knot_work, nc_work, this%degree)
                ktol = knot_tolerance(knot_work, this%degree+1, nc_work+1)
                if (Xth(i) < lo - ktol .or. Xth(i) > hi + ktol) then
                    call this%err%set(&
                        code       = 100,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Knot insertion value is outside the active parameter interval.",&
                        location   = "insert_knots",&
                        suggestion = "Insert knots only inside [knot_start, knot_end].")
                    return
                end if
                s = compute_multiplicity(knot_work, Xth(i))
                if (s + r(i) > this%degree + merge(0,1,periodic)) then
                    call this%err%set(&
                        code       = 100,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Knot insertion would exceed the supported knot multiplicity.",&
                        location   = "insert_knots",&
                        suggestion = "Keep the final multiplicity within the current topology's limit.")
                    return
                end if
                k = findspan(nc_work-1, this%degree, Xth(i), knot_work)

                if (periodic) then
                    if (present(B) .or. present(Bs)) then
                        call insert_knot_periodic_A_5_1(this%degree, knot_work, Xc, Xth(i), k, s, r(i), &
                            n_new, knot_new, Xcw_new, A_re_loc)
                    else
                        call insert_knot_periodic_A_5_1(this%degree, knot_work, Xc, Xth(i), k, s, r(i), &
                            n_new, knot_new, Xcw_new)
                    end if
                else
                    if (present(B) .or. present(Bs)) then
                        call insert_knot_A_5_1(this%degree, knot_work, Xc, Xth(i), k, s, r(i), &
                            n_new, knot_new, Xcw_new, A_re_loc)
                    else
                        call insert_knot_A_5_1(this%degree, knot_work, Xc, Xth(i), k, s, r(i), &
                            n_new, knot_new, Xcw_new)
                    end if
                end if
                if (n_new < 0) then
                    call this%err%set(&
                        code       = 108,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Knot insertion failed to preserve the spline topology.",&
                        location   = "insert_knots",&
                        suggestion = "Check the knot value, multiplicity, and periodic representation.")
                    return
                end if
                if (present(B) .or. present(Bs)) then
                    call sparse_left_matmul(A_re_loc, A1, A_work)
                    call move_alloc(A_work, A1)
                end if

                call move_alloc(Xcw_new, Xc)
                call move_alloc(knot_new, knot_work)
                nc_work = n_new + 1
            end do
            call this%set(knot=knot_work, Xc=Xc, wrap_parameters=this%wrap_parameters)
        end if
        if (present(B) .or. present(Bs)) then
            mS = size(this%Xc,1)
            nS = size(A1,2)
            if (present(B) .and. .not. present(Bs)) then
                allocate(B(mS*ncoord, nS*ncoord), source=0.0_rk)
                if (this%is_rational()) then
                    do concurrent (j = 1:nS, i = 1:mS, c = 1:ncoord, &
                        A1(i,j) < 0.0_rk .or. A1(i,j) > 0.0_rk .or. ieee_is_nan(A1(i,j)))
                        B((i-1)*ncoord + c, (j-1)*ncoord + c) = A1(i,j) * Wc_old(j) / this%Wc(i)
                    end do
                else
                    do concurrent (j = 1:nS, i = 1:mS, c = 1:ncoord, &
                        A1(i,j) < 0.0_rk .or. A1(i,j) > 0.0_rk .or. ieee_is_nan(A1(i,j)))
                        B((i-1)*ncoord + c, (j-1)*ncoord + c) = A1(i,j)
                    end do
                end if
            else
                allocate(S_loc(mS,nS), source=0.0_rk)

                if (this%is_rational()) then
                    do concurrent (j = 1:nS, i = 1:mS, &
                        A1(i,j) < 0.0_rk .or. A1(i,j) > 0.0_rk .or. ieee_is_nan(A1(i,j)))
                        S_loc(i,j) = A1(i,j) * Wc_old(j) / this%Wc(i)
                    end do
                else
                    S_loc = A1
                end if

                if (present(B)) then
                    allocate(B(mS*ncoord, nS*ncoord), source=0.0_rk)
                    do concurrent (j = 1:nS, i = 1:mS, c = 1:ncoord, &
                        S_loc(i,j) < 0.0_rk .or. S_loc(i,j) > 0.0_rk .or. ieee_is_nan(S_loc(i,j)))
                        B((i-1)*ncoord + c, (j-1)*ncoord + c) = S_loc(i,j)
                    end do
                end if
                if (present(Bs)) then
                    call move_alloc(S_loc, Bs)
                else
                    deallocate(S_loc)
                end if
            end if
            if (allocated(Wc_old)) deallocate(Wc_old)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Raise the polynomial degree by `t` while preserving active-domain geometry.
    !! Knot multiplicities are elevated consistently and rational data are
    !! transformed homogeneously. For `t>0`, a verified periodic representation
    !! retains its cyclic knot extension and repeated control rows. Other
    !! non-clamped input is converted to an equivalent open-clamped
    !! representation on the active interval.
    !! `t=0` is an identity operation. Optional `Bs` and `B` have the shapes and
    !! Euclidean fixed-weight mapping semantics documented by
    !! [[nurbs_curve:insert_knots]]; for rational data they are not homogeneous
    !! elevation operators.
    !!
    !! **Algorithm:** H. Prautzsch, "Degree Elevation of B-Spline Curves,"
    !! *Computer Aided Geometric Design* 1 (1984), 193--198.
    !! [doi:10.1016/0167-8396(84)90031-1](https://doi.org/10.1016/0167-8396(84)90031-1).
    !! Implemented in the form of Piegl and Tiller, Algorithm A5.9.
    pure subroutine elevate_degree(this, t, B, Bs)
        class(nurbs_curve),               intent(inout) :: this
        integer,                          intent(in)    :: t
            !! Nonnegative degree increment.
        real(rk), allocatable, optional,  intent(out)   :: B(:,:)
            !! Optional coordinate-expanded elevation operator.
        real(rk), allocatable, optional,  intent(out)   :: Bs(:,:)
            !! Optional scalar elevation operator.

        integer :: ncoord, mS, nS, c, i, j, nc_new, weight_exponent
        real(rk), allocatable :: knot_new(:)
        real(rk), allocatable :: Xc(:,:), Xcw(:,:), Xcw_new(:,:)
        real(rk), allocatable :: Tdir(:,:), Wc_old(:), S_loc(:,:)
        logical :: success

        if (.not. this%err%ok) return
        if (t < 0) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Invalid degree elevation request.",&
                location   = "elevate_degree",&
                suggestion = "Use a nonnegative degree elevation count.")
            return
        end if

        ncoord = size(this%Xc,2)

        if (present(B) .or. present(Bs)) then
            if (this%is_rational()) then
                allocate(Wc_old(size(this%Xc,1)))
                Wc_old = this%Wc
            end if
        end if

        if (this%is_rational()) then
            weight_exponent = exponent(maxval(this%Wc))
            allocate(Xcw(size(this%Xc,1), ncoord+1))
            do concurrent (i = 1:size(this%Xc,1))
                Xcw(i,1:ncoord) = this%Xc(i,1:ncoord)*scale(this%Wc(i), -weight_exponent)
            end do
            Xcw(:,ncoord+1) = scale(this%Wc(:), -weight_exponent)

            if (present(B) .or. present(Bs)) then
                call elevate_degree_A_5_9(&
                    t        = t,&
                    knot     = this%knot,&
                    degree   = this%degree,&
                    Xcw      = Xcw,&
                    nc_new   = nc_new,&
                    knot_new = knot_new,&
                    Xcw_new  = Xcw_new,&
                    Tmap     = Tdir,&
                    success  = success)
            else
                call elevate_degree_A_5_9(&
                    t        = t,&
                    knot     = this%knot,&
                    degree   = this%degree,&
                    Xcw      = Xcw,&
                    nc_new   = nc_new,&
                    knot_new = knot_new,&
                    Xcw_new  = Xcw_new,&
                    success  = success)
            end if

            if (.not. success) then
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Degree elevation failed to construct a valid spline-space transformation.",&
                    location   = "elevate_degree",&
                    suggestion = "Check the knot vector, degree, and active spline domain.")
                return
            end if

            associate (C => Xcw_new(:,1:ncoord), W => Xcw_new(:,ncoord+1))
                do i = 1, ncoord
                    C(:,i) = C(:,i) / W(:)
                end do
                W = scale(W, weight_exponent)
                call this%set(knot=knot_new, Xc=C, Wc=W, wrap_parameters=this%wrap_parameters)
            end associate
            deallocate(Xcw, Xcw_new)

        else
            Xc = this%Xc
            if (present(B) .or. present(Bs)) then
                call elevate_degree_A_5_9(&
                    t        = t,&
                    knot     = this%knot,&
                    degree   = this%degree,&
                    Xcw      = Xc,&
                    nc_new   = nc_new,&
                    knot_new = knot_new,&
                    Xcw_new  = Xcw_new,&
                    Tmap     = Tdir,&
                    success  = success)
            else
                call elevate_degree_A_5_9(&
                    t        = t,&
                    knot     = this%knot,&
                    degree   = this%degree,&
                    Xcw      = Xc,&
                    nc_new   = nc_new,&
                    knot_new = knot_new,&
                    Xcw_new  = Xcw_new,&
                    success  = success)
            end if
            if (.not. success) then
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Degree elevation failed to construct a valid spline-space transformation.",&
                    location   = "elevate_degree",&
                    suggestion = "Check the knot vector, degree, and active spline domain.")
                return
            end if
            call this%set(knot=knot_new, Xc=Xcw_new, wrap_parameters=this%wrap_parameters)
            deallocate(Xc, Xcw_new)
        end if
        if (present(B) .or. present(Bs)) then
            mS = size(this%Xc,1)
            nS = size(Tdir,2)
            if (present(B) .and. .not. present(Bs)) then
                allocate(B(mS*ncoord, nS*ncoord), source=0.0_rk)
                if (this%is_rational()) then
                    do concurrent (j = 1:nS, i = 1:mS, c = 1:ncoord, &
                        Tdir(i,j) < 0.0_rk .or. Tdir(i,j) > 0.0_rk .or. ieee_is_nan(Tdir(i,j)))
                        B((i-1)*ncoord + c, (j-1)*ncoord + c) = Tdir(i,j) * Wc_old(j) / this%Wc(i)
                    end do
                else
                    do concurrent (j = 1:nS, i = 1:mS, c = 1:ncoord, &
                        Tdir(i,j) < 0.0_rk .or. Tdir(i,j) > 0.0_rk .or. ieee_is_nan(Tdir(i,j)))
                        B((i-1)*ncoord + c, (j-1)*ncoord + c) = Tdir(i,j)
                    end do
                end if
            else
                allocate(S_loc(mS,nS), source=0.0_rk)

                if (this%is_rational()) then
                    do concurrent (j = 1:nS, i = 1:mS, &
                        Tdir(i,j) < 0.0_rk .or. Tdir(i,j) > 0.0_rk .or. ieee_is_nan(Tdir(i,j)))
                        S_loc(i,j) = Tdir(i,j) * Wc_old(j) / this%Wc(i)
                    end do
                else
                    S_loc = Tdir
                end if

                if (present(B)) then
                    allocate(B(mS*ncoord, nS*ncoord), source=0.0_rk)
                    do concurrent (j = 1:nS, i = 1:mS, c = 1:ncoord, &
                        S_loc(i,j) < 0.0_rk .or. S_loc(i,j) > 0.0_rk .or. ieee_is_nan(S_loc(i,j)))
                        B((i-1)*ncoord + c, (j-1)*ncoord + c) = S_loc(i,j)
                    end do
                end if
                if (present(Bs)) then
                    call move_alloc(S_loc, Bs)
                else
                    deallocate(S_loc)
                end if
            end if

            if (allocated(Wc_old)) deallocate(Wc_old)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Validate the state required by all curve derivative evaluators.
    pure subroutine check_derivative_state(this, location, valid)
        class(nurbs_curve), intent(inout) :: this
        character(len=*), intent(in) :: location
        logical, intent(out) :: valid

        valid = .false.
        if (.not. this%err%ok) return
        if (.not. allocated(this%knot)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Curve derivatives require a valid knot vector and control count.",&
                location   = location,&
                suggestion = "Set the curve before evaluating basis derivatives.")
            return
        end if
        valid = this%degree >= 0 .and. this%nc > this%degree .and. &
            size(this%knot) == this%nc + this%degree + 1
        if (valid .and. allocated(this%Wc)) then
            valid = size(this%Wc) == this%nc
        end if
        if (.not. valid) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Curve derivative state is inconsistent.",&
                location   = location,&
                suggestion = "Use matching finite knots, degree, control count, and positive weights.")
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate one arbitrary-order curve basis derivative without global work arrays.
    pure subroutine compute_derivative_order_1d_local(Xt, knot, degree, nc, order, wrap_parameters, values, first_active, elem, Wc)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree, nc, order
        logical, intent(in) :: wrap_parameters
        real(rk), intent(out) :: values(:)
        integer, intent(out), optional :: first_active
        integer, intent(in), contiguous, optional :: elem(:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk) :: ders(0:max(0,order),0:degree), identity(0:0,0:0), local_values(0:degree)
        integer :: first, i, local

        values = 0.0_rk
        identity(0,0) = 1.0_rk
        call basis_bspline_der_all_active(&
            Xt     = map_parameter(Xt, knot, nc, degree, wrap_parameters),&
            knot   = knot,&
            nc     = nc,&
            degree = degree,&
            nder   = order,&
            first  = first,&
            ders   = ders)
        if (present(Wc)) then
            call tensor_basis_derivative_local(&
                first1 = first,&
                nc1    = nc,&
                ders1  = ders,&
                first2 = 1,&
                nc2    = 1,&
                ders2  = identity,&
                first3 = 1,&
                nc3    = 1,&
                ders3  = identity,&
                values = local_values,&
                Wc     = Wc)
        else
            call tensor_basis_derivative_local(&
                first1 = first,&
                nc1    = nc,&
                ders1  = ders,&
                first2 = 1,&
                nc2    = 1,&
                ders2  = identity,&
                first3 = 1,&
                nc3    = 1,&
                ders3  = identity,&
                values = local_values)
        end if

        if (present(first_active)) then
            first_active = first
            do i = 0, min(degree,size(values)-1)
                values(i+1) = local_values(i)
            end do
            return
        end if

        if (present(elem)) then
            do i = 1, min(size(values),size(elem))
                local = elem(i) - first
                if (local >= 0 .and. local <= degree) values(i) = local_values(local)
            end do
        else
            do local = 0, degree
                if (first + local >= 1 .and. first + local <= min(nc,size(values))) then
                    values(first+local) = local_values(local)
                end if
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute active arbitrary-order derivatives at explicit curve parameters.
    !! For local offset \(a=0,\ldots,p\) and parameter row \(g\),
    !!
    !! \[
    !! \mathtt{Dgc(a+1,g)}=
    !! \frac{\mathrm d^{k}R_{\mathtt{first\_active(g)}+a}}
    !!      {\mathrm du^{k}},\qquad k=\mathtt{order}.
    !! \]
    !!
    !! Parameter wrapping is applied before support selection. An unwrapped
    !! parameter outside the active interval returns a zero column and
    !! `first_active=1`; that index does not denote a nonzero support there.
    pure subroutine derivative_order_active_vector(this, Xt, order, first_active, Dgc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:)
            !! Finite parameters `[npoint]`.
        integer, intent(in) :: order
            !! Nonnegative derivative order (k); zero returns basis values.
        integer, allocatable, intent(out) :: first_active(:)
            !! First active control-point index for every parameter.
        real(rk), allocatable, intent(out) :: Dgc(:,:)
            !! Local derivatives `[degree+1,npoint]`.
        integer :: point
        logical :: valid

        call check_derivative_state(this, "derivative_order_active_vector", valid)
        if (.not. valid) then
            allocate(first_active(0))
            allocate(Dgc(0,0))
            return
        end if
        if (order < 0 .or. .not. all(ieee_is_finite(Xt))) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Active curve derivative inputs must be finite with nonnegative order.",&
                location   = "derivative_order_active_vector",&
                suggestion = "Use finite parameters and order >= 0.")
            allocate(first_active(0))
            allocate(Dgc(0,0))
            return
        end if

        allocate(Dgc(this%degree+1,size(Xt)), source=0.0_rk)
        allocate(first_active(size(Xt)), source=0)
        if (this%is_rational()) then
#if defined(__NVCOMPILER)
            do point = 1, size(Xt)
#else
            do concurrent (point = 1:size(Xt))
#endif
                call compute_derivative_order_1d_local(&
                    Xt              = Xt(point),&
                    knot            = this%knot,&
                    degree          = this%degree,&
                    nc              = this%nc,&
                    order           = order,&
                    wrap_parameters = this%wrap_parameters,&
                    values          = Dgc(:,point),&
                    first_active    = first_active(point),&
                    Wc              = this%Wc)
            end do
        else
#if defined(__NVCOMPILER)
            do point = 1, size(Xt)
#else
            do concurrent (point = 1:size(Xt))
#endif
                call compute_derivative_order_1d_local(&
                    Xt              = Xt(point),&
                    knot            = this%knot,&
                    degree          = this%degree,&
                    nc              = this%nc,&
                    order           = order,&
                    wrap_parameters = this%wrap_parameters,&
                    values          = Dgc(:,point),&
                    first_active    = first_active(point))
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one active arbitrary-order curve derivative.
    !! Entry `Dgc(a+1)` is the derivative of basis
    !! `first_active+a`, for `a=0:degree`. Wrapping and outside-domain behavior
    !! are identical to the vector overload.
    pure subroutine derivative_order_active_scalar(this, Xt, order, first_active, Dgc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xt
            !! Finite parameter value.
        integer, intent(in) :: order
            !! Nonnegative derivative order.
        integer, intent(out) :: first_active
            !! First active control-point index.
        real(rk), allocatable, intent(out) :: Dgc(:)
            !! Local derivative values `[degree+1]`.
        logical :: valid

        first_active = 0
        call check_derivative_state(this, "derivative_order_active_scalar", valid)
        if (.not. valid) then
            allocate(Dgc(0))
            return
        end if
        if (order < 0 .or. .not. ieee_is_finite(Xt)) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Active curve derivative input must be finite with nonnegative order.",&
                location   = "derivative_order_active_scalar",&
                suggestion = "Use a finite parameter and order >= 0.")
            allocate(Dgc(0))
            return
        end if

        allocate(Dgc(this%degree+1), source=0.0_rk)
        if (this%is_rational()) then
            call compute_derivative_order_1d_local(&
                Xt              = Xt,&
                knot            = this%knot,&
                degree          = this%degree,&
                nc              = this%nc,&
                order           = order,&
                wrap_parameters = this%wrap_parameters,&
                values          = Dgc,&
                first_active    = first_active,&
                Wc              = this%Wc)
        else
            call compute_derivative_order_1d_local(&
                Xt              = Xt,&
                knot            = this%knot,&
                degree          = this%degree,&
                nc              = this%nc,&
                order           = order,&
                wrap_parameters = this%wrap_parameters,&
                values          = Dgc,&
                first_active    = first_active)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one requested derivative order at a vector of curve parameters.
    !! The dense result satisfies
    !! `Dgc(g,i)=d**order R_i(Xt(g))/du**order`. Explicit `Xt` replaces the
    !! parameter cache; otherwise positive `res` generates an inclusive active-
    !! interval grid, and with neither argument the existing cache is reused.
    !! The method updates `ng` and the parameter cache but does not recompute
    !! cached geometry `Xg`. Nonfinite explicit parameters are not diagnosed by
    !! this overload and produce zero basis rows when they cannot select a span.
    pure subroutine derivative_order_vector(this, order, Dgc, res, Xt)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in) :: order
            !! Nonnegative derivative order; zero returns the dense basis.
        real(rk), allocatable, intent(out) :: Dgc(:,:)
            !! Dense derivative matrix `[ng,nc]`.
        integer, intent(in), optional :: res
            !! Optional positive generated parameter count.
        real(rk), intent(in), contiguous, optional :: Xt(:)
            !! Optional explicit parameter vector; takes precedence over `res`.
        integer :: point
        logical :: valid

        call check_derivative_state(this, "derivative_order_vector", valid)
        if (.not. valid) then
            allocate(Dgc(0,0))
            return
        end if
        if (order < 0) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Derivative order must be nonnegative.",&
                location   = "derivative_order_vector",&
                suggestion = "Use order >= 0; order zero returns the basis functions.")
            allocate(Dgc(0,0))
            return
        end if

        if (present(Xt)) then
            this%Xt = Xt
        else if (present(res)) then
            if (res < 1) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Curve derivative resolution must be positive.",&
                    location   = "derivative_order_vector",&
                    suggestion = "Use res >= 1 or supply a nonempty Xt array.")
                allocate(Dgc(0,0))
                return
            end if
            if (.not. allocated(this%Xt)) allocate(this%Xt(res))
            if (size(this%Xt) /= res) then
                deallocate(this%Xt)
                allocate(this%Xt(res))
            end if
            call fill_uniform(&
                a = knot_start(this%knot, this%nc, this%degree),&
                b = knot_end(this%knot, this%nc, this%degree),&
                x = this%Xt)
        else if (.not. allocated(this%Xt)) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "No curve parameter values are available.",&
                location   = "derivative_order_vector",&
                suggestion = "Supply Xt or res before requesting vector derivatives.")
            allocate(Dgc(0,0))
            return
        end if

        this%ng = size(this%Xt)
        allocate(Dgc(this%ng,this%nc), source=0.0_rk)
        if (this%is_rational()) then
#if defined(__NVCOMPILER)
            ! NVFORTRAN 26.3 cannot lower the pure local kernel for stdpar GPU execution.
            do point = 1, this%ng
#else
            do concurrent (point = 1:this%ng)
#endif
                call compute_derivative_order_1d_local(&
                    Xt     = this%Xt(point),&
                    knot   = this%knot,&
                    degree = this%degree,&
                    nc     = this%nc,&
                    order  = order,&
                    wrap_parameters = this%wrap_parameters,&
                    values = Dgc(point,:),&
                    Wc     = this%Wc)
            end do
        else
#if defined(__NVCOMPILER)
            do point = 1, this%ng
#else
            do concurrent (point = 1:this%ng)
#endif
                call compute_derivative_order_1d_local(&
                    Xt     = this%Xt(point),&
                    knot   = this%knot,&
                    degree = this%degree,&
                    nc     = this%nc,&
                    order  = order,&
                    wrap_parameters = this%wrap_parameters,&
                    values = Dgc(point,:))
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one requested derivative order at one curve parameter.
    !! Without `elem`, `Dgc(i)=d**order R_i(Xt)/du**order` for all controls.
    !! With `elem`, entries follow the supplied one-based control-index order;
    !! indices outside `1:nc` return zero. A finite `Xt` is a caller
    !! precondition; a nonfinite value is not diagnosed by this overload.
    pure subroutine derivative_order_scalar(this, Xt, order, Dgc, elem)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xt
            !! Parameter value, mapped according to the wrapping policy.
        integer, intent(in) :: order
            !! Nonnegative derivative order.
        real(rk), allocatable, intent(out) :: Dgc(:)
            !! Dense `[nc]` or selected `[size(elem)]` derivative vector.
        integer, intent(in), contiguous, optional :: elem(:)
            !! Optional ordered one-based control indices.
        logical :: valid

        call check_derivative_state(this, "derivative_order_scalar", valid)
        if (.not. valid) then
            allocate(Dgc(0))
            return
        end if
        if (order < 0) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Derivative order must be nonnegative.",&
                location   = "derivative_order_scalar",&
                suggestion = "Use order >= 0; order zero returns the basis functions.")
            allocate(Dgc(0))
            return
        end if

        if (present(elem)) then
            allocate(Dgc(size(elem)), source=0.0_rk)
        else
            allocate(Dgc(this%nc), source=0.0_rk)
        end if
        if (this%is_rational()) then
            call compute_derivative_order_1d_local(&
                Xt     = Xt,&
                knot   = this%knot,&
                degree = this%degree,&
                nc     = this%nc,&
                order  = order,&
                wrap_parameters = this%wrap_parameters,&
                values = Dgc,&
                elem   = elem,&
                Wc     = this%Wc)
        else
            call compute_derivative_order_1d_local(&
                Xt     = Xt,&
                knot   = this%knot,&
                degree = this%degree,&
                nc     = this%nc,&
                order  = order,&
                wrap_parameters = this%wrap_parameters,&
                values = Dgc,&
                elem   = elem)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate first parametric derivatives at a vector of parameters.
    !!
    !! Explicit `Xt` takes precedence over `res`; otherwise `res` creates an
    !! inclusive uniform grid on the active interval. With neither argument,
    !! the existing parameter cache is reused. Output rows correspond to
    !! parameters and columns to all control variables. This convenience
    !! overload updates the parameter cache and `ng`, but not cached geometry.
    !! It assumes a valid spline state and either explicit parameters, a
    !! positive resolution, or an allocated cache; those preconditions are not
    !! diagnosed here.
    pure subroutine derivative_vector(this, res, Xt, dTgc, Tgc)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: res
            !! Optional positive uniform-grid size.
        real(rk), intent(in), contiguous, optional :: Xt(:)
            !! Optional parameter vector `[ng]`.
        real(rk), allocatable, intent(out) :: dTgc(:,:)
            !! First derivatives \(\partial R_i/\partial u\), shape `[ng,nc]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
            !! Optional basis values computed on the same grid, shape `[ng,nc]`.
        real(rk), allocatable :: Tgc_(:,:)

        if (.not. this%err%ok) return

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) then
                if (size(Xt,1) /= size(this%Xt,1)) deallocate(this%Xt)
            end if
            this%Xt = Xt
        else if (present(res)) then
            if (allocated(this%Xt)) then
                if (size(this%Xt,1) /= res) then
                    deallocate(this%Xt)
                    allocate(this%Xt(res))
                end if
            else
                allocate(this%Xt(res))
            end if
            call fill_uniform(knot_start(this%knot, this%nc, this%degree), knot_end(this%knot, this%nc, this%degree), this%Xt)
        end if

        ! Set number of geometry points
        this%ng = size(this%Xt,1)

        if (this%is_rational()) then ! NURBS
            if (present(Tgc)) then
                call compute_dTgc(this%Xt, this%knot, this%degree, this%nc, this%ng, this%Wc, dTgc, Tgc, this%wrap_parameters)
            else
                call compute_dTgc(this%Xt, this%knot, this%degree, this%nc, this%ng, this%Wc, dTgc, Tgc_, this%wrap_parameters)
            end if
        else ! B-Spline
            if (present(Tgc)) then
                call compute_dTgc(this%Xt, this%knot, this%degree, this%nc, this%ng, dTgc, Tgc, this%wrap_parameters)
            else
                call compute_dTgc(this%Xt, this%knot, this%degree, this%nc, this%ng, dTgc, Tgc_, this%wrap_parameters)
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate first parametric derivatives at one parameter.
    !! Optional `elem` selects and orders a local subset of control variables.
    !! An initialized spline state is an unchecked precondition of this
    !! convenience overload.
    pure subroutine derivative_scalar(this, Xt, dTgc, Tgc, elem)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xt
            !! Parameter value, mapped according to the wrapping policy.
        integer, intent(in), optional :: elem(:)
            !! Optional one-based control indices.
        real(rk), allocatable, intent(out) :: dTgc(:)
            !! Selected derivatives \(\partial R_i/\partial u\).
        real(rk), allocatable, intent(out), optional :: Tgc(:)
            !! Optional selected basis values.
        real(rk), allocatable :: Tgc_(:)

        if (.not. this%err%ok) return

        if (this%is_rational()) then ! NURBS
            if (present(elem)) then
                if (present(Tgc)) then
                    call compute_dTgc(Xt, this%knot, this%degree, this%nc, this%Wc, dTgc, Tgc, elem, this%wrap_parameters)
                else
                    call compute_dTgc(Xt, this%knot, this%degree, this%nc, this%Wc, dTgc, Tgc_, elem, this%wrap_parameters)
                end if
            else
                if (present(Tgc)) then
                    call compute_dTgc(&
                        Xt              = Xt,&
                        knot            = this%knot,&
                        degree          = this%degree,&
                        nc              = this%nc,&
                        Wc              = this%Wc,&
                        dTgc            = dTgc,&
                        Tgc             = Tgc,&
                        wrap_parameters = this%wrap_parameters)
                else
                    call compute_dTgc(&
                        Xt              = Xt,&
                        knot            = this%knot,&
                        degree          = this%degree,&
                        nc              = this%nc,&
                        Wc              = this%Wc,&
                        dTgc            = dTgc,&
                        Tgc             = Tgc_,&
                        wrap_parameters = this%wrap_parameters)
                end if
            end if
        else ! B-Spline
            if (present(Tgc)) then
                call compute_dTgc(Xt, this%knot, this%degree, this%nc, dTgc, Tgc, elem, this%wrap_parameters)
            else
                call compute_dTgc(Xt, this%knot, this%degree, this%nc, dTgc, Tgc_, elem, this%wrap_parameters)
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate second parametric derivatives at a vector of parameters.
    !! Optional lower derivatives and values are returned from the same
    !! rational or polynomial evaluation path. Parameter selection, cache side
    !! effects, and unchecked convenience-overload preconditions are the same
    !! as for the vector [[nurbs_curve:derivative]] overload.
    pure subroutine derivative2_vector(this, res, Xt, d2Tgc, dTgc, Tgc)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: res
            !! Optional positive uniform-grid size.
        real(rk), intent(in), contiguous, optional :: Xt(:)
            !! Optional parameter vector `[ng]`; takes precedence over `res`.
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
            !! Second derivatives \(\partial^2R_i/\partial u^2\), shape `[ng,nc]`.
        real(rk), allocatable, intent(out), optional :: dTgc(:,:)
            !! Optional first derivatives, shape `[ng,nc]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
            !! Optional basis values, shape `[ng,nc]`.
        real(rk), allocatable :: dTgc_(:,:), Tgc_(:,:)

        if (.not. this%err%ok) return

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) then
                if (size(Xt,1) /= size(this%Xt,1)) deallocate(this%Xt)
            end if
            this%Xt = Xt
        else if (present(res)) then
            if (allocated(this%Xt)) then
                if (size(this%Xt,1) /= res) then
                    deallocate(this%Xt)
                    allocate(this%Xt(res))
                end if
            else
                allocate(this%Xt(res))
            end if
            call fill_uniform(knot_start(this%knot, this%nc, this%degree), knot_end(this%knot, this%nc, this%degree), this%Xt)
        end if

        ! Set number of geometry points
        this%ng = size(this%Xt,1)

        if (this%is_rational()) then ! NURBS
            if (present(dTgc)) then
                if (present(Tgc)) then
                    call compute_d2Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, this%Wc, &
                        d2Tgc, dTgc, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, this%Wc, &
                        d2Tgc, dTgc, Tgc_, this%wrap_parameters)
                end if
            else
                if (present(Tgc)) then
                    call compute_d2Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, this%Wc, &
                        d2Tgc, dTgc_, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, this%Wc, &
                        d2Tgc, dTgc_, Tgc_, this%wrap_parameters)
                end if
            end if
        else ! B-Spline
            if (present(dTgc)) then
                if (present(Tgc)) then
                    call compute_d2Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, &
                        d2Tgc, dTgc, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, &
                        d2Tgc, dTgc, Tgc_, this%wrap_parameters)
                end if
            else
                if (present(Tgc)) then
                    call compute_d2Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, &
                        d2Tgc, dTgc_, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, &
                        d2Tgc, dTgc_, Tgc_, this%wrap_parameters)
                end if
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate second parametric derivatives at one parameter.
    !! An initialized spline state is an unchecked precondition of this
    !! convenience overload.
    pure subroutine derivative2_scalar(this, Xt, d2Tgc, dTgc, Tgc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xt
            !! Parameter value, mapped according to the wrapping policy.
        real(rk), allocatable, intent(out) :: d2Tgc(:)
            !! Dense second-derivative vector `[nc]`.
        real(rk), allocatable, intent(out), optional :: dTgc(:)
            !! Optional dense first-derivative vector `[nc]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:)
            !! Optional dense basis vector `[nc]`.
        real(rk), allocatable :: dTgc_(:), Tgc_(:)

        if (.not. this%err%ok) return

        if (this%is_rational()) then ! NURBS
            if (present(dTgc)) then
                if (present(Tgc)) then
                    call compute_d2Tgc(Xt, this%knot, this%degree, this%nc, this%Wc, d2Tgc, dTgc, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(Xt, this%knot, this%degree, this%nc, this%Wc, d2Tgc, dTgc, Tgc_, this%wrap_parameters)
                end if
            else
                if (present(Tgc)) then
                    call compute_d2Tgc(Xt, this%knot, this%degree, this%nc, this%Wc, d2Tgc, dTgc_, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(Xt, this%knot, this%degree, this%nc, this%Wc, d2Tgc, dTgc_, Tgc_, this%wrap_parameters)
                end if
            end if
        else ! B-Spline
            if (present(dTgc)) then
                if (present(Tgc)) then
                    call compute_d2Tgc(Xt, this%knot, this%degree, this%nc, d2Tgc, dTgc, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(Xt, this%knot, this%degree, this%nc, d2Tgc, dTgc, Tgc_, this%wrap_parameters)
                end if
            else
                if (present(Tgc)) then
                    call compute_d2Tgc(Xt, this%knot, this%degree, this%nc, d2Tgc, dTgc_, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(Xt, this%knot, this%degree, this%nc, d2Tgc, dTgc_, Tgc_, this%wrap_parameters)
                end if
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the dense basis at a vector of parameters.
    !!
    !! Explicit `Xt` takes precedence over `res`; otherwise `res` creates an
    !! inclusive uniform grid. With neither argument, the current parameter
    !! cache is reused. This convenience overload updates the parameter cache
    !! and `ng`, but not cached geometry. It assumes a valid spline state and
    !! available explicit, generated, or cached parameters; that precondition
    !! is not diagnosed here.
    pure subroutine basis_vector(this, res, Xt, Tgc)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: res
            !! Optional positive uniform-grid size.
        real(rk), intent(in), contiguous, optional :: Xt(:)
            !! Optional parameter vector `[ng]`.
        real(rk), allocatable, intent(out) :: Tgc(:,:)
            !! Basis matrix `[ng,nc]`; each valid row sums to one.
        if (.not. this%err%ok) return

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) then
                if (size(Xt,1) /= size(this%Xt,1)) deallocate(this%Xt)
            end if
            this%Xt = Xt
        else if (present(res)) then
            if (allocated(this%Xt)) then
                if (size(this%Xt,1) /= res) then
                    deallocate(this%Xt)
                    allocate(this%Xt(res))
                end if
            else
                allocate(this%Xt(res))
            end if
            call fill_uniform(knot_start(this%knot, this%nc, this%degree), knot_end(this%knot, this%nc, this%degree), this%Xt)
        end if

        ! Set number of geometry points
        this%ng = size(this%Xt,1)

        if (this%is_rational()) then ! NURBS
            Tgc = compute_Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, this%Wc, this%wrap_parameters)
        else ! B-Spline
            Tgc = compute_Tgc(this%Xt, this%knot, this%degree, this%nc, this%ng, this%wrap_parameters)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the basis at one parameter.
    !! Optional `elem` selects and orders a local subset of control variables.
    !! An initialized spline state is an unchecked precondition of this
    !! convenience overload.
    pure subroutine basis_scalar(this, Xt, Tgc, elem)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xt
            !! Parameter value, mapped according to the wrapping policy.
        integer, intent(in), optional :: elem(:)
            !! Optional one-based control indices.
        real(rk), allocatable, intent(out) :: Tgc(:)
            !! Dense `[nc]` or selected `[size(elem)]` basis vector.

        if (.not. this%err%ok) return

        if (this%is_rational()) then ! NURBS
            if (present(elem)) then
                Tgc = compute_Tgc(Xt, this%knot, this%degree, this%nc, this%Wc, elem, this%wrap_parameters)
            else
                Tgc = compute_Tgc(Xt, this%knot, this%degree, this%nc, this%Wc, wrap_parameters=this%wrap_parameters)
            end if
        else ! B-Spline
            if (present(elem)) then
                Tgc = compute_Tgc(Xt, this%knot, this%degree, this%nc, elem, this%wrap_parameters)
            else
                Tgc = compute_Tgc(Xt, this%knot, this%degree, this%nc, wrap_parameters=this%wrap_parameters)
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report whether nonuniform explicit weights change the basis.
    !! With positive stored weights, the result is true precisely when
    !!
    !! \[
    !! \max_i|w_i-w_1|>
    !! 32\,\epsilon_{rk}\max_iw_i.
    !! \]
    !!
    !! Weights within that threshold are evaluated through the polynomial path,
    !! not merely reported as approximately uniform. An active diagnostic or
    !! absent explicit weights returns false.
    pure elemental function is_rational(this) result(r)
        class(nurbs_curve), intent(in) :: this
        logical :: r

        r = .false.
        if (.not. this%err%ok) return
        r = this%rational .and. allocated(this%Wc)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report whether modulo parameter mapping is enabled.
    pure elemental function get_parameter_wrapping(this) result(wrap_parameters)
        class(nurbs_curve), intent(in) :: this
        logical :: wrap_parameters

        wrap_parameters = .false.
        if (.not. this%err%ok) return
        wrap_parameters = this%wrap_parameters
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Classify the stored parameter topology.
    !!
    !! The result is `"periodic"` only when the knot extension, repeated
    !! control points, and any stored weights satisfy the cyclic identities in
    !! [[periodic_topology]]. Every other valid curve is `"bounded"`; an
    !! incomplete or inconsistent object is `"invalid"`. Parameter wrapping is
    !! an independent evaluation policy and does not affect this result.
    pure function get_parameter_topology(this) result(topology)
        class(nurbs_curve), intent(in) :: this
        character(len=8) :: topology

        topology = "invalid"
        if (.not. this%is_valid()) return

        topology = "bounded"
        if (allocated(this%Wc)) then
            if (periodic_topology(this%knot, this%degree, this%Xc, this%Wc)) topology = "periodic"
        else
            if (periodic_topology(this%knot, this%degree, this%Xc)) topology = "periodic"
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Replace one-based control-polygon visualization connectivity.
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
    !> Replace one-based sampled-curve visualization connectivity.
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
    !> Replace one-based IGA element connectivity after shape validation.
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
    !> Return a copy of control-polygon visualization connectivity.
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
    !> Return a copy of sampled-curve visualization connectivity.
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
    !> Return a copy of active-span IGA element connectivity.
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
    !> Attempt geometry-preserving removal of interior knots.
    !!
    !! Up to `r(i)` copies of `Xth(i)` are attempted in array order. `Xth(i)`
    !! must equal an existing stored knot exactly and must lie more than the
    !! current `knot_tolerance` from either active endpoint. The algorithm
    !! accepts the longest removable prefix of each request and stops that
    !! request at the first failed componentwise, roundoff-scaled control test.
    !! Accepted removals update knots, controls, and weights together. For
    !! rational curves the test is performed on homogeneous controls and is
    !! not a certified Euclidean geometry-error bound.
    !!
    !! Verified periodic representations are rejected when any removal is
    !! requested. Algorithm A5.8 is a bounded-curve algorithm and does not
    !! update the cyclic knot extension and repeated control rows; applying it
    !! directly would silently destroy the periodic topology.
    !!
    !! **Algorithm:** W. Tiller, "Knot-removal Algorithms for NURBS Curves and
    !! Surfaces," *Computer-Aided Design* 24 (1992), 445--453.
    !! [doi:10.1016/0010-4485(92)90012-Y](https://doi.org/10.1016/0010-4485(92)90012-Y).
    !! Implemented in the form of Piegl and Tiller, Algorithm A5.8.
    pure subroutine remove_knots(this,Xth,r)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xth(:)
            !! Exact stored values of interior knots outside the endpoint-tolerance bands.
        integer, intent(in), contiguous :: r(:)
            !! Nonnegative maximum removal counts not exceeding current multiplicity.
        integer :: k, i, s, d, j, nc_new, t, nc_work, weight_exponent
        real(rk) :: lo, hi, ktol
        real(rk), allocatable :: Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:), knot_work(:)
        logical :: changed

        if (.not. this%err%ok) return
        if (size(Xth) /= size(r)) then
            call this%err%set(&
                code       = 101,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Knot removal input sizes do not match.",&
                location   = "remove_knots",&
                suggestion = "Pass Xth and r arrays with the same size.")
            return
        end if
        if (any(.not. ieee_is_finite(Xth)) .or. any(r < 0)) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Invalid knot removal request.",&
                location   = "remove_knots",&
                suggestion = "Use finite knot values and nonnegative removal counts.")
            return
        end if
        if (any(r > 0) .and. this%get_parameter_topology() == "periodic") then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Periodic knot removal is not supported without preserving cyclic topology.",&
                location   = "remove_knots",&
                suggestion = "Retain the periodic knots or first construct an explicit bounded representation.")
            return
        end if

        if (this%is_rational()) then ! NURBS

            d = size(this%Xc,2)
            weight_exponent = exponent(maxval(this%Wc))
            allocate(Xcw(size(this%Xc,1),d+1))
            do j = 1, size(this%Xc,1)
                Xcw(j,1:d) = this%Xc(j,1:d)*scale(this%Wc(j), -weight_exponent)
            end do
            Xcw(:,d+1) = scale(this%Wc(:), -weight_exponent)
            knot_work = this%knot
            nc_work = size(Xcw,1)
            changed = .false.

            do i = 1, size(Xth)
                if (r(i) == 0) cycle
                lo = knot_start(knot_work, nc_work, this%degree)
                hi = knot_end(knot_work, nc_work, this%degree)
                ktol = knot_tolerance(knot_work, this%degree+1, nc_work+1)
                if (Xth(i) < lo - ktol .or. Xth(i) > hi + ktol) then
                    call this%err%set(&
                        code       = 100,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Knot removal value is outside the active parameter interval.",&
                        location   = "remove_knots",&
                        suggestion = "Remove knots only inside [knot_start, knot_end].")
                    return
                end if
                if (Xth(i) <= lo + ktol .or. Xth(i) >= hi - ktol) then
                    call this%err%set(&
                        code       = 100,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Active-boundary knots cannot be removed by geometry-preserving knot removal.",&
                        location   = "remove_knots",&
                        suggestion = "Remove only knots strictly inside the active parameter interval.")
                    return
                end if
                s = compute_multiplicity(knot_work,Xth(i))
                if (s == 0 .or. r(i) > s) then
                    call this%err%set(&
                        code       = 100,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Requested knot removal multiplicity is not available.",&
                        location   = "remove_knots",&
                        suggestion = "Remove only knots present in the knot vector and use r <= current multiplicity.")
                    return
                end if
                k = findspan(nc_work-1,this%degree,Xth(i),knot_work)
                k = k + 1

                call remove_knots_A_5_8(&
                    p        = this%degree,&
                    knot     = knot_work,&
                    Pw       = Xcw,&
                    u        = Xth(i),&
                    r        = k,&
                    s        = s,&
                    num      = r(i),&
                    t        = t,&
                    knot_new = knot_new,&
                    Pw_new   = Xcw_new)

                if (t > 0) then
                    call move_alloc(Xcw_new, Xcw)
                    call move_alloc(knot_new, knot_work)
                    nc_work = size(Xcw,1)
                    changed = .true.
                end if
            end do

            if (changed) then
                nc_new = size(Xcw,1)
                allocate(Xc_new(nc_new,d), Wc_new(nc_new))
                do concurrent (j = 1:d, i = 1:nc_new)
                    Xc_new(i,j) = Xcw(i,j)/Xcw(i,d+1)
                end do
                Wc_new(:) = scale(Xcw(:,d+1), weight_exponent)
                call this%set(knot=knot_work, Xc=Xc_new, Wc=Wc_new, wrap_parameters=this%wrap_parameters)
            end if

        else ! B-Spline

            Xc_new = this%Xc
            knot_work = this%knot
            nc_work = size(Xc_new,1)
            changed = .false.
            do i = 1, size(Xth)
                if (r(i) == 0) cycle
                lo = knot_start(knot_work, nc_work, this%degree)
                hi = knot_end(knot_work, nc_work, this%degree)
                ktol = knot_tolerance(knot_work, this%degree+1, nc_work+1)
                if (Xth(i) < lo - ktol .or. Xth(i) > hi + ktol) then
                    call this%err%set(&
                        code       = 100,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Knot removal value is outside the active parameter interval.",&
                        location   = "remove_knots",&
                        suggestion = "Remove knots only inside [knot_start, knot_end].")
                    return
                end if
                if (Xth(i) <= lo + ktol .or. Xth(i) >= hi - ktol) then
                    call this%err%set(&
                        code       = 100,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Active-boundary knots cannot be removed by geometry-preserving knot removal.",&
                        location   = "remove_knots",&
                        suggestion = "Remove only knots strictly inside the active parameter interval.")
                    return
                end if
                s = compute_multiplicity(knot_work,Xth(i))
                if (s == 0 .or. r(i) > s) then
                    call this%err%set(&
                        code       = 100,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Requested knot removal multiplicity is not available.",&
                        location   = "remove_knots",&
                        suggestion = "Remove only knots present in the knot vector and use r <= current multiplicity.")
                    return
                end if
                k = findspan(nc_work-1,this%degree,Xth(i),knot_work)
                k = k + 1

                call remove_knots_A_5_8(&
                    p        = this%degree,&
                    knot     = knot_work,&
                    Pw       = Xc_new,&
                    u        = Xth(i),&
                    r        = k,&
                    s        = s,&
                    num      = r(i),&
                    t        = t,&
                    knot_new = knot_new,&
                    Pw_new   = Xcw_new)

                if (t > 0) then
                    call move_alloc(Xcw_new, Xc_new)
                    call move_alloc(knot_new, knot_work)
                    nc_work = size(Xc_new,1)
                    changed = .true.
                end if
            end do

            if (changed) call this%set(knot=knot_work, Xc=Xc_new, wrap_parameters=this%wrap_parameters)

        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct an exact planar rational circle.
    !!
    !! The circle is represented by three quadratic rational arcs over
    !! `[0,1]`, with seven control points and positive weights. The result is
    !! embedded in the xy plane at `center(3)`. The geometric radius is
    !! `abs(radius)`; a negative value rotates the parameterized template by
    !! \(\pi\), while zero collapses it to the center.
    !!
    !! The exact rational construction follows the conic representation in
    !! Piegl and Tiller, *The NURBS Book*, 2nd ed., Springer, 1997, Chapter 7.
    pure subroutine set_circle(this, center, radius)
        class(nurbs_curve), intent(inout) :: this
            !! Curve replaced by the constructed circle.
        real(rk), intent(in), contiguous :: center(:)
            !! Circle center; at least three finite components are required.
        real(rk), intent(in) :: radius
            !! Finite signed scale; use a nonzero value for a nondegenerate circle.
        real(rk), allocatable :: Xc(:,:), Wc(:), knot(:)
        integer :: i, d

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
        do concurrent (i = 1:size(Xc, 1), d = 1:size(Xc, 2))
            Xc(i,d) = center(d) + Xc(i,d)*radius
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
    !> Construct an exact planar rational C-shaped circular arc.
    !!
    !! The open curve consists of two quadratic rational arcs spanning
    !! 240 degrees. It is embedded in the xy plane at `center(3)`. Its geometric
    !! radius is `abs(radius)`; a negative value rotates the template by
    !! \(\pi\), while zero collapses it to the center.
    !!
    pure subroutine set_C(this, center, radius)
        class(nurbs_curve), intent(inout) :: this
            !! Curve replaced by the constructed arc.
        real(rk), intent(in), contiguous :: center(:)
            !! Circle center; at least three finite components are required.
        real(rk), intent(in) :: radius
            !! Finite signed scale; use a nonzero value for a nondegenerate arc.
        real(rk), allocatable :: Xc(:,:), Wc(:), knot(:)
        integer :: i, d

        if (.not. this%err%ok) return

        ! Define control points for C-shape
        allocate(Xc(5, 3))
        Xc(1,:)= [ 1.0_rk,  0.0_rk,              0.0_rk]
        Xc(2,:)= [ 1.0_rk,  sqrt(3.0_rk),        0.0_rk]
        Xc(3,:)= [-0.5_rk,  sqrt(3.0_rk)/2.0_rk, 0.0_rk]
        Xc(4,:)= [-2.0_rk,  0.0_rk,              0.0_rk]
        Xc(5,:)= [-0.5_rk, -sqrt(3.0_rk)/2.0_rk, 0.0_rk]

        ! Scale and translate the control points
        do concurrent (i = 1:size(Xc, 1), d = 1:size(Xc, 2))
            Xc(i,d) = center(d) + Xc(i,d)*radius
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
    !> Build IGA connectivity for every positive-length active knot span.
    !! Each row contains `degree+1` one-based stored-control indices. For a
    !! verified periodic topology, apply [[nurbs_curve:cmp_dof_map]] to this
    !! array before assembling a system with independent cyclic degrees of
    !! freedom.
    pure function cmp_elem(this) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        call elemConn_Cn(this%nc, this%degree, active_knots(this%knot, this%nc, this%degree), &
            active_knot_multiplicity(this%knot, this%nc, this%degree),&
            elemConn)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Map stored controls to independent cyclic IGA degrees of freedom.
    !!
    !! Bounded curves return the identity map. A verified degree-\(p\)
    !! periodic representation stores \(p\) repeated trailing controls, so its
    !! map has `nc` entries but only `nc-degree` distinct identifiers. Apply
    !! the returned one-based map to [[nurbs_curve:cmp_elem]] connectivity.
    pure function cmp_dof_map(this) result(map)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: map(:)
        integer :: i, independent_count

        if (.not. this%is_valid()) return

        allocate(map(this%nc))
        independent_count = this%nc
        if (this%get_parameter_topology() == "periodic") independent_count = this%nc - this%degree
        do concurrent (i = 1:this%nc)
            map(i) = modulo(i-1, independent_count) + 1
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Rotate the control polygon using active extrinsic x-y-z Euler rotations.
    !! Angles `alpha`, `beta`, and `theta` are in degrees and apply
    !! \(R_z(\theta)R_y(\beta)R_x(\alpha)\). Cached physical samples are not
    !! changed. See [[rotation]] for the behavior of embeddings with fewer than
    !! three coordinates.
    pure subroutine rotate_Xc(this, alpha, beta, theta)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: alpha, beta, theta

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) return

        call rotate_points(this%Xc, this%nc, alpha, beta, theta)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Rotate cached geometry using active extrinsic x-y-z Euler rotations.
    !! Angles are in degrees; control points are not changed.
    pure subroutine rotate_Xg(this, alpha, beta, theta)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: alpha, beta, theta

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xg)) return

        call rotate_points(this%Xg, this%ng, alpha, beta, theta)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Translate every control point by `vec`.
    !! `vec` must have at least the stored coordinate dimension; extra entries
    !! are ignored and finiteness is not checked. Cached samples are unchanged.
    pure subroutine translate_Xc(this, vec)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: vec(:)
        integer :: i, d

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) return

        do concurrent (i = 1:this%nc, d = 1:size(this%Xc,2))
            this%Xc(i,d) = this%Xc(i,d) + vec(d)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Translate every cached geometry point by `vec`.
    !! `vec` must have at least the cached coordinate dimension; extra entries
    !! are ignored and finiteness is not checked. Control points are unchanged.
    pure subroutine translate_Xg(this, vec)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: vec(:)
        integer :: i, d

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xg)) return

        do concurrent (i = 1:this%ng, d = 1:size(this%Xg,2))
            this%Xg(i,d) = this%Xg(i,d) + vec(d)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Display previously exported curve representations with PyVista.
    !!
    !! The viewer discovers scalar point fields stored in the sampled `Xg`
    !! file and exposes field and colormap selectors in its window. This
    !! procedure never exports files.
    impure subroutine show(this, vtkfile_Xc, vtkfile_Xg, vtkfile_Xth_in_Xg)
        class(nurbs_curve), intent(inout) :: this
        character(len=*), intent(in) :: vtkfile_Xc
            !! Existing control-geometry VTK path.
        character(len=*), intent(in) :: vtkfile_Xg
            !! Existing sampled-geometry VTK path containing optional point fields.
        character(len=*), intent(in), optional :: vtkfile_Xth_in_Xg
            !! Optional existing parameter-line VTK path.

        if (.not. this%err%ok) return

#ifndef NOSHOW_PYVISTA
        call show_pyvista_singlepatch(&
            vtkfile_Xc        = vtkfile_Xc,&
            vtkfile_Xg        = vtkfile_Xg,&
            vtkfile_Xth_in_Xg = vtkfile_Xth_in_Xg,&
            rank_name         = "curve")
#endif
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct an exact planar rational semicircle.
    !!
    !! Two quadratic rational quarter-arcs span the upper half-plane over
    !! `[0,1]`; the result is embedded at `center(3)`. A negative template scale
    !! rotates the semicircle by \(\pi\); zero collapses it to the center.
    !!
    pure subroutine set_half_circle(this, center, radius)
        class(nurbs_curve), intent(inout) :: this
            !! Curve replaced by the constructed semicircle.
        real(rk), intent(in), contiguous :: center(:)
            !! Circle center; at least three finite components are required.
        real(rk), intent(in) :: radius
            !! Finite legacy template scale. The represented radius is `abs(radius)/2`; use a nonzero value for a
            !! nondegenerate curve.
        real(rk), allocatable :: Xc(:,:), Wc(:), knot(:)
        integer :: i, d

        if (.not. this%err%ok) return

        ! Define control points for half circle
        allocate(Xc(5, 3))
        Xc(1,:) = [ 0.5_rk, 0.0_rk, 0.0_rk]
        Xc(2,:) = [ 0.5_rk, 0.5_rk, 0.0_rk]
        Xc(3,:) = [ 0.0_rk, 0.5_rk, 0.0_rk]
        Xc(4,:) = [-0.5_rk, 0.5_rk, 0.0_rk]
        Xc(5,:) = [-0.5_rk, 0.0_rk, 0.0_rk]

        ! Scale and translate the control points
        do concurrent (i = 1:size(Xc, 1), d = 1:size(Xc, 2))
            Xc(i,d) = center(d) + Xc(i,d)*radius
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
    !> Find the closest point among cached curve samples.
    !!
    !! This is a discrete search, not continuous projection. Call `create`
    !! first; accuracy is limited by the cached sampling. Squared Euclidean
    !! distance uses the first
    !! `min(size(point_Xg),physical_dimension)` coordinates. Extra query
    !! coordinates are ignored, while missing physical coordinates are omitted
    !! from the metric. A returned physical vector longer than the embedding is
    !! padded with zeros. Equal-distance ties select the earliest cached sample.
    !! The curve must have a valid allocated cache and the query must be
    !! nonempty; these preconditions are not revalidated here. If `err` is
    !! already active, the routine returns before defining optional outputs.
    pure subroutine nearest_point(this, point_Xg, nearest_Xg, nearest_Xt, id)
        class(nurbs_curve), intent(in) :: this
        real(rk), intent(in), contiguous :: point_Xg(:)
            !! Query point.
        real(rk), intent(out), optional :: nearest_Xg(size(point_Xg))
            !! Optional closest cached point.
        real(rk), intent(out), optional :: nearest_Xt
            !! Optional parameter of the selected sample.
        integer, intent(out), optional :: id
            !! Optional one-based cached sample identifier.
        integer :: id_, i, d, dim_, ncopy, npoint
        real(rk) :: dist, best_dist, dx

        if (.not. this%err%ok) return

        id_ = 1
        npoint = min(size(this%Xg,1), size(this%Xt))
        dim_ = min(size(point_Xg), size(this%Xg,2))
        ncopy = min(size(point_Xg), size(this%Xg,2))
        best_dist = huge(1.0_rk)
        select case (dim_)
        case (1)
            do i = 1, npoint
                dx = this%Xg(i,1) - point_Xg(1)
                dist = dx*dx
                if (dist < best_dist) then
                    best_dist = dist
                    id_ = i
                    if (best_dist <= 0.0_rk) exit
                end if
            end do
        case (2)
            do i = 1, npoint
                dx = this%Xg(i,1) - point_Xg(1)
                dist = dx*dx
                dx = this%Xg(i,2) - point_Xg(2)
                dist = dist + dx*dx
                if (dist < best_dist) then
                    best_dist = dist
                    id_ = i
                    if (best_dist <= 0.0_rk) exit
                end if
            end do
        case (3)
            do i = 1, npoint
                dx = this%Xg(i,1) - point_Xg(1)
                dist = dx*dx
                dx = this%Xg(i,2) - point_Xg(2)
                dist = dist + dx*dx
                dx = this%Xg(i,3) - point_Xg(3)
                dist = dist + dx*dx
                if (dist < best_dist) then
                    best_dist = dist
                    id_ = i
                    if (best_dist <= 0.0_rk) exit
                end if
            end do
        case default
            do i = 1, npoint
                dist = 0.0_rk
                do d = 1, dim_
                    dx = this%Xg(i,d) - point_Xg(d)
                    dist = dist + dx*dx
                end do
                if (dist < best_dist) then
                    best_dist = dist
                    id_ = i
                    if (best_dist <= 0.0_rk) exit
                end if
            end do
        end select
        if (present(id)) id = id_
        if (present(nearest_Xg)) then
            do concurrent (d = 1:ncopy)
                nearest_Xg(d) = this%Xg(id_,d)
            end do
            do concurrent (d = ncopy+1:size(nearest_Xg))
                nearest_Xg(d) = 0.0_rk
            end do
        end if
        if (present(nearest_Xt)) nearest_Xt = this%Xt(id_)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute a continuous nearest-point candidate on the active curve domain.
    !!
    !! Ten streamed parameter seeds choose an initial point. Projected Newton
    !! iterations minimize
    !! \(f(u)=\tfrac12\|\mathbf{C}(u)-\mathbf{x}\|^2\) on the active interval
    !! \(\Omega=[u_{\min},u_{\max}]\). Convergence is tested with the gradient
    !! mapping
    !!
    !! \[
    !! g_P(u)=u-\Pi_\Omega\!\left(u-f'(u)\right),
    !! \]
    !!
    !! where \(\Pi_\Omega\) is Euclidean projection onto the interval.
    !! \(g_P=0\) is equivalent to first-order KKT stationarity: \(f'=0\) in
    !! the interior, \(f'\ge0\) at the lower bound, and \(f'\le0\) at the upper
    !! bound. The Newton proposal is projected onto \(\Omega\); a non-descent
    !! proposal is replaced by the projected negative-gradient displacement.
    !! A trial changes the iterate only after satisfying the Armijo
    !! sufficient-decrease inequality. Failure to accept one of 50 reductions
    !! retains the current iterate and sets code 109. Exhausting `maxit` also
    !! sets code 109; a rejected scalar Hessian sets code 108. On successful
    !! return, the stationary point is local and need not be the global nearest
    !! point.
    !!
    !! Line-search reference: L. Armijo, "Minimization of Functions Having
    !! Lipschitz Continuous First Partial Derivatives," *Pacific Journal of
    !! Mathematics* 16 (1966), 1--3.
    !! [doi:10.2140/pjm.1966.16.1](https://doi.org/10.2140/pjm.1966.16.1).
    !! Bound projection follows D. P. Bertsekas, "Projected Newton Methods for
    !! Optimization Problems with Simple Constraints," *SIAM Journal on
    !! Control and Optimization* 20 (1982), 221--246.
    !! [doi:10.1137/0320018](https://doi.org/10.1137/0320018).
    impure subroutine nearest_point2(this, point_Xg, tol, maxit, nearest_Xt, nearest_Xg)

        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: point_Xg(:)
            !! Finite query with at least the curve's physical dimension; extra components are ignored.
        real(rk), intent(in) :: tol
            !! Nonnegative absolute gradient-mapping tolerance; not internally validated.
        integer, intent(in) :: maxit
            !! Positive maximum iteration count; not internally validated.
        real(rk), intent(out) :: nearest_Xt
            !! Parameter of the returned candidate.
        real(rk), intent(out), optional :: nearest_Xg(size(this%Xc,2))
            !! Optional physical point at `nearest_Xt`.

        real(rk) :: xk, obj, obj_trial, grad, projected_grad, hess, dk, alphak
        real(rk) :: tau, beta, lower_bounds, upper_bounds, xt, grad_step
        real(rk) :: Xg(size(this%Xc,2))
        real(rk) :: residual(size(this%Xc,2)), dXg(size(this%Xc,2)), d2Xg(size(this%Xc,2))
        real(rk), allocatable :: dTgc(:), d2Tgc(:)
        integer :: k, l, iseed, nseed, i, d, dim_
        logical :: convergenz, line_search_accepted
        real(rk) :: dist, best_dist

        nearest_Xt = 0.0_rk
        if (present(nearest_Xg)) nearest_Xg = 0.0_rk
        if (.not. this%err%ok) return

        dk  = 0.0_rk
        k   = 0
        dim_ = size(this%Xc,2)

        ! bounds
        lower_bounds = knot_start(this%knot, this%nc, this%degree)
        upper_bounds = knot_end(this%knot, this%nc, this%degree)

        ! initial guess (streamed coarse search)
        nseed = 10
        xk = lower_bounds
        best_dist = huge(1.0_rk)
        do iseed = 1, nseed
            if (nseed > 1) then
                xt = lower_bounds + (upper_bounds - lower_bounds)*real(iseed-1, rk)/real(nseed-1, rk)
            else
                xt = lower_bounds
            end if
            Xg = this%cmp_Xg(xt)
            dist = 0.0_rk
            do d = 1, dim_
                dist = dist + (Xg(d) - point_Xg(d))**2
            end do
            if (dist < best_dist) then
                best_dist = dist
                xk = xt
            end if
        end do

        ! clamp initial guess to bounds
        xk = max(min(xk, upper_bounds), lower_bounds)
        nearest_Xt = xk
        if (present(nearest_Xg)) nearest_Xg = this%cmp_Xg(xk)

        convergenz = .false.

        do while (.not. convergenz .and. k < maxit)

            ! objective, gradient, hessian
            Xg = this%cmp_Xg(xk)
            call this%derivative2(Xt=xk, d2Tgc=d2Tgc, dTgc=dTgc)

            obj = 0.0_rk
            dXg = 0.0_rk
            d2Xg = 0.0_rk
            do d = 1, dim_
                residual(d) = Xg(d) - point_Xg(d)
                obj = obj + residual(d)*residual(d)
            end do
            obj = 0.5_rk*obj
            do i = 1, this%nc
                do d = 1, dim_
                    dXg(d) = dXg(d) + dTgc(i)*this%Xc(i,d)
                    d2Xg(d) = d2Xg(d) + d2Tgc(i)*this%Xc(i,d)
                end do
            end do
            grad = 0.0_rk
            do d = 1, dim_
                grad = grad + residual(d)*dXg(d)
            end do
            hess = 0.0_rk
            do d = 1, dim_
                hess = hess + dXg(d)*dXg(d) + residual(d)*d2Xg(d)
            end do
            projected_grad = xk - max(min(xk-grad, upper_bounds), lower_bounds)

            if (abs(projected_grad) <= tol) then
                convergenz = .true.
                nearest_Xt = xk
                if (present(nearest_Xg)) nearest_Xg = Xg
            else
                ! Newton step
                if (abs(hess) <= 32.0_rk*epsilon(1.0_rk)*max(1.0_rk, abs(hess))) then
                    call this%err%set(&
                        code       = 108,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Singular Hessian in nearest-point iteration.",&
                        location   = "nearest_point2",&
                        suggestion = "Check the curve geometry, tolerance, or initial sampling interval.")
                    nearest_Xt = xk
                    if (present(nearest_Xg)) nearest_Xg = Xg
                    return
                end if
                dk = - grad / hess
                dk = max(min(xk+dk, upper_bounds), lower_bounds) - xk
                grad_step = grad*dk
                if (grad_step >= 0.0_rk) then
                    dk = max(min(xk-grad, upper_bounds), lower_bounds) - xk
                    grad_step = grad*dk
                end if
                if (.not. ieee_is_finite(grad_step) .or. grad_step >= 0.0_rk) then
                    nearest_Xt = xk
                    if (present(nearest_Xg)) nearest_Xg = Xg
                    call this%err%set(&
                        code       = 109,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "No feasible descent direction in nearest-point iteration.",&
                        location   = "nearest_point2",&
                        suggestion = "Relax tolerance or check geometry smoothness and conditioning.")
                    return
                end if

                tau  = 0.5_rk
                beta = 1.0e-4_rk
                alphak = 1.0_rk
                line_search_accepted = .false.
                do l = 1, 50
                    xt = xk + alphak*dk
                    Xg = this%cmp_Xg(xt)
                    obj_trial = 0.0_rk
                    do d = 1, dim_
                        obj_trial = obj_trial + (Xg(d) - point_Xg(d))**2
                    end do
                    obj_trial = 0.5_rk*obj_trial
                    if (ieee_is_finite(obj_trial) .and. &
                        obj_trial <= obj + alphak*beta*grad_step) then
                        line_search_accepted = .true.
                        exit
                    end if
                    alphak = tau*alphak
                end do

                if (.not. line_search_accepted) then
                    nearest_Xt = xk
                    if (present(nearest_Xg)) nearest_Xg = this%cmp_Xg(xk)
                    call this%err%set(&
                        code       = 109,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Armijo line search failed in nearest-point iteration.",&
                        location   = "nearest_point2",&
                        suggestion = "Relax tolerance or check geometry smoothness and conditioning.")
                    return
                end if

                xk = max(min(xt, upper_bounds), lower_bounds)
                k = k + 1
            end if
        end do

        if (.not. convergenz) then
            nearest_Xt = xk
            if (present(nearest_Xg)) nearest_Xg = this%cmp_Xg(xk)
            call this%err%set(&
                code       = 109,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Nearest-point iteration did not converge.",&
                location   = "nearest_point2",&
                suggestion = "Increase maxit, relax tolerance, or check geometry conditioning.")
        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute element-local IGA basis data at one quadrature point.
    !!
    !! `ie` addresses nonzero active spans and `ig` addresses the selected
    !! Gauss-Legendre rule. The local basis order agrees with `cmp_elem`.
    !! `dTgc_dXg(a,:)` is the tangential physical gradient embedded in Cartesian
    !! space. `dL` already includes quadrature weight and both parameter-to-span
    !! and curve metric factors. Supplying `d2Tgc_dXg2` computes the second
    !! derivative with respect to arc length. By default, a zero tangent returns
    !! zero physical derivatives and zero `dL`, which preserves nonnegative
    !! geometric-measure behavior. With `strict=.true.`, a numerically singular
    !! tangent sets diagnostic code 108 instead. For a one-coordinate physical
    !! map, strict mode also requires a positive derivative; an embedded curve
    !! has no canonical signed orientation without an external reference.
    !! The check is local to the requested quadrature point and therefore does
    !! not prove global injectivity between quadrature points.
    !!
    !! With \(\mathbf t=\partial\mathbf C/\partial u\), the first physical
    !! derivative returned for basis \(R_a\) is
    !!
    !! \[
    !! \nabla_\Gamma R_a
    !!   =\frac{\partial R_a/\partial u}{\mathbf t\mathbin{\cdot}\mathbf t}
    !!    \,\mathbf t ,
    !! \]
    !!
    !! while the optional scalar second derivative is the arc-length derivative
    !!
    !! \[
    !! \frac{\mathrm d^2R_a}{\mathrm ds^2}
    !!   =\frac{R_{a,uu}(\mathbf t\mathbin{\cdot}\mathbf t)
    !!           -R_{a,u}(\mathbf t\mathbin{\cdot}\mathbf C_{,uu})}
    !!          {(\mathbf t\mathbin{\cdot}\mathbf t)^2}.
    !! \]
    !!
    !! and the weighted quadrature contribution is
    !!
    !! \[
    !! dL=\omega_g\,(u_{e+1}-u_e)\,\lVert\mathbf t\rVert .
    !! \]
    !!
    !! IGA formulation reference: T. J. R. Hughes, J. A. Cottrell, and
    !! Y. Bazilevs, *CMAME* 194 (2005), 4135--4195.
    !! [doi:10.1016/j.cma.2004.10.008](https://doi.org/10.1016/j.cma.2004.10.008).
    !!
    !! Elementwise analysis assumes an invertible geometry map with a smooth
    !! inverse: Y. Bazilevs et al., *MMMAS* 16 (2006), 1031--1090.
    !! [doi:10.1142/S0218202506001455](https://doi.org/10.1142/S0218202506001455).
    pure subroutine ansatz(this, ie, ig, Tgc, dTgc_dXg, dL, ngauss, d2Tgc_dXg2, strict)
        class(nurbs_curve), intent(inout) :: this

        integer, intent(in) :: ie
            !! One-based active element index.
        integer, intent(in) :: ig
            !! One-based quadrature-point index within the element rule.
        real(rk), intent(out) :: dL
            !! Weighted differential length contribution.
        real(rk), allocatable, intent(out) :: Tgc(:)
            !! Local basis values, size `degree+1`.
        real(rk), allocatable, intent(out) :: dTgc_dXg(:,:)
            !! Local physical gradients `[degree+1,ncoord]`.
        integer, intent(in), optional :: ngauss
            !! Optional positive quadrature-point count; default `degree+1`.
        real(rk), allocatable, intent(out), optional :: d2Tgc_dXg2(:)
            !! Optional local second arc-length derivatives.
        logical, intent(in), optional :: strict
            !! Reject singular tangents and reversed one-coordinate mappings.
        real(rk), allocatable :: Xth(:), Xksi(:), Wksi(:)
        integer, allocatable :: elem_c(:,:)
        real(rk), allocatable :: dTgc_dXt(:), gradient(:,:), hessian(:,:,:)
        real(rk) :: Xt, dXt_dXksi, norm_dXg, metric, tangent_curvature
        real(rk) :: control_scale, tol
        real(rk) :: dXg_dXt(size(this%Xc, 2)), d2Xg_dXt2(size(this%Xc, 2))
        real(rk) :: ders(0:2,0:this%degree), identity(0:2,0:0)
        integer :: a, d, dim_, nloc, degree_, first
        logical :: strict_

        if (.not. this%err%ok) return

        strict_ = .false.
        if (present(strict)) strict_ = strict

        if (present(ngauss)) then
            degree_ = ngauss - 1
        else
            degree_ = this%degree
        end if
        if (degree_ < 0) then
            allocate(Tgc(0), dTgc_dXg(0,0))
            if (present(d2Tgc_dXg2)) allocate(d2Tgc_dXg2(0))
            dL = 0.0_rk
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "The quadrature count must be positive.",&
                location   = "ansatz",&
                suggestion = "Set ngauss to at least one.")
            return
        end if
        call gauss_leg([0.0_rk, 1.0_rk], degree_, Xksi, Wksi)

        Xth = active_knots(this%knot, this%nc, this%degree)
        elem_c = this%cmp_elem()
        if (ie < 1 .or. ie >= size(Xth) .or. ig < 1 .or. ig > size(Wksi)) then
            allocate(Tgc(0), dTgc_dXg(0,0))
            if (present(d2Tgc_dXg2)) allocate(d2Tgc_dXg2(0))
            dL = 0.0_rk
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Invalid ansatz element or quadrature index.",&
                location   = "ansatz",&
                suggestion = "Use an active element and a valid positive quadrature count.")
            return
        end if

        dXt_dXksi = Xth(ie+1) - Xth(ie)
        Xt = Xth(ie) + Xksi(ig)*dXt_dXksi

        if (present(d2Tgc_dXg2)) then
            nloc = size(elem_c,2)
            allocate(Tgc(nloc), gradient(nloc,1), hessian(nloc,1,1))
            identity = 0.0_rk
            identity(0,0) = 1.0_rk
            call basis_bspline_der_all_active(&
                Xt     = Xt,&
                knot   = this%knot,&
                nc     = this%nc,&
                degree = this%degree,&
                nder   = 2,&
                first  = first,&
                ders   = ders)
            if (this%is_rational()) then
                call tensor_basis_derivatives2_local(&
                    first1   = first,&
                    nc1      = this%nc,&
                    ders1    = ders,&
                    first2   = 1,&
                    nc2      = 1,&
                    ders2    = identity,&
                    first3   = 1,&
                    nc3      = 1,&
                    ders3    = identity,&
                    basis    = Tgc,&
                    gradient = gradient,&
                    hessian  = hessian,&
                    Wc       = this%Wc)
            else
                call tensor_basis_derivatives2_local(&
                    first1   = first,&
                    nc1      = this%nc,&
                    ders1    = ders,&
                    first2   = 1,&
                    nc2      = 1,&
                    ders2    = identity,&
                    first3   = 1,&
                    nc3      = 1,&
                    ders3    = identity,&
                    basis    = Tgc,&
                    gradient = gradient,&
                    hessian  = hessian)
            end if
        else if (this%is_rational()) then
            call compute_dTgc(Xt, this%knot, this%degree, this%nc, this%Wc, dTgc_dXt, Tgc, elem_c(ie,:), this%wrap_parameters)
        else
            call compute_dTgc(Xt, this%knot, this%degree, this%nc, dTgc_dXt, Tgc, elem_c(ie,:), this%wrap_parameters)
        end if

        nloc = size(Tgc)
        dim_ = size(this%Xc, 2)
        dXg_dXt = 0.0_rk
        d2Xg_dXt2 = 0.0_rk
        if (present(d2Tgc_dXg2)) then
            do a = 1, nloc
                do d = 1, dim_
                    dXg_dXt(d) = dXg_dXt(d) + this%Xc(elem_c(ie,a),d)*gradient(a,1)
                    d2Xg_dXt2(d) = d2Xg_dXt2(d) + this%Xc(elem_c(ie,a),d)*hessian(a,1,1)
                end do
            end do
        else
            do a = 1, nloc
                do d = 1, dim_
                    dXg_dXt(d) = dXg_dXt(d) + this%Xc(elem_c(ie,a),d)*dTgc_dXt(a)
                end do
            end do
        end if

        norm_dXg = norm2(dXg_dXt)
        allocate(dTgc_dXg(nloc, dim_), source=0.0_rk)
        if (present(d2Tgc_dXg2)) allocate(d2Tgc_dXg2(nloc), source=0.0_rk)
        if (strict_) then
            control_scale = 0.0_rk
            do a = 2, nloc
                do d = 1, dim_
                    control_scale = max(&
                        control_scale,&
                        abs(this%Xc(elem_c(ie,a),d) - this%Xc(elem_c(ie,1),d)))
                end do
            end do
            tol = 64.0_rk*epsilon(1.0_rk)*control_scale/abs(dXt_dXksi)
            if (norm_dXg <= tol) then
                dL = 0.0_rk
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Singular curve tangent in strict analysis mode.",&
                    location   = "ansatz",&
                    suggestion = "Use a regular geometry map at every analysis quadrature point.")
                return
            end if
            if (dim_ == 1 .and. dXg_dXt(1) < 0.0_rk) then
                dL = 0.0_rk
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Inverted curve Jacobian in strict analysis mode.",&
                    location   = "ansatz",&
                    suggestion = "Orient the one-dimensional geometry map consistently.")
                return
            end if
        end if
        if (norm_dXg > 0.0_rk) then
            if (present(d2Tgc_dXg2)) then
                do d = 1, dim_
                    do a = 1, nloc
                        dTgc_dXg(a,d) = gradient(a,1)*dXg_dXt(d)/(norm_dXg*norm_dXg)
                    end do
                end do
                metric = norm_dXg*norm_dXg
                tangent_curvature = dot_product(dXg_dXt, d2Xg_dXt2)
                do a = 1, nloc
                    d2Tgc_dXg2(a) = (hessian(a,1,1)*metric-gradient(a,1)*tangent_curvature)/(metric*metric)
                end do
            else
                do d = 1, dim_
                    do a = 1, nloc
                        dTgc_dXg(a,d) = dTgc_dXt(a)*dXg_dXt(d)/(norm_dXg*norm_dXg)
                    end do
                end do
            end if
        end if

        dL = norm_dXg*abs(dXt_dXksi)*Wksi(ig)
    end subroutine
    !===============================================================================

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Integrate arc length over all nonzero active knot spans.
    !!
    !! The computed approximation is
    !!
    !! \[L_h=\sum_e\sum_{g=1}^{n_q}
    !!   \omega_g\,\|\mathbf C'(u_{eg})\|\,|u_{e+1}-u_e|.\]
    !!
    !! Each nonzero active span uses an independent Gauss-Legendre rule, so
    !! repeated knots do not bridge elements. Degenerate tangents contribute
    !! zero measure at the affected quadrature points. Although an \(n_q\)-point
    !! Gauss rule is exact for polynomial integrands through degree \(2n_q-1\),
    !! the arc-length integrand is generally nonpolynomial, including for a
    !! polynomial spline curve; this method therefore provides numerical
    !! quadrature rather than an exact-length guarantee.
    !! The rational tangent is evaluated after exact radix-power weight
    !! normalization as \((S'-\mathbf C W')/W\); the implementation never forms
    !! \(W^2\), whose underflow or overflow would violate projective-scale
    !! invariance.
    !!
    !! A valid initialized curve and `ngauss>=1` are caller preconditions; the
    !! latter is not diagnosed here. If an existing diagnostic is active, the
    !! routine returns `length=0` without integration.
    pure subroutine cmp_length(this, length, ngauss)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(out) :: length
            !! Nonnegative integrated curve length.
        integer, intent(in), optional :: ngauss
            !! Optional positive points per span; default `degree+1`.
        real(rk), allocatable :: Xth(:), Xksi(:), Wksi(:)
        real(rk) :: N(0:this%degree), dN(0:this%degree)
        real(rk) :: Xt, dXt_dXksi, denom, ddenom, wi
        real(rk) :: dL, norm2_dXg, Sw_d, dSw_d, dXg_dXt_d
        integer :: ie, ig, j, d, idx, first, dim_, j0, j1, ngauss_, nelem, weight_exponent
        logical :: rational

        length = 0.0_rk
        if (.not. this%err%ok) return

        if (present(ngauss)) then
            ngauss_ = ngauss
        else
            ngauss_ = this%degree + 1
        end if

        call gauss_leg([0.0_rk, 1.0_rk], ngauss_-1, Xksi, Wksi)
        Xth = active_knots(this%knot, this%nc, this%degree)
        nelem = size(Xth) - 1
        rational = this%is_rational()
        dim_ = size(this%Xc, 2)
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do ie = 1, nelem
#else
        do concurrent (ie = 1: nelem) &
            local(N, dN, Xt, dXt_dXksi, denom, ddenom, wi, dL, norm2_dXg, Sw_d, dSw_d, &
                dXg_dXt_d, ig, j, d, idx, first, j0, j1, weight_exponent) &
            reduce(+:length)
#endif
            dL = 0.0_rk
            do ig = 1, ngauss_
                dXt_dXksi = Xth(ie+1) - Xth(ie)
                Xt = Xth(ie) + Xksi(ig)*dXt_dXksi
                norm2_dXg = 0.0_rk

                call basis_bspline_der(Xt, this%knot, this%nc, this%degree, 1, first, N, dN)

                j0 = max(0, 1 - first)
                j1 = min(this%degree, this%nc - first)
                if (rational) then
                    weight_exponent = exponent(this%Wc(first+j0))
                    do j = j0 + 1, j1
                        weight_exponent = max(weight_exponent, exponent(this%Wc(first+j)))
                    end do
                    denom = 0.0_rk
                    ddenom = 0.0_rk
                    do j = j0, j1
                        idx = first + j
                        wi = scale(this%Wc(idx), -weight_exponent)
                        denom = denom + N(j)*wi
                        ddenom = ddenom + dN(j)*wi
                    end do

                    if (denom > 0.0_rk) then
                        do d = 1, dim_
                            Sw_d = 0.0_rk
                            dSw_d = 0.0_rk
                            do j = j0, j1
                                idx = first + j
                                wi = scale(this%Wc(idx), -weight_exponent)
                                Sw_d = Sw_d + N(j)*wi*this%Xc(idx,d)
                                dSw_d = dSw_d + dN(j)*wi*this%Xc(idx,d)
                            end do
                            Sw_d = Sw_d/denom
                            dXg_dXt_d = (dSw_d - Sw_d*ddenom)/denom
                            norm2_dXg = norm2_dXg + dXg_dXt_d*dXg_dXt_d
                        end do
                    end if
                else
                    do d = 1, dim_
                        dXg_dXt_d = 0.0_rk
                        do j = j0, j1
                            idx = first + j
                            dXg_dXt_d = dXg_dXt_d + dN(j)*this%Xc(idx,d)
                        end do
                        norm2_dXg = norm2_dXg + dXg_dXt_d*dXg_dXt_d
                    end do
                end if

                dL = dL + sqrt(norm2_dXg)*abs(dXt_dXksi)*Wksi(ig)
            end do
            length = length + dL
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate a rational curve at a vector of parameters.
    !!
    !! The result has shape `[ng_,size(Xc,2)]`, where `ng_` is `ng` when
    !! present and otherwise `size(Xt)`. Each row uses only the `degree+1`
    !! locally active controls and is divided by
    !! \(W(u)=\sum_iN_{i,p}(u)w_i\). A zero denominator leaves that row zero.
    pure function compute_Xg_nurbs_1d(Xt, knot, degree, nc, ng, Xc, Wc, wrap_parameters) result(Xg)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in), optional :: ng
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous :: Wc(:)
        logical, intent(in), optional :: wrap_parameters
        real(rk), allocatable :: Xg(:,:)
        real(rk) :: N(0:degree), bw, wi, wsum, wtol
        integer :: i, j, d, idx, first, ng_, dim_, j0, j1, weight_exponent
        logical :: wrap_parameters_

        if (present(ng)) then
            ng_ = ng
        else
            ng_ = size(Xt, 1)
        end if

        dim_ = size(Xc,2)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
        allocate(Xg(ng_, dim_))
        wtol = 0.0_rk
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N, wsum, bw, wi, j, d, idx, first, j0, j1, weight_exponent)
#endif
            do d = 1, dim_
                Xg(i,d) = 0.0_rk
            end do
            call basis_bspline_der(map_parameter(Xt(i), knot, nc, degree, wrap_parameters_), &
                knot, nc, degree, 0, first, N)
            j0 = max(0, 1 - first)
            j1 = min(degree, nc - first)
            weight_exponent = exponent(Wc(first+j0))
            do j = j0 + 1, j1
                weight_exponent = max(weight_exponent, exponent(Wc(first+j)))
            end do
            wsum = 0.0_rk
            do j = j0, j1
                idx = first + j
                wi = scale(Wc(idx), -weight_exponent)
                bw = N(j)*wi
                wsum = wsum + bw
                do d = 1, dim_
                    Xg(i,d) = Xg(i,d) + bw*Xc(idx,d)
                end do
            end do
            if (abs(wsum) > wtol) then
                do d = 1, dim_
                    Xg(i,d) = Xg(i,d)/wsum
                end do
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate a rational curve at one parameter.
    !!
    !! The fixed-size result contains one value per physical coordinate. The
    !! parameter is used directly; parameter wrapping is handled by the object
    !! API before this scalar implementation is selected.
    pure function compute_Xg_nurbs_1d_1point(Xt, knot, degree, nc, Xc, Wc) result(Xg)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk) :: Xg(size(Xc,2))
        real(rk) :: N(0:degree), bw, wi, wsum, wtol
        integer :: j, d, idx, first, dim_, j0, j1, weight_exponent

        Xg = 0.0_rk
        dim_ = size(Xc,2)
        wtol = 0.0_rk
        call basis_bspline_der(Xt, knot, nc, degree, 0, first, N)
        j0 = max(0, 1 - first)
        j1 = min(degree, nc - first)
        weight_exponent = exponent(Wc(first+j0))
        do j = j0 + 1, j1
            weight_exponent = max(weight_exponent, exponent(Wc(first+j)))
        end do
        wsum = 0.0_rk
        do j = j0, j1
            idx = first + j
            wi = scale(Wc(idx), -weight_exponent)
            bw = N(j)*wi
            wsum = wsum + bw
            do d = 1, dim_
                Xg(d) = Xg(d) + bw*Xc(idx,d)
            end do
        end do
        if (abs(wsum) > wtol) then
            do d = 1, dim_
                Xg(d) = Xg(d)/wsum
            end do
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate a polynomial B-spline curve at a vector of parameters.
    !!
    !! The result has shape `[ng_,size(Xc,2)]`; `ng_` defaults to `size(Xt)`.
    !! Local-support accumulation costs \(O(ng_\,p\,n_{coord})\) and avoids a
    !! dense basis matrix.
    pure function compute_Xg_bspline_1d(Xt, knot, degree, nc, ng, Xc, wrap_parameters) result(Xg)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in), optional :: ng
        real(rk), intent(in), contiguous :: Xc(:,:)
        logical, intent(in), optional :: wrap_parameters
        real(rk), allocatable :: Xg(:,:)
        real(rk) :: N(0:degree)
        integer :: i, j, d, idx, first, ng_, dim_
        logical :: wrap_parameters_

        if (present(ng)) then
            ng_ = ng
        else
            ng_ = size(Xt, 1)
        end if

        dim_ = size(Xc,2)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
        allocate(Xg(ng_, dim_))
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N, j, d, idx, first)
#endif
            do d = 1, dim_
                Xg(i,d) = 0.0_rk
            end do
            call basis_bspline_der(map_parameter(Xt(i), knot, nc, degree, wrap_parameters_), &
                knot, nc, degree, 0, first, N)
            do j = max(0, 1 - first), min(degree, nc - first)
                idx = first + j
                do d = 1, dim_
                    Xg(i,d) = Xg(i,d) + N(j)*Xc(idx,d)
                end do
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate a polynomial B-spline curve at one parameter.
    !!
    !! The result has length `size(Xc,2)` and is accumulated from only the
    !! locally active basis functions. The parameter is not wrapped here.
    pure function compute_Xg_bspline_1d_1point(Xt, knot, degree, nc, Xc) result(Xg)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk) :: Xg(size(Xc,2))
        real(rk) :: N(0:degree)
        integer :: j, d, idx, first, dim_

        Xg = 0.0_rk
        dim_ = size(Xc,2)
        call basis_bspline_der(Xt, knot, nc, degree, 0, first, N)
        do j = max(0, 1 - first), min(degree, nc - first)
            idx = first + j
            do d = 1, dim_
                Xg(d) = Xg(d) + N(j)*Xc(idx,d)
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate rational curve basis values and derivatives through order two.
    !!
    !! Outputs have shape `[ng_,nc]`, with parameters in rows and global basis
    !! indices in columns. Quotient-rule normalization is applied in one pass;
    !! rows remain zero when the weight denominator vanishes.
    pure subroutine compute_d2Tgc_nurbs_1d_vector(Xt, knot, degree, nc, ng, Wc, d2Tgc, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in), optional :: ng
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        logical, intent(in), optional :: wrap_parameters
        real(rk) :: N(0:degree), dN(0:degree), d2N(0:degree)
        real(rk) :: denom, ddenom, d2denom, wi, wtol
        integer :: i, j, idx, first, ng_, j0, j1, weight_exponent
        logical :: wrap_parameters_

        if (present(ng)) then
            ng_ = ng
        else
            ng_ = size(Xt, 1)
        end if

        allocate(d2Tgc(ng_, nc), dTgc(ng_, nc), Tgc(ng_, nc), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
        wtol = 0.0_rk

#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) &
            local(N, dN, d2N, denom, ddenom, d2denom, wi, j, idx, first, j0, j1, weight_exponent)
#endif
            call basis_bspline_der(map_parameter(Xt(i), knot, nc, degree, wrap_parameters_), &
                knot, nc, degree, 2, first, N, dN, d2N)
            j0 = max(0, 1 - first)
            j1 = min(degree, nc - first)
            weight_exponent = exponent(Wc(first+j0))
            do j = j0 + 1, j1
                weight_exponent = max(weight_exponent, exponent(Wc(first+j)))
            end do
            denom = 0.0_rk
            ddenom = 0.0_rk
            d2denom = 0.0_rk
            do j = j0, j1
                idx = first + j
                wi = scale(Wc(idx), -weight_exponent)
                denom = denom + N(j)*wi
                ddenom = ddenom + dN(j)*wi
                d2denom = d2denom + d2N(j)*wi
            end do
            if (abs(denom) > wtol) then
                do j = j0, j1
                    idx = first + j
                    wi = scale(Wc(idx), -weight_exponent)
                    Tgc(i,idx) = N(j)*wi/denom
                    dTgc(i,idx) = (dN(j)*wi - Tgc(i,idx)*ddenom)/denom
                    d2Tgc(i,idx) = (d2N(j)*wi - 2.0_rk*dTgc(i,idx)*ddenom - Tgc(i,idx)*d2denom)/denom
                end do
            end if
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate dense rational curve basis values and first two derivatives.
    !!
    !! Each allocated output has length `nc`. Entries outside the local support
    !! are zero, and optional wrapping maps `Xt` into the active interval before
    !! basis evaluation.
    pure subroutine compute_d2Tgc_nurbs_1d_scalar(Xt, knot, degree, nc, Wc, d2Tgc, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:)
        real(rk), allocatable, intent(out) :: dTgc(:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        logical, intent(in), optional :: wrap_parameters
        real(rk) :: N(0:degree), dN(0:degree), d2N(0:degree)
        real(rk) :: denom, ddenom, d2denom, wi, wtol
        integer :: j, idx, first, j0, j1, weight_exponent
        logical :: wrap_parameters_

        allocate(d2Tgc(nc), dTgc(nc), Tgc(nc), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
        wtol = 0.0_rk
        call basis_bspline_der(map_parameter(Xt, knot, nc, degree, wrap_parameters_), &
            knot, nc, degree, 2, first, N, dN, d2N)
        j0 = max(0, 1 - first)
        j1 = min(degree, nc - first)
        weight_exponent = exponent(Wc(first+j0))
        do j = j0 + 1, j1
            weight_exponent = max(weight_exponent, exponent(Wc(first+j)))
        end do
        denom = 0.0_rk
        ddenom = 0.0_rk
        d2denom = 0.0_rk
        do j = j0, j1
            idx = first + j
            wi = scale(Wc(idx), -weight_exponent)
            denom = denom + N(j)*wi
            ddenom = ddenom + dN(j)*wi
            d2denom = d2denom + d2N(j)*wi
        end do
        if (abs(denom) > wtol) then
            do j = j0, j1
                idx = first + j
                wi = scale(Wc(idx), -weight_exponent)
                Tgc(idx) = N(j)*wi/denom
                dTgc(idx) = (dN(j)*wi - Tgc(idx)*ddenom)/denom
                d2Tgc(idx) = (d2N(j)*wi - 2.0_rk*dTgc(idx)*ddenom - Tgc(idx)*d2denom)/denom
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate polynomial curve basis values and derivatives through order two.
    !!
    !! Outputs have shape `[ng_,nc]`. Only local-support columns are populated
    !! for each parameter row; derivatives above the spline degree are zero.
    pure subroutine compute_d2Tgc_bspline_1d_vector(Xt, knot, degree, nc, ng, d2Tgc, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in), optional :: ng
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        logical, intent(in), optional :: wrap_parameters
        integer :: i, ng_
        real(rk) :: N(0:degree), dN(0:degree), d2N(0:degree)
        integer :: j, idx, first
        logical :: wrap_parameters_

        if (present(ng)) then
            ng_ = ng
        else
            ng_ = size(Xt, 1)
        end if

        allocate(d2Tgc(ng_, nc), dTgc(ng_, nc), Tgc(ng_, nc), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N, dN, d2N, j, idx, first)
#endif
            call basis_bspline_der(map_parameter(Xt(i), knot, nc, degree, wrap_parameters_), &
                knot, nc, degree, 2, first, N, dN, d2N)
            do j = max(0, 1 - first), min(degree, nc - first)
                idx = first + j
                Tgc(i,idx) = N(j)
                dTgc(i,idx) = dN(j)
                d2Tgc(i,idx) = d2N(j)
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate dense polynomial curve basis values and first two derivatives.
    !!
    !! Each allocated output has length `nc`; entries outside the active support
    !! are zero. Optional wrapping is applied before span selection.
    pure subroutine compute_d2Tgc_bspline_1d_scalar(Xt, knot, degree, nc, d2Tgc, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        real(rk), allocatable, intent(out) :: d2Tgc(:)
        real(rk), allocatable, intent(out) :: dTgc(:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        logical, intent(in), optional :: wrap_parameters
        real(rk) :: N(0:degree), dN(0:degree), d2N(0:degree)
        integer :: j, idx, first
        logical :: wrap_parameters_

        allocate(d2Tgc(nc), dTgc(nc), Tgc(nc), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
        call basis_bspline_der(map_parameter(Xt, knot, nc, degree, wrap_parameters_), &
            knot, nc, degree, 2, first, N, dN, d2N)
        do j = max(0, 1 - first), min(degree, nc - first)
            idx = first + j
            Tgc(idx) = N(j)
            dTgc(idx) = dN(j)
            d2Tgc(idx) = d2N(j)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate rational curve basis values and first derivatives at many parameters.
    !!
    !! `Tgc` and `dTgc` have shape `[ng_,nc]`. The quotient rule uses the same
    !! local polynomial basis pass for values and derivatives, while optional
    !! parameter wrapping is applied independently to every row.
    pure subroutine compute_dTgc_nurbs_1d_vector(Xt, knot, degree, nc, ng, Wc, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in), optional :: ng
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        logical, intent(in), optional :: wrap_parameters
        real(rk) :: N(0:degree), dN(0:degree)
        real(rk) :: denom, ddenom, wi, wtol
        integer :: i, j, idx, first, ng_, j0, j1, weight_exponent
        logical :: wrap_parameters_

        if (present(ng)) then
            ng_ = ng
        else
            ng_ = size(Xt, 1)
        end if

        allocate(dTgc(ng_, nc), Tgc(ng_, nc), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
        wtol = 0.0_rk
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) &
            local(N, dN, denom, ddenom, wi, j, idx, first, j0, j1, weight_exponent)
#endif
            call basis_bspline_der(map_parameter(Xt(i), knot, nc, degree, wrap_parameters_), &
                knot, nc, degree, 1, first, N, dN)
            j0 = max(0, 1 - first)
            j1 = min(degree, nc - first)
            weight_exponent = exponent(Wc(first+j0))
            do j = j0 + 1, j1
                weight_exponent = max(weight_exponent, exponent(Wc(first+j)))
            end do
            denom = 0.0_rk
            ddenom = 0.0_rk
            do j = j0, j1
                idx = first + j
                wi = scale(Wc(idx), -weight_exponent)
                denom = denom + N(j)*wi
                ddenom = ddenom + dN(j)*wi
            end do
            if (abs(denom) > wtol) then
                do j = j0, j1
                    idx = first + j
                    wi = scale(Wc(idx), -weight_exponent)
                    Tgc(i,idx) = N(j)*wi/denom
                    dTgc(i,idx) = (dN(j)*wi - Tgc(i,idx)*ddenom)/denom
                end do
            end if
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate rational curve basis values and first derivatives at one parameter.
    !!
    !! Without `elem`, dense outputs of length `nc` are returned. With `elem`,
    !! output order follows that one-based control-index list and `Wc` may be
    !! either the full global weight vector or an element-local vector of the
    !! same length as `elem`.
    pure subroutine compute_dTgc_nurbs_1d_scalar(Xt, knot, degree, nc, Wc, dTgc, Tgc, elem, wrap_parameters)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        real(rk), intent(in), contiguous :: Wc(:)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), optional :: wrap_parameters
        real(rk), allocatable, intent(out) :: dTgc(:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        real(rk) :: N(0:degree), dN(0:degree)
        real(rk) :: denom, ddenom
        integer :: i, idx, nloc, first, ia, i0, i1, weight_exponent
        real(rk) :: wi, wtol
        logical :: local_weights
        logical :: wrap_parameters_

        wtol = 0.0_rk
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters

        if (.not. present(elem)) then
            allocate(dTgc(nc), Tgc(nc), source=0.0_rk)
            call basis_bspline_der(map_parameter(Xt, knot, nc, degree, wrap_parameters_), &
                knot, nc, degree, 1, first, N, dN)
            i0 = max(0, 1 - first)
            i1 = min(degree, nc - first)
            weight_exponent = exponent(Wc(first+i0))
            do i = i0 + 1, i1
                weight_exponent = max(weight_exponent, exponent(Wc(first+i)))
            end do
            denom = 0.0_rk
            ddenom = 0.0_rk
            do i = i0, i1
                idx = first + i
                wi = scale(Wc(idx), -weight_exponent)
                denom = denom + N(i)*wi
                ddenom = ddenom + dN(i)*wi
            end do
            if (abs(denom) > wtol) then
                do i = i0, i1
                    idx = first + i
                    wi = scale(Wc(idx), -weight_exponent)
                    Tgc(idx) = N(i)*wi/denom
                    dTgc(idx) = (dN(i)*wi - Tgc(idx)*ddenom)/denom
                end do
            end if
        else
            nloc = size(elem)
            allocate(dTgc(nloc), Tgc(nloc), source=0.0_rk)
            call basis_bspline_der(map_parameter(Xt, knot, nc, degree, wrap_parameters_), &
                knot, nc, degree, 1, first, N, dN)
            denom = 0.0_rk
            ddenom = 0.0_rk
            local_weights = size(Wc) == nloc .and. size(Wc) /= nc
            if (local_weights) then
                weight_exponent = exponent(maxval(Wc))
                do i = 1, nloc
                    idx = elem(i)
                    if (idx < 1 .or. idx > nc) cycle
                    ia = idx - first
                    if (ia < 0 .or. ia > degree) cycle
                    wi = scale(Wc(i), -weight_exponent)
                    denom = denom + N(ia)*wi
                    ddenom = ddenom + dN(ia)*wi
                end do
            else
                i0 = max(0, 1 - first)
                i1 = min(degree, nc - first)
                weight_exponent = exponent(Wc(first+i0))
                do i = i0 + 1, i1
                    weight_exponent = max(weight_exponent, exponent(Wc(first+i)))
                end do
                do i = i0, i1
                    idx = first + i
                    wi = scale(Wc(idx), -weight_exponent)
                    denom = denom + N(i)*wi
                    ddenom = ddenom + dN(i)*wi
                end do
            end if
            if (abs(denom) > wtol) then
                do i = 1, nloc
                    idx = elem(i)
                    if (idx < 1 .or. idx > nc) cycle
                    ia = idx - first
                    if (local_weights) then
                        wi = scale(Wc(i), -weight_exponent)
                    else
                        wi = scale(Wc(idx), -weight_exponent)
                    end if
                    if (ia >= 0 .and. ia <= degree) then
                        Tgc(i) = N(ia)*wi/denom
                        dTgc(i) = (dN(ia)*wi - Tgc(i)*ddenom)/denom
                    end if
                end do
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate polynomial curve basis values and first derivatives at many parameters.
    !!
    !! Outputs have shape `[ng_,nc]`; each row contains at most `degree+1`
    !! nonzero basis entries. `ng_` defaults to `size(Xt)`.
    pure subroutine compute_dTgc_bspline_1d_vector(Xt, knot, degree, nc, ng, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in), optional :: ng
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        logical, intent(in), optional :: wrap_parameters
        integer :: i, j, idx, first, ng_
        real(rk) :: N(0:degree), dN(0:degree)
        logical :: wrap_parameters_

        if (present(ng)) then
            ng_ = ng
        else
            ng_ = size(Xt, 1)
        end if

        allocate(dTgc(ng_, nc), Tgc(ng_, nc), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N, dN, j, idx, first)
#endif
            call basis_bspline_der(map_parameter(Xt(i), knot, nc, degree, wrap_parameters_), &
                knot, nc, degree, 1, first, N, dN)
            do j = max(0, 1 - first), min(degree, nc - first)
                idx = first + j
                Tgc(i,idx) = N(j)
                dTgc(i,idx) = dN(j)
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate polynomial curve basis values and first derivatives at one parameter.
    !!
    !! Without `elem`, outputs use global basis order and length `nc`. With
    !! `elem`, they use the supplied one-based control-index order and have
    !! length `size(elem)`; inactive or invalid listed indices remain zero.
    pure subroutine compute_dTgc_bspline_1d_scalar(Xt, knot, degree, nc, dTgc, Tgc, elem, wrap_parameters)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in), optional :: elem(:)
        logical, intent(in), optional :: wrap_parameters
        real(rk), allocatable, intent(out) :: dTgc(:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        real(rk) :: N(0:degree), dN(0:degree)
        integer :: i, idx, nloc, first, ia
        logical :: wrap_parameters_

        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters

        if (.not. present(elem)) then
            allocate(dTgc(nc), Tgc(nc), source=0.0_rk)
            call basis_bspline_der(map_parameter(Xt, knot, nc, degree, wrap_parameters_), &
                knot, nc, degree, 1, first, N, dN)
            do i = max(0, 1 - first), min(degree, nc - first)
                idx = first + i
                Tgc(idx) = N(i)
                dTgc(idx) = dN(i)
            end do
        else
            nloc = size(elem)
            allocate(dTgc(nloc), Tgc(nloc), source=0.0_rk)
            call basis_bspline_der(map_parameter(Xt, knot, nc, degree, wrap_parameters_), &
                knot, nc, degree, 1, first, N, dN)
            do i = 1, nloc
                idx = elem(i)
                if (idx < 1 .or. idx > nc) cycle
                ia = idx - first
                if (ia < 0 .or. ia > degree) cycle
                Tgc(i) = N(ia)
                dTgc(i) = dN(ia)
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the dense rational curve basis at many parameters.
    !!
    !! The allocated matrix has shape `[ng_,nc]` and each row is normalized to
    !! partition of unity when its positive-weight denominator is nonzero.
    !! `ng_` defaults to `size(Xt)`.
    pure function compute_Tgc_nurbs_1d_vector(Xt, knot, degree, nc, ng, Wc, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in), optional :: ng
        real(rk), intent(in), contiguous :: Wc(:)
        logical, intent(in), optional :: wrap_parameters
        real(rk), allocatable :: Tgc(:,:)
        real(rk) :: N(0:degree), wsum, wi, wtol
        integer :: i, j, idx, first, ng_, j0, j1, weight_exponent
        logical :: wrap_parameters_

        if (present(ng)) then
            ng_ = ng
        else
            ng_ = size(Xt, 1)
        end if

        allocate(Tgc(ng_, nc), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
        wtol = 0.0_rk
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N, wsum, wi, j, idx, first, j0, j1, weight_exponent)
#endif
            call basis_bspline_der(map_parameter(Xt(i), knot, nc, degree, wrap_parameters_), &
                knot, nc, degree, 0, first, N)
            j0 = max(0, 1 - first)
            j1 = min(degree, nc - first)
            weight_exponent = exponent(Wc(first+j0))
            do j = j0 + 1, j1
                weight_exponent = max(weight_exponent, exponent(Wc(first+j)))
            end do
            wsum = 0.0_rk
            do j = j0, j1
                idx = first + j
                wi = scale(Wc(idx), -weight_exponent)
                wsum = wsum + N(j)*wi
            end do
            if (abs(wsum) > wtol) then
                do j = j0, j1
                    idx = first + j
                    wi = scale(Wc(idx), -weight_exponent)
                    Tgc(i,idx) = N(j)*wi/wsum
                end do
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the rational curve basis at one parameter.
    !!
    !! Without `elem`, the result is dense with length `nc`. With `elem`, its
    !! order follows the supplied global control indices; `Wc` may be global or
    !! element-local. Invalid and inactive listed indices remain zero.
    pure function compute_Tgc_nurbs_1d_scalar(Xt, knot, degree, nc, Wc, elem, wrap_parameters) result(Tgc)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        real(rk), intent(in), contiguous :: Wc(:)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), optional :: wrap_parameters
        real(rk), allocatable :: Tgc(:)
        real(rk) :: N(0:degree), wsum, wi, wtol
        integer :: j, idx, ia, first, nloc, j0, j1, weight_exponent
        logical :: local_weights
        logical :: wrap_parameters_

        if (present(elem)) then
            nloc = size(elem)
            allocate(Tgc(nloc), source=0.0_rk)
        else
            allocate(Tgc(nc), source=0.0_rk)
        end if
        wtol = 0.0_rk
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
        call basis_bspline_der(map_parameter(Xt, knot, nc, degree, wrap_parameters_), &
            knot, nc, degree, 0, first, N)
        wsum = 0.0_rk
        if (present(elem)) then
            local_weights = size(Wc) == nloc .and. size(Wc) /= nc
            if (local_weights) then
                weight_exponent = exponent(maxval(Wc))
                do j = 1, nloc
                    idx = elem(j)
                    if (idx < 1 .or. idx > nc) cycle
                    ia = idx - first
                    if (ia < 0 .or. ia > degree) cycle
                    wi = scale(Wc(j), -weight_exponent)
                    wsum = wsum + N(ia)*wi
                end do
            else
                j0 = max(0, 1 - first)
                j1 = min(degree, nc - first)
                weight_exponent = exponent(Wc(first+j0))
                do j = j0 + 1, j1
                    weight_exponent = max(weight_exponent, exponent(Wc(first+j)))
                end do
                do j = j0, j1
                    idx = first + j
                    wi = scale(Wc(idx), -weight_exponent)
                    wsum = wsum + N(j)*wi
                end do
            end if
            if (abs(wsum) > wtol) then
                do j = 1, nloc
                    idx = elem(j)
                    if (idx < 1 .or. idx > nc) cycle
                    ia = idx - first
                    if (ia < 0 .or. ia > degree) cycle
                    if (local_weights) then
                        wi = scale(Wc(j), -weight_exponent)
                    else
                        wi = scale(Wc(idx), -weight_exponent)
                    end if
                    Tgc(j) = N(ia)*wi/wsum
                end do
            end if
        else
            j0 = max(0, 1 - first)
            j1 = min(degree, nc - first)
            weight_exponent = exponent(Wc(first+j0))
            do j = j0 + 1, j1
                weight_exponent = max(weight_exponent, exponent(Wc(first+j)))
            end do
            do j = j0, j1
                idx = first + j
                wi = scale(Wc(idx), -weight_exponent)
                wsum = wsum + N(j)*wi
            end do
            if (abs(wsum) > wtol) then
                do j = j0, j1
                    idx = first + j
                    wi = scale(Wc(idx), -weight_exponent)
                    Tgc(idx) = N(j)*wi/wsum
                end do
            end if
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the dense polynomial curve basis at many parameters.
    !!
    !! The allocated result has shape `[ng_,nc]`, with at most `degree+1`
    !! nonzero entries per row. Optional wrapping maps each parameter into the
    !! active interval before span selection.
    pure function compute_Tgc_bspline_1d_vector(Xt, knot, degree, nc, ng, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in), optional :: ng
        logical, intent(in), optional :: wrap_parameters
        real(rk), allocatable :: Tgc(:,:)
        real(rk) :: N(0:degree)
        integer :: i, j, idx, first, ng_
        logical :: wrap_parameters_

        if (present(ng)) then
            ng_ = ng
        else
            ng_ = size(Xt, 1)
        end if

        allocate(Tgc(ng_, nc), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N, j, idx, first)
#endif
            call basis_bspline_der(map_parameter(Xt(i), knot, nc, degree, wrap_parameters_), &
                knot, nc, degree, 0, first, N)
            do j = max(0, 1 - first), min(degree, nc - first)
                idx = first + j
                Tgc(i,idx) = N(j)
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the polynomial curve basis at one parameter.
    !!
    !! The result is dense with length `nc` unless `elem` requests an ordered
    !! element-local subset. Listed controls outside the active support, and
    !! invalid listed indices, retain zero values.
    pure function compute_Tgc_bspline_1d_scalar(Xt, knot, degree, nc, elem, wrap_parameters) result(Tgc)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in), optional :: elem(:)
        logical, intent(in), optional :: wrap_parameters
        real(rk), allocatable :: Tgc(:)
        real(rk) :: N(0:degree)
        integer :: j, idx, ia, first, nloc
        logical :: wrap_parameters_

        if (present(elem)) then
            nloc = size(elem)
            allocate(Tgc(nloc), source=0.0_rk)
        else
            allocate(Tgc(nc), source=0.0_rk)
        end if
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) wrap_parameters_ = wrap_parameters
        call basis_bspline_der(map_parameter(Xt, knot, nc, degree, wrap_parameters_), &
            knot, nc, degree, 0, first, N)
        if (present(elem)) then
            do j = 1, nloc
                idx = elem(j)
                if (idx < 1 .or. idx > nc) cycle
                ia = idx - first
                if (ia < 0 .or. ia > degree) cycle
                Tgc(j) = N(ia)
            end do
        else
            do j = max(0, 1 - first), min(degree, nc - first)
                idx = first + j
                Tgc(idx) = N(j)
            end do
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Fit polynomial B-spline control points to parameterized data.
    !!
    !! The existing degree, knot vector, and control count define the trial
    !! space. The first `ndata` rows minimize the unweighted coordinatewise
    !! residual sum. Local support is accumulated directly into symmetric
    !! banded normal equations and solved by banded Cholesky.
    !!
    !! This routine uses the polynomial basis even if explicit weights are
    !! stored; it neither reads nor clears those weights. Consequently
    !! `this%is_rational()` must be false if the fitted controls are to describe
    !! the same objective through subsequent object evaluation. `Xt` and
    !! `Xdata` must contain at least `ndata` finite rows, `Xdata` must have at
    !! least one coordinate, and the parameters must give full-rank active-
    !! domain coverage. These array and finiteness preconditions are not
    !! completely diagnosed. Only `Xc` is replaced; cached `Xg` is not
    !! recomputed.
    !!
    !! @warning Normal equations square the condition number. This routine is
    !! intended for well-parameterized, full-rank spline fits; rank-deficient or
    !! poorly scaled data produce diagnostic code 108 rather than a QR/SVD fit.
    !! @endwarning
    pure subroutine lsq_fit_bspline(this, Xt, Xdata, ndata)
        use forcad_utils, only: solve_spd_banded
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:)
            !! Data parameters with at least `ndata` entries.
        real(rk), intent(in), contiguous :: Xdata(:,:)
            !! Data coordinates `[at least ndata,ncoord]`.
        integer, intent(in) :: ndata
            !! Number of data rows used; must satisfy `ndata>=nc`.
        real(rk), allocatable :: TtT(:,:), TtX(:,:), Xsol(:,:)
        real(rk) :: N(0:this%degree), ba
        integer :: i, a, b, k, ia, ib, first, dim_

        if (.not. this%err%ok) return

        if (this%nc > ndata) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Too few data points for the requested number of control points.",&
                location   = "lsq_fit_bspline",&
                suggestion = "Use nc <= ndata: reduce nc or increase the number of data points.")
            return
        end if

        dim_ = size(Xdata, 2)
        allocate(TtT(0:this%degree, this%nc), TtX(this%nc, dim_), source=0.0_rk)

        do i = 1, ndata
            call basis_bspline_der(&
                map_parameter(Xt(i), this%knot, this%nc, this%degree, this%wrap_parameters),&
                this%knot, this%nc, this%degree, 0, first, N)
            do a = 0, this%degree
                ia = first + a
                if (ia < 1 .or. ia > this%nc) cycle
                ba = N(a)
                do k = 1, dim_
                    TtX(ia,k) = TtX(ia,k) + ba*Xdata(i,k)
                end do
                do b = 0, this%degree
                    ib = first + b
                    if (ib < 1 .or. ib > this%nc) cycle
                    if (ia >= ib) TtT(ia-ib,ib) = TtT(ia-ib,ib) + ba*N(b)
                end do
            end do
        end do
        Xsol = solve_spd_banded(TtT, TtX)
        if (size(Xsol,1) /= this%nc .or. size(Xsol,2) /= dim_) then
            call this%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Banded least-squares solve failed.",&
                location   = "lsq_fit_bspline",&
                suggestion = "Check parameter coverage or reduce the number of control points.")
            return
        end if
        this%Xc = Xsol
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Alternately fit NURBS control points and positive rational weights.
    !!
    !! For fixed weights, control coordinates are solved from banded normal
    !! equations. Log-weights are then updated by a damped Gauss-Newton normal
    !! equation. Subtracting their mean fixes the projective scale at geometric
    !! mean weight one. Each trial weight vector is used to recompute its
    !! least-squares controls before the residual is tested. A trial is accepted
    !! only when its residual norm does not increase; otherwise controls and
    !! weights remain at the last accepted iterate and the damping is increased.
    !! Accepted trials reduce the damping. Log-weight steps outside the finite
    !! exponential range are rejected before exponentiation.
    !!
    !! With residuals
    !! \(r_{ik}=C_k(u_i)-X_{ik}\), the monitored quantity is
    !!
    !! \[
    !! \rho=\frac{\sqrt{\sum_{i,k}r_{ik}^2}}
    !!            {n_{data}n_{coord}},
    !! \]
    !!
    !! which is a scaled Euclidean norm, not the root-mean-square residual.
    !! After an accepted trial, iteration stops when
    !! \(0\le\rho_{old}-\rho\le\mathtt{tol}\max(1,\rho_{old})\), or when
    !! \(\lVert r\rVert_2^2\le\epsilon_{\mathrm{mach}}
    !! \max(1,\lVert X\rVert_2^2)\). The latter recognizes a fit at the
    !! attainable accuracy of the normal-equation solves without accepting an
    !! increasing trial. The default tolerance is
    !! \(\sqrt{\epsilon_{\mathrm{mach}}}\). Exhausting `maxit` sets diagnostic
    !! code 109 and retains the last accepted controls and weights. Only `Xc`,
    !! `Wc`, and rational-state metadata are updated; cached `Xg` is not
    !! recomputed.
    !!
    !! @warning This nonlinear problem is nonconvex and uses normal equations;
    !! the result depends on the initial weights and may be a local fit. Inputs
    !! must provide at least `ndata` finite parameter/data rows and at least one
    !! coordinate, with `ndata>=nc`, `maxit>0`, `tol>=0`, and finite,
    !! nonnegative regularization/damping values.
    !! @endwarning
    !!
    !! Damped least-squares references: K. Levenberg, "A Method for the Solution
    !! of Certain Non-Linear Problems in Least Squares," *Quarterly of Applied
    !! Mathematics* 2 (1944), 164--168
    !! ([doi:10.1090/qam/10666](https://doi.org/10.1090/qam/10666)); D. W.
    !! Marquardt, "An Algorithm for Least-Squares Estimation of Nonlinear
    !! Parameters," *Journal of the Society for Industrial and Applied
    !! Mathematics* 11 (1963), 431--441
    !! ([doi:10.1137/0111030](https://doi.org/10.1137/0111030)).
    pure subroutine lsq_fit_nurbs(this, Xt, Xdata, ndata, maxit, tol, lambda_xc, mu0, reg_logw)
        use forcad_utils, only: solve_spd_banded
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous  :: Xt(:)
            !! Data parameters with at least `ndata` entries.
        real(rk), intent(in), contiguous  :: Xdata(:,:)
            !! Data coordinates `[at least ndata,ncoord]`.
        integer,  intent(in)              :: ndata
            !! Number of rows used; must satisfy `ndata>=nc`.
        integer,  intent(in),  optional   :: maxit
            !! Optional positive iteration limit; default 30.
        real(rk), intent(in),  optional   :: tol
            !! Optional nonnegative relative decrease tolerance; default `sqrt(epsilon(rk))`.
        real(rk), intent(in),  optional   :: lambda_xc
            !! Optional nonnegative diagonal regularization for controls.
        real(rk), intent(in),  optional   :: mu0
            !! Optional nonnegative initial damping of the log-weight update.
        real(rk), intent(in),  optional   :: reg_logw
            !! Optional nonnegative additional log-weight diagonal regularization.
        real(rk), parameter :: damping_max = sqrt(huge(1.0_rk))
        real(rk), parameter :: log_weight_max = log(huge(1.0_rk)) - 1.0_rk
        real(rk), parameter :: log_weight_min = log(tiny(1.0_rk)) + 1.0_rk
        real(rk), allocatable :: TtT(:,:), TtX(:,:), v(:), Xsol(:,:)
        real(rk), allocatable :: Jband(:,:), Jrhs(:,:), delta_u(:,:)
        real(rk) :: N(0:this%degree)
        real(rk) :: Cg(size(Xdata, 2))
        real(rk) :: cost_prev, cost_now, cost_sum, cost_reduction, damp, data_scale, denom, epss, ha, hb, resid, ta, tb
        real(rk) :: tol_, lamx_, mu, regw
        integer :: dim_, it, maxit_, i, j, k, a, b, ia, ib, first
        logical :: converged, trial_pending

        if (.not. this%err%ok) return

        if (this%nc > ndata) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Too few data points for the requested number of control points.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Use nc <= ndata: reduce nc or increase number of data points.")
            return
        end if

        if (this%nc < 2) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Too few control points for NURBS least-squares fitting.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Use nc >= 2 before calling lsq_fit_nurbs.")
            return
        end if

        dim_   = size(Xdata, 2)
        maxit_ = 30
        tol_   = sqrt(epsilon(0.0_rk))
        lamx_  = 0.0_rk
        mu     = sqrt(epsilon(0.0_rk))
        regw   = sqrt(epsilon(0.0_rk))
        epss   = 10.0_rk*epsilon(0.0_rk)
        if (present(maxit))     maxit_ = maxit
        if (present(tol))       tol_   = tol
        if (present(lambda_xc)) lamx_  = lambda_xc
        if (present(mu0))       mu     = mu0
        if (present(reg_logw))  regw   = reg_logw

        if (maxit_ < 1 .or. .not. ieee_is_finite(tol_) .or. tol_ < 0.0_rk .or. &
            .not. ieee_is_finite(lamx_) .or. lamx_ < 0.0_rk .or. &
            .not. ieee_is_finite(mu) .or. mu < 0.0_rk .or. &
            .not. ieee_is_finite(regw) .or. regw < 0.0_rk) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Invalid nonlinear fitting controls.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Use maxit > 0 and finite, nonnegative tolerances and damping.")
            return
        end if
        mu = min(mu, damping_max)
        regw = min(regw, damping_max)

        if (allocated(this%Wc)) then
            if (size(this%Wc) /= this%nc) deallocate(this%Wc)
        end if
        if (.not. allocated(this%Wc)) allocate(this%Wc(this%nc), source=1.0_rk)

        v = log(max(this%Wc, epss))
        v = v - sum(v)/real(this%nc, rk)
        this%Wc = exp(v)
        this%rational = max(maxval(this%Wc)-this%Wc(1),this%Wc(1)-minval(this%Wc)) > &
            32.0_rk*epsilon(1.0_rk)*maxval(this%Wc)

        allocate(TtT(0:this%degree, this%nc), TtX(this%nc, dim_))
        allocate(Jband(0:this%degree, this%nc), Jrhs(this%nc,1))

        data_scale = 0.0_rk
        do i = 1, ndata
            do k = 1, dim_
                data_scale = data_scale + Xdata(i,k)*Xdata(i,k)
            end do
        end do
        data_scale = max(1.0_rk, data_scale)
        cost_prev = huge(1.0_rk)
        converged = .false.
        trial_pending = .false.

        do it = 1, maxit_
            TtT = 0.0_rk
            TtX = 0.0_rk
            do i = 1, ndata
                call basis_bspline_der(&
                    map_parameter(Xt(i), this%knot, this%nc, this%degree, this%wrap_parameters),&
                    this%knot, this%nc, this%degree, 0, first, N)
                denom = 0.0_rk
                do a = 0, this%degree
                    ia = first + a
                    if (ia < 1 .or. ia > this%nc) cycle
                    denom = denom + N(a)*this%Wc(ia)
                end do
                if (abs(denom) < epss) denom = sign(epss, denom)
                do a = 0, this%degree
                    ia = first + a
                    if (ia < 1 .or. ia > this%nc) cycle
                    ta = N(a)*this%Wc(ia)/denom
                    do k = 1, dim_
                        TtX(ia,k) = TtX(ia,k) + ta*Xdata(i,k)
                    end do
                    do b = 0, this%degree
                        ib = first + b
                        if (ib < 1 .or. ib > this%nc) cycle
                        if (ia >= ib) then
                            tb = N(b)*this%Wc(ib)/denom
                            TtT(ia-ib,ib) = TtT(ia-ib,ib) + ta*tb
                        end if
                    end do
                end do
            end do
            if (lamx_ > 0.0_rk) then
                do concurrent (j = 1:this%nc)
                    TtT(0,j) = TtT(0,j) + lamx_
                end do
            end if
            Xsol = solve_spd_banded(TtT, TtX)
            if (size(Xsol,1) /= this%nc .or. size(Xsol,2) /= dim_ .or. &
                any(.not. ieee_is_finite(Xsol))) then
                if (trial_pending) then
                    do concurrent (j = 1:this%nc)
                        this%Wc(j) = exp(v(j))
                    end do
                    mu = min(damping_max, 10.0_rk*max(mu, sqrt(epsilon(1.0_rk))))
                    trial_pending = .false.
                    cycle
                end if
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Banded rational control-point solve failed.",&
                    location   = "lsq_fit_nurbs",&
                    suggestion = "Check parameter coverage or increase lambda_xc.")
                return
            end if

            cost_sum = 0.0_rk
            Jband = 0.0_rk
            Jrhs = 0.0_rk
            do i = 1, ndata
                call basis_bspline_der(&
                    map_parameter(Xt(i), this%knot, this%nc, this%degree, this%wrap_parameters),&
                    this%knot, this%nc, this%degree, 0, first, N)
                denom = 0.0_rk
                do a = 0, this%degree
                    ia = first + a
                    if (ia < 1 .or. ia > this%nc) cycle
                    denom = denom + N(a)*this%Wc(ia)
                end do
                if (abs(denom) < epss) denom = sign(epss, denom)
                Cg = 0.0_rk
                do a = 0, this%degree
                    ia = first + a
                    if (ia < 1 .or. ia > this%nc) cycle
                    ta = N(a)*this%Wc(ia)/denom
                    do k = 1, dim_
                        Cg(k) = Cg(k) + ta*Xsol(ia,k)
                    end do
                end do
                do k = 1, dim_
                    resid = Cg(k) - Xdata(i,k)
                    cost_sum = cost_sum + resid*resid
                    do a = 0, this%degree
                        ia = first + a
                        if (ia < 1 .or. ia > this%nc) cycle
                        ha = N(a)*this%Wc(ia)*(Xsol(ia,k) - Cg(k))/denom
                        Jrhs(ia,1) = Jrhs(ia,1) + ha*resid
                        do b = 0, this%degree
                            ib = first + b
                            if (ib < 1 .or. ib > this%nc) cycle
                            if (ia >= ib) then
                                hb = N(b)*this%Wc(ib)*(Xsol(ib,k) - Cg(k))/denom
                                Jband(ia-ib,ib) = Jband(ia-ib,ib) + ha*hb
                            end if
                        end do
                    end do
                end do
            end do
            cost_now = sqrt(cost_sum) / real(ndata*dim_, rk)

            if (trial_pending) then
                if (.not. ieee_is_finite(cost_now) .or. cost_now > cost_prev) then
                    do concurrent (j = 1:this%nc)
                        this%Wc(j) = exp(v(j))
                    end do
                    mu = min(damping_max, 10.0_rk*max(mu, sqrt(epsilon(1.0_rk))))
                    trial_pending = .false.
                    cycle
                end if
                cost_reduction = cost_prev - cost_now
                this%Xc = Xsol
                v = log(this%Wc)
                v = v - sum(v)/real(this%nc, rk)
                this%Wc = exp(v)
                this%rational = max(maxval(this%Wc)-this%Wc(1),this%Wc(1)-minval(this%Wc)) > &
                    32.0_rk*epsilon(1.0_rk)*maxval(this%Wc)
                mu = max(sqrt(epsilon(1.0_rk)), 0.3_rk*mu)
                trial_pending = .false.
                if (cost_reduction <= tol_*max(1.0_rk, cost_prev)) then
                    converged = .true.
                    exit
                end if
                cost_prev = cost_now
            else
                if (.not. ieee_is_finite(cost_now)) then
                    call this%err%set(&
                        code       = 108,&
                        severity   = 1,&
                        category   = "forcad_nurbs_curve",&
                        message    = "Nonlinear fitting produced a nonfinite residual.",&
                        location   = "lsq_fit_nurbs",&
                        suggestion = "Check the data, weights, and parameter coverage.")
                    return
                end if
                this%Xc = Xsol
                cost_prev = cost_now
            end if

            if (cost_sum <= epsilon(1.0_rk)*data_scale) then
                converged = .true.
                exit
            end if

            damp = max(mu + regw, sqrt(epsilon(1.0_rk)))
            do concurrent (j = 1:this%nc)
                Jband(0,j) = Jband(0,j) + damp
            end do

            delta_u = solve_spd_banded(Jband, Jrhs)
            if (size(delta_u,1) /= this%nc .or. size(delta_u,2) /= 1 .or. &
                any(.not. ieee_is_finite(delta_u))) then
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_curve",&
                    message    = "Sparse rational weight solve failed.",&
                    location   = "lsq_fit_nurbs",&
                    suggestion = "Increase mu0 or reg_logw for the rational weight update.")
                return
            end if
            delta_u(:,1) = v - delta_u(:,1)
            delta_u(:,1) = delta_u(:,1) - sum(delta_u(:,1))/real(this%nc, rk)
            if (any(.not. ieee_is_finite(delta_u(:,1))) .or. &
                maxval(delta_u(:,1)) > log_weight_max .or. &
                minval(delta_u(:,1)) < log_weight_min) then
                mu = min(damping_max, 10.0_rk*max(mu, sqrt(epsilon(1.0_rk))))
            else
                do concurrent (j = 1:this%nc)
                    this%Wc(j) = exp(delta_u(j,1))
                end do
                trial_pending = .true.
            end if
        end do

        if (.not. converged) then
            if (trial_pending) then
                do concurrent (j = 1:this%nc)
                    this%Wc(j) = exp(v(j))
                end do
            end if
            this%rational = max(maxval(this%Wc)-this%Wc(1),this%Wc(1)-minval(this%Wc)) > &
                32.0_rk*epsilon(1.0_rk)*maxval(this%Wc)
            call this%err%set(&
                code       = 109,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Nonlinear NURBS fitting did not converge.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Increase maxit or mu0, or relax tol for this data set.")
        end if
    end subroutine
    !===============================================================================
end module forcad_nurbs_curve
