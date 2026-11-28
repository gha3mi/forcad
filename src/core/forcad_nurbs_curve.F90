! author: Seyed Ali Ghasemi
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
!! \(w_i\mapsto\alpha w_i\), \(\alpha>0\). Rational values scale nonzero
!! weighted basis contributions; derivative recurrences share a weight gauge.
!! Homogeneous refinement selects a common representable gauge for weights and
!! weighted controls and restores the original weights afterward.
!! Rational geometric derivatives use the nested quotient form
!! \((S'-\mathbf C W')/W\), avoiding an explicit \(W^2\).
!!
!! Omitting weights or using exactly equal weights selects the polynomial
!! B-spline form; every nonconstant positive weight list uses rational
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
!! with multiplicity \(p+1\). Geometry and basis queries accept `side=-1` or
!! `side=1` to evaluate the adjacent span at the knot. The default takes the
!! right trace internally. At a closed endpoint both sides use the inward trace.
!!
!! Refinement operates on homogeneous controls
!! \(\mathbf H_i=(w_i\mathbf P_i,w_i)\) and preserves
!! \(\mathbf C_{\mathrm{old}}(u)=\mathbf C_{\mathrm{new}}(u)\) to working
!! precision on the active domain. Knot insertion and positive degree elevation
!! preserve a verified periodic topology, including its exterior extension and
!! repeated control rows. Other non-clamped representations are converted by
!! degree elevation to an equivalent active-domain open-clamped form. Interior
!! knot removal applies the bounded removal test to a guarded cyclic extension
!! for periodic data and commits only complete periodic removals. For integration,
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
    use forcad_utils, only: valid_trace_side, tensor_ansatz, tensor_ansatz_element, curve_frame
    use forcad_utils, only: elemConn_C0, ndgrid, compute_multiplicity, compute_knot_vector, basis_bspline_der, &
        geometry_derivative, coarsen_control_lines, reduce_control_lines, refine_control_net, get_array, &
        spline_connectivity, spline_span_plan, bezier_element_data, repelem, unique, rotate_points, transform_geometry, &
        gauss_leg, export_vtk_legacy, valid_knot_vector, knot_start, knot_end, active_knots, &
        active_span_count, tensor_basis_derivative_at, &
        basis_bspline_der_all_active, tensor_basis_derivatives2_local, map_parameter, periodic_topology, valid_weights, &
        nonuniform_weights, query_parameter_axis, normalize_rational_basis, restrict_control_net, map_parameter_values, &
        reparameterize_knots, reorient_patch, nearest_spline_point, nearest_sample_index
    use fordebug, only: debug

    implicit none

    private
    public nurbs_curve
#ifndef FORCAD_PUBLIC_DOCS
    public compute_Tgc, compute_dTgc
#endif

    !===============================================================================
    ! author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
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
    !! every constructor validates these requirements before replacing geometry.
    !!
    !! Diagnostic codes use the following categories: 100 invalid argument or
    !! index; 101 shape mismatch; 102--105 missing geometry state; 106
    !! underdetermined fitting problem; 107 invalid parameters; 108 singular
    !! system or mapping; and 109 iterative nonconvergence.
    !!
    !! @note Query methods return allocatable copies of owned arrays. Private
    !! storage prevents direct external mutation; mutating methods validate the
    !! core state they replace. Successful setters discard sampled points and
    !! connectivity while retaining parameters for create() to refresh.
    !! @endnote
    !! @warning Fitting and the separate
    !! `rotate_Xc` or `translate_Xc` operations do not automatically recompute
    !! every sampled or connectivity cache. [[nurbs_curve:transform]] instead
    !! keeps an existing `Xg` cache affine-consistent with `Xc`. Rebuild
    !! connectivity after a spline-space change. The fixed-shape
    !! modify_Xc/modify_Wc operations discard sampled points and display
    !! connectivity while keeping parameters for create() to refresh.
    !! @endwarning
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
        procedure, private :: cmp_Xg_scalar, cmp_Xg_vector
        generic :: cmp_Xg => cmp_Xg_scalar, cmp_Xg_vector !!> Evaluate points without changing sampled geometry.
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
        procedure :: cmp_Xth_in_Xg         !!> Build the embedded active parameter line.
        procedure :: export_Xth_in_Xg      !!> Export parameter-grid lines embedded in physical space.
        procedure :: export_iges           !!> Export the curve as IGES entity 126.
        procedure, private :: modify_Xc_scalar, modify_Xc_all
        generic :: modify_Xc => modify_Xc_scalar, modify_Xc_all !!> Replace finite controls, retaining the sampling parameters.
        procedure, private :: modify_Wc_scalar, modify_Wc_all
        generic :: modify_Wc => modify_Wc_scalar, modify_Wc_all !!> Replace positive weights, retaining the sampling parameters.
        procedure :: get_multiplicity      !!> Return multiplicities of all distinct stored knots.
        procedure :: get_continuity        !!> Return \(p-s\) metadata for all distinct stored knots.
        procedure :: cmp_nc                !!> Compute the control-point count implied by knots and degree.
        procedure, private :: get_nc_dir   !!> Return the control-point count after validating direction 1.
        procedure, private :: get_nc_all   !!> Return the control-point count.
        generic :: get_nc => get_nc_all, get_nc_dir !!> Query the control-point count.
        procedure :: insert_knots          !!> Insert knots while preserving active-domain geometry.
        procedure :: elevate_degree        !!> Elevate degree while preserving the active-domain geometry.
        procedure :: reduce_degree         !!> Reduce degree under one physical geometry-error budget.
        procedure, private :: bezier_element, bezier_all
        generic :: bezier => bezier_element, bezier_all !!> Extract one exact Bezier element or decompose the patch.
        procedure :: restrict              !!> Extract a bounded parameter interval without reparameterization.
        procedure :: reverse               !!> Reverse the parameter direction on the original domain.
        procedure, private :: reparameterize_interval, reparameterize_axes
        generic :: reparameterize => reparameterize_interval, reparameterize_axes !!> Rescale or reorder parameter axes.
        procedure :: split                 !!> Split into two bounded curves at an interior parameter.
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
        procedure :: transform             !!> Apply one affine map to control points and cached geometry.
        procedure :: rotate_Xc             !!> Rotate control points in three-dimensional embedding space.
        procedure :: rotate_Xg             !!> Rotate cached geometry points in three-dimensional embedding space.
        procedure :: translate_Xc          !!> Translate control points by a physical-space vector.
        procedure :: translate_Xg          !!> Translate cached geometry points by a physical-space vector.
        procedure :: nearest_point         !!> Select the nearest cached geometry sample.
        procedure :: nearest_point2        !!> Compute a bounded continuous nearest-point candidate with optional warm start.
        procedure :: ansatz                !!> Evaluate element-local physical basis data and differential length.
        procedure :: ansatz_element        !!> Evaluate all quadrature data for one analysis element.
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
        module procedure compute_Xg_vector
        module procedure compute_Xg_scalar
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
    !! The `Tgc` output returns basis values computed in the same pass.
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

    !> Evaluate only the requested dense basis data using local-support scratch.
    interface compute_basis_data
        module procedure compute_basis_data_vector
        module procedure compute_basis_data_scalar
    end interface


contains

    !> Physical coordinate count without evaluating an unavailable allocatable.
    pure integer function curve_ncoord(this) result(ncoord)
        class(nurbs_curve), intent(in) :: this
        ncoord = 0
        if (allocated(this%Xc)) ncoord = size(this%Xc,2)
    end function curve_ncoord


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Apply one affine transformation to the curve in place.
    !!
    !! Control points and an existing sampled-geometry cache are transformed by
    !!
    !! \[
    !! \widetilde{\mathbf x}
    !! =\mathbf c+\mathbf A(\mathbf x-\mathbf c)+\mathbf b ,
    !! \]
    !!
    !! with `linear` \(=\mathbf A\), optional `center` \(=\mathbf c\), and
    !! optional `translation` \(=\mathbf b\). Because the NURBS basis is a
    !! partition of unity,
    !!
    !! \[
    !! \sum_i R_i\widetilde{\mathbf P}_i
    !! =\mathbf c+\mathbf A
    !!   \left(\sum_iR_i\mathbf P_i-\mathbf c\right)+\mathbf b ,
    !! \]
    !!
    !! so transforming the control polygon gives the exact affine image of the
    !! curve. Knots, weights, parameters, and connectivity remain unchanged.
    !! A singular `linear` map is mathematically valid but may reduce geometric
    !! rank and therefore invalidate an IGA mapping.
    pure subroutine transform(this, linear, translation, center)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: linear(:,:)
            !! Square linear map `[ncoord,ncoord]`.
        real(rk), intent(in), contiguous, optional :: translation(:)
            !! Translation vector `[ncoord]`; zero when absent.
        real(rk), intent(in), contiguous, optional :: center(:)
            !! Transformation center `[ncoord]`; the origin when absent.

        call transform_geometry(this%Xc,this%Xg,linear,translation,center,this%err,"forcad_nurbs_curve")
    end subroutine
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
            if (this%rational .neqv. nonuniform_weights(this%Wc)) return
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
    !! this size relation. `Xc` must contain at least one control row and one
    !! finite coordinate component. The operation validates coordinates, knots,
    !! and weights before replacing any core state.
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

        if (size(Xc,1) < 1 .or. size(Xc,2) < 1) then
            call this%err%set(&
                code       = 101,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Control points need at least one row and one coordinate.",&
                location   = "set1",&
                suggestion = "Provide Xc with shape [nc,ncoord], where nc and ncoord are positive.")
            return
        end if
        if (.not. all(ieee_is_finite(Xc))) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Control-point coordinates must all be finite.",&
                location   = "set1",&
                suggestion = "Remove NaN and infinite values from Xc before calling set().")
            return
        end if

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

        call invalidate_samples(this,space_changed=.true.)
        this%Xc = Xc
        this%nc = nc_
        this%knot = knot
        this%degree = degree_
        this%wrap_parameters = wrap_parameters_
        if (present(Wc)) then
            this%Wc = Wc
            this%rational = nonuniform_weights(this%Wc)
        else
            if (allocated(this%Wc)) deallocate(this%Wc)
            this%rational = .false.
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate all quadrature points of one curve analysis element.
    !!
    !! This is the batched counterpart of [[nurbs_curve:ansatz]]. It evaluates
    !! the same rational basis, tangential physical gradients, optional
    !! arc-length second derivatives, and weighted differential lengths, while
    !! constructing the quadrature rule and locating the knot span only once.
    !! The local basis index is the first array dimension so every quadrature
    !! column is contiguous:
    !!
    !! - `Tgc(nlocal,nquadrature)`;
    !! - `dTgc_dXg(nlocal,ncoord,nquadrature)`;
    !! - `dL(nquadrature)`;
    !! - optional `Xg(ncoord,nquadrature)` and `dXg_dXt(ncoord,1,nquadrature)`; and
    !! - optional `d2Tgc_dXg2(nlocal,nquadrature)`.
    !!
    !! Allocatable outputs have `intent(inout)`: storage with the required
    !! shape is reused, while a shape change triggers reallocation. This avoids
    !! one allocation per quadrature point in element assembly. If evaluation
    !! fails at any point, every returned entry is reset to zero and the
    !! diagnostic identifies the failure. Optional `err` uses caller-owned
    !! diagnostics and leaves `this%err` unchanged for parallel element assembly;
    !! without it, the existing sticky object diagnostic is retained.
    !!
    !! The formulas, strict regularity policy, and references are those of
    !! [[nurbs_curve:ansatz]].
    pure subroutine ansatz_element(this, ie, Tgc, dTgc_dXg, dL, ngauss, Xg, d2Tgc_dXg2, strict, err, dXg_dXt, &
        dmeasure_dXc, dmeasure_dlogWc, dgradient_dXc, dgradient_dlogWc, djacobian_dXc, djacobian_dlogWc)
        class(nurbs_curve), target, intent(inout) :: this
        integer, intent(in) :: ie
            !! One-based active element index.
        real(rk), allocatable, intent(inout) :: Tgc(:,:), dTgc_dXg(:,:,:), dL(:)
            !! Basis `[nlocal,nq]`, physical gradients `[nlocal,ncoord,nq]`, weighted lengths `[nq]`.
        integer, intent(in), optional :: ngauss
            !! Positive quadrature count; default degree+1.
        real(rk), allocatable, target, intent(inout), optional :: Xg(:,:), d2Tgc_dXg2(:,:)
            !! Points `[ncoord,nq]` and second arc-length derivatives `[nlocal,nq]`.
        logical, intent(in), optional :: strict
            !! Reject singular tangents and reversed one-coordinate mappings.
        type(debug), target, intent(out), optional :: err
            !! Caller-owned diagnostic; leaves the object diagnostic unchanged.
        real(rk), allocatable, target, intent(inout), optional :: dXg_dXt(:,:,:)
            !! Raw parametric Jacobians `[ncoord,1,nq]`, before tangent normalization.
        real(rk), allocatable, target, intent(inout), optional :: dmeasure_dXc(:,:,:), dmeasure_dlogWc(:,:)
            !! Derivatives of weighted length with respect to active Euclidean controls and log weights.
        real(rk), allocatable, target, intent(inout), optional :: dgradient_dXc(:,:,:,:,:), dgradient_dlogWc(:,:,:,:)
            !! Physical gradient derivatives indexed by basis, component, active control, and (for coordinates) design component.
        real(rk), allocatable, intent(inout), optional :: djacobian_dXc(:,:,:), djacobian_dlogWc(:,:,:,:)
            !! Raw Jacobian derivatives for controls (A,alpha,q) and log weights (i,alpha,A,q).
        integer :: counts(1)
        counts=this%degree+1
        if (present(ngauss)) counts=ngauss
        call tensor_ansatz_element(this%Xc,this%Wc,this%rational,[this%degree],[this%nc],this%knot, &
            object_err=this%err,err=err,ie=ie,basis=Tgc,gradient=dTgc_dXg,measure=dL,ngauss=counts,strict=strict, &
            Xg=Xg,dXg_dXt=dXg_dXt,curve_second=d2Tgc_dXg2,dmeasure_dXc=dmeasure_dXc,dmeasure_dlogWc=dmeasure_dlogWc, &
            dgradient_dXc=dgradient_dXc,dgradient_dlogWc=dgradient_dlogWc, &
            djacobian_dXc=djacobian_dXc,djacobian_dlogWc=djacobian_dlogWc)
    end subroutine ansatz_element
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
        real(rk), allocatable :: Xg(:,:)
        integer, allocatable :: elemConn(:,:)

        call this%cmp_Xth_in_Xg(&
            Xg       = Xg,&
            elemConn = elemConn,&
            res      = res)
        if (.not. this%err%ok) return

        call export_vtk_legacy(&
            filename    = filename,&
            points      = Xg,&
            elemConn    = elemConn,&
            vtkCellType = 4,&
            encoding    = encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build the active parameter line embedded in physical space.
    !!
    !! `Xg` contains the evaluated line points and the single row of `elemConn`
    !! is VTK polyline connectivity. No file is written.
    impure subroutine cmp_Xth_in_Xg(this, Xg, elemConn, res)
        class(nurbs_curve), intent(inout) :: this
        real(rk), allocatable, intent(out) :: Xg(:,:)
        integer, allocatable, intent(out) :: elemConn(:,:)
        integer, intent(in), optional :: res
        real(rk), allocatable :: U(:), Ur(:), Ur_eval(:)
        integer :: nspan, res_min, n, s, r, o, ncoord

        if (.not. this%err%ok) return

        if (.not. allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Control points are not set.",&
                location   = "cmp_Xth_in_Xg",&
                suggestion = "Call set(...) first before exporting.")
            return
        end if

        if (.not. allocated(this%knot)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Knot vector is not set.",&
                location   = "cmp_Xth_in_Xg",&
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
                location   = "cmp_Xth_in_Xg",&
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
                location   = "cmp_Xth_in_Xg",&
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
            Xg = compute_Xg(Ur_eval, this%knot, this%degree, this%nc, n, this%Xc, wrap_parameters=this%wrap_parameters)
        end if

        do concurrent (r = 1:n)
            elemConn(1,r) = r
        end do
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
        real(rk), allocatable :: controls(:,:)

        if (.not. this%err%ok) return
        allocate(controls(size(Xc),3), source=0.0_rk)
        controls(:,1) = Xc
        call this%set1(knot, controls, Wc, degree, wrap_parameters)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct a spline space from breakpoints and requested continuities.
    !!
    !! `continuity(i)` is converted to knot multiplicity `degree-continuity(i)`
    !! for each parameter node according to `compute_knot_vector`. Optional
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
            call this%set1(knot_, Xc, Wc, degree=degree)
            return
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

        call invalidate_samples(this,space_changed=.true.)
        this%knot = knot_
        this%degree = degree
        this%nc = nc_
        if (allocated(this%Xc)) deallocate(this%Xc)

        if (present(Wc)) then
            this%Wc = Wc
            this%rational = nonuniform_weights(this%Wc)
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
        real(rk), allocatable :: knots(:)
        integer :: nc

        if (.not. this%err%ok) return
        nc = size(Xc,1)
        allocate(knots(2*nc))
        knots(:nc) = 0.0_rk
        knots(nc+1:) = 1.0_rk
        call this%set1(knots, Xc, Wc, degree=nc-1)
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
        real(rk), allocatable :: knots(:)
        integer :: m, i

        if (.not. this%err%ok) return
        if (degree < 0 .or. nc <= degree .or. size(Xc,1) /= nc) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Invalid curve dimensions for generated knot-vector construction.",&
                location   = "set4",&
                suggestion = "Require 0 <= degree < nc and size(Xc,1) == nc.")
            return
        end if
        m = nc + degree + 1
        allocate(knots(m))
        knots(:degree+1) = 0.0_rk
        knots(degree+2:nc) = [(real(i,rk)/real(nc-degree,rk), i=1,nc-degree-1)]
        knots(nc+1:) = 1.0_rk
        call this%set1(knots, Xc, Wc, degree)
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
    !! Parameters, geometry samples, and display connectivity are committed
    !! together after validation; a failed request preserves the previous cache.
    pure subroutine create(this, res, Xt)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: res
            !! Number of generated samples when `Xt` is absent.
        real(rk), intent(in), contiguous, optional :: Xt(:)
            !! Optional explicit parameter vector.
        real(rk), allocatable :: parameters(:), points(:,:)
        integer, allocatable :: connectivity(:,:)
        logical :: valid

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Control points are not set.",&
                location   = "create",&
                suggestion = "Call set(...) before create().")
            return
        end if
        call query_parameters(this, "create", parameters, valid, res, Xt)
        if (.not. valid) return

        if (this%is_rational()) then
            points = compute_Xg(parameters, this%knot, this%degree, this%nc, size(parameters), &
                this%Xc, this%Wc, this%wrap_parameters)
        else
            points = compute_Xg(parameters, this%knot, this%degree, this%nc, size(parameters), &
                this%Xc, wrap_parameters=this%wrap_parameters)
        end if
        if (.not. all(ieee_is_finite(points))) then
            call this%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Sampled coordinates exceed the finite numerical range.",&
                location   = "create",&
                suggestion = "Rescale the control coordinates or use representable parameter values.")
            return
        end if
        connectivity = elemConn_C0(size(parameters), 1)
        this%ng = size(parameters)
        call move_alloc(parameters, this%Xt)
        call move_alloc(points, this%Xg)
        call move_alloc(connectivity, this%elemConn_Xg_vis)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate one curve point or parametric derivative.
    !! Returns a vector of physical coordinate dimension. Outside-domain values
    !! evaluate to the zero basis unless parameter wrapping is enabled.
    !! A nonfinite parameter or unavailable geometry returns an empty vector.
    pure function cmp_Xg_scalar(this, Xt, order, side) result(Xg)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in), optional :: side !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in) :: Xt
        integer, intent(in), optional :: order !! Nonnegative parametric derivative order; omitted means position.
        real(rk), allocatable :: Xg(:), rows(:,:)

        rows=cmp_Xg_vector(this,[Xt],order,side)
        Xg=reshape(rows,[size(rows)])
    end function cmp_Xg_scalar


    !> Evaluate arbitrary parameter rows without modifying caches or diagnostics.
    !! Strided inputs are accepted. A valid empty batch returns [0,ncoord];
    !! nonfinite parameters or missing geometry return [0,0].
    pure function cmp_Xg_vector(this, Xt, order, side) result(Xg)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in), optional :: side !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: order !! Nonnegative order; omitted means position.
        real(rk), allocatable :: Xg(:,:)
        integer :: selected(1), trace(1)

        trace=curve_trace(side)
        selected=0
        if (present(order)) selected(1)=order
        if (any(selected<0) .or. .not. curve_geometry_ready(this) .or. &
            .not. valid_trace_side(trace,1) .or. .not. all(ieee_is_finite(Xt))) then
            allocate(Xg(0,0))
            return
        end if
        if (any(selected>0)) then
            call geometry_derivative(Xt,[this%nc],[this%degree],this%knot,this%Xc,this%Wc, &
                selected,[this%wrap_parameters],Xg,side=trace)
        else if (this%is_rational()) then
            Xg = compute_Xg(Xt,this%knot,this%degree,this%nc,size(Xt),this%Xc,this%Wc,this%wrap_parameters,side=side)
        else
            Xg = compute_Xg(Xt,this%knot,this%degree,this%nc,size(Xt),this%Xc,wrap_parameters=this%wrap_parameters,side=side)
        end if
    end function cmp_Xg_vector


    !> Constant-work validation of stored geometry shapes before read-only queries.
    pure logical function curve_geometry_ready(this) result(ready)
        class(nurbs_curve), intent(in) :: this
        ready = .false.
        if (.not. this%err%ok .or. .not. allocated(this%Xc) .or. .not. allocated(this%knot)) return
        if (this%degree < 0 .or. this%nc <= this%degree) return
        if (size(this%Xc,1) /= this%nc .or. size(this%Xc,2) < 1) return
        if (size(this%knot) /= this%nc+this%degree+1) return
        if (allocated(this%Wc)) then
            if (size(this%Wc) /= this%nc) return
        end if
        ready = .true.
    end function curve_geometry_ready
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
        real(rk) :: min_X1, max_X1

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

        min_X1 = minval(X(:,1))
        max_X1 = maxval(X(:,1))
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

        call map_parameter_values(X(:,1),min_X1,max_X1,knot_start(this%knot,this%nc,this%degree), &
            knot_end(this%knot,this%nc,this%degree),Xt,this%err)
        if (.not.this%err%ok) return
        Xg_new = this%cmp_Xg(Xt)
        if (size(Xg_new,1)/=size(X,1) .or. .not.all(ieee_is_finite(Xg_new))) then
            call this%err%set(code=108,severity=1,category="forcad_nurbs_curve",location="put_to_nurbs", &
                message="Mapped curve points must be finite.", &
                suggestion="Use geometry and ordering data with representable mapped points.")
            return
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

        call get_array(this%Xc,Xc,this%err%ok)
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
    !! Zero denotes an invalid index or unavailable control data.
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

        call get_array(this%Xg,Xg,this%err%ok)
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
    !! Zero denotes an invalid index or unavailable sampled data.
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

        call get_array(this%Wc,Wc,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return explicit control weight `n`.
    !! Zero denotes an invalid index or unavailable weight data.
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
    !! The result has the retained grid size, or zero when no grid exists.
    !! It survives modify_Xc/Wc even while get_ng() reports no sampled points.
    pure function get_Xt(this) result(Xt)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: Xt(:)

        call get_array(this%Xt,Xt,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the number of cached geometry points.
    !! Zero denotes unavailable sampled points, including after modify_Xc/Wc.
    !! Their retained parameters can still be reused by create().
    pure elemental function get_ng(this) result(ng)
        class(nurbs_curve), intent(in) :: this
        integer :: ng

        ng = 0
        if (.not. this%err%ok .or. .not. allocated(this%Xg)) return

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

        call get_array(this%knot,knot,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return knot `i`.
    !! Zero denotes an invalid index or object state.
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
    !! Optional `p` defaults to one. Missing samples return an empty array.
    pure function cmp_elem_Xg_vis(this, p) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), optional :: p

        if (.not. this%err%ok .or. .not. allocated(this%Xg)) then
            allocate(elemConn(0,0))
            return
        end if

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
    !! including its exact uniform-weight classification.
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
    !> Replace one finite control coordinate, preserving the sampling grid.
    !! Sampled points and their display connectivity are invalidated. Call
    !! [[nurbs_curve:create]] with no arguments to refresh the retained grid.
    pure subroutine modify_Xc_scalar(this, X, num, dir)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: X
        integer, intent(in) :: num, dir

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) then
            call this%err%set(code=102, severity=1, category="forcad_nurbs_curve", &
                message="Control points are not set.", location="modify_Xc", &
                suggestion="Call set(...) before modifying controls.")
        else if (num < 1 .or. num > this%nc .or. dir < 1 .or. dir > size(this%Xc,2)) then
            call this%err%set(code=100, severity=1, category="forcad_nurbs_curve", &
                message="Invalid control point index or direction.", location="modify_Xc", &
                suggestion="Use 1 <= num <= nc and a valid coordinate direction.")
        else if (.not. ieee_is_finite(X)) then
            call this%err%set(code=100, severity=1, category="forcad_nurbs_curve", &
                message="Control coordinates must be finite.", location="modify_Xc", &
                suggestion="Provide a finite coordinate value.")
        else
            this%Xc(num,dir) = X
            call invalidate_samples(this)
        end if
    end subroutine modify_Xc_scalar


    !> Atomically replace all control coordinates with the same shape.
    !! Failed validation preserves controls and caches. Success invalidates only
    !! sampled points and display connectivity, retaining the parameter grid.
    pure subroutine modify_Xc_all(this, Xc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xc(:,:)

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) then
            call this%err%set(code=102, severity=1, category="forcad_nurbs_curve", &
                message="Control points are not set.", location="modify_Xc", &
                suggestion="Call set(...) before modifying controls.")
        else if (any(shape(Xc) /= shape(this%Xc)) .or. .not. all(ieee_is_finite(Xc))) then
            call this%err%set(code=101, severity=1, category="forcad_nurbs_curve", &
                message="Replacement controls must have the existing shape and finite coordinates.", &
                location="modify_Xc", suggestion="Use shape(Xc) == shape(get_Xc()) and finite coordinates.")
        else
            this%Xc = Xc
            call invalidate_samples(this)
        end if
    end subroutine modify_Xc_all


    !> Replace one positive weight and recompute the rational state.
    !! Polynomial spaces begin with implicit unit weights.
    !! Sampling parameters survive; sampled points must be refreshed by create().
    pure subroutine modify_Wc_scalar(this, W, num)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: W
        integer, intent(in) :: num

        if (.not. this%err%ok) return
        if (num < 1 .or. num > this%nc) then
            call this%err%set(code=100, severity=1, category="forcad_nurbs_curve", &
                message="Invalid weight index.", location="modify_Wc", &
                suggestion="Use 1 <= num <= nc.")
        else if (.not. ieee_is_finite(W) .or. W <= 0.0_rk) then
            call this%err%set(code=105, severity=1, category="forcad_nurbs_curve", &
                message="Weights must be finite and positive.", location="modify_Wc", &
                suggestion="Provide a strictly positive finite weight.")
        else
            if (.not. allocated(this%Wc)) allocate(this%Wc(this%nc), source=1.0_rk)
            this%Wc(num) = W
            this%rational = nonuniform_weights(this%Wc)
            call invalidate_samples(this)
        end if
    end subroutine modify_Wc_scalar


    !> Atomically replace or initialize the complete positive weight vector.
    !! The control count stays fixed. Failed validation preserves weights and
    !! caches; success retains parameters for a subsequent create() refresh.
    pure subroutine modify_Wc_all(this, Wc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Wc(:)

        if (.not. this%err%ok) return
        if (size(Wc) /= this%nc .or. .not. valid_weights(Wc)) then
            call this%err%set(code=105, severity=1, category="forcad_nurbs_curve", &
                message="Replacement weights must match nc and be finite and positive.", &
                location="modify_Wc", suggestion="Supply nc strictly positive finite weights.")
        else
            this%Wc = Wc
            this%rational = nonuniform_weights(this%Wc)
            call invalidate_samples(this)
        end if
    end subroutine modify_Wc_all


    !> Drop derived geometry, optionally including topology, while keeping parameters.
    pure subroutine invalidate_samples(this,space_changed)
        class(nurbs_curve), intent(inout) :: this
        logical, intent(in), optional :: space_changed
        if (allocated(this%Xg)) deallocate(this%Xg)
        if (allocated(this%elemConn_Xg_vis)) deallocate(this%elemConn_Xg_vis)
        if (present(space_changed)) then
            if (space_changed) then
                if (allocated(this%elemConn)) deallocate(this%elemConn)
                if (allocated(this%elemConn_Xc_vis)) deallocate(this%elemConn_Xc_vis)
            end if
        end if
    end subroutine invalidate_samples


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

        c=this%get_multiplicity()
        c=this%get_degree()-c
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


    !> Extract one exact element; optional polynomial extraction satisfies N=C*B.
    !! Element and support IDs use direction-1-fastest stored-control ordering.
    !! Bounds have shape [2,1]. The piece retains literal span coordinates,
    !! with clamped knots, right/left endpoint traces and disabled wrapping.
    !! Caller err bypasses the source diagnostic; source geometry and caches are untouched.
    pure subroutine bezier_element(this, ie, piece, C, ids, bounds, err)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: ie
        type(nurbs_curve), intent(out), optional :: piece
        real(rk), allocatable, intent(out), optional :: C(:,:), bounds(:,:)
        integer, allocatable, intent(out), optional :: ids(:)
        type(debug), intent(out), optional :: err

        call bezier_extract(this,ie,piece,C,ids,bounds,err)
    end subroutine bezier_element


    !> Share extraction between direct queries and a decomposition's precomputed span plan.
    pure subroutine bezier_extract(this, ie, piece, C, ids, bounds, err, spans)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: ie
        integer, intent(in), optional :: spans(:)
        type(nurbs_curve), intent(out), optional :: piece
        real(rk), allocatable, intent(out), optional :: C(:,:), bounds(:,:)
        integer, allocatable, intent(out), optional :: ids(:)
        type(debug), intent(out), optional :: err
        real(rk), allocatable :: Xc(:,:), Wc(:)
        integer, allocatable :: nodes(:)
        real(rk) :: limits(2,1)
        integer :: p(1)
        type(debug) :: status

        call bezier_check(this,status,present(err))
        if (status%ok) then
            p=[this%degree]
            call bezier_element_data(this%knot,[this%nc],p,ie,this%Xc,this%Wc, &
                present(piece),Xc,Wc,nodes,limits,status,C=C,rational=this%rational,spans=spans)
            if (status%ok .and. present(piece)) then
                call piece%set(knot=repelem(limits(:,1),[p(1)+1,p(1)+1]), &
                    Xc=Xc,Wc=Wc,degree=p(1))
                status=piece%err
            end if
        end if
        if (status%ok) then
            if (present(ids)) ids=nodes
            if (present(bounds)) bounds=limits
        else
            if (present(C)) then
                if (allocated(C)) deallocate(C)
                allocate(C(0,0))
            end if
            if (present(ids)) allocate(ids(0))
            if (present(bounds)) allocate(bounds(0,0))
        end if
        if (present(piece)) piece%err=status
        if (present(err)) err=status
    end subroutine bezier_extract


    !> Decompose all positive active spans in cmp_elem ordering, without source sampling.
    !! Failure returns an empty array; request err to inspect the diagnostic.
    pure subroutine bezier_all(this, pieces, err)
        class(nurbs_curve), intent(in) :: this
        type(nurbs_curve), allocatable, intent(out) :: pieces(:)
        type(debug), intent(out), optional :: err
        type(debug) :: status
        integer, allocatable :: plan(:,:)
        integer :: ie, ne(1)

        call bezier_check(this,status,present(err))
        if (status%ok) then
            plan=spline_span_plan([this%nc],[this%degree],this%knot)
            ne=count(plan>0,dim=1)
            allocate(pieces(product(ne)))
            do ie=1,size(pieces)
                call bezier_extract(this,ie,pieces(ie),err=status,spans=plan(ie,:))
                if (.not.status%ok) exit
            end do
        end if
        if (.not.status%ok) then
            if (allocated(pieces)) deallocate(pieces)
            allocate(pieces(0))
        end if
        if (present(err)) err=status
    end subroutine bezier_all


    !> Validate owned geometry presence independently of optional caller diagnostics.
    pure subroutine bezier_check(this, err, bypass)
        class(nurbs_curve), intent(in) :: this
        type(debug), intent(out) :: err
        logical, intent(in) :: bypass

        err=this%err
        if (bypass) call err%reset()
        if (.not.err%ok) return
        if (.not.allocated(this%Xc) .or. &
            .not.allocated(this%knot)) call err%set( &
            code=103,severity=1,category="forcad_nurbs_curve",location="bezier", &
            message="Spline geometry is not initialized.",suggestion="Set knots and controls before Bezier extraction.")
    end subroutine bezier_check


    !> Extract the curve on the literal active-domain interval `[a,b]`.
    !! Original parameter values and degree are retained, with clamped endpoints
    !! and wrapping disabled. The lower endpoint uses its right trace and the
    !! upper endpoint its left trace, including at a discontinuous interior knot.
    !! The source remains unchanged; the result has no sample or connectivity
    !! caches. Invalid bounds or source diagnostics belong to the returned curve.
    pure function restrict(this,a,b) result(patch)
        class(nurbs_curve), intent(in) :: this
        real(rk), intent(in) :: a,b
        type(nurbs_curve) :: patch
        real(rk), allocatable :: knots(:), Xc(:,:), Wc(:)

        call prepare_result(this,patch,"restrict")
        if (.not.patch%err%ok) return
        call restrict_control_net(this%Xc,this%Wc,[this%nc],1,this%knot,this%degree,a,b,knots,Xc,Wc,patch%err)
        if (patch%err%ok) call patch%set(knot=knots,Xc=Xc,Wc=Wc,degree=this%degree)
    end function restrict


    !> Split at a literal parameter strictly inside the active interval.
    !! Each result follows restrict(): fresh caches, original parameters and no
    !! wrapping. At a discontinuity the left result ends with the left trace and
    !! the right result starts with the right trace; these points need not coincide.
    !! Failure leaves both results empty with the same diagnostic.
    pure subroutine split(this,Xt,left,right)
        class(nurbs_curve), intent(in) :: this
        real(rk), intent(in) :: Xt
        type(nurbs_curve), intent(out) :: left,right
        real(rk) :: lo,hi

        call prepare_result(this,left,"split")
        right%err=left%err
        if (.not.left%err%ok) return
        lo=knot_start(this%knot,this%nc,this%degree)
        hi=knot_end(this%knot,this%nc,this%degree)
        if (.not.ieee_is_finite(Xt) .or. Xt<=lo .or. Xt>=hi) then
            call left%err%set(code=107,severity=1,category="forcad_nurbs_curve",location="split", &
                message="The split parameter must be finite and strictly inside the active interval.", &
                suggestion="Use an interior parameter without wrapping or snapping.")
        else
            left=this%restrict(lo,Xt)
            right=this%restrict(Xt,hi)
        end if
        if (.not.left%err%ok) right%err=left%err
        if (.not.right%err%ok) then
            left%err=right%err
            call left%finalize()
            call right%finalize()
        end if
    end subroutine split


    !> Reverse the complete knot vector and control order on the original domain.
    !! Degree, explicit weights and wrapping are retained, with fresh caches.
    !! Away from breaks, the kth parameter derivative changes by (-1)**k.
    !! The usual right interior trace becomes the source left trace at a break;
    !! wrapping still maps the upper endpoint to the lower endpoint.
    pure function reverse(this) result(patch)
        class(nurbs_curve), intent(in) :: this
        type(nurbs_curve) :: patch

        call prepare_result(this,patch,"reverse")
        if (patch%err%ok) call orient(this,[-1],patch)
    end function reverse

    !> New axis d takes source axis abs(axes(d)); a negative entry reverses its original interval.
    pure function reparameterize_axes(this,axes) result(patch)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: axes(:)
        type(nurbs_curve) :: patch

        call prepare_result(this,patch,"reparameterize")
        if (patch%err%ok) call orient(this,axes,patch)
    end function reparameterize_axes

    pure subroutine orient(this,axes,patch)
        class(nurbs_curve), intent(in) :: this
        integer, intent(in) :: axes(:)
        type(nurbs_curve), intent(inout) :: patch

        call reorient_patch([this%nc],[this%degree],axes,this%knot,this%Xc,this%Wc,patch%knot,patch%Xc,patch%Wc,patch%err)
        if (.not.patch%err%ok) return
        patch%degree=this%degree
        patch%nc=this%nc
        patch%wrap_parameters=this%wrap_parameters
        patch%rational=this%rational
    end subroutine orient


    !> Affinely map the complete knot vector to the finite active interval [a,b].
    !! Controls, explicit weights, degree and wrapping are retained without samples.
    !! A kth derivative scales by ((old_hi-old_lo)/(b-a))**k. An unrepresentable
    !! knot or a collapsed distinct knot span returns a diagnostic in the result.
    pure function reparameterize_interval(this,a,b) result(patch)
        class(nurbs_curve), intent(in) :: this
        real(rk), intent(in) :: a,b
        type(nurbs_curve) :: patch
        real(rk), allocatable :: knots(:)

        call prepare_result(this,patch,"reparameterize")
        if (.not.patch%err%ok) return
        call reparameterize_knots(this%knot,this%nc,this%degree,a,b,knots,patch%err, &
            category="forcad_nurbs_curve")
        if (patch%err%ok) call patch%set(knot=knots,Xc=this%Xc,Wc=this%Wc,degree=this%degree, &
            wrap_parameters=this%wrap_parameters)
    end function reparameterize_interval


    !> Validate a geometry-only result without copying source geometry or caches.
    pure subroutine prepare_result(this,patch,location)
        class(nurbs_curve), intent(in) :: this
        type(nurbs_curve), intent(out) :: patch
        character(*), intent(in) :: location

        patch%err=this%err
        if (.not.patch%err%ok) return
        if (.not.this%is_valid()) call patch%err%set(code=102,severity=1,category="forcad_nurbs_curve",location=location, &
            message="Complete valid curve geometry is required.", &
            suggestion="Initialize knots, controls and any weights before this operation.")
    end subroutine prepare_result


    !> Decode the packed workspace and commit through the validated constructor.
    pure subroutine commit_refinement_controls(this, knots, controls, weight_exponent)
        class(nurbs_curve), intent(inout) :: this
        real(rk), contiguous, intent(in) :: knots(:)
        real(rk), allocatable, intent(inout) :: controls(:,:)
        integer, intent(in) :: weight_exponent
        integer :: d, ncoord

        ! A curve needs no tensor permutation; decode its workspace in place.
        ncoord = size(this%Xc,2)
        if (allocated(this%Wc)) then
            do d = 1, ncoord
                controls(:,d) = controls(:,d)/controls(:,ncoord+1)
            end do
            controls(:,ncoord+1) = scale(controls(:,ncoord+1),weight_exponent)
            call this%set(knot=knots,Xc=controls(:,:ncoord),Wc=controls(:,ncoord+1), &
                wrap_parameters=this%wrap_parameters)
        else
            call this%set(knot=knots,Xc=controls,wrap_parameters=this%wrap_parameters)
        end if
    end subroutine commit_refinement_controls


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
        real(rk), allocatable :: knots(:), Xc(:,:), Wc(:)

        if (.not.this%err%ok) return
        if (.not.this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_curve",location="insert_knots", &
                message="Refinement requires valid geometry and nonnegative counts.", &
                suggestion="Set complete spline geometry and check the refinement request.")
            return
        end if
        call refine_control_net(this%knot,this%degree,this%Xc,this%Wc,[this%nc],1, &
            knots,Xc,Wc,this%err,B,Bs,Xth=Xth,r=r)
        if (.not.this%err%ok .or. .not.allocated(Xc)) return
        call this%set(knot=knots,Xc=Xc,Wc=Wc,wrap_parameters=this%get_parameter_wrapping())
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
        real(rk), allocatable :: knots(:), Xc(:,:), Wc(:)

        if (.not.this%err%ok) return
        if (.not.this%is_valid() .or. t<0) then
            call this%err%set(code=merge(100,103,t<0),severity=1,category="forcad_nurbs_curve",location="elevate_degree", &
                message="Refinement requires valid geometry and nonnegative counts.", &
                suggestion="Set complete spline geometry and check the refinement request.")
            return
        end if
        call refine_control_net(this%knot,this%degree,this%Xc,this%Wc,[this%nc],1, &
            knots,Xc,Wc,this%err,B,Bs,t=t)
        if (.not.this%err%ok .or. .not.allocated(Xc)) return
        call this%set(knot=knots,Xc=Xc,Wc=Wc,wrap_parameters=this%get_parameter_wrapping())
    end subroutine
    !===============================================================================


    !> Reduce degree under one Euclidean error budget for the complete call.
    !! Each candidate is compared with the call-entry geometry. Return the longest
    !! accepted prefix in reduced; over-budget proposals are normal rejections.
    !! Commit once, clearing derived caches only when at least one step succeeds.
    pure subroutine reduce_degree(this, t, tol, reduced, error_bound)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in) :: t
            !! Requested decrement, 0 <= t <= the selected degree.
        real(rk), intent(in) :: tol
            !! Finite nonnegative physical error budget, including every accepted step.
        integer, optional, intent(out) :: reduced
            !! Number of committed decrements; zero on failure or complete rejection.
        real(rk), optional, intent(out) :: error_bound
            !! Conservative whole-call bound for the committed geometry.
        real(rk), allocatable :: knots(:),controls(:,:)
        real(rk) :: bound
        integer :: accepted,weight_exponent

        if (present(reduced)) reduced=0
        if (present(error_bound)) error_bound=0.0_rk
        if (.not.this%err%ok) return
        if (.not.this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_curve",location="reduce_degree", &
                message="Degree reduction requires valid geometry.",suggestion="Initialize valid knots, controls and weights.")
            return
        end if
        knots=this%knot
        call reduce_control_lines(knots,this%degree,this%Xc,this%Wc,[this%nc],1,t,tol, &
            this%err,controls,weight_exponent,accepted,bound)
        if (.not.this%err%ok .or. accepted==0) return
        call commit_refinement_controls(this,knots,controls,weight_exponent)
        if (.not.this%err%ok) return
        if (present(reduced)) reduced=accepted
        if (present(error_bound)) error_bound=bound
    end subroutine reduce_degree

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Validate the state required by all curve derivative evaluators.
    pure subroutine check_derivative_state(this, location, valid, side)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: side
        character(len=*), intent(in) :: location
        logical, intent(out) :: valid

        valid = .false.
        if (.not. this%err%ok) return
        if (.not. valid_trace_side(side)) then
            call this%err%set(code=100,severity=1,category="forcad_nurbs_curve",location=location, &
                message="Invalid knot trace selector.",suggestion="Use -1, 0, or 1 for each parameter direction.")
            return
        end if
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


    !> Select finite basis-query parameters without changing sampled geometry.
    pure subroutine query_parameters(this, location, parameters, valid, res, Xt, side)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: side
        character(len=*), intent(in) :: location
        real(rk), allocatable, intent(out) :: parameters(:)
        logical, intent(out) :: valid
        integer, intent(in), optional :: res
        real(rk), intent(in), optional :: Xt(:)

        call check_derivative_state(this,location,valid,side)
        if (.not. valid) then
            allocate(parameters(0))
            return
        end if
        call query_parameter_axis(this%knot, this%nc, this%degree, this%Xt, parameters, valid, res, Xt)
        if (.not. valid) then
            if (allocated(parameters)) deallocate(parameters)
            allocate(parameters(0))
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Curve queries require nonempty finite parameter values.",&
                location   = location,&
                suggestion = "Supply finite Xt, res >= 1, or create a parameter cache first.")
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> Validate a single public parameter before entering evaluation kernels.
    pure subroutine check_query_parameter(this, location, Xt, valid, side)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: side
        character(len=*), intent(in) :: location
        real(rk), intent(in) :: Xt
        logical, intent(out) :: valid

        call check_derivative_state(this,location,valid,side)
        if (.not. valid) return
        valid = ieee_is_finite(Xt)
        if (.not. valid) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Curve query parameters must be finite.",&
                location   = location,&
                suggestion = "Use a finite parameter value.")
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure integer function curve_trace(side) result(trace)
        integer, intent(in), optional :: side
        trace=0
        if (present(side)) trace=side
    end function curve_trace

    !> Evaluate one arbitrary-order curve basis derivative without global work arrays.
    pure subroutine evaluate_basis_order(this,Xt,order,values,Wc,side,first_active,elem)
        class(nurbs_curve), intent(in) :: this
        real(rk), intent(in) :: Xt
        integer, intent(in) :: order
        real(rk), intent(out) :: values(:)
        real(rk), contiguous, intent(in), optional :: Wc(:)
        integer, intent(in) :: side
        integer, intent(out), optional :: first_active
        integer, contiguous, intent(in), optional :: elem(:)
        integer :: first(1)

        if (present(first_active)) then
            call tensor_basis_derivative_at([Xt],[this%nc],[this%degree],this%knot,[order],[this%wrap_parameters], &
                values,first_active=first,Wc=Wc,side=[side])
            first_active=first(1)
        else
            call tensor_basis_derivative_at([Xt],[this%nc],[this%degree],this%knot,[order],[this%wrap_parameters], &
                values,elem=elem,Wc=Wc,side=[side])
        end if
    end subroutine evaluate_basis_order
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
    pure subroutine derivative_order_active_vector(this, Xt, order, first_active, Dgc, side)
        class(nurbs_curve), target, intent(inout) :: this
        integer, intent(in), optional :: side !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in), contiguous :: Xt(:)
            !! Finite parameters `[npoint]`.
        integer, intent(in) :: order
            !! Nonnegative derivative order (k); zero returns basis values.
        integer, allocatable, intent(out) :: first_active(:)
            !! First active control-point index for every parameter.
        real(rk), allocatable, intent(out) :: Dgc(:,:)
            !! Local derivatives `[degree+1,npoint]`.
        integer :: point, trace
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call check_derivative_state(this, "derivative_order_active_vector", valid, side)
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
        trace=curve_trace(side)
        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
#if defined(__NVCOMPILER)
        do point = 1, size(Xt)
#else
        do concurrent (point = 1:size(Xt))
#endif
            call evaluate_basis_order(this,Xt(point),order,Dgc(:,point),weights,trace, &
                first_active=first_active(point))
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one active arbitrary-order curve derivative.
    !! Entry `Dgc(a+1)` is the derivative of basis
    !! `first_active+a`, for `a=0:degree`. Wrapping and outside-domain behavior
    !! are identical to the vector overload.
    pure subroutine derivative_order_active_scalar(this, Xt, order, first_active, Dgc, side)
        class(nurbs_curve), target, intent(inout) :: this
        integer, intent(in), optional :: side !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in) :: Xt
            !! Finite parameter value.
        integer, intent(in) :: order
            !! Nonnegative derivative order.
        integer, intent(out) :: first_active
            !! First active control-point index.
        real(rk), allocatable, intent(out) :: Dgc(:)
            !! Local derivative values `[degree+1]`.
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        first_active = 0
        call check_derivative_state(this, "derivative_order_active_scalar", valid, side)
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
        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
        call evaluate_basis_order(this,Xt,order,Dgc,weights,curve_trace(side), &
            first_active=first_active)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one requested derivative order at a vector of curve parameters.
    !! The dense result satisfies
    !! `Dgc(g,i)=d**order R_i(Xt(g))/du**order`. Explicit `Xt` takes
    !! precedence over positive `res`; with neither argument the sampling
    !! parameters are reused. Queries do not change sampled geometry or its
    !! matching parameters. Invalid inputs return an empty result and diagnostic.
    pure subroutine derivative_order_vector(this, order, Dgc, res, Xt, side)
        class(nurbs_curve), target, intent(inout) :: this
        integer, intent(in), optional :: side !! Trace per direction: -1 left, 0 default, +1 right.
        integer, intent(in) :: order
            !! Nonnegative derivative order; zero returns the dense basis.
        real(rk), allocatable, intent(out) :: Dgc(:,:)
            !! Dense derivative matrix `[ng,nc]`.
        integer, intent(in), optional :: res
            !! Optional positive generated parameter count.
        real(rk), intent(in), contiguous, optional :: Xt(:)
            !! Optional explicit parameter vector; takes precedence over `res`.
        integer :: point, trace
        real(rk), allocatable :: parameters(:)
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call check_derivative_state(this, "derivative_order_vector", valid, side)
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

        call query_parameters(this, "derivative_order_vector", parameters, valid, res, Xt)
        if (.not. valid) then
            allocate(Dgc(0,0))
            return
        end if

        allocate(Dgc(size(parameters),this%nc), source=0.0_rk)
        trace=curve_trace(side)
        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
#if defined(__NVCOMPILER)
        ! NVFORTRAN 26.3 cannot lower the pure local kernel for stdpar GPU execution.
        do point = 1, size(parameters)
#else
        do concurrent (point = 1:size(parameters))
#endif
            call evaluate_basis_order(this,parameters(point),order,Dgc(point,:),weights,trace)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one requested derivative order at one curve parameter.
    !! Without `elem`, `Dgc(i)=d**order R_i(Xt)/du**order` for all controls.
    !! With `elem`, entries follow the supplied one-based control-index order;
    !! indices outside `1:nc` return zero. Nonfinite parameters return an
    !! empty result and diagnostic.
    pure subroutine derivative_order_scalar(this, Xt, order, Dgc, elem, side)
        class(nurbs_curve), target, intent(inout) :: this
        integer, intent(in), optional :: side !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in) :: Xt
            !! Parameter value, mapped according to the wrapping policy.
        integer, intent(in) :: order
            !! Nonnegative derivative order.
        real(rk), allocatable, intent(out) :: Dgc(:)
            !! Dense `[nc]` or selected `[size(elem)]` derivative vector.
        integer, intent(in), contiguous, optional :: elem(:)
            !! Optional ordered one-based control indices.
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call check_query_parameter(this, "derivative_order_scalar", Xt, valid, side)
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
        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
        call evaluate_basis_order(this,Xt,order,Dgc,weights,curve_trace(side),elem=elem)
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
    !! parameters and columns to all control variables. The sampled geometry
    !! and its matching parameters remain unchanged. Inputs must be finite;
    !! invalid state or parameters return empty outputs and a diagnostic.
    pure subroutine derivative_vector(this, res, Xt, dTgc, Tgc, side)
        class(nurbs_curve), target, intent(inout) :: this
        integer, intent(in), optional :: side
        integer, intent(in), optional :: res
            !! Optional positive uniform-grid size.
        real(rk), intent(in), contiguous, optional :: Xt(:)
            !! Optional parameter vector `[ng]`.
        real(rk), allocatable, intent(out) :: dTgc(:,:)
            !! First derivatives \(\partial R_i/\partial u\), shape `[ng,nc]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
            !! Optional basis values computed on the same grid, shape `[ng,nc]`.
        call query_basis_grid(this, "derivative_vector", res, Xt, side, Tgc=Tgc, dTgc=dTgc)
    end subroutine derivative_vector
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate first parametric derivatives at one parameter.
    !! Optional `elem` selects and orders a local subset of control variables.
    !! Invalid spline state or a nonfinite parameter returns empty outputs
    !! and a recoverable diagnostic.
    pure subroutine derivative_scalar(this, Xt, dTgc, Tgc, elem, side)
        class(nurbs_curve), target, intent(inout) :: this
        integer, intent(in), optional :: side
        real(rk), intent(in) :: Xt
            !! Parameter value, mapped according to the wrapping policy.
        integer, intent(in), optional :: elem(:)
            !! Optional one-based control indices.
        real(rk), allocatable, intent(out) :: dTgc(:)
            !! Selected derivatives \(\partial R_i/\partial u\).
        real(rk), allocatable, intent(out), optional :: Tgc(:)
            !! Optional selected basis values.
        call query_basis_point(this, "derivative_scalar", Xt, elem, side, Tgc=Tgc, dTgc=dTgc)
    end subroutine derivative_scalar
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate second parametric derivatives at a vector of parameters.
    !! Optional lower derivatives and values are returned from the same
    !! rational or polynomial evaluation path. Parameter selection and
    !! validation are the same as for the vector [[nurbs_curve:derivative]]
    !! overload; sampled geometry and its matching parameters remain unchanged.
    pure subroutine derivative2_vector(this, res, Xt, d2Tgc, dTgc, Tgc, side)
        class(nurbs_curve), target, intent(inout) :: this
        integer, intent(in), optional :: side
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
        call query_basis_grid(this, "derivative2_vector", res, Xt, side, Tgc=Tgc, dTgc=dTgc, d2Tgc=d2Tgc)
    end subroutine derivative2_vector
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate second parametric derivatives at one parameter.
    !! Optional elem selects global control indices without changing rational
    !! normalization. Repeated indices retain their order; invalid indices give
    !! zero entries. Invalid spline state or a nonfinite parameter returns empty outputs
    !! and a recoverable diagnostic.
    pure subroutine derivative2_scalar(this, Xt, d2Tgc, dTgc, Tgc, elem, side)
        class(nurbs_curve), target, intent(inout) :: this
        integer, intent(in), optional :: side
        real(rk), intent(in) :: Xt
            !! Parameter value, mapped according to the wrapping policy.
        real(rk), allocatable, intent(out) :: d2Tgc(:)
            !! Dense `[nc]` or selected `[size(elem)]` second derivatives.
        real(rk), allocatable, intent(out), optional :: dTgc(:)
            !! Optional first derivatives in the same output order.
        real(rk), allocatable, intent(out), optional :: Tgc(:)
            !! Optional basis values in the same output order.
        integer, intent(in), optional :: elem(:)
            !! Optional ordered global indices; outputs have length size(elem).
        call query_basis_point(this, "derivative2_scalar", Xt, elem, side, Tgc=Tgc, dTgc=dTgc, d2Tgc=d2Tgc)
    end subroutine derivative2_scalar
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the dense basis at a vector of parameters.
    !!
    !! Explicit `Xt` takes precedence over `res`; otherwise `res` creates an
    !! inclusive uniform grid. With neither argument, the current parameter
    !! cache is reused. The sampled geometry and its matching parameters
    !! remain unchanged. Invalid state, nonfinite or empty parameters, and
    !! nonpositive resolution return an empty result and diagnostic.
    pure subroutine basis_vector(this, res, Xt, Tgc, side)
        class(nurbs_curve), target, intent(inout) :: this
        integer, intent(in), optional :: side
        integer, intent(in), optional :: res
            !! Optional positive uniform-grid size.
        real(rk), intent(in), contiguous, optional :: Xt(:)
            !! Optional parameter vector `[ng]`.
        real(rk), allocatable, intent(out) :: Tgc(:,:)
            !! Basis matrix `[ng,nc]`; each valid row sums to one.
        call query_basis_grid(this, "basis_vector", res, Xt, side, Tgc=Tgc)
    end subroutine basis_vector

    !> Evaluate the requested basis values, gradients, and Hessians on a grid.
    pure subroutine query_basis_grid(this, location, res, Xt, side, Tgc, dTgc, d2Tgc)
        class(nurbs_curve), target, intent(inout) :: this
        character(*), intent(in) :: location
        integer, intent(in), optional :: res
        real(rk), intent(in), contiguous, optional :: Xt(:)
        integer, intent(in), optional :: side
        real(rk), allocatable, intent(out), optional :: Tgc(:,:), dTgc(:,:), d2Tgc(:,:)
        real(rk), allocatable :: parameters(:)
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call query_parameters(this, location, parameters, valid, res, Xt,side=side)
        if (.not. valid) then
            if (present(d2Tgc)) allocate(d2Tgc(0,0))
            if (present(dTgc)) allocate(dTgc(0,0))
            if (present(Tgc)) allocate(Tgc(0,0))
            return
        end if

        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
        call compute_basis_data(parameters, this%knot, this%degree, this%nc, Tgc=Tgc, dTgc=dTgc, d2Tgc=d2Tgc, &
            wrap_parameters=this%wrap_parameters, Wc=weights,side=side)
    end subroutine query_basis_grid
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the basis at one parameter.
    !! Optional `elem` selects and orders a local subset of control variables.
    !! Invalid spline state or a nonfinite parameter returns empty outputs
    !! and a recoverable diagnostic.
    pure subroutine basis_scalar(this, Xt, Tgc, elem, side)
        class(nurbs_curve), target, intent(inout) :: this
        integer, intent(in), optional :: side
        real(rk), intent(in) :: Xt
            !! Parameter value, mapped according to the wrapping policy.
        integer, intent(in), optional :: elem(:)
            !! Optional one-based control indices.
        real(rk), allocatable, intent(out) :: Tgc(:)
            !! Dense `[nc]` or selected `[size(elem)]` basis vector.

        call query_basis_point(this, "basis_scalar", Xt, elem, side, Tgc=Tgc)
    end subroutine basis_scalar

    !> Evaluate the requested basis values, gradients, and Hessians at one point.
    pure subroutine query_basis_point(this, location, Xt, elem, side, Tgc, dTgc, d2Tgc)
        class(nurbs_curve), target, intent(inout) :: this
        character(*), intent(in) :: location
        real(rk), intent(in) :: Xt
        integer, intent(in), optional :: elem(:)
        integer, intent(in), optional :: side
        real(rk), allocatable, intent(out), optional :: Tgc(:), dTgc(:), d2Tgc(:)
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call check_query_parameter(this, location, Xt, valid,side=side)
        if (.not. valid) then
            if (present(d2Tgc)) allocate(d2Tgc(0))
            if (present(dTgc)) allocate(dTgc(0))
            if (present(Tgc)) allocate(Tgc(0))
            return
        end if

        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
        call compute_basis_data(Xt, this%knot, this%degree, this%nc, Tgc=Tgc, dTgc=dTgc, d2Tgc=d2Tgc, &
            elem=elem, wrap_parameters=this%wrap_parameters, Wc=weights,side=side)
    end subroutine query_basis_point
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report whether nonuniform explicit weights change the basis.
    !! Exactly equal weights select the polynomial basis. Every nonconstant
    !! positive weight list selects rational evaluation; no tolerance discards
    !! weight variation. An active diagnostic or absent weights returns false.
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
    !! `periodic_topology`. Every other valid curve is `"bounded"`; an
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

        call get_array(this%elemConn_Xc_vis,elemConn,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return a copy of sampled-curve visualization connectivity.
    pure function get_elem_Xg_vis(this) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        call get_array(this%elemConn_Xg_vis,elemConn,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return a copy of cached active-span IGA element connectivity.
    !! Returns `[0,0]` when no connectivity has been cached; `cmp_elem()`
    !! computes connectivity independently and does not populate this cache.
    pure function get_elem(this) result(elemConn)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        call get_array(this%elemConn,elemConn,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Remove interior knots with optional physical geometry-error control.
    !! Requests use exact stored knots and are processed in array order. Each
    !! request accepts its longest removable prefix; a rejected candidate is
    !! not an error. Without `tol`, the original roundoff-scaled homogeneous
    !! control test applies. With `tol`, every candidate is compared with the
    !! call-entry geometry in a common basis, bounding Euclidean error over the
    !! complete active domain. The bound includes all accepted requests in this
    !! call; it is conservative, subject to floating-point arithmetic, and is
    !! not a sampled or measured maximum. Periodic removals preserve complete
    !! cyclic knot orbits. Invalid requests leave geometry unchanged and return
    !! zero counts and bound.
    !!
    !! **Algorithm:** Piegl and Tiller, Algorithm A5.8; W. Tiller,
    !! [doi:10.1016/0010-4485(92)90012-Y](https://doi.org/10.1016/0010-4485(92)90012-Y).
    pure subroutine remove_knots(this, Xth, r, removed, tol, error_bound)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xth(:)
            !! Exact stored knot values strictly inside the active interval.
        integer, intent(in), contiguous :: r(:)
            !! Nonnegative requested counts, not exceeding current multiplicity.
        integer, allocatable, intent(out), optional :: removed(:)
            !! Committed counts per request; zero on failure or rejection.
        real(rk), intent(in), optional :: tol
            !! Finite nonnegative bound in physical units for this complete call.
        real(rk), intent(out), optional :: error_bound
            !! Conservative Euclidean bound for the committed update; zero on failure.
        integer, allocatable :: counts(:)
        real(rk), allocatable :: controls(:,:), knots(:)
        integer :: weight_exponent

        if (present(removed)) allocate(removed(size(r)),source=0)
        if (present(error_bound)) error_bound=0.0_rk
        if (.not.this%err%ok) return
        if (.not.this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_curve",location="remove_knots", &
                message="Knot removal requires valid geometry.",suggestion="Initialize valid controls, knots and weights.")
            return
        end if
        knots=this%knot
        call coarsen_control_lines(knots,this%degree,this%Xc,this%Wc,[this%nc],1,Xth,r, &
            this%err,controls,weight_exponent,counts,tol,error_bound)
        if (.not.this%err%ok .or. .not.any(counts>0)) return
        call commit_refinement_controls(this,knots,controls,weight_exponent)
        if (.not.this%err%ok) then
            if (present(error_bound)) error_bound=0.0_rk
            return
        end if
        if (present(removed)) removed=counts
    end subroutine remove_knots
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

        if (.not.this%err%ok .or. &
            .not.allocated(this%knot)) then
            allocate(elemConn(0,0))
        else
            elemConn=spline_connectivity([this%nc],[this%degree],this%knot)
        end if
    end function cmp_elem
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
    !! Unavailable or empty samples, an empty query, or an active diagnostic
    !! return zero outputs and sample identifier zero. Nonfinite metric inputs
    !! likewise return zero; distances are compared without squared-distance overflow.
    !! After modify_Xc/Wc,
    !! call create() to refresh the retained parameter grid before searching.
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
        integer :: id_, ncopy, npoint

        if (present(nearest_Xg)) nearest_Xg = 0.0_rk
        if (present(nearest_Xt)) nearest_Xt = 0.0_rk
        if (present(id)) id = 0
        if (.not.this%err%ok .or. .not.allocated(this%Xg) .or. .not.allocated(this%Xt)) return
        npoint=min(size(this%Xg,1),size(this%Xt))
        id_=nearest_sample_index(this%Xg(:npoint,:),point_Xg)
        if (id_==0) return
        if (present(id)) id=id_
        if (present(nearest_Xt)) nearest_Xt=this%Xt(id_)
        if (present(nearest_Xg)) then
            ncopy=min(size(point_Xg),size(this%Xg,2))
            nearest_Xg(:ncopy)=this%Xg(id_,:ncopy)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute a bounded continuous nearest-point candidate without sampling the curve.
    !! Active spans supply streamed seeds, including narrow and degree-zero pieces.
    !! Local basis jets drive safeguarded Newton, Gauss-Newton and descent steps.
    !! This numerical candidate has no global closest-point certificate.
    !!
    !! Optional initial_Xt supplies a finite warm start inside the active interval.
    !! Query coordinates must cover the physical dimension; extras are ignored.
    !! An optional nearest_Xg must have exactly the physical coordinate count.
    !! Tolerance must be finite and nonnegative, and maxit nonnegative; maxit=0
    !! retains the seeded candidate with a nonconvergence diagnostic. Validation
    !! failures zero outputs; iteration failures retain the best bounded candidate.
    !! Optional err resets and bypasses the object's sticky diagnostic. Geometry,
    !! sampled parameters, points and connectivity remain unchanged.
    pure subroutine nearest_point2(this, point_Xg, tol, maxit, nearest_Xt, nearest_Xg, initial_Xt, err)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: point_Xg(:), tol
        integer, intent(in) :: maxit
        real(rk), intent(out) :: nearest_Xt
        real(rk), intent(out), optional :: nearest_Xg(:)
        real(rk), intent(in), optional :: initial_Xt
        type(debug), intent(out), optional :: err
        real(rk), allocatable :: initial(:)
        real(rk) :: parameters(1)
        type(debug) :: status

        nearest_Xt=0.0_rk
        if (present(nearest_Xg)) nearest_Xg=0.0_rk
        if (present(err)) then
            call err%reset()
        else
            if (.not.this%err%ok) return
        end if
        parameters=0.0_rk
        if (.not.allocated(this%Xc) .or. .not.allocated(this%knot)) then
            call status%set(code=103,severity=1,category="forcad_nurbs_curve",location="nearest_point2", &
                message="Curve geometry is not initialized.", &
                suggestion="Set knots and controls before nearest-point projection.")
        else
            if (present(initial_Xt)) initial=[initial_Xt]
            call nearest_spline_point(this%Xc,this%Wc,[this%nc],[this%degree],this%knot,point_Xg,tol,maxit, &
                parameters,nearest_Xg,status,initial_Xt=initial,wrap=[this%wrap_parameters])
        end if
        nearest_Xt=parameters(1)
        if (present(err)) then
            err=status
        else
            this%err=status
        end if
    end subroutine nearest_point2
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
    !! has no canonical signed orientation without an external reference. At a zero
    !! tangent either measure-sensitivity output also requires diagnostic 108.
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
    pure subroutine ansatz(this, ie, ig, Tgc, dTgc_dXg, dL, ngauss, d2Tgc_dXg2, strict, err, Xg, dXg_dXt, &
        dmeasure_dXc, dmeasure_dlogWc, dgradient_dXc, dgradient_dlogWc, djacobian_dXc, djacobian_dlogWc)
        class(nurbs_curve), target, intent(inout) :: this
        integer, intent(in) :: ie, ig
            !! One-based active element and quadrature-point indices.
        real(rk), allocatable, intent(out) :: Tgc(:), dTgc_dXg(:,:)
            !! Local basis `[degree+1]` and physical gradients `[degree+1,ncoord]`.
        real(rk), intent(out) :: dL
            !! Complete quadrature-weighted differential length.
        integer, intent(in), optional :: ngauss
            !! Positive quadrature count; default degree+1.
        real(rk), allocatable, target, intent(out), optional :: d2Tgc_dXg2(:)
            !! Local second arc-length derivatives.
        logical, intent(in), optional :: strict
            !! Reject singular tangents and reversed one-coordinate mappings.
        type(debug), target, intent(out), optional :: err
            !! Caller-owned diagnostic; leaves the object diagnostic unchanged.
        real(rk), allocatable, target, intent(out), optional :: Xg(:), dXg_dXt(:,:)
            !! Physical point `[ncoord]` and raw parametric Jacobian `[ncoord,1]`.
        real(rk), allocatable, target, intent(out), optional :: dmeasure_dXc(:,:), dmeasure_dlogWc(:)
            !! Derivatives of weighted length with respect to active Euclidean controls and log weights.
        real(rk), allocatable, target, intent(out), optional :: dgradient_dXc(:,:,:,:), dgradient_dlogWc(:,:,:)
            !! Physical gradient derivatives indexed by basis, component, active control, and (for coordinates) design component.
        real(rk), allocatable, intent(out), optional :: djacobian_dXc(:,:), djacobian_dlogWc(:,:,:)
            !! Raw Jacobian derivatives for controls (A,alpha) and log weights (i,alpha,A).
        integer :: counts(1)
        counts=this%degree+1
        if (present(ngauss)) counts=ngauss
        call tensor_ansatz(this%Xc,this%Wc,this%rational,[this%degree],[this%nc],this%knot, &
            object_err=this%err,err=err,ie=ie,ig=ig,basis=Tgc,gradient=dTgc_dXg,measure=dL,ngauss=counts,strict=strict, &
            Xg=Xg,dXg_dXt=dXg_dXt,curve_second=d2Tgc_dXg2,dmeasure_dXc=dmeasure_dXc,dmeasure_dlogWc=dmeasure_dlogWc, &
            dgradient_dXc=dgradient_dXc,dgradient_dlogWc=dgradient_dlogWc, &
            djacobian_dXc=djacobian_dXc,djacobian_dlogWc=djacobian_dlogWc)
    end subroutine ansatz


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
    !! Tangents contract the same local, weight-scaled basis derivatives used
    !! by the IGA queries, preserving invariance under a common weight scale.
    !!
    !! Geometry and quadrature counts are validated before integration. Strict
    !! mode uses the tangent and orientation criterion of ansatz. Caller-owned
    !! err bypasses the sticky object diagnostic; valid quadrature contributions
    !! remain in length if a numerical or strict-mode check fails.
    pure subroutine cmp_length(this, length, ngauss, strict, err)
        class(nurbs_curve), target, intent(inout) :: this
        real(rk), intent(out) :: length
            !! Sum of valid nonnegative quadrature contributions, including on failure.
        integer, intent(in), optional :: ngauss
            !! Positive quadrature count per span; default degree+1.
        logical, intent(in), optional :: strict
            !! Reject singular tangents and inverted one-coordinate mappings.
        type(debug), target, intent(out), optional :: err
            !! Caller-owned diagnostic; resets and leaves the object diagnostic unchanged.
        type(debug), pointer :: diagnostic
        real(rk), allocatable :: Xksi(:), Wksi(:)
        real(rk), pointer, contiguous :: weights(:)
        real(rk) :: basis(0:this%degree), gradient(0:this%degree,1), unused_hessian(0,1,1)
        real(rk) :: tangent(curve_ncoord(this)), parameter, width, speed, contribution, partial, cell_sum
        integer :: span, point, first, nq, invalid_points
        logical :: ready, strict_, valid

        length = 0.0_rk
        if (present(err)) then
            diagnostic => err
            call diagnostic%reset()
        else
            diagnostic => this%err
        end if
        if (.not. diagnostic%ok) return
        nq = this%degree+1
        if (present(ngauss)) nq = ngauss
        strict_ = .false.
        if (present(strict)) strict_ = strict
        ready = .false.
        if (allocated(this%Xc) .and. allocated(this%knot)) then
            ready = valid_knot_vector(this%knot,this%nc,this%degree) .and. size(this%Xc,1)==this%nc .and. &
                size(this%Xc,2)>0 .and. all(ieee_is_finite(this%Xc))
            if (allocated(this%Wc)) ready = ready .and. size(this%Wc)==this%nc .and. valid_weights(this%Wc)
        end if
        if (.not. ready .or. nq<1) then
            call diagnostic%set(code=merge(100,103,ready),severity=1,category="forcad_nurbs_curve", &
                message="Length integration requires valid geometry and a positive quadrature count.",location="cmp_length", &
                suggestion="Set finite controls, valid knots and positive weights; use ngauss >= 1.")
            return
        end if
        call gauss_leg([0.0_rk,1.0_rk],nq-1,Xksi,Wksi)
        nullify(weights)
        if (this%rational .and. allocated(this%Wc)) weights => this%Wc
        partial = 0.0_rk
        invalid_points = 0
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do span = this%degree+1, this%nc
#else
        do concurrent (span = this%degree+1:this%nc) &
            local(basis,gradient,unused_hessian,tangent,parameter,width,speed,contribution,valid,cell_sum,point,first) &
            reduce(+:length) reduce(max:invalid_points) reduce(max:partial)
#endif
            width = this%knot(span+1)-this%knot(span)
            if (width<=0.0_rk) cycle
            cell_sum = 0.0_rk
            do point = 1, nq
                parameter = this%knot(span)+Xksi(point)*width
                if (.not. ieee_is_finite(parameter) .or. parameter<=this%knot(span) .or. &
                    parameter>=this%knot(span+1)) then
                    invalid_points = max(invalid_points,107)
                    cycle
                end if
                call compute_basis_jet(parameter,this%knot,this%degree,this%nc,1,first,basis,gradient,unused_hessian,Wc=weights)
                call curve_frame(this%Xc(first:first+this%degree,:),gradient,width,Wksi(point),strict_, &
                    tangent,speed,contribution,valid)
                if (.not.valid) then
                    invalid_points = max(invalid_points,108)
                    cycle
                end if
                if (cell_sum>huge(cell_sum)-contribution) then
                    invalid_points = max(invalid_points,108)
                    cycle
                end if
                cell_sum = cell_sum+contribution
            end do
            length = length+cell_sum
            partial = max(partial,cell_sum)
        end do
        if (.not. ieee_is_finite(length)) then
            length = partial
            invalid_points = max(invalid_points,108)
        end if
        if (invalid_points>0) call diagnostic%set(code=invalid_points,severity=1,category="forcad_nurbs_curve", &
            message="Some curve quadrature points have unrepresentable, singular or inverted mappings.",location="cmp_length", &
            suggestion="Use regular geometry and a scale that represents its derivatives and measure.")
    end subroutine cmp_length
    !===============================================================================




    !> Evaluate one local basis jet, computing only the requested derivative order.
    pure subroutine compute_basis_jet(Xt, knot, degree, nc, order, first, basis, gradient, hessian, wrap_parameters, Wc, side)
        real(rk), intent(in) :: Xt
        integer, intent(in), optional :: side
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree, nc, order
        integer, intent(out) :: first
        real(rk), intent(out), contiguous :: basis(0:), gradient(0:,:)
        real(rk), intent(out), contiguous, optional :: hessian(0:,:,:)
        logical, intent(in), optional :: wrap_parameters
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk) :: ders(0:order,0:degree), identity(0:2,0:0)

        call basis_bspline_der_all_active(Xt,knot,nc,degree,order,first,ders,side,wrap_parameters)
        identity = 0.0_rk
        identity(0,0) = 1.0_rk
        select case (order)
        case (0)
            basis=ders(0,:)
            if (.not.present(Wc)) return
            if (size(Wc)<nc .or. first+degree>size(Wc)) then
                basis=0.0_rk
            else
                call normalize_rational_basis(basis,Wc(first:first+degree))
            end if
        case (1)
            call tensor_basis_derivatives2_local(first,nc,ders,1,1,identity,1,1,identity,basis,gradient,Wc=Wc)
        case default
            call tensor_basis_derivatives2_local(first,nc,ders,1,1,identity,1,1,identity,basis,gradient,hessian,Wc)
        end select
    end subroutine compute_basis_jet

    !> Scatter the requested local basis data directly into dense parameter rows.
    pure subroutine compute_basis_data_vector(Xt, knot, degree, nc, Tgc, dTgc, d2Tgc, ng, wrap_parameters, Wc, side)
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: side
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree, nc
        real(rk), allocatable, intent(out), optional :: Tgc(:,:), dTgc(:,:), d2Tgc(:,:)
        integer, intent(in), optional :: ng
        logical, intent(in), optional :: wrap_parameters
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk) :: basis(0:degree), gradient(0:degree,1), hessian(0:degree,1,1)
        logical :: wrap
        integer :: trace, point, first, local, global, order, npoint

        trace=curve_trace(side)
        wrap=.false.
        if (present(wrap_parameters)) wrap=wrap_parameters
        order = 0
        if (present(dTgc)) order = 1
        if (present(d2Tgc)) order = 2
        npoint = size(Xt)
        if (present(ng)) npoint = ng
        if (present(Tgc)) allocate(Tgc(npoint,nc), source=0.0_rk)
        if (present(dTgc)) allocate(dTgc(npoint,nc), source=0.0_rk)
        if (present(d2Tgc)) allocate(d2Tgc(npoint,nc), source=0.0_rk)
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do point = 1, npoint
#else
        do concurrent (point = 1:npoint) local(basis,gradient,hessian,first,local,global)
#endif
            call compute_basis_jet(Xt(point),knot,degree,nc,order,first,basis,gradient,hessian,wrap,Wc,trace)
            do local = max(0,1-first), min(degree,nc-first)
                global = first+local
                if (present(Tgc)) Tgc(point,global) = basis(local)
                if (present(dTgc)) dTgc(point,global) = gradient(local,1)
                if (present(d2Tgc)) d2Tgc(point,global) = hessian(local,1,1)
            end do
        end do
    end subroutine compute_basis_data_vector


    !> Select local basis data after normalization over the complete active support.
    pure subroutine compute_basis_data_scalar(Xt, knot, degree, nc, Tgc, dTgc, d2Tgc, elem, wrap_parameters, Wc, side)
        real(rk), intent(in) :: Xt
        integer, intent(in), optional :: side
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree, nc
        real(rk), allocatable, intent(out), optional :: Tgc(:), dTgc(:), d2Tgc(:)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), optional :: wrap_parameters
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk) :: basis(0:degree), gradient(0:degree,1), hessian(0:degree,1,1)
        integer :: first, local, output, noutput, nwork, global, order, target

        order = 0
        if (present(dTgc)) order = 1
        if (present(d2Tgc)) order = 2
        noutput = nc
        nwork = degree+1
        if (present(elem)) then
            noutput = size(elem)
            nwork = noutput
        end if
        if (present(Tgc)) allocate(Tgc(noutput), source=0.0_rk)
        if (present(dTgc)) allocate(dTgc(noutput), source=0.0_rk)
        if (present(d2Tgc)) allocate(d2Tgc(noutput), source=0.0_rk)
        call compute_basis_jet(Xt,knot,degree,nc,order,first,basis,gradient,hessian,wrap_parameters,Wc,side=side)
        do output = 1, nwork
            local = output-1
            global = first+local
            if (present(elem)) then
                global = elem(output)
                local = global-first
            end if
            if (global < 1 .or. global > nc .or. local < 0 .or. local > degree) cycle
            target = global
            if (present(elem)) target = output
            if (present(Tgc)) Tgc(target) = basis(local)
            if (present(dTgc)) dTgc(target) = gradient(local,1)
            if (present(d2Tgc)) d2Tgc(target) = hessian(local,1,1)
        end do
    end subroutine compute_basis_data_scalar


    !> Evaluate physical points by contracting only active basis values and controls.
    pure function compute_Xg_vector(Xt, knot, degree, nc, ng, Xc, Wc, wrap_parameters, side) result(Xg)
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: side
        real(rk), intent(in), contiguous :: knot(:), Xc(:,:)
        integer, intent(in) :: degree, nc, ng
        real(rk), intent(in), contiguous, optional :: Wc(:)
        logical, intent(in), optional :: wrap_parameters
        real(rk), allocatable :: Xg(:,:)
        logical :: wrap
        integer :: trace, point

        trace=curve_trace(side)
        wrap=.false.
        if (present(wrap_parameters)) wrap=wrap_parameters
        allocate(Xg(ng,size(Xc,2)),source=0.0_rk)
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do point=1,ng
#else
        do concurrent (point=1:ng)
#endif
            Xg(point,:)=compute_Xg_scalar(Xt(point),knot,degree,nc,Xc,Wc,wrap,trace)
        end do
    end function compute_Xg_vector


    !> Evaluate one physical point with the same local basis kernel as batch queries.
    pure function compute_Xg_scalar(Xt, knot, degree, nc, Xc, Wc, wrap_parameters, side) result(Xg)
        real(rk), intent(in) :: Xt
        integer, intent(in), optional :: side
        real(rk), intent(in), contiguous :: knot(:), Xc(:,:)
        logical, intent(in), optional :: wrap_parameters
        integer, intent(in) :: degree, nc
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk) :: Xg(size(Xc,2)), basis(0:degree), gradient(0:degree,1), hessian(0:degree,1,1)
        integer :: first, local, coordinate

        Xg = 0.0_rk
        call compute_basis_jet(Xt,knot,degree,nc,0,first,basis,gradient,hessian,wrap_parameters,Wc,side)
        do coordinate = 1, size(Xc,2)
            do local = max(0,1-first), min(degree,nc-first)
                Xg(coordinate) = Xg(coordinate)+basis(local)*Xc(first+local,coordinate)
            end do
        end do
    end function compute_Xg_scalar
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

        call compute_basis_data(Xt,knot,degree,nc,Tgc=Tgc, dTgc=dTgc,ng=ng,wrap_parameters=wrap_parameters,Wc=Wc)
    end subroutine compute_dTgc_nurbs_1d_vector


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
        logical :: local_weights

        local_weights = .false.
        if (present(elem)) local_weights = size(Wc)==size(elem) .and. size(Wc)/=nc
        if (local_weights) then
            call compute_basis_data(Xt,knot,degree,nc,Tgc=Tgc, dTgc=dTgc,elem=elem,wrap_parameters=wrap_parameters)
            call normalize_rational_basis(Tgc,Wc,dTgc)
        else
            call compute_basis_data(Xt,knot,degree,nc,Tgc=Tgc, dTgc=dTgc,elem=elem, &
                wrap_parameters=wrap_parameters,Wc=Wc)
        end if
    end subroutine compute_dTgc_nurbs_1d_scalar


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

        call compute_basis_data(Xt,knot,degree,nc,Tgc=Tgc, dTgc=dTgc,ng=ng,wrap_parameters=wrap_parameters)
    end subroutine compute_dTgc_bspline_1d_vector


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

        call compute_basis_data(Xt,knot,degree,nc,Tgc=Tgc, dTgc=dTgc,elem=elem,wrap_parameters=wrap_parameters)
    end subroutine compute_dTgc_bspline_1d_scalar


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

        call compute_basis_data(Xt,knot,degree,nc,Tgc=Tgc,ng=ng,wrap_parameters=wrap_parameters,Wc=Wc)
    end function compute_Tgc_nurbs_1d_vector


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
        logical :: local_weights

        local_weights = .false.
        if (present(elem)) local_weights = size(Wc)==size(elem) .and. size(Wc)/=nc
        if (local_weights) then
            call compute_basis_data(Xt,knot,degree,nc,Tgc=Tgc,elem=elem,wrap_parameters=wrap_parameters)
            call normalize_rational_basis(Tgc,Wc)
        else
            call compute_basis_data(Xt,knot,degree,nc,Tgc=Tgc,elem=elem, &
                wrap_parameters=wrap_parameters,Wc=Wc)
        end if
    end function compute_Tgc_nurbs_1d_scalar


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

        call compute_basis_data(Xt,knot,degree,nc,Tgc=Tgc,ng=ng,wrap_parameters=wrap_parameters)
    end function compute_Tgc_bspline_1d_vector


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

        call compute_basis_data(Xt,knot,degree,nc,Tgc=Tgc,elem=elem,wrap_parameters=wrap_parameters)
    end function compute_Tgc_bspline_1d_scalar


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Fit polynomial B-spline control points to parameterized data.
    !!
    !! The existing degree, knot vector, and control count define the trial
    !! space. The first `ndata` rows minimize the unweighted coordinatewise
    !! residual sum. Local support is accumulated directly into symmetric
    !! banded normal equations and solved by banded Cholesky. If that
    !! factorization rejects a numerically full-rank design, the original
    !! collocation matrix is solved by Householder QR with column pivoting.
    !!
    !! This routine fits a polynomial B-spline model. A curve with nonuniform
    !! weights is rejected without changing its controls; use
    !! [[nurbs_curve:lsq_fit_nurbs]] for a rational model. Uniform weights are
    !! harmless because they produce the same polynomial basis. `Xt` and
    !! `Xdata` must contain at least `ndata` finite rows, `Xdata` must have at
    !! least one coordinate, and the parameters must give full-rank active-
    !! domain coverage. These array and finiteness preconditions are not
    !! completely diagnosed. Only `Xc` is replaced; cached `Xg` is not
    !! recomputed.
    !!
    !! **Numerical method.** The banded normal equations are the low-memory
    !! primary path. `solve_least_squares_qr` is allocated only after their
    !! failure and rejects numerically rank-deficient data with code 108.
    pure subroutine lsq_fit_bspline(this, Xt, Xdata, ndata)
        use forcad_utils, only: solve_least_squares_qr, solve_spd_banded
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:)
            !! Data parameters with at least `ndata` entries.
        real(rk), intent(in), contiguous :: Xdata(:,:)
            !! Data coordinates `[at least ndata,ncoord]`.
        integer, intent(in) :: ndata
            !! Number of data rows used; must satisfy `ndata>=nc`.
        real(rk), allocatable :: TtT(:,:), TtX(:,:), Xsol(:,:), design(:,:)
        real(rk) :: N(0:this%degree), ba
        integer :: i, a, b, k, ia, ib, first, dim_

        if (.not. this%err%ok) return

        if (this%is_rational()) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "Polynomial fitting cannot retain nonuniform rational weights.",&
                location   = "lsq_fit_bspline",&
                suggestion = "Use lsq_fit_nurbs() or reset the curve with uniform weights.")
            return
        end if

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
            allocate(design(ndata,this%nc), source=0.0_rk)
            do i = 1, ndata
                call basis_bspline_der(&
                    map_parameter(Xt(i), this%knot, this%nc, this%degree, this%wrap_parameters),&
                    this%knot, this%nc, this%degree, 0, first, N)
                do a = 0, this%degree
                    ia = first + a
                    if (ia < 1 .or. ia > this%nc) cycle
                    design(i,ia) = N(a)
                end do
            end do
            Xsol = solve_least_squares_qr(design, Xdata(1:ndata,:))
        end if
        if (size(Xsol,1) /= this%nc .or. size(Xsol,2) /= dim_) then
            call this%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_nurbs_curve",&
                message    = "The least-squares design is numerically rank deficient.",&
                location   = "lsq_fit_bspline",&
                suggestion = "Use distinct parameters with full spline-space coverage or reduce nc.")
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
        this%rational = nonuniform_weights(this%Wc)

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
                this%rational = nonuniform_weights(this%Wc)
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
            this%rational = nonuniform_weights(this%Wc)
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
