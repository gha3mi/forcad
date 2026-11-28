! author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Tensor-product B-spline and NURBS surfaces.
!!
!! For degrees \(p\) and \(q\), the rational tensor-product surface is
!!
!! \[\mathbf{S}(u,v)=
!!   \frac{\sum_{j=1}^{n_2}\sum_{i=1}^{n_1}
!!   N_{i,p}(u)M_{j,q}(v)w_{ij}\mathbf{P}_{ij}}
!!   {\sum_{j=1}^{n_2}\sum_{i=1}^{n_1}
!!   N_{i,p}(u)M_{j,q}(v)w_{ij}}.\]
!!
!! Control subscripts in this module's formulas are one-based to match Fortran
!! storage. This is equivalent to conventional zero-based notation after
!! subtracting one from each control index.
!!
!! Let \(B_{ij}=N_{i,p}M_{j,q}\) and
!! \(W=\sum_{ij}B_{ij}w_{ij}\). The rational tensor basis
!! \(R_{ij}=B_{ij}w_{ij}/W\) is nonnegative and forms a partition of unity on
!! the active domain because the weights are positive. At parameter values
!! where the requested mixed derivative exists,
!! \(\sum_{ij}\partial^{a+b}R_{ij}/(\partial u^a\partial v^b)=0\) for
!! every positive total derivative order \(a+b\).
!! Polynomial tensor-basis derivatives vanish when an order exceeds the degree
!! in its direction, but rational quotient derivatives generally need not.
!! A common positive scaling of all weights is only a projective gauge and
!! cannot change the surface. Rational values scale nonzero weighted basis
!! contributions; derivative recurrences share a weight gauge. Homogeneous
!! refinement selects a common representable gauge for weights and weighted
!! controls and restores the original weights afterward. Geometric derivatives
!! use nested division and do not explicitly form \(W^2\).
!!
!! Control points are flattened with direction 1 varying fastest:
!! `linear(i,j)=i+(j-1)*nc(1)`. The physical coordinate is the second array
!! dimension. Omitting weights selects the polynomial tensor-product
!! B-spline form. Exactly equal weights cancel; every nonconstant positive
!! weight list uses rational normalization.
!!
!! Each direction accepts any valid finite, nondecreasing knot vector,
!! including uniform or nonuniform, clamped, one-sided clamped, unclamped, and
!! periodic-form vectors. A complete cyclic knot extension together with
!! repeated control and weight layers is verified directionally by
!! [[nurbs_surface:get_parameter_topology]]; malformed periodic-form data remain
!! valid bounded splines. Its active interval is
!! `[knot(degree+1), knot(nc+1)]`. Optional parameter wrapping performs modulo
!! evaluation independently by direction and does not impose geometric closure.
!! At an interior knot of multiplicity \(s_d\), the spline space has guaranteed
!! directional continuity \(C^{p_d-s_d}\), although special control data can be
!! smoother. Geometry and basis queries accept one trace selector per direction
!! (−1 left, 0 default, +1 right). The default takes the right trace internally.
!! At a closed endpoint both sides use the inward trace.
!!
!! Directional refinement applies an exact one-dimensional operator to every
!! homogeneous control line, preserving
!! \(\mathbf S_{\mathrm{old}}(u,v)=\mathbf S_{\mathrm{new}}(u,v)\) to working
!! precision on the active domain. Insertion and degree elevation preserve a
!! verified periodic direction. Other non-clamped directional representations
!! are converted by degree elevation to an equivalent active-domain
!! open-clamped form. Interior knot removal preserves periodic directions by
!! applying the bounded removal test to a guarded cyclic extension.
!! For an embedded surface with Jacobian
!! \(\mathbf J=[\mathbf S_{,u}\ \mathbf S_{,v}]\), metric
!! \(\mathbf G=\mathbf J^T\mathbf J\), and parameter area \(\mathrm du\mathrm dv\),
!!
!! \[
!! \mathrm dA=\sqrt{\det\mathbf G}\,\mathrm du\,\mathrm dv .
!! \]
!!
!! The [[nurbs_surface]] API provides evaluation, mixed basis derivatives, IGA
!! element data, geometry-preserving refinement, surface quadrature, fitting,
!! discrete and local continuous nearest-point searches, transformations, and
!! VTK/IGES export. Recoverable failures are recorded in
!! [[nurbs_surface:err]].
!!
!! **Primary reference:** L. Piegl and W. Tiller, *The NURBS Book*, 2nd ed.,
!! Springer, 1997.
!! [doi:10.1007/978-3-642-59223-2](https://doi.org/10.1007/978-3-642-59223-2).
!! **Finite-range reference:** IEEE, *IEEE Standard for Floating-Point
!! Arithmetic*, IEEE Std 754-2019, 2019.
!! [IEEE 754-2019](https://standards.ieee.org/ieee/754/6210/).
module forcad_nurbs_surface

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    use forcad_kinds, only: rk
    use forcad_nurbs_curve, only: nurbs_curve
    use forcad_utils, only: valid_trace_side, tensor_ansatz, tensor_ansatz_element, tensor_basis_jet, tensor_measure
    use forcad_utils, only: elemConn_C0, ndgrid, compute_multiplicity, compute_knot_vector, basis_bspline_der, &
        geometry_derivative, tensor_basis_derivative_at, coarsen_control_lines, reduce_control_lines, refine_control_net, &
        get_array, tetragon_Xc, spline_connectivity, spline_span_plan, bezier_element_data, repelem, unique, rotate_points, &
        transform_geometry, export_vtk_legacy, valid_knot_vector, knot_start, knot_end, &
        knot_tolerance, active_knots, active_span, active_span_count, infer_knot_shape, &
        map_parameter, &
        periodic_topology, valid_weights, nonuniform_weights, query_parameter_axis, &
        normalize_rational_basis, unpack_refinement_controls, slice_control_net, &
        isoparametric_grid, restrict_control_net, reparameterize_knots, reorient_patch, map_parameter_values, &
        nearest_spline_point, nearest_sample_index
    use fordebug, only: debug

    implicit none

    private
    public nurbs_surface
#ifndef FORCAD_PUBLIC_DOCS
    public compute_Tgc, compute_dTgc, surface_geometry_derivative
#endif

    !===============================================================================
    ! author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Mutable tensor-product B-spline or NURBS surface.
    !!
    !! Initialize the object with [[nurbs_surface:set]]. The two entries of
    !! `degree`, `nc`, `ng`, and wrapping options correspond to parameter
    !! directions 1 and 2. `degree(d)` is the polynomial degree, while spline
    !! order is `degree(d)+1`. Weights must be finite and strictly positive.
    !! Control coordinates are required to be finite and have at least one
    !! component; coordinate finiteness is a caller precondition in setter
    !! overloads that do not explicitly diagnose it.
    !!
    !! [[nurbs_surface:create]] stores a Cartesian product of parameter samples.
    !! Sample and control arrays follow direction-1-fastest ordering. Exact
    !! insertion and elevation operate directionally in homogeneous coordinates;
    !! removal is accepted only when an internal roundoff-scaled geometry test passes.
    !!
    !! Successful setters discard sampled points and connectivity while keeping
    !! parameters for create() to refresh. Fixed-shape modify_Xc/modify_Wc
    !! operations discard samples and retain control and IGA connectivity.
    !!
    !! @warning Fitting and the separate
    !! `rotate_Xc` or `translate_Xc` operations do not automatically recompute
    !! every sampled or connectivity cache. [[nurbs_surface:transform]] instead
    !! keeps an existing `Xg` cache affine-consistent with `Xc`. Rebuild
    !! connectivity after a control-count or spline-space change.
    !! @endwarning
    !!
    !! Diagnostic codes are: 100 invalid argument or index; 101 shape mismatch;
    !! 102--105 missing geometry state; 106 underdetermined fitting problem; 107
    !! invalid parameters; 108 singular system, Jacobian, or metric; and 109
    !! iterative nonconvergence.
    !!
    !! @note Getter generics return allocatable copies of object-owned arrays.
    !! @endnote
    type nurbs_surface
        real(rk), allocatable, private :: Xc(:,:)  !! Flattened control net `[product(nc),ncoord]`.
        real(rk), allocatable, private :: Xg(:,:)  !! Cached sampled points `[product(ng),ncoord]`.
        real(rk), allocatable, private :: Wc(:)    !! Positive flattened control weights `[product(nc)]`.
        real(rk), allocatable, private :: Xt1(:)   !! Cached parameters in direction 1 `[ng(1)]`.
        real(rk), allocatable, private :: Xt2(:)   !! Cached parameters in direction 2 `[ng(2)]`.
        real(rk), allocatable, private :: Xt(:,:)  !! Flattened parameter grid `[product(ng),2]`.
        real(rk), allocatable, private :: knot1(:) !! Knot vector in direction 1.
        real(rk), allocatable, private :: knot2(:) !! Knot vector in direction 2.
        integer, private :: degree(2) = 0          !! Polynomial degrees \([p,q]\).
        integer, private :: nc(2) = 0              !! Control-point counts by direction.
        integer, private :: ng(2) = 0              !! Cached sample counts by direction.
        logical, private :: rational = .false.     !! True when nonuniform weights affect the basis.
        logical, private :: wrap_parameters(2) = .false. !! Modulo evaluation policy by direction.
        integer, allocatable, private :: elemConn_Xc_vis(:,:) !! Control-net VTK connectivity.
        integer, allocatable, private :: elemConn_Xg_vis(:,:) !! Sampled-surface VTK connectivity.
        integer, allocatable, private :: elemConn(:,:)        !! IGA element-to-control-point map.

        type(debug) :: err !! Recoverable diagnostic state for the most recent failed operation.
    contains
        procedure, private :: set1         !!> Set explicit directional knot vectors, control points, and weights.
        procedure, private :: set2         !!> Build knot vectors from breakpoints, degrees, and continuities.
        procedure, private :: set3         !!> Construct a Bezier or rational Bezier surface.
        procedure, private :: set4         !!> Construct an open uniform surface from degrees and sizes.
        generic :: set => set1, set2, set3, set4 !!> Replace surface geometry or initialize a fitting space.
        procedure :: create                !!> Sample a tensor grid and cache parameters, points, and connectivity.
        procedure, private :: cmp_Xg_scalar, cmp_Xg_vector
        generic :: cmp_Xg => cmp_Xg_scalar, cmp_Xg_vector !!> Evaluate a point or paired batch without changing samples.
        procedure :: slice                 !!> Extract an exact isoparametric curve by fixing one parameter.
        procedure, private :: bezier_element, bezier_all
        generic :: bezier => bezier_element, bezier_all !!> Extract one exact Bezier element or decompose the patch.
        procedure :: restrict              !!> Extract a bounded patch on a literal parameter subinterval.
        procedure :: split                 !!> Split a patch at an interior parameter, preserving both one-sided traces.
        procedure :: reverse               !!> Reverse one parameter direction and its control/weight ordering.
        procedure, private :: reparameterize_interval, reparameterize_axes
        generic :: reparameterize => reparameterize_interval, reparameterize_axes !!> Rescale or reorder parameter axes.
        procedure, private :: get_Xc_all   !!> Return the flattened control net.
        procedure, private :: get_Xci      !!> Return one flattened control point.
        procedure, private :: get_Xcid     !!> Return one coordinate of one control point.
        generic :: get_Xc => get_Xc_all, get_Xci, get_Xcid !!> Query control-point data.
        procedure, private :: get_Xg_all   !!> Return all cached geometry points.
        procedure, private :: get_Xgi      !!> Return one cached geometry point.
        procedure, private :: get_Xgid     !!> Return one coordinate of one cached point.
        generic :: get_Xg => get_Xg_all, get_Xgi, get_Xgid !!> Query cached geometry points.
        procedure, private :: get_Wc_all   !!> Return all flattened control weights.
        procedure, private :: get_Wci      !!> Return one control weight.
        generic :: get_Wc => get_Wc_all, get_Wci !!> Query control weights.
        procedure :: get_Xt                !!> Return cached parameter samples in one direction.
        procedure, private :: get_knot_all !!> Return the knot vector selected by direction.
        procedure, private :: get_knoti    !!> Return one knot value selected by direction and index.
        generic :: get_knot => get_knoti, get_knot_all !!> Query knot-vector data.
        procedure :: get_ng                !!> Return cached sample counts by direction.
        procedure, private :: get_nc_dir   !!> Return the control-point count in one direction.
        procedure, private :: get_nc_all   !!> Return control-point counts in both directions.
        generic :: get_nc => get_nc_all, get_nc_dir !!> Query directional control-point counts.
        procedure :: cmp_degree            !!> Infer both degrees from knots and control-net dimensions.
        procedure, private :: get_degree_all !!> Return both polynomial degrees.
        procedure, private :: get_degree_dir !!> Return the degree in one direction.
        generic :: get_degree => get_degree_all, get_degree_dir !!> Query polynomial degrees.
        procedure :: finalize              !!> Release all owned and derived surface data while preserving diagnostics.
        procedure :: cmp_elem_Xc_vis       !!> Build quadrilateral connectivity for the control net.
        procedure :: cmp_elem_Xg_vis       !!> Build quadrilateral connectivity for cached samples.
        procedure :: cmp_elem_Xth          !!> Build quadrilateral connectivity for the parameter grid.
        procedure :: cmp_elem              !!> Build nonzero-span tensor-product IGA connectivity.
        procedure :: cmp_dof_map           !!> Map stored periodic controls to independent cyclic degrees of freedom.
        procedure :: get_elem_Xc_vis       !!> Return control-net visualization connectivity.
        procedure :: get_elem_Xg_vis       !!> Return sampled-surface visualization connectivity.
        procedure :: get_elem              !!> Return IGA element-to-control-point connectivity.
        procedure :: set_elem_Xc_vis       !!> Replace control-net visualization connectivity.
        procedure :: set_elem_Xg_vis       !!> Replace sampled-surface visualization connectivity.
        procedure :: set_elem              !!> Replace stored IGA connectivity.
        procedure :: export_Xc             !!> Export the control net in legacy VTK format.
        procedure :: export_Xg             !!> Export cached surface samples in legacy VTK format.
        procedure :: export_Xth            !!> Export the sampled parameter surface in legacy VTK format.
        procedure :: cmp_Xth_in_Xg         !!> Build embedded active isoparametric grid lines.
        procedure :: export_Xth_in_Xg      !!> Export parameter-grid lines embedded in physical space.
        procedure :: export_iges           !!> Export the surface as IGES entity 128.
        procedure, private :: modify_Xc_scalar, modify_Xc_all
        generic :: modify_Xc => modify_Xc_scalar, modify_Xc_all !!> Replace one coordinate or a fixed-shape control net.
        procedure, private :: modify_Wc_scalar, modify_Wc_all
        generic :: modify_Wc => modify_Wc_scalar, modify_Wc_all !!> Replace one weight or all positive control weights.
        procedure :: get_multiplicity      !!> Return multiplicities of all distinct stored knots in one direction.
        procedure :: get_continuity        !!> Return \(p_d-s\) metadata for each distinct stored knot.
        procedure :: cmp_nc                !!> Compute control-point counts implied by knots and degrees.
        procedure, private :: basis_vector !!> Evaluate all tensor-product basis functions at many points.
        procedure, private :: basis_scalar !!> Evaluate all tensor-product basis functions at one point.
        generic :: basis => basis_vector, basis_scalar !!> Evaluate the dense rational or polynomial basis.
        procedure, private :: derivative_vector !!> Evaluate parametric basis gradients at many points.
        procedure, private :: derivative_scalar !!> Evaluate the parametric basis gradient at one point.
        generic :: derivative => derivative_vector, derivative_scalar !!> Evaluate first parametric derivatives.
        procedure, private :: derivative2_vector !!> Evaluate packed basis Hessians at many parameter pairs.
        procedure, private :: derivative2_scalar !!> Evaluate one packed basis Hessian.
        generic :: derivative2 => derivative2_vector, derivative2_scalar !!> Evaluate complete parametric Hessians.
        procedure, private :: derivative_order_vector !!> Evaluate a requested mixed derivative at many points.
        procedure, private :: derivative_order_scalar !!> Evaluate a requested mixed derivative at one point.
        generic :: derivative_order => derivative_order_vector, derivative_order_scalar !!> Evaluate dense mixed derivatives.
        procedure, private :: derivative_order_active_vector !!> Evaluate locally active mixed derivatives at many points.
        procedure, private :: derivative_order_active_scalar !!> Evaluate locally active mixed derivatives at one point.
        !> Evaluate active mixed derivatives.
        generic :: derivative_order_active => &
            derivative_order_active_vector, derivative_order_active_scalar
        procedure :: insert_knots          !!> Insert knots while preserving active-domain geometry.
        procedure :: elevate_degree        !!> Elevate degree while preserving active-domain geometry in one direction.
        procedure :: reduce_degree         !!> Reduce degree under one physical geometry-error budget.
        procedure :: is_valid              !!> Report whether complete, internally consistent surface geometry is stored.
        procedure :: is_rational           !!> Report whether nonuniform weights affect the current basis.
        procedure :: get_parameter_wrapping !!> Report modulo parameter mapping in one direction.
        procedure :: get_parameter_topology !!> Classify one stored parameter topology.
        procedure :: put_to_nurbs          !!> Map structured points to the current parameter domain and cache them.
        procedure :: remove_knots          !!> Remove directional knots that pass an internal geometry test.
        procedure :: transform             !!> Apply one affine map to control points and cached geometry.
        procedure :: rotate_Xc             !!> Rotate control points in three-dimensional embedding space.
        procedure :: rotate_Xg             !!> Rotate cached geometry points in three-dimensional embedding space.
        procedure :: translate_Xc          !!> Translate the control net by a physical-space vector.
        procedure :: translate_Xg          !!> Translate cached geometry points by a physical-space vector.
        procedure :: nearest_point         !!> Select the nearest cached geometry sample.
        procedure :: nearest_point2        !!> Minimize squared distance with bounded Newton iterations and backtracking.
        procedure :: ansatz                !!> Evaluate element-local physical basis data and differential area.
        procedure :: ansatz_element        !!> Evaluate all quadrature data for one analysis element.
        procedure :: cmp_area              !!> Integrate surface area with elementwise tensor Gauss quadrature.
        procedure :: lsq_fit_bspline       !!> Fit polynomial control points by structured least squares.
        procedure :: lsq_fit_nurbs         !!> Fit control points and positive weights by damped nonlinear least squares.

        ! Shapes
        procedure :: set_tetragon          !!> Construct a bilinear tetragonal surface.
        procedure :: set_ring              !!> Construct an exact rational annular surface.
        procedure :: set_half_ring         !!> Construct an exact rational half-annulus.
        procedure :: set_C                 !!> Construct the predefined planar C-shaped surface.
    end type
    !===============================================================================

    !> Evaluate polynomial or rational surface points from explicit spline data.
    !! This implementation generic is private; user code should normally call
    !! [[nurbs_surface:cmp_Xg]].
    interface compute_Xg
        module procedure compute_Xg_vector
        module procedure compute_Xg_scalar
    end interface

    !> Evaluate dense tensor-product surface basis values.
    !! Scalar overloads return `product(nc)` values; vector overloads return one
    !! row per parameter pair. Rational overloads normalize with `Wc`.
    !! Scalar `elem` selections preserve ordering and repeated indices. When
    !! size(Wc)==size(elem) differs from product(nc), weights are local to that
    !! selection; otherwise normalization uses the complete active support.
    !! For \(B_{ij}=N_{i,p}(u)M_{j,q}(v)\),
    !! \(R_{ij}=B_{ij}w_{ij}/\sum_{ab}B_{ab}w_{ab}\).
    interface compute_Tgc
        module procedure compute_Tgc_nurbs_2d_vector
        module procedure compute_Tgc_bspline_2d_vector
        module procedure compute_Tgc_nurbs_2d_scalar
        module procedure compute_Tgc_bspline_2d_scalar
    end interface

    !> Evaluate dense first parametric derivatives of the surface basis.
    !! Derivative component 1 is with respect to the first parameter and
    !! component 2 to the second. Values and gradients share one local pass.
    !! Scalar `elem` selections support either full weights or local weights
    !! when size(Wc)==size(elem) differs from product(nc).
    !! For \(\alpha\in\{u,v\}\), rational derivatives use
    !!
    !! \[
    !! R_{ij,\alpha}
    !! =\frac{w_{ij}(B_{ij,\alpha}W-B_{ij}W_{,\alpha})}{W^2}.
    !! \]
    interface compute_dTgc
        module procedure compute_dTgc_nurbs_2d_vector
        module procedure compute_dTgc_bspline_2d_vector
        module procedure compute_dTgc_nurbs_2d_scalar
        module procedure compute_dTgc_bspline_2d_scalar
    end interface

    !> Evaluate requested basis values, gradients, and packed Hessians through one local tensor calculation.
    interface compute_basis_data
        module procedure compute_basis_data_vector
        module procedure compute_basis_data_scalar
    end interface

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Apply one affine transformation to the surface in place.
    !!
    !! Control points and an existing sampled-geometry cache are transformed by
    !!
    !! \[
    !! \widetilde{\mathbf x}
    !! =\mathbf c+\mathbf A(\mathbf x-\mathbf c)+\mathbf b .
    !! \]
    !!
    !! The optional vectors `center` and `translation` represent \(\mathbf c\)
    !! and \(\mathbf b\). Affine covariance follows from partition of unity, so
    !! the transformed control net represents the exact affine image of the
    !! original NURBS surface. Knots, weights, parameters, and connectivity are
    !! unchanged. A singular matrix remains an exact map but can reduce surface
    !! rank and invalidate an IGA parameterization.
    pure subroutine transform(this, linear, translation, center)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: linear(:,:)
            !! Square linear map `[ncoord,ncoord]`.
        real(rk), intent(in), contiguous, optional :: translation(:)
            !! Translation vector `[ncoord]`; zero when absent.
        real(rk), intent(in), contiguous, optional :: center(:)
            !! Transformation center `[ncoord]`; the origin when absent.

        call transform_geometry(this%Xc,this%Xg,linear,translation,center,this%err,"forcad_nurbs_surface")
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report whether the defining surface geometry is complete and consistent.
    !!
    !! For directional degrees \(p_d\) and control counts \(n_d\), both knot
    !! vectors must satisfy
    !! \(\operatorname{size}(U_d)=n_d+p_d+1\). The flattened control net must
    !! contain \(n_1n_2\) finite points with at least one physical coordinate.
    !! Explicit weights, when stored, must comprise the same number of finite
    !! positive values and agree with the cached rational classification.
    !! Sample and connectivity caches are derived state and are not inspected.
    pure function is_valid(this) result(ok)
        class(nurbs_surface), intent(in) :: this
            !! Surface geometry to inspect.
        logical :: ok

        ok = .false.
        if (.not. this%err%ok) return
        if (.not. allocated(this%knot1) .or. .not. allocated(this%knot2) .or. &
            .not. allocated(this%Xc)) return
        if (.not. valid_knot_vector(this%knot1, this%nc(1), this%degree(1))) return
        if (.not. valid_knot_vector(this%knot2, this%nc(2), this%degree(2))) return
        if (size(this%Xc,1) /= product(this%nc) .or. size(this%Xc,2) < 1) return
        if (.not. all(ieee_is_finite(this%Xc))) return

        if (allocated(this%Wc)) then
            if (size(this%Wc) /= product(this%nc) .or. .not. valid_weights(this%Wc)) return
            if (this%rational .neqv. &
                nonuniform_weights(this%Wc)) return
        else if (this%rational) then
            return
        end if

        ok = .true.
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report modulo parameter mapping in one surface direction.
    !! Invalid directions or an active diagnostic return `.false.`.
    pure elemental function get_parameter_wrapping(this, dir) result(wrap_parameters)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        logical :: wrap_parameters

        wrap_parameters = .false.
        if (.not. this%err%ok .or. dir < 1 .or. dir > 2) return
        wrap_parameters = this%wrap_parameters(dir)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Classify one stored parameter topology.
    !!
    !! The result is `"periodic"` only when the selected knot extension,
    !! every repeated control layer, and any stored weight layer satisfy the
    !! cyclic identities in `periodic_topology`. Every other valid direction
    !! is `"bounded"`. Invalid directions and incomplete or inconsistent
    !! objects return `"invalid"`. Parameter wrapping is independent.
    pure function get_parameter_topology(this, dir) result(topology)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        character(len=8) :: topology

        topology = "invalid"
        if (.not. this%is_valid() .or. dir < 1 .or. dir > 2) return

        topology = "bounded"
        if (dir==1) then
            if (periodic_topology(this%knot1,this%degree(1),this%Xc,Wc=this%Wc,nc=this%nc,dir=1)) topology="periodic"
        else
            if (periodic_topology(this%knot2,this%degree(2),this%Xc,Wc=this%Wc,nc=this%nc,dir=2)) topology="periodic"
        end if
    end function
    !===============================================================================




    !===============================================================================
    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate all quadrature points of one surface analysis element.
    !!
    !! This batched form of [[nurbs_surface:ansatz]] constructs the tensor
    !! Gauss rule and locates the two knot spans once per element. It returns
    !! contiguous quadrature columns with shapes
    !! `Tgc(nlocal,nquadrature)`,
    !! `dTgc_dXg(nlocal,ncoord,nquadrature)`, and `dA(nquadrature)`.
    !! Optional physical points use `Xg(ncoord,nquadrature)` and optional
    !! physical Hessians use
    !! `d2Tgc_dXg2(nlocal,ncoord,ncoord,nquadrature)`.
    !!
    !! Allocatable outputs have `intent(inout)` and retain their allocation when
    !! the required shape is unchanged. Rational values and derivatives use the
    !! same exact radix-power weight normalization and quotient recurrence as
    !! [[nurbs_surface:ansatz]]. A failure at any quadrature point zeros all
    !! element outputs before setting `err`, so callers cannot assemble a
    !! partially evaluated element.
    !!
    !! The planar inverse-Jacobian and embedded contravariant transformations,
    !! strict orientation policy, and mathematical references are documented by
    !! [[nurbs_surface:ansatz]].
    pure subroutine ansatz_element(this, ie, Tgc, dTgc_dXg, dA, ngauss, Xg, d2Tgc_dXg2, strict, err, dXg_dXt, &
        dmeasure_dXc, dmeasure_dlogWc, dgradient_dXc, dgradient_dlogWc, djacobian_dXc, djacobian_dlogWc)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in) :: ie
            !! One-based direction-1-fastest active element index.
        real(rk), allocatable, intent(inout) :: Tgc(:,:)
            !! Basis values `[nlocal,nquadrature]`.
        real(rk), allocatable, intent(inout) :: dTgc_dXg(:,:,:)
            !! Physical gradients `[nlocal,ncoord,nquadrature]`.
        real(rk), allocatable, intent(inout) :: dA(:)
            !! Weighted differential areas `[nquadrature]`.
        integer, intent(in), optional :: ngauss(2)
            !! Optional positive directional quadrature counts; default `degree+1`.
        real(rk), allocatable, target, intent(inout), optional :: Xg(:,:)
            !! Optional physical points `[ncoord,nquadrature]`.
        real(rk), allocatable, target, intent(inout), optional :: d2Tgc_dXg2(:,:,:,:)
            !! Optional physical Hessians `[nlocal,ncoord,ncoord,nquadrature]`.
        logical, intent(in), optional :: strict
            !! Reject planar orientation reversal in addition to rank deficiency.
        type(debug), intent(out), optional :: err
            !! Caller-owned diagnostic; when present, object diagnostics are neither read nor changed.
        real(rk), allocatable, target, intent(inout), optional :: dXg_dXt(:,:,:)
            !! Raw geometry Jacobians `[ncoord,2,nquadrature]`.
        real(rk), allocatable, target, intent(inout), optional :: dmeasure_dXc(:,:,:), dmeasure_dlogWc(:,:)
            !! Derivatives of each weighted area w.r.t. local controls and log weights, with Gauss point last.
        real(rk), allocatable, intent(inout), optional :: dgradient_dXc(:,:,:,:,:), dgradient_dlogWc(:,:,:,:)
            !! Physical basis gradient derivatives with respect to controls `(A,i,B,j,q)` and log weights `(A,i,B,q)`.
        real(rk), allocatable, intent(inout), optional :: djacobian_dXc(:,:,:), djacobian_dlogWc(:,:,:,:)
            !! Raw Jacobian derivatives for controls (A,alpha,q) and log weights (i,alpha,A,q).

        call tensor_ansatz_element(this%Xc,this%Wc,this%rational,this%degree,this%nc,this%knot1,this%knot2, &
            this%err,err,ie,Tgc,dTgc_dXg,dA,ngauss,strict,Xg,dXg_dXt,d2Tgc_dXg2, &
            dmeasure_dXc,dmeasure_dlogWc,dgradient_dXc,dgradient_dlogWc, &
            djacobian_dXc=djacobian_dXc,djacobian_dlogWc=djacobian_dlogWc)
    end subroutine ansatz_element
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Parameterize rows by two coordinates and cache the mapped surface.
    !!
    !! Columns 1 and 2 of `X` are independently mapped from their data ranges
    !! to the active intervals \([u_a,u_b]\) and \([v_a,v_b]\):
    !!
    !! \[
    !! u_i=u_a+\frac{X_{i1}-x_{\min}}{x_{\max}-x_{\min}}(u_b-u_a),
    !! \qquad
    !! v_i=v_a+\frac{X_{i2}-y_{\min}}{y_{\max}-y_{\min}}(v_b-v_a).
    !! \]
    !!
    !! Remaining columns do not affect the parameterization. The routine
    !! atomically replaces the flattened parameter-point pairs
    !! \(((u_i,v_i),\mathbf S(u_i,v_i))\) and stores their count as
    !! `ng=[size(X,1),1]`. Because arbitrary input rows need not form a Cartesian
    !! product, directional parameter vectors are cleared. When `elemConn` is
    !! absent, sampled visualization connectivity is reset to an empty
    !! quadrilateral topology.
    pure subroutine put_to_nurbs(this, X, elemConn)
        class(nurbs_surface), intent(inout) :: this
            !! Initialized surface whose sampled geometry cache is replaced.
        real(rk), intent(in), contiguous :: X(:,:)
            !! Ordering data `[npoint,ncolumn]`; columns 1 and 2 need nonzero ranges.
        integer, intent(in), contiguous, optional :: elemConn(:,:)
            !! Optional one-based visualization connectivity for the cached points.
        real(rk), allocatable :: Xt(:,:), Xg_new(:,:), parameters(:), knot(:)
        real(rk) :: lo, hi
        integer :: d

        if (.not. this%err%ok) return

        if (.not. allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Control points are not set.",&
                location   = "put_to_nurbs",&
                suggestion = "Call set(...) first before mapping points to the NURBS surface.")
            return
        end if
        if (.not. allocated(this%knot1) .or. .not. allocated(this%knot2)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Knot vector is not set.",&
                location   = "put_to_nurbs",&
                suggestion = "Call set(...) first before mapping points to the NURBS surface.")
            return
        end if
        if (size(X,1) < 1 .or. size(X,2) < 2) then
            call this%err%set(&
                code       = 101,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Input points must have at least one point and two coordinates.",&
                location   = "put_to_nurbs",&
                suggestion = "Provide X(:,:) with size(X,1) >= 1 and size(X,2) >= 2.")
            return
        end if

        allocate(Xt(size(X,1),2))
        do d=1,2
            lo=minval(X(:,d))
            hi=maxval(X(:,d))
            if (hi<=lo) then
                call this%err%set(code=100,severity=1,category="forcad_nurbs_surface",location="put_to_nurbs", &
                    message="Input points have zero extent in a parameter direction.", &
                    suggestion="Provide points with distinct first and second coordinates.")
                return
            end if
            knot=this%get_knot(d)
            call map_parameter_values(X(:,d),lo,hi,knot_start(knot,this%nc(d),this%degree(d)), &
                knot_end(knot,this%nc(d),this%degree(d)),parameters,this%err)
            if (.not.this%err%ok) return
            Xt(:,d)=parameters
        end do
        Xg_new=this%cmp_Xg(Xt)
        if (size(Xg_new,1)/=size(X,1) .or. .not.all(ieee_is_finite(Xg_new))) then
            call this%err%set(code=108,severity=1,category="forcad_nurbs_surface",location="put_to_nurbs", &
                message="Mapped geometry is not finite or complete.",suggestion="Check the spline geometry and parameter ranges.")
            return
        end if
        call move_alloc(Xt, this%Xt)
        call move_alloc(Xg_new, this%Xg)
        this%ng = [size(X,1),1]
        if (allocated(this%Xt1)) deallocate(this%Xt1)
        if (allocated(this%Xt2)) deallocate(this%Xt2)

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
    !> Replace the surface with explicit directional knots and control data.
    !!
    !! `Xc` is direction-1-fastest. If `degree` is absent, an open-clamped
    !! candidate inferred from the initial knot runs takes precedence. If that
    !! candidate is invalid, exhaustive degree/count search must find exactly
    !! one valid pair. Pass explicit degrees whenever this convention is
    !! ambiguous or should not determine the interpretation of non-open data.
    !! `Xc` must contain at least one flattened control row and one finite
    !! coordinate component. Coordinate, knot, and positive-weight invariants
    !! are validated before any core state is replaced.
    !!
    pure subroutine set1(this, knot1, knot2, Xc, Wc, degree, wrap_parameters)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: knot1(:)
            !! Complete knot vector in direction 1.
        real(rk), intent(in), contiguous :: knot2(:)
            !! Complete knot vector in direction 2.
        real(rk), intent(in), contiguous :: Xc(:,:)
            !! Flattened control net `[nc(1)*nc(2),ncoord]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional positive flattened weights `[product(nc)]`.
        integer, intent(in), contiguous, optional :: degree(:)
            !! Optional polynomial degrees `[2]`.
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
            !! Optional modulo evaluation flags `[2]`.
        integer :: degree_(2), nc_(2), ncp
        logical :: ok
        logical :: wrap_parameters_(2)

        if (.not. this%err%ok) return

        if (size(Xc,1) < 1 .or. size(Xc,2) < 1) then
            call this%err%set(&
                code       = 101,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Control points need at least one row and one coordinate.",&
                location   = "set1",&
                suggestion = "Provide Xc with shape [product(nc),ncoord], with positive dimensions.")
            return
        end if
        if (.not. all(ieee_is_finite(Xc))) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Control-point coordinates must all be finite.",&
                location   = "set1",&
                suggestion = "Remove NaN and infinite values from Xc before calling set().")
            return
        end if

        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            if (size(wrap_parameters) < 2) then
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = "forcad_nurbs_surface",&
                    message    = "Surface parameter wrapping needs one value per parametric direction.",&
                    location   = "set1",&
                    suggestion = "Pass one logical value per parametric direction.")
                return
            end if
            wrap_parameters_ = wrap_parameters(1:2)
        end if

        if (present(degree)) then
            if (size(degree) < 2) then
                call this%err%set(&
                    code       = 103,&
                    severity   = 1,&
                    category   = "forcad_nurbs_surface",&
                    message    = "Explicit surface degree must provide one entry per parametric direction.",&
                    location   = "set1",&
                    suggestion = "Pass degree=[p1,p2].")
                return
            end if
            degree_ = degree(1:2)
            nc_ = [size(knot1) - degree_(1) - 1, size(knot2) - degree_(2) - 1]
            ok = valid_knot_vector(knot1, nc_(1), degree_(1)) .and. &
                valid_knot_vector(knot2, nc_(2), degree_(2))
        else
            call infer_knot_shape(knot1, knot2, size(Xc,1), degree_, nc_, ok)
        end if
        if (.not. ok .or. product(nc_) /= size(Xc,1)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Invalid or ambiguous knot vectors for the supplied control points.",&
                location   = "set1",&
                suggestion = "Use valid knot vectors and pass degree=[p1,p2] for non-open or ambiguous knot vectors.")
            return
        end if
        ncp = product(nc_)

        if (present(Wc)) then
            if (size(Wc) /= ncp) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_surface",&
                    message    = "Weights length mismatch: size(Wc) must equal number of control points.",&
                    location   = "set1",&
                    suggestion = "Provide Wc with size(Wc) == nc(1)*nc(2).")
                return
            else if (.not. valid_weights(Wc)) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_surface",&
                    message    = "Invalid rational weights: all weights must be finite and positive.",&
                    location   = "set1",&
                    suggestion = "Use strictly positive finite Wc values for rational NURBS geometry.")
                return
            end if
        end if

        call invalidate_samples(this,space_changed=.true.)
        this%knot1 = knot1
        this%knot2 = knot2
        this%degree = degree_
        this%nc = nc_
        this%Xc = Xc
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
    !> Construct both spline directions from breakpoints and continuities.
    !! `continuity1/2` are converted to knot multiplicities by
    !! `compute_knot_vector`. Optional control data must match the resulting
    !! tensor-product count. In direction \(d\), breakpoints must be finite and
    !! strictly increasing. Values in \([-1,p_d-1]\) retain a breakpoint with
    !! multiplicity \(p_d-c_i\); `-1` gives full multiplicity \(p_d+1\), while
    !! `c_i=p_d` omits that breakpoint. Both generated knot vectors are
    !! validated before object state changes. If `Xc` is absent, any previous
    !! control net is released and the result is a tensor-product spline space
    !! for subsequent fitting, not complete surface geometry. Optional weights
    !! are retained only when explicitly supplied.
    pure subroutine set2(this, Xth_dir1, Xth_dir2, degree, continuity1, continuity2, Xc, Wc)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xth_dir1(:)
            !! Finite, strictly increasing breakpoints in direction 1.
        real(rk), intent(in), contiguous :: Xth_dir2(:)
            !! Finite, strictly increasing breakpoints in direction 2.
        integer, intent(in), contiguous :: degree(:)
            !! Polynomial degrees `[2]`.
        integer, intent(in), contiguous :: continuity1(:)
            !! One value per direction-1 breakpoint in `-1:degree(1)`; the upper value omits it.
        integer, intent(in), contiguous :: continuity2(:)
            !! One value per direction-2 breakpoint in `-1:degree(2)`; the upper value omits it.
        real(rk), intent(in), contiguous, optional :: Xc(:,:)
            !! Optional flattened control net `[product(computed_nc),ncoord]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional positive weights `[product(computed_nc)]`.
        real(rk), allocatable :: knot1_(:), knot2_(:)
        integer :: degree_(2), nc_(2), ncp

        if (.not. this%err%ok) return

        if (size(degree) < 2) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Surface degree must provide one entry per parametric direction.",&
                location   = "set2",&
                suggestion = "Pass degree=[p1,p2].")
            return
        end if

        degree_ = degree(1:2)
        knot1_ = compute_knot_vector(Xth_dir1, degree_(1), continuity1)
        knot2_ = compute_knot_vector(Xth_dir2, degree_(2), continuity2)
        nc_ = [size(knot1_) - degree_(1) - 1, size(knot2_) - degree_(2) - 1]
        if (.not. valid_knot_vector(knot1_, nc_(1), degree_(1)) .or. &
            .not. valid_knot_vector(knot2_, nc_(2), degree_(2))) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Invalid knot vector generated from parameter nodes, degree, and continuity.",&
                location   = "set2",&
                suggestion = "Use nondecreasing parameter nodes and knot multiplicity <= degree+1.")
            return
        end if
        ncp = product(nc_)

        if (present(Xc)) then
            call this%set1(knot1_, knot2_, Xc, Wc, degree=degree_)
            return
        end if

        if (present(Wc)) then
            if (size(Wc) /= ncp) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_surface", &
                    message    = "Weights size mismatch in set2",&
                    location   = "set2", &
                    suggestion = "size(Wc) must equal nc(1)*nc(2)." )
                return
            else if (.not. valid_weights(Wc)) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_surface",&
                    message    = "Invalid rational weights: all weights must be finite and positive.",&
                    location   = "set2",&
                    suggestion = "Use strictly positive finite Wc values for rational NURBS geometry.")
                return
            end if
        end if

        call invalidate_samples(this,space_changed=.true.)
        this%knot1 = knot1_
        this%knot2 = knot2_
        this%degree = degree_
        this%nc = nc_
        this%wrap_parameters = .false.
        if (allocated(this%Xc)) deallocate(this%Xc)
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
    !> Construct a tensor-product Bezier or rational Bezier surface on \([0,1]^2\).
    !! Directional degrees are `nc(1:2)-1`; all endpoint knots are fully clamped.
    pure subroutine set3(this, nc, Xc, Wc)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), contiguous :: nc(:)
            !! Control-point counts `[2]`.
        real(rk), intent(in), contiguous :: Xc(:,:)
            !! Direction-1-fastest control net `[product(nc),ncoord]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional positive weights `[product(nc)]`.
        integer :: nc_(2), ncp
        real(rk), allocatable :: knot1_(:), knot2_(:)

        if (.not. this%err%ok) return

        if (size(nc) < 2) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Bezier surface control-point count must provide one entry per parametric direction.",&
                location   = "set3",&
                suggestion = "Pass nc=[nc1,nc2].")
            return
        end if

        nc_ = nc(1:2)
        ncp = product(nc_)
        if (any(nc_ <= 0) .or. size(Xc,1) /= ncp .or. size(Xc,2) < 1) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Invalid Bezier surface dimensions.",&
                location   = "set3",&
                suggestion = "Require nc > 0 per direction, size(Xc,1) == product(nc), and at least one coordinate dimension.")
            return
        end if

        allocate(knot1_(2*nc_(1)), knot2_(2*nc_(2)), source=0.0_rk)
        knot1_(nc_(1)+1:) = 1.0_rk
        knot2_(nc_(2)+1:) = 1.0_rk
        call this%set1(knot1_, knot2_, Xc, Wc=Wc, degree=nc_-1)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct a clamped open-uniform surface on \([0,1]^2\).
    !! Interior knots are uniformly spaced and simple in each direction.
    pure subroutine set4(this, degree, nc, Xc, Wc)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), contiguous :: degree(:)
            !! Polynomial degrees `[2]`, each smaller than `nc`.
        integer, intent(in), contiguous :: nc(:)
            !! Directional control-point counts `[2]`.
        real(rk), intent(in), contiguous :: Xc(:,:)
            !! Direction-1-fastest control net `[product(nc),ncoord]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional positive weights `[product(nc)]`.
        integer :: m(2), i
        real(rk), allocatable :: knot1_(:), knot2_(:)

        if (.not. this%err%ok) return

        if (size(degree) < 2 .or. size(nc) < 2) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Invalid surface rank for generated knot-vector construction.",&
                location   = "set4",&
                suggestion = "Provide degree and nc arrays with at least two entries.")
            return
        end if

        if (any(degree(1:2) < 0) .or. any(nc(1:2) <= 0) .or. any(degree(1:2) >= nc(1:2)) .or. &
            size(Xc, 1) /= nc(1)*nc(2) .or. size(Xc, 2) < 1) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Invalid surface dimensions for generated knot-vector construction.",&
                location   = "set4",&
                suggestion = "Require degree >= 0, nc > degree per direction, size(Xc,1) == product(nc), and at least one coordinate dimension.")
            return
        end if

        m = nc(1:2) + degree(1:2) + 1
        allocate(knot1_(m(1)), knot2_(m(2)), source=0.0_rk)
        knot1_(degree(1)+2:m(1)-degree(1)-1) = &
            [(real(i,rk)/(nc(1)-degree(1)), i=1,nc(1)-degree(1)-1)]
        knot2_(degree(2)+2:m(2)-degree(2)-1) = &
            [(real(i,rk)/(nc(2)-degree(2)), i=1,nc(2)-degree(2)-1)]
        knot1_(m(1)-degree(1):) = 1.0_rk
        knot2_(m(2)-degree(2):) = 1.0_rk
        call this%set1(knot1_, knot2_, Xc, Wc=Wc, degree=degree(1:2))
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Sample and cache surface geometry.
    !!
    !! Directional `Xt1/Xt2` replace cached vectors; `res1/res2` generate
    !! uniformly spaced vectors when the corresponding explicit vector is
    !! absent. Existing vectors are reused when neither is supplied. Unless
    !! explicit point rows `Xt(:,1:2)` are supplied, `ndgrid(Xt1,Xt2)` forms a
    !! direction-1-fastest Cartesian grid. Paired rows ignore directional inputs
    !! and clear directional caches. With no inputs, create() reuses all retained
    !! parameter rows. Validation and evaluation finish before any cache changes.
    pure subroutine create(this, res1, res2, Xt1, Xt2, Xt)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: res1, res2
            !! Positive directional grid counts; explicit parameter vectors take precedence.
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:)
            !! Directional vectors defining a Cartesian product.
        real(rk), intent(in), contiguous, optional :: Xt(:,:)
            !! Paired rows [npoint,at least 2], taking precedence over directional inputs.
        real(rk), allocatable :: parameters(:,:), axis1(:), axis2(:), points(:,:)
        integer, allocatable :: connectivity(:,:)
        integer :: counts(2)
        logical :: valid1, valid2, grid_request

        if (.not. this%err%ok) return
        if (.not. this%is_valid()) then
            call this%err%set(code=103, severity=1, category="forcad_nurbs_surface", &
                message="Sampling requires valid surface controls, knots, and weights.", location="create", &
                suggestion="Initialize a complete finite surface before sampling.")
            return
        end if
        grid_request = present(res1) .or. present(res2) .or. present(Xt1) .or. present(Xt2)
        if (present(Xt)) then
            if (size(Xt,1)<1 .or. size(Xt,2)<2) then
                call this%err%set(code=107, severity=1, category="forcad_nurbs_surface", &
                    message="Parameter rows require at least one row and two columns.", location="create", &
                    suggestion="Supply Xt(npoint,2) or additional ignored columns.")
                return
            end if
            if (.not. all(ieee_is_finite(Xt(:,1:2)))) then
                call this%err%set(code=107, severity=1, category="forcad_nurbs_surface", &
                    message="Surface parameters must be finite.", location="create", &
                    suggestion="Remove nonfinite values from the first two parameter columns.")
                return
            end if
            parameters = Xt
            counts = [size(Xt,1),1]
        else if (.not. grid_request .and. allocated(this%Xt)) then
            parameters = this%Xt
            counts = this%ng
        else
            call query_parameter_axis(this%knot1,this%nc(1),this%degree(1),this%Xt1,axis1,valid1,res1,Xt1)
            call query_parameter_axis(this%knot2,this%nc(2),this%degree(2),this%Xt2,axis2,valid2,res2,Xt2)
            if (.not. valid1 .or. .not. valid2) then
                call this%err%set(code=107, severity=1, category="forcad_nurbs_surface", &
                    message="Sampling requires finite directional parameters or positive resolutions.", location="create", &
                    suggestion="Provide both grid directions, or reuse an existing parameter grid.")
                return
            end if
            counts = [size(axis1),size(axis2)]
            call ndgrid(axis1,axis2,parameters)
        end if
        if (this%is_rational()) then
            points = compute_Xg(parameters,this%knot1,this%knot2,this%degree,this%nc,counts, &
                this%Xc,Wc=this%Wc,wrap_parameters=this%wrap_parameters)
        else
            points = compute_Xg(parameters,this%knot1,this%knot2,this%degree,this%nc,counts, &
                this%Xc,wrap_parameters=this%wrap_parameters)
        end if
        if (.not. all(ieee_is_finite(points))) then
            call this%err%set(code=108, severity=1, category="forcad_nurbs_surface", &
                message="Sampling produced nonfinite geometry.", location="create", &
                suggestion="Check the coordinate scales and parameter domain.")
            return
        end if
        connectivity = elemConn_C0(counts(1),counts(2),1,1)
        call move_alloc(parameters,this%Xt)
        call move_alloc(points,this%Xg)
        call move_alloc(connectivity,this%elemConn_Xg_vis)
        this%ng = counts
        if (allocated(axis1)) call move_alloc(axis1,this%Xt1)
        if (allocated(axis2)) call move_alloc(axis2,this%Xt2)
        if (present(Xt)) then
            if (allocated(this%Xt1)) deallocate(this%Xt1)
            if (allocated(this%Xt2)) deallocate(this%Xt2)
        end if
    end subroutine create

    !===============================================================================


    !===============================================================================
    !> Inspect geometry dimensions without scanning the control net on every point query.
    pure logical function geometry_ready(this) result(ready)
        class(nurbs_surface), intent(in) :: this
        ready = .false.
        if (.not. this%err%ok .or. .not. allocated(this%Xc) .or. &
            .not. allocated(this%knot1) .or. .not. allocated(this%knot2)) return
        if (any(this%degree<0) .or. any(this%nc<=this%degree)) return
        if (size(this%knot1)/=this%nc(1)+this%degree(1)+1 .or. &
            size(this%knot2)/=this%nc(2)+this%degree(2)+1) return
        if (size(this%Xc,1)/=product(this%nc) .or. size(this%Xc,2)<1) return
        if (allocated(this%Wc)) then
            if (size(this%Wc)/=product(this%nc)) return
        else if (this%rational) then
            return
        end if
        ready = .true.
    end function geometry_ready
    !===============================================================================


    !===============================================================================
    !> Evaluate a point or mixed parametric derivative without changing caches or diagnostics.
    !! At least two finite parameters are required; later entries are ignored.
    pure function cmp_Xg_scalar(this, Xt, order, side) result(Xg)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in), optional :: side(:) !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: order(:) !! Nonnegative parametric derivative order; omitted means position.
        real(rk), allocatable :: Xg(:), rows(:,:)

        rows=cmp_Xg_vector(this,reshape(Xt,[1,size(Xt)]),order,side)
        Xg=reshape(rows,[size(rows)])
    end function cmp_Xg_scalar
    !===============================================================================


    !===============================================================================
    !> Extract an exact isoparametric curve by fixing parameter direction `dir`.
    !! The remaining direction keeps its knots, degree, control ordering and wrapping.
    !! Rational controls are contracted homogeneously; no geometry is sampled.
    !! Interior knots use the right trace; an unwrapped upper endpoint uses the left trace.
    !! Source wrapping is applied, mapping the upper endpoint to the lower endpoint when enabled.
    !! Errors are returned on the curve, including an existing source diagnostic.
    !! Source geometry, parameters and samples remain unchanged.
    pure function slice(this, dir, Xt) result(curve)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), intent(in) :: Xt
        type(nurbs_curve) :: curve
        real(rk), allocatable :: Xc(:,:), Wc(:), knot(:)

        call check_direction(this,dir,"slice",curve%err)
        if (.not.curve%err%ok) return
        knot=this%get_knot(dir)
        call slice_control_net(this%Xc,this%Wc,this%nc,dir,knot,this%degree(dir),Xt,this%wrap_parameters(dir),Xc,Wc,curve%err)
        if (curve%err%ok) call curve%set(knot=this%get_knot(3-dir),Xc=Xc,Wc=Wc, &
            degree=this%degree(3-dir),wrap_parameters=this%wrap_parameters(3-dir))
    end function slice
    !===============================================================================


    !> Extract one exact element; optional polynomial extraction satisfies N=C*B.
    !! Element and support IDs use direction-1-fastest stored-control ordering.
    !! Bounds have shape [2,2]. The piece retains literal span coordinates,
    !! with clamped knots, right/left endpoint traces and disabled wrapping.
    !! Caller err bypasses the source diagnostic; source geometry and caches are untouched.
    pure subroutine bezier_element(this, ie, piece, C, ids, bounds, err)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: ie
        type(nurbs_surface), intent(out), optional :: piece
        real(rk), allocatable, intent(out), optional :: C(:,:), bounds(:,:)
        integer, allocatable, intent(out), optional :: ids(:)
        type(debug), intent(out), optional :: err

        call bezier_extract(this,ie,piece,C,ids,bounds,err)
    end subroutine bezier_element


    !> Share extraction between direct queries and a decomposition's precomputed span plan.
    pure subroutine bezier_extract(this, ie, piece, C, ids, bounds, err, spans)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: ie
        integer, intent(in), optional :: spans(:)
        type(nurbs_surface), intent(out), optional :: piece
        real(rk), allocatable, intent(out), optional :: C(:,:), bounds(:,:)
        integer, allocatable, intent(out), optional :: ids(:)
        type(debug), intent(out), optional :: err
        real(rk), allocatable :: Xc(:,:), Wc(:)
        integer, allocatable :: nodes(:)
        real(rk) :: limits(2,2)
        integer :: p(2)
        type(debug) :: status

        call bezier_check(this,status,present(err))
        if (status%ok) then
            p=this%degree
            call bezier_element_data(this%knot1,this%nc,p,ie,this%Xc,this%Wc, &
                present(piece),Xc,Wc,nodes,limits,status,C=C,knot2=this%knot2,rational=this%rational,spans=spans)
            if (status%ok .and. present(piece)) then
                call piece%set(knot1=repelem(limits(:,1),[p(1)+1,p(1)+1]), &
                    knot2=repelem(limits(:,2),[p(2)+1,p(2)+1]), &
                    Xc=Xc,Wc=Wc,degree=p)
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
        class(nurbs_surface), intent(in) :: this
        type(nurbs_surface), allocatable, intent(out) :: pieces(:)
        type(debug), intent(out), optional :: err
        type(debug) :: status
        integer, allocatable :: plan(:,:)
        integer :: ie, d, element, ne(2), selected(2)

        call bezier_check(this,status,present(err))
        if (status%ok) then
            plan=spline_span_plan(this%nc,this%degree,this%knot1,knot2=this%knot2)
            ne=count(plan>0,dim=1)
            allocate(pieces(product(ne)))
            do ie=1,size(pieces)
                element=ie-1
                do d=1,2
                    selected(d)=plan(mod(element,ne(d))+1,d)
                    element=element/ne(d)
                end do
                call bezier_extract(this,ie,pieces(ie),err=status,spans=selected)
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
        class(nurbs_surface), intent(in) :: this
        type(debug), intent(out) :: err
        logical, intent(in) :: bypass

        err=this%err
        if (bypass) call err%reset()
        if (.not.err%ok) return
        if (.not.allocated(this%Xc) .or. &
            .not.allocated(this%knot1) .or. &
            .not.allocated(this%knot2)) call err%set( &
            code=103,severity=1,category="forcad_nurbs_surface",location="bezier", &
            message="Spline geometry is not initialized.",suggestion="Set knots and controls before Bezier extraction.")
    end subroutine bezier_check


    !> Restrict one parameter to a literal interval inside its active domain.
    !! The selected direction is clamped on [a,b] without rescaling parameters;
    !! its wrapping is disabled and the other direction is unchanged. At a
    !! discontinuity, the lower endpoint takes the right trace and the upper
    !! endpoint the left trace. The result owns fresh geometry and diagnostics.
    pure function restrict(this, dir, a, b) result(patch)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), intent(in) :: a, b
        type(nurbs_surface) :: patch
        real(rk), allocatable :: knot(:), knots(:), Xc(:,:), Wc(:)
        logical :: wrapping(2)

        call check_direction(this,dir,"restrict",patch%err)
        if (.not.patch%err%ok) return
        knot=this%get_knot(dir)
        call restrict_control_net(this%Xc,this%Wc,this%nc,dir,knot,this%degree(dir),a,b,knots,Xc,Wc,patch%err)
        if (.not.patch%err%ok) return
        wrapping=this%wrap_parameters
        wrapping(dir)=.false.
        call set_direction(patch,dir,knots,this%get_knot(3-dir),Xc,Wc,this%degree,wrapping)
    end function restrict


    !> Reverse one parameter while retaining the full knot vector and explicit weights.
    !! On smooth spans, Q(u)=S(lo+hi-u) in the selected direction. At a break,
    !! the new right trace is the old left trace; the endpoint traces are exchanged.
    !! Wrapping is preserved, including this choice of trace at a nonclosed seam.
    !! The result owns geometry only; source geometry and caches remain unchanged.
    pure function reverse(this, dir) result(patch)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        type(nurbs_surface) :: patch
        integer :: axes(2)

        call check_direction(this,dir,"reverse",patch%err)
        if (.not.patch%err%ok) return
        axes=[1,2]
        axes(dir)=-dir
        call orient(this,axes,patch)
    end function reverse

    !> New axis d takes source axis abs(axes(d)); a negative entry reverses its original interval.
    pure function reparameterize_axes(this,axes) result(patch)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: axes(:)
        type(nurbs_surface) :: patch

        call check_direction(this,1,"reparameterize",patch%err)
        if (patch%err%ok) call orient(this,axes,patch)
    end function reparameterize_axes

    pure subroutine orient(this,axes,patch)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: axes(:)
        type(nurbs_surface), intent(inout) :: patch
        integer :: source(2)

        call reorient_patch(this%nc,this%degree,axes,this%knot1,this%Xc,this%Wc,patch%knot1,patch%Xc,patch%Wc,patch%err, &
            this%knot2,patch%knot2)
        if (.not.patch%err%ok) return
        source=abs(axes)
        patch%degree=this%degree(source)
        patch%nc=this%nc(source)
        patch%wrap_parameters=this%wrap_parameters(source)
        patch%rational=this%rational
    end subroutine orient


    !> Affinely map one active domain to a finite increasing interval, preserving exterior knots.
    !! Parameters retain their relative position; geometry, degree, control order,
    !! explicit weights and wrapping are preserved. Mapped derivatives scale by
    !! ((old_hi-old_lo)/(b-a))**order in the selected direction.
    pure function reparameterize_interval(this, dir, a, b) result(patch)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), intent(in) :: a, b
        type(nurbs_surface) :: patch
        real(rk), allocatable :: knot(:), knots(:)

        call check_direction(this,dir,"reparameterize",patch%err)
        if (.not.patch%err%ok) return
        knot=this%get_knot(dir)
        call reparameterize_knots(knot,this%nc(dir),this%degree(dir),a,b,knots,patch%err, &
            category="forcad_nurbs_surface")
        if (patch%err%ok) call set_direction(patch,dir,knots,this%get_knot(3-dir), &
            this%Xc,this%Wc,this%degree,this%wrap_parameters)
    end function reparameterize_interval


    !> Commit one changed knot vector; callers snapshot the unchanged vector before mutation.
    pure subroutine set_direction(this, dir, knot, other_knot, Xc, Wc, degree, wrapping)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in) :: dir, degree(2)
        real(rk), intent(in), contiguous :: knot(:), other_knot(:), Xc(:,:)
        real(rk), allocatable, intent(in) :: Wc(:)
        logical, intent(in) :: wrapping(2)

        if (dir==1) then
            call this%set(knot1=knot,knot2=other_knot,Xc=Xc,Wc=Wc,degree=degree,wrap_parameters=wrapping)
        else
            call this%set(knot1=other_knot,knot2=knot,Xc=Xc,Wc=Wc,degree=degree,wrap_parameters=wrapping)
        end if
    end subroutine set_direction


    !> Split at a strictly interior parameter, retaining the original coordinates.
    !! Each output is a restriction with the appropriate one-sided boundary trace.
    !! Failure clears both outputs and reports the same diagnostic on each.
    pure subroutine split(this, dir, Xt, left, right)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), intent(in) :: Xt
        type(nurbs_surface), intent(out) :: left, right
        real(rk), allocatable :: knot(:)
        type(debug) :: failure

        call check_direction(this,dir,"split",left%err)
        if (left%err%ok) then
            knot=this%get_knot(dir)
            if (.not.ieee_is_finite(Xt) .or. Xt<=knot(this%degree(dir)+1) .or. Xt>=knot(this%nc(dir)+1)) then
                call left%err%set(code=107,severity=1,category="forcad_nurbs_surface",location="split", &
                    message="The split parameter must be finite and strictly inside the active interval.", &
                    suggestion="Use an interior parameter without modulo wrapping.")
            else
                left=this%restrict(dir,knot(this%degree(dir)+1),Xt)
                if (left%err%ok) right=this%restrict(dir,Xt,knot(this%nc(dir)+1))
            end if
        end if
        if (left%err%ok .and. right%err%ok) return
        failure=left%err
        if (left%err%ok) failure=right%err
        call left%finalize()
        call right%finalize()
        left%err=failure
        right%err=failure
    end subroutine split


    !> Inspect directional knot sizes without allocating getter copies.
    pure function knot_sizes(this) result(lengths)
        class(nurbs_surface), intent(in) :: this
        integer :: lengths(2)
        lengths=0
        if (allocated(this%knot1)) lengths(1)=size(this%knot1)
        if (allocated(this%knot2)) lengths(2)=size(this%knot2)
    end function knot_sizes


    !> Validate source geometry and direction for geometry-only extraction results.
    pure subroutine check_direction(this, dir, location, err)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        character(len=*), intent(in) :: location
        type(debug), intent(out) :: err

        err=this%err
        if (.not.err%ok) return
        if (dir<1 .or. dir>2) then
            call err%set(code=100,severity=1,category="forcad_nurbs_surface",location=location, &
                message="Invalid parameter direction.",suggestion="Use dir=1 or dir=2.")
        else if (.not.this%is_valid()) then
            call err%set(code=102,severity=1,category="forcad_nurbs_surface",location=location, &
                message="Surface geometry is not initialized or is invalid.",suggestion="Set complete surface geometry first.")
        end if
    end subroutine check_direction


    !===============================================================================
    !> Evaluate paired parameter rows without changing caches or diagnostics.
    !! Only the first two columns are used. A valid empty batch returns [0,ncoord];
    !! invalid parameters or missing geometry return [0,0]. Strided inputs are accepted.
    pure function cmp_Xg_vector(this, Xt, order, side) result(Xg)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in), optional :: side(:) !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in) :: Xt(:,:)
        integer, intent(in), optional :: order(:) !! Nonnegative order; omitted means position.
        real(rk), allocatable :: Xg(:,:)

        call surface_geometry_derivative(this,Xt,order,Xg,side)
    end function cmp_Xg_vector

    !> Evaluate patch geometry, optionally retaining a common binary exponent per derivative row.
    pure subroutine surface_geometry_derivative(patch,Xt,order,Xg,side,row_power)
        class(nurbs_surface), intent(in) :: patch
        real(rk), intent(in) :: Xt(:,:)
        integer, optional, intent(in) :: order(:),side(:)
        real(rk), allocatable, intent(out) :: Xg(:,:)
        integer, optional, allocatable, intent(out) :: row_power(:)
        integer :: selected(2)

        selected=0
        if (present(order)) then
            selected=-1
            if (size(order)==2) selected=order
        end if
        if (any(selected<0) .or. .not. geometry_ready(patch) .or. size(Xt,2)<2 .or. .not. valid_trace_side(side,2)) then
            allocate(Xg(0,0))
            if (present(row_power)) allocate(row_power(0))
            return
        end if
        if (.not. all(ieee_is_finite(Xt(:,1:2)))) then
            allocate(Xg(0,0))
            if (present(row_power)) allocate(row_power(0))
            return
        end if
        if (any(selected>0)) then
            call geometry_derivative(Xt(:,1),patch%nc,patch%degree,patch%knot1,patch%Xc,patch%Wc, &
                selected,patch%wrap_parameters,Xg,knot2=patch%knot2,Xt2=Xt(:,2),side=side,row_power=row_power)
        else if (patch%is_rational()) then
            Xg = compute_Xg(Xt,patch%knot1,patch%knot2,patch%degree,patch%nc,Xc=patch%Xc, &
                wrap_parameters=patch%wrap_parameters,Wc=patch%Wc,side=side)
        else
            Xg = compute_Xg(Xt,patch%knot1,patch%knot2,patch%degree,patch%nc,Xc=patch%Xc, &
                wrap_parameters=patch%wrap_parameters,side=side)
        end if
        if (present(row_power).and.all(selected==0)) allocate(row_power(size(Xg,1)),source=0)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocated copy of the flattened control net.
    !!
    !! The result has shape `[product(nc),ncoord]` and uses
    !! direction-1-fastest ordering. An active diagnostic or unset control net
    !! produces an allocated `[0,0]` result.
    pure function get_Xc_all(this) result(Xc)
        class(nurbs_surface), intent(in) :: this
        real(rk), allocatable :: Xc(:,:)

        call get_array(this%Xc,Xc,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return one control point from the flattened control net.
    !!
    !! `n` is the one-based, direction-1-fastest control index. An invalid
    !! index or unavailable control net produces an allocated empty result.
    pure function get_Xci(this, n) result(Xc)
        class(nurbs_surface), intent(in) :: this
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
    !> Return one Cartesian component of one flattened control point.
    !!
    !! The elemental query returns zero when `n` or `dir` is outside
    !! the stored control-net bounds, or when the object has an active
    !! diagnostic.
    pure elemental function get_Xcid(this, n, dir) result(Xc)
        class(nurbs_surface), intent(in) :: this
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
    !> Return an allocated copy of all cached physical samples.
    !!
    !! The result has shape `[product(ng),ncoord]` and direction 1 varies
    !! fastest. An active diagnostic or absent cache produces `[0,0]`.
    pure function get_Xg_all(this) result(Xg)
        class(nurbs_surface), intent(in) :: this
        real(rk), allocatable :: Xg(:,:)

        call get_array(this%Xg,Xg,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return one cached physical sample.
    !!
    !! `n` is the one-based flattened sample index. Invalid or unavailable
    !! data produces an allocated empty result.
    pure function get_Xgi(this, n) result(Xg)
        class(nurbs_surface), intent(in) :: this
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
    !> Return one Cartesian component of one cached physical sample.
    !!
    !! The elemental query returns zero for an invalid sample or coordinate
    !! index and while the object has an active diagnostic.
    pure elemental function get_Xgid(this, n, dir) result(Xg)
        class(nurbs_surface), intent(in) :: this
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
    !> Return an allocated copy of the flattened control weights.
    !!
    !! The result has size `product(nc)`. Polynomial surfaces, unset
    !! weights, and active diagnostics produce an allocated empty result.
    pure function get_Wc_all(this) result(Wc)
        class(nurbs_surface), intent(in) :: this
        real(rk), allocatable :: Wc(:)

        call get_array(this%Wc,Wc,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return one flattened control weight.
    !!
    !! The elemental query returns zero when no weights are stored, `n` is
    !! invalid, or the object has an active diagnostic.
    pure elemental function get_Wci(this, n) result(Wc)
        class(nurbs_surface), intent(in) :: this
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
    !> Return cached parameter samples in direction `dir`.
    !!
    !! Valid directions are 1 and 2. For a Cartesian cache created from
    !! directional vectors, the result contains the corresponding vector.
    !! A row-wise cache created by `put_to_nurbs`, an invalid direction, absent
    !! directional data, or an active diagnostic produces an allocated empty
    !! result.
    pure function get_Xt(this, dir) result(Xt)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), allocatable :: Xt(:)

        select case(dir)
        case(1)
            call get_array(this%Xt1,Xt,this%err%ok)
        case(2)
            call get_array(this%Xt2,Xt,this%err%ok)
        case default
            allocate(Xt(0))
        end select
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the two directional sample counts.
    !!
    !! A Cartesian cache reports its two directional counts. A row-wise cache
    !! reports `[npoint,1]`. The result is zero before sampling and while a
    !! diagnostic is active. The parameter grid is retained after control or
    !! weight edits, but counts remain zero until create() regenerates geometry.
    pure function get_ng(this) result(ng)
        class(nurbs_surface), intent(in) :: this
        integer :: ng(2)

        ng = 0
        if (.not. this%err%ok .or. .not. allocated(this%Xg)) return

        ng = this%ng
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Infer one or both polynomial degrees from the stored spline dimensions.
    !!
    !! When a control count is known, direction \(d\) uses
    !! \(p_d=m_d-n_d-1\), where \(m_d\) is the knot-vector size and \(n_d\)
    !! is the number of control points. Otherwise the leading knot
    !! multiplicity is interpreted as \(p_d+1\). Optional `dir` restricts
    !! the update to direction 1 or 2.
    pure subroutine cmp_degree(this, dir)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: dir
        integer :: first, last, d, lengths(2), values(2)

        if (.not.this%err%ok) return
        first=1
        last=2
        if (present(dir)) then
            if (dir<1 .or. dir>2) then
                call this%err%set(code=100,severity=1,category="forcad_nurbs_surface",location="cmp_degree", &
                    message="Invalid parameter direction.",suggestion="Use dir=1 or dir=2.")
                return
            end if
            first=dir
            last=dir
        end if
        lengths=knot_sizes(this)
        values=this%degree
        do d=first,last
            if (lengths(d)==0) then
                call this%err%set(code=103,severity=1,category="forcad_nurbs_surface",location="cmp_degree", &
                    message="Knot vector is not set.",suggestion="Call set(...) before inferring spline dimensions.")
                return
            end if
            if (this%nc(d)>0) then
                values(d)=lengths(d)-this%nc(d)-1
            else if (d==1) then
                values(d)=compute_multiplicity(this%knot1,this%knot1(1))-1
            else
                values(d)=compute_multiplicity(this%knot2,this%knot2(1))-1
            end if
        end do
        this%degree=values
    end subroutine cmp_degree
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return both polynomial degrees \([p,q]\).
    !!
    !! An active diagnostic produces zeros.
    pure function get_degree_all(this) result(degree)
        class(nurbs_surface), intent(in) :: this
        integer :: degree(2)

        degree = 0
        if (.not. this%err%ok) return

        degree = this%degree
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the polynomial degree in one parameter direction.
    !!
    !! Valid directions are 1 and 2. Invalid input or an active diagnostic
    !! returns zero.
    pure elemental function get_degree_dir(this,dir) result(degree)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        integer :: degree

        degree=0
        if (this%err%ok .and. dir>=1 .and. dir<=size(this%degree)) degree=this%degree(dir)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocated copy of one directional knot vector.
    !!
    !! Valid directions are 1 and 2. Invalid input, unavailable knots, or an
    !! active diagnostic produces an allocated empty result.
    pure function get_knot_all(this, dir) result(knot)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), allocatable :: knot(:)

        select case(dir)
        case(1)
            call get_array(this%knot1,knot,this%err%ok)
        case(2)
            call get_array(this%knot2,knot,this%err%ok)
        case default
            allocate(knot(0))
        end select
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return one knot from a selected parameter direction.
    !!
    !! The elemental query returns zero for an invalid direction or index,
    !! unavailable knots, or an active diagnostic.
    pure elemental function get_knoti(this, dir, i) result(knot)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        integer, intent(in) :: i
        real(rk) :: knot

        knot = 0.0_rk
        if (.not. this%err%ok) return

        if (dir == 1) then
            if (allocated(this%knot1)) then
                if (i >= 1 .and. i <= size(this%knot1)) knot = this%knot1(i)
            end if
        else if (dir == 2) then
            if (allocated(this%knot2)) then
                if (i >= 1 .and. i <= size(this%knot2)) knot = this%knot2(i)
            end if
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Release all owned surface arrays and reset derived state.
    !!
    !! The diagnostic object is deliberately preserved so callers can inspect
    !! an error after cleanup.
    pure elemental subroutine finalize(this)
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
        this%rational = .false.
        this%wrap_parameters = .false.
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build structured quadrilateral connectivity for the control net.
    !!
    !! Optional `p=[p_1,p_2]` selects the visualization cell order;
    !! the default is bilinear. Numbering follows direction-1-fastest storage.
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
    !> Build structured quadrilateral connectivity for cached geometry samples.
    !!
    !! Optional `p=[p_1,p_2]` selects the visualization cell order;
    !! the default is bilinear. Numbering follows direction-1-fastest storage.
    !! An absent physical sample cache returns allocated empty connectivity.
    pure function cmp_elem_Xg_vis(this, p) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), contiguous, optional :: p(:)

        if (.not. this%err%ok .or. .not. allocated(this%Xg)) then
            allocate(elemConn(0,0))
            return
        end if

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
    !> Build quadrilateral connectivity for the active parameter-span grid.
    !!
    !! Grid vertices are the distinct active knots in both directions.
    !! Optional `p` selects cell order and defaults to `[1,1]`.
    pure function cmp_elem_Xth(this, p) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), contiguous, optional :: p(:)

        if (.not. this%err%ok) return

        if (present(p)) then
            elemConn = elemConn_C0(active_span_count(this%knot1, this%nc(1), this%degree(1)) + 1, &
                active_span_count(this%knot2, this%nc(2), this%degree(2)) + 1, p(1), p(2))
        else
            elemConn = elemConn_C0(active_span_count(this%knot1, this%nc(1), this%degree(1)) + 1, &
                active_span_count(this%knot2, this%nc(2), this%degree(2)) + 1, 1, 1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export the control net as legacy VTK quadrilateral cells.
    !!
    !! Stored visualization connectivity is used when available; otherwise
    !! bilinear connectivity is generated. Optional `point_data` is
    !! point-associated data with matching `field_names`. `encoding`
    !! is forwarded to the VTK writer.
    impure subroutine export_Xc(this, filename, point_data, field_names, encoding)
        class(nurbs_surface), intent(inout) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), contiguous, optional :: point_data(:,:)
        character(len=*), intent(in), contiguous, optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        if (.not.allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
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

        call export_vtk_legacy(filename=filename, points=this%Xc, elemConn=elemConn, vtkCellType=9, &
                               point_data=point_data, field_names=field_names, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export cached physical samples as legacy VTK quadrilateral cells.
    !!
    !! Call [[nurbs_surface:create]] before this operation. Stored visualization
    !! connectivity is used when available; otherwise bilinear connectivity is
    !! generated.
    impure subroutine export_Xg(this, filename, point_data, field_names, encoding)
        class(nurbs_surface), intent(inout) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), contiguous, optional :: point_data(:,:)
        character(len=*), intent(in), contiguous, optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        if (.not.allocated(this%Xg)) then
            call this%err%set(&
                code       = 104,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
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

        call export_vtk_legacy(filename=filename, points=this%Xg, elemConn=elemConn, vtkCellType=9, &
                               point_data=point_data, field_names=field_names, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export the active tensor-product parameter grid in parameter space.
    !!
    !! Vertices are pairs of distinct active knots and quadrilateral cells span
    !! adjacent values. The written coordinates are \((u,v)\), not physical
    !! surface points.
    impure subroutine export_Xth(this, filename, point_data, field_names, encoding)
        class(nurbs_surface), intent(in) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), contiguous, optional :: point_data(:,:)
        character(len=*), intent(in), contiguous, optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Xth(:,:), Xth1(:), Xth2(:)

        if (.not. this%err%ok) return

        elemConn = this%cmp_elem_Xth()
        Xth1 = active_knots(this%knot1, this%nc(1), this%degree(1))
        Xth2 = active_knots(this%knot2, this%nc(2), this%degree(2))
        call ndgrid(Xth1, Xth2, Xth)

        call export_vtk_legacy(filename=filename, points=Xth, elemConn=elemConn, vtkCellType=9, &
                               point_data=point_data, field_names=field_names, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export active isoparametric grid lines embedded in physical space.
    !!
    !! Every distinct active knot in one direction defines a line sampled over
    !! all active spans of the other direction. Optional `res` sets the
    !! minimum number of points per span and defaults to 10.
    impure subroutine export_Xth_in_Xg(this, filename, res, encoding)
        class(nurbs_surface), intent(inout) :: this
        character(len=*),     intent(in) :: filename
        integer, intent(in),  optional   :: res
        character(len=*), intent(in), optional :: encoding

        real(rk), allocatable :: Xg_all(:,:)
        integer, allocatable :: elemConn(:,:)

        call this%cmp_Xth_in_Xg(&
            Xg       = Xg_all,&
            elemConn = elemConn,&
            res      = res)
        if (.not. this%err%ok) return

        call export_vtk_legacy(&
            filename    = filename,&
            points      = Xg_all,&
            elemConn    = elemConn,&
            vtkCellType = 4,&
            encoding    = encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build active isoparametric grid lines embedded in physical space.
    !!
    !! The returned rows of `elemConn` are VTK polyline connectivities into
    !! `Xg`. No file is written.
    impure subroutine cmp_Xth_in_Xg(this, Xg, elemConn, res)
        class(nurbs_surface), intent(inout) :: this
        real(rk), allocatable, intent(out) :: Xg(:,:)
        integer, allocatable, intent(out) :: elemConn(:,:)
        integer, intent(in), optional :: res
        real(rk), allocatable :: Xt_all(:,:)

        if (.not.this%err%ok) return
        if (.not.allocated(this%Xc)) then
            call this%err%set(code=102,severity=1,category="forcad_nurbs_surface",location="cmp_Xth_in_Xg", &
                message="Control points are not set.",suggestion="Call set(...) first before exporting.")
            return
        end if
        if (.not.allocated(this%knot1) .or. .not.allocated(this%knot2)) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_surface",location="cmp_Xth_in_Xg", &
                message="Knot vector is not set.",suggestion="Call set(...) first before exporting.")
            return
        end if
        if (size(this%Xc,2)<2 .or. size(this%Xc,2)>3) then
            call this%err%set(code=100,severity=1,category="forcad_nurbs_surface",location="cmp_Xth_in_Xg", &
                message="Invalid geometry dimension.",suggestion="Use two or three physical coordinates for grid export.")
            return
        end if
        if (.not.this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_surface",location="cmp_Xth_in_Xg", &
                message="Invalid spline geometry.",suggestion="Set valid knots, control points and positive weights first.")
            return
        end if
        call isoparametric_grid(this%knot1,this%knot2,this%nc,this%degree,Xt_all,elemConn,res)
        Xg=this%cmp_Xg(Xt_all)
    end subroutine cmp_Xth_in_Xg
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export the current surface as IGES entity 128.
    !!
    !! Both knot vectors, degrees, control-net ordering, weights, and active
    !! parameter ranges are transferred to a Rational B-Spline Surface entity.
    !! One- and two-dimensional embeddings are padded with zero coordinates;
    !! for embeddings above three dimensions, only the first three coordinates
    !! are written because entity 128 stores three Cartesian coordinates.
    !! Closed and periodic properties are classified independently of parameter
    !! wrapping. A direction is closed only when its two boundary traces agree
    !! throughout every positive tangential knot span. On each span, the
    !! comparison uses \(2p_t+1\) interior points: after clearing the positive
    !! rational denominators, each coordinate difference has degree at most
    !! \(2p_t\), so these values determine equality of the two trace pieces.
    !! A direction is periodic only when it is closed and its final \(p\)
    !! control layers repeat the first \(p\) layers while its knot sequence
    !! repeats after \(n_c-p\) intervals. Rational periodicity also requires
    !! repeated weights. The polynomial/rational property follows
    !! [[nurbs_surface:is_rational]], including its roundoff-scale uniform-weight
    !! classification.
    !!
    !! The routine requires initialized surface data. File-system and writer
    !! failures are delegated to `forIGES` and are not converted to
    !! [[nurbs_surface:err]].
    !!
    !! **Format reference:** U.S. National Bureau of Standards, *Initial
    !! Graphics Exchange Specification (IGES), Version 4.0*, NBSIR 88-3813
    !! (1988), Rational B-Spline Surface Entity 128.
    !! [NIST publication](https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir88-3813.pdf).
    impure subroutine export_iges(this, filename)
        use forIGES, only: Gsection_t, Dentry_t, entity128_t, DElist_t, PElist_t,&
                           makeSsection, makeGsection, makeDPsections, writeIGESfile, wp

        class(nurbs_surface), intent(inout) :: this
            !! Valid initialized surface. Existing diagnostic errors suppress export.
        character(len=*),     intent(in)    :: filename
            !! Nonempty destination path for the IGES text file.

        type(Gsection_t)  :: G
        type(Dentry_t)    :: D
        type(entity128_t) :: surf128
        type(DElist_t)    :: Dlist
        type(PElist_t)    :: Plist
        character(80), allocatable :: Ssection(:), Gsection(:), Dsection(:), Psection(:), Ssec_out(:)
        real(wp), allocatable :: X(:,:), Y(:,:), Z(:,:), W(:,:), S(:), T(:)
        integer :: i, j, idx, dim_, dim_export
        integer :: sample, span, sample_count
        integer :: K1, K2, M1, M2, N1, N2, prop3
        real(wp) :: U(0:1), V(0:1)
        real(rk), allocatable :: Xg_start(:), Xg_end(:)
        real(rk) :: limits(2,2), Xt_boundary(2)
        real(rk) :: geometry_tol
        real(rk) :: span_midpoint, span_radius, pi_value
        logical :: closed(2), periodic(2)

        if (.not. this%err%ok) return

        ! Parameters consistent with the IGES definition
        K1 = this%nc(1) - 1
        K2 = this%nc(2) - 1
        M1 = this%degree(1)
        M2 = this%degree(2)
        dim_ = size(this%Xc, 2)
        dim_export = min(dim_, 3)

        ! Compute required N1 and N2 based on IGES standard
        N1 = 1 + K1 - M1
        N2 = 1 + K2 - M2

        allocate(S(-M1:K1+1), T(-M2:K2+1))
        allocate(X(0:K1,0:K2), Y(0:K1,0:K2), Z(0:K1,0:K2), W(0:K1,0:K2))

        ! Copy knots explicitly, matching IGES indexing exactly
        do i = -M1, K1 + 1
            S(i) = real(this%knot1(i + M1 + 1), kind=wp)
        end do

        do i = -M2, K2 + 1
            T(i) = real(this%knot2(i + M2 + 1), kind=wp)
        end do

        ! Correctly map control points and weights
        X = real(0.0_rk, kind=wp)
        Y = real(0.0_rk, kind=wp)
        Z = real(0.0_rk, kind=wp)
        if (this%is_rational()) then
            do j = 0, K2
                do i = 0, K1
                    idx = j * this%nc(1) + i + 1
                    if (dim_ >= 1) X(i,j) = real(this%Xc(idx,1), kind=wp)
                    if (dim_ >= 2) Y(i,j) = real(this%Xc(idx,2), kind=wp)
                    if (dim_ >= 3) Z(i,j) = real(this%Xc(idx,3), kind=wp)
                    W(i,j) = real(this%Wc(idx), kind=wp)
                end do
            end do
            prop3 = 0
        else
            do j = 0, K2
                do i = 0, K1
                    idx = j * this%nc(1) + i + 1
                    if (dim_ >= 1) X(i,j) = real(this%Xc(idx,1), kind=wp)
                    if (dim_ >= 2) Y(i,j) = real(this%Xc(idx,2), kind=wp)
                    if (dim_ >= 3) Z(i,j) = real(this%Xc(idx,3), kind=wp)
                    W(i,j) = real(1.0_rk, kind=wp)
                end do
            end do
            prop3 = 1
        end if

        limits(1,:) = [&
            knot_start(this%knot1, this%nc(1), this%degree(1)),&
            knot_end(this%knot1, this%nc(1), this%degree(1))]
        limits(2,:) = [&
            knot_start(this%knot2, this%nc(2), this%degree(2)),&
            knot_end(this%knot2, this%nc(2), this%degree(2))]
        U = real(limits(1,:), kind=wp)
        V = real(limits(2,:), kind=wp)

        allocate(Xg_start(dim_), Xg_end(dim_))
        geometry_tol = 64.0_rk*real(M1+1, rk)*real(M2+1, rk)*epsilon(1.0_rk)*&
            max(1.0_rk, maxval(abs(this%Xc(:,1:dim_export))))
        pi_value = acos(-1.0_rk)
        closed = .true.

        sample_count = 2*M2 + 1
        do span = M2 + 1, this%nc(2)
            if (this%knot2(span+1) <= this%knot2(span)) cycle
            span_midpoint = 0.5_rk*(this%knot2(span) + this%knot2(span+1))
            span_radius = 0.5_rk*(this%knot2(span+1) - this%knot2(span))
            do sample = 1, sample_count
                Xt_boundary = [&
                    limits(1,1),&
                    span_midpoint + span_radius*cos(&
                        real(2*sample-1, rk)*pi_value/real(2*sample_count, rk))]
                if (prop3 == 0) then
                    Xg_start = compute_Xg(&
                        Xt_boundary,&
                        this%knot1,&
                        this%knot2,&
                        this%degree,&
                        this%nc,&
                        this%Xc,&
                        Wc=this%Wc)
                else
                    Xg_start = compute_Xg(&
                        Xt_boundary,&
                        this%knot1,&
                        this%knot2,&
                        this%degree,&
                        this%nc,&
                        this%Xc)
                end if
                Xt_boundary(1) = limits(1,2)
                if (prop3 == 0) then
                    Xg_end = compute_Xg(&
                        Xt_boundary,&
                        this%knot1,&
                        this%knot2,&
                        this%degree,&
                        this%nc,&
                        this%Xc,&
                        Wc=this%Wc)
                else
                    Xg_end = compute_Xg(&
                        Xt_boundary,&
                        this%knot1,&
                        this%knot2,&
                        this%degree,&
                        this%nc,&
                        this%Xc)
                end if
                if (any(abs(Xg_start(1:dim_export) - Xg_end(1:dim_export)) > geometry_tol)) then
                    closed(1) = .false.
                    exit
                end if
            end do
            if (.not. closed(1)) exit
        end do

        sample_count = 2*M1 + 1
        do span = M1 + 1, this%nc(1)
            if (this%knot1(span+1) <= this%knot1(span)) cycle
            span_midpoint = 0.5_rk*(this%knot1(span) + this%knot1(span+1))
            span_radius = 0.5_rk*(this%knot1(span+1) - this%knot1(span))
            do sample = 1, sample_count
                Xt_boundary = [&
                    span_midpoint + span_radius*cos(&
                        real(2*sample-1, rk)*pi_value/real(2*sample_count, rk)),&
                    limits(2,1)]
                if (prop3 == 0) then
                    Xg_start = compute_Xg(&
                        Xt_boundary,&
                        this%knot1,&
                        this%knot2,&
                        this%degree,&
                        this%nc,&
                        this%Xc,&
                        Wc=this%Wc)
                else
                    Xg_start = compute_Xg(&
                        Xt_boundary,&
                        this%knot1,&
                        this%knot2,&
                        this%degree,&
                        this%nc,&
                        this%Xc)
                end if
                Xt_boundary(2) = limits(2,2)
                if (prop3 == 0) then
                    Xg_end = compute_Xg(&
                        Xt_boundary,&
                        this%knot1,&
                        this%knot2,&
                        this%degree,&
                        this%nc,&
                        this%Xc,&
                        Wc=this%Wc)
                else
                    Xg_end = compute_Xg(&
                        Xt_boundary,&
                        this%knot1,&
                        this%knot2,&
                        this%degree,&
                        this%nc,&
                        this%Xc)
                end if
                if (any(abs(Xg_start(1:dim_export) - Xg_end(1:dim_export)) > geometry_tol)) then
                    closed(2) = .false.
                    exit
                end if
            end do
            if (.not. closed(2)) exit
        end do

        if (prop3 == 0) then
            periodic(1) = closed(1) .and. periodic_topology(&
                this%knot1, this%degree(1), this%Xc, this%Wc, this%nc, 1)
            periodic(2) = closed(2) .and. periodic_topology(&
                this%knot2, this%degree(2), this%Xc, this%Wc, this%nc, 2)
        else
            periodic(1) = closed(1) .and. periodic_topology(&
                this%knot1, this%degree(1), this%Xc, nc=this%nc, dir=1)
            periodic(2) = closed(2) .and. periodic_topology(&
                this%knot2, this%degree(2), this%Xc, nc=this%nc, dir=2)
        end if

        ! Initialize IGES entity 128 (Rational B-spline Surface).
        surf128%entity_type = 128
        surf128%DEP         = 1
        surf128%form        = 0
        surf128%K1          = K1
        surf128%K2          = K2
        surf128%M1          = M1
        surf128%M2          = M2
        surf128%PROP1       = merge(1, 0, closed(1))
        surf128%PROP2       = merge(1, 0, closed(2))
        surf128%PROP3       = prop3
        surf128%PROP4       = merge(1, 0, periodic(1))
        surf128%PROP5       = merge(1, 0, periodic(2))
        allocate(surf128%S(lbound(S, 1):ubound(S, 1)))
        allocate(surf128%T(lbound(T, 1):ubound(T, 1)))
        allocate(surf128%W(0:K1,0:K2), surf128%X(0:K1,0:K2), surf128%Y(0:K1,0:K2), surf128%Z(0:K1,0:K2))
        surf128%S = S
        surf128%T = T
        surf128%W = W
        surf128%X = X
        surf128%Y = Y
        surf128%Z = Z
        surf128%U = U
        surf128%V = V

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
    !> Drop derived geometry, optionally including topology, while keeping parameters.
    pure subroutine invalidate_samples(this,space_changed)
        class(nurbs_surface), intent(inout) :: this
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

    !> Replace one finite control coordinate; num and dir are one-based.
    !! Physical samples are cleared; create() reuses the retained parameter grid.
    pure subroutine modify_Xc_scalar(this, X, num, dir)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: X
        integer, intent(in) :: num, dir

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) then
            call this%err%set(code=102, severity=1, category="forcad_nurbs_surface", &
                message="Control points are not set.", location="modify_Xc", &
                suggestion="Call set(...) before modifying control coordinates.")
            return
        end if
        if (num<1 .or. num>size(this%Xc,1) .or. dir<1 .or. dir>size(this%Xc,2) .or. &
            .not. ieee_is_finite(X)) then
            call this%err%set(code=100, severity=1, category="forcad_nurbs_surface", &
                message="Invalid control index, coordinate direction, or nonfinite coordinate.", location="modify_Xc", &
                suggestion="Use valid one-based indices and a finite coordinate.")
            return
        end if
        this%Xc(num,dir) = X
        call invalidate_samples(this)
    end subroutine modify_Xc_scalar

    !> Atomically replace a finite control net of the same shape.
    !! Physical samples are cleared; knots and control connectivity remain valid.
    !! create() reuses the retained parameter grid.
    pure subroutine modify_Xc_all(this, Xc)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: Xc(:,:)

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) then
            call this%err%set(code=102, severity=1, category="forcad_nurbs_surface", &
                message="Control points are not set.", location="modify_Xc", &
                suggestion="Call set(...) before replacing control coordinates.")
            return
        end if
        if (any(shape(Xc)/=shape(this%Xc)) .or. .not. all(ieee_is_finite(Xc))) then
            call this%err%set(code=101, severity=1, category="forcad_nurbs_surface", &
                message="Replacement controls must retain their shape and be finite.", location="modify_Xc", &
                suggestion="Supply a finite Xc array with the current control-net shape.")
            return
        end if
        this%Xc = Xc
        call invalidate_samples(this)
    end subroutine modify_Xc_all

    !> Replace one positive weight; polynomial surfaces initially have unit weights.
    !! Physical samples are cleared; create() reuses the retained parameter grid.
    pure subroutine modify_Wc_scalar(this, W, num)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: W
        integer, intent(in) :: num

        if (.not. this%err%ok) return
        if (num<1 .or. num>product(this%nc)) then
            call this%err%set(code=100, severity=1, category="forcad_nurbs_surface", &
                message="Invalid weight index.", location="modify_Wc", &
                suggestion="Use 1 <= num <= product(nc).")
            return
        end if
        if (.not. ieee_is_finite(W) .or. W<=0.0_rk) then
            call this%err%set(code=105, severity=1, category="forcad_nurbs_surface", &
                message="Weights must be finite and positive.", location="modify_Wc", &
                suggestion="Provide a strictly positive finite weight.")
            return
        end if
        if (.not. allocated(this%Wc)) allocate(this%Wc(product(this%nc)),source=1.0_rk)
        this%Wc(num) = W
        this%rational = nonuniform_weights(this%Wc)
        call invalidate_samples(this)
    end subroutine modify_Wc_scalar

    !> Atomically replace all positive weights at fixed control counts.
    !! Physical samples are cleared; create() reuses the retained parameter grid.
    pure subroutine modify_Wc_all(this, Wc)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: Wc(:)

        if (.not. this%err%ok) return
        if (size(Wc)/=product(this%nc) .or. .not. valid_weights(Wc)) then
            call this%err%set(code=105, severity=1, category="forcad_nurbs_surface", &
                message="Replacement weights must match the controls and be finite and positive.", location="modify_Wc", &
                suggestion="Supply one strictly positive finite weight for each stored control.")
            return
        end if
        this%Wc = Wc
        this%rational = nonuniform_weights(this%Wc)
        call invalidate_samples(this)
    end subroutine modify_Wc_all


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return multiplicities of the distinct knots in one direction.
    !!
    !! Result entries correspond, in nondecreasing order, to the unique values
    !! in the complete selected knot vector. Invalid input or unavailable state
    !! produces an allocated empty result.
    pure function get_multiplicity(this, dir) result(m)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        integer, allocatable :: m(:)

        if (this%err%ok) then
            if (dir==1 .and. allocated(this%knot1)) m=compute_multiplicity(this%knot1)
            if (dir==2 .and. allocated(this%knot2)) m=compute_multiplicity(this%knot2)
        end if
        if (.not.allocated(m)) allocate(m(0))
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return \(p_d-s_i\) at every distinct stored directional knot.
    !!
    !! For degree \(p_d\) and knot multiplicity \(s_i\), entry \(i\) is
    !! \(c_i=p_d-s_i\). At interior active knots this is guaranteed spline-space
    !! continuity; a particular surface can be smoother. Endpoint entries and
    !! entries outside the active interval describe representation multiplicity,
    !! not continuity between two active spans.
    pure function get_continuity(this, dir) result(c)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        integer, allocatable :: c(:)

        c=this%get_multiplicity(dir)
        c=this%get_degree(dir)-c
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Infer one or both directional control counts from knots and degree.
    !!
    !! Each selected direction uses \(n_d=m_d-p_d-1\), where \(m_d\) is the
    !! knot-vector size. Optional `dir` must be 1 or 2.
    pure subroutine cmp_nc(this, dir)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: dir
        integer :: first, last, d, lengths(2), values(2)

        if (.not.this%err%ok) return
        first=1
        last=2
        if (present(dir)) then
            if (dir<1 .or. dir>2) then
                call this%err%set(code=103,severity=1,category="forcad_nurbs_surface",location="cmp_nc", &
                    message="Invalid parameter direction.",suggestion="Use dir=1 or dir=2.")
                return
            end if
            first=dir
            last=dir
        end if
        lengths=knot_sizes(this)
        values=this%nc
        do d=first,last
            if (lengths(d)==0) then
                call this%err%set(code=103,severity=1,category="forcad_nurbs_surface",location="cmp_nc", &
                    message="Knot vector is not set.",suggestion="Call set(...) before inferring spline dimensions.")
                return
            end if
            values(d)=lengths(d)-this%degree(d)-1
        end do
        this%nc=values
    end subroutine cmp_nc
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return control-point counts \([n_1,n_2]\).
    !!
    !! An active diagnostic produces zeros.
    pure function get_nc_all(this) result(nc)
        class(nurbs_surface), intent(in) :: this
        integer :: nc(2)

        nc = 0
        if (.not. this%err%ok) return

        nc = this%nc
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the control-point count in one direction.
    !!
    !! Invalid input, unavailable knots, or an active diagnostic returns zero.
    pure elemental function get_nc_dir(this, dir) result(nc)
        class(nurbs_surface), intent(in) :: this
        integer, intent(in) :: dir
        integer :: nc

        nc = 0
        if (.not. this%err%ok) return

        if (dir == 1) then

            if (allocated(this%knot1)) nc = this%nc(1)

        else if (dir == 2) then

            if (allocated(this%knot2)) nc = this%nc(2)
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Validate the state required by all surface derivative evaluators.
    pure subroutine check_derivative_state(this, location, valid, side)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: side(:)
        character(len=*), intent(in) :: location
        logical, intent(out) :: valid

        valid = .false.
        if (.not. this%err%ok) return
        if (.not. valid_trace_side(side,2)) then
            call this%err%set(code=100,severity=1,category="forcad_nurbs_surface",location=location, &
                message="Invalid knot trace selector.",suggestion="Use -1, 0, or 1 for each parameter direction.")
            return
        end if
        if (.not. allocated(this%knot1) .or. .not. allocated(this%knot2)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Surface derivatives require two valid knot vectors and control counts.",&
                location   = location,&
                suggestion = "Set the surface before evaluating basis derivatives.")
            return
        end if
        valid = all(this%degree >= 0) .and. all(this%nc > this%degree) .and. &
            size(this%knot1) == this%nc(1) + this%degree(1) + 1 .and. &
            size(this%knot2) == this%nc(2) + this%degree(2) + 1
        if (valid .and. allocated(this%Wc)) then
            valid = size(this%Wc) == product(this%nc)
        end if
        if (.not. valid) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Surface derivative state is inconsistent.",&
                location   = location,&
                suggestion = "Use matching finite knots, degrees, control counts, and positive weights.")
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> Validate stored basis state and the two active scalar parameter entries.
    pure subroutine check_scalar_query(this, Xt, location, valid, side)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: side(:)
        real(rk), intent(in) :: Xt(:)
        character(len=*), intent(in) :: location
        logical, intent(out) :: valid

        call check_derivative_state(this,location,valid,side)
        if (.not. valid) return
        valid = size(Xt)>=2
        if (valid) valid = all(ieee_is_finite(Xt(1:2)))
        if (.not. valid) call this%err%set(code=107,severity=1,category="forcad_nurbs_surface", &
            message="Surface queries require at least two finite parameter entries.",location=location, &
            suggestion="Pass a finite pair Xt=[u,v]; later entries are ignored.")
    end subroutine check_scalar_query
    !===============================================================================


    !===============================================================================
    !> Select paired cached rows or a directional grid without changing samples.
    pure subroutine prepare_surface_query(this, location, parameters, query_ng, &
        valid, res1, res2, Xt1, Xt2, side)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), optional :: side(:)
        character(len=*), intent(in) :: location
        real(rk), allocatable, intent(out) :: parameters(:,:)
        integer, intent(out) :: query_ng(2)
        logical, intent(out) :: valid
        integer, intent(in), optional :: res1, res2
        real(rk), intent(in), optional :: Xt1(:), Xt2(:)
        real(rk), allocatable :: axis1(:), axis2(:)
        logical :: valid1, valid2

        query_ng = 0
        call check_derivative_state(this,location,valid,side)
        if (.not. valid) return
        if (.not. present(res1) .and. .not. present(res2) .and. &
            .not. present(Xt1) .and. .not. present(Xt2) .and. allocated(this%Xt)) then
            parameters = this%Xt
            query_ng = this%ng
            return
        end if
        call query_parameter_axis(this%knot1, this%nc(1), this%degree(1), this%Xt1, &
            axis1, valid1, res1, Xt1)
        call query_parameter_axis(this%knot2, this%nc(2), this%degree(2), this%Xt2, &
            axis2, valid2, res2, Xt2)
        valid = valid1 .and. valid2
        if (.not. valid) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Surface queries require finite parameters or positive resolutions in both directions.",&
                location   = location,&
                suggestion = "Supply Xt1/Xt2 or positive res1/res2, or create a sample before omitting parameters.")
            return
        end if
        query_ng = [size(axis1),size(axis2)]
        call ndgrid(axis1, axis2, parameters)
    end subroutine prepare_surface_query




    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate one arbitrary mixed surface basis derivative on its local support.
    pure subroutine evaluate_basis_order(this,Xt,order,values,Wc,side,first_active,elem)
        class(nurbs_surface), intent(in) :: this
        real(rk), intent(in) :: Xt(:)
        integer, intent(in) :: order(:)
        real(rk), intent(out) :: values(:)
        real(rk), contiguous, intent(in), optional :: Wc(:)
        integer, intent(in), optional :: side(:)
        integer, intent(out), optional :: first_active(:)
        integer, contiguous, intent(in), optional :: elem(:)

        call tensor_basis_derivative_at(Xt,this%nc,this%degree,this%knot1,order,this%wrap_parameters,values, &
            first_active,elem,Wc,knot2=this%knot2,side=side)
    end subroutine evaluate_basis_order
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute active mixed derivatives on an explicit surface parameter grid.
    !! Grid points use
    !! \(g=i_1+(i_2-1)\operatorname{size}(Xt1)\). For local offsets
    !! \(a_d=0,\ldots,p_d\), define
    !!
    !! \[
    !! \ell=1+a_1+(p_1+1)a_2,\qquad
    !! A=(f_2+a_2-1)n_{c,1}+f_1+a_1,
    !! \]
    !!
    !! where \(f_d=\mathtt{first\_active(d,g)}\). Then
    !! `Dgc(ell,g)` is
    !! \(\partial^{k_1+k_2}R_A/(\partial u^{k_1}\partial v^{k_2})\)
    !! for `order=[k1,k2]`. Wrapping is applied directionally. If an unwrapped
    !! coordinate is outside its active interval, the derivative column is zero
    !! and the corresponding first-active index is 1.
    pure subroutine derivative_order_active_vector(this, Xt1, Xt2, order, first_active, Dgc, side)
        class(nurbs_surface), target, intent(inout) :: this
        integer, intent(in), optional :: side(:) !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in), contiguous :: Xt1(:), Xt2(:)
            !! Finite directional parameter vectors defining a Cartesian grid.
        integer, intent(in) :: order(2)
            !! Nonnegative mixed-derivative orders `[du,dv]`.
        integer, allocatable, intent(out) :: first_active(:,:)
            !! First active control-point indices `[2,npoint]`.
        real(rk), allocatable, intent(out) :: Dgc(:,:)
            !! Local derivatives `[product(degree+1),size(Xt1)*size(Xt2)]`.
        integer :: point1, point2
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call check_derivative_state(this, "derivative_order_active_vector", valid, side)
        if (.not. valid) then
            allocate(first_active(0,0))
            allocate(Dgc(0,0))
            return
        end if
        if (any(order < 0) .or. .not. all(ieee_is_finite(Xt1)) .or. .not. all(ieee_is_finite(Xt2))) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Active surface derivative inputs must be finite with nonnegative orders.",&
                location   = "derivative_order_active_vector",&
                suggestion = "Use finite parameter grids and order=[du,dv] with nonnegative entries.")
            allocate(first_active(0,0))
            allocate(Dgc(0,0))
            return
        end if

        allocate(Dgc(product(this%degree+1),size(Xt1)*size(Xt2)), source=0.0_rk)
        allocate(first_active(2,size(Xt1)*size(Xt2)), source=0)
        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
#if defined(__NVCOMPILER)
        do point2 = 1, size(Xt2)
            do point1 = 1, size(Xt1)
#else
        do concurrent (point2 = 1:size(Xt2), point1 = 1:size(Xt1))
#endif
            call evaluate_basis_order(this,[Xt1(point1),Xt2(point2)],order,Dgc(:,point1+(point2-1)*size(Xt1)),weights,side, &
                first_active=first_active(:,point1+(point2-1)*size(Xt1)))
#if defined(__NVCOMPILER)
            end do
#endif
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one active mixed surface derivative.
    !! Local and global control indices use the same direction-1-fastest
    !! formulas as the vector overload, with one parameter pair.
    pure subroutine derivative_order_active_scalar(this, Xt, order, first_active, Dgc, side)
        class(nurbs_surface), target, intent(inout) :: this
        integer, intent(in), optional :: side(:) !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in), contiguous :: Xt(:)
            !! Finite parameter pair `[u,v]`.
        integer, intent(in) :: order(2)
            !! Nonnegative mixed-derivative orders `[du,dv]`.
        integer, intent(out) :: first_active(2)
            !! First active control-point index in each direction.
        real(rk), allocatable, intent(out) :: Dgc(:)
            !! Local derivatives `[product(degree+1)]`.
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        first_active = 0
        call check_derivative_state(this, "derivative_order_active_scalar", valid, side)
        if (.not. valid) then
            allocate(Dgc(0))
            return
        end if
        if (size(Xt) /= 2 .or. any(order < 0)) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Active surface derivative input requires two finite parameters and nonnegative orders.",&
                location   = "derivative_order_active_scalar",&
                suggestion = "Use Xt=[u,v] and order=[du,dv] with nonnegative entries.")
            allocate(Dgc(0))
            return
        end if
        if (.not. all(ieee_is_finite(Xt))) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Active surface derivative parameters must be finite.",&
                location   = "derivative_order_active_scalar",&
                suggestion = "Use finite values for Xt=[u,v].")
            allocate(Dgc(0))
            return
        end if

        allocate(Dgc(product(this%degree+1)), source=0.0_rk)
        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
        call evaluate_basis_order(this,Xt,order,Dgc,weights,side, &
            first_active=first_active)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one requested mixed derivative on a tensor grid of surface parameters.
    !! For `order=[k1,k2]`, the dense matrix stores
    !!
    !! \[
    !! \mathtt{Dgc(g,A)}=
    !! \frac{\partial^{k_1+k_2}R_A}
    !!      {\partial u^{k_1}\partial v^{k_2}},\qquad
    !! g=i_1+(i_2-1)n_{g,1}.
    !! \]
    !!
    !! Explicit directional vectors take precedence over positive resolutions;
    !! omitted directions reuse cached sample parameters. Query parameters and
    !! counts are local, so stored samples remain unchanged. With no arguments,
    !! cached paired rows are reused. Nonfinite explicit
    !! parameters and unavailable omitted directions produce a diagnostic.
    pure subroutine derivative_order_vector(this, order, Dgc, res1, res2, Xt1, Xt2, side)
        class(nurbs_surface), target, intent(inout) :: this
        integer, intent(in), optional :: side(:) !! Trace per direction: -1 left, 0 default, +1 right.
        integer, intent(in) :: order(2)
            !! Nonnegative mixed-derivative orders `[du,dv]`.
        real(rk), allocatable, intent(out) :: Dgc(:,:)
            !! Dense derivative matrix `[product(ng),product(nc)]`.
        integer, intent(in), optional :: res1, res2
            !! Optional positive generated counts by direction.
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:)
            !! Optional explicit directional parameter vectors.
        real(rk), allocatable :: parameters(:,:)
        integer :: point, query_ng(2)
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call check_derivative_state(this, "derivative_order_vector", valid, side)
        if (.not. valid) then
            allocate(Dgc(0,0))
            return
        end if
        if (any(order < 0)) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Surface derivative orders must be nonnegative.",&
                location   = "derivative_order_vector",&
                suggestion = "Use order=[du,dv] with du >= 0 and dv >= 0.")
            allocate(Dgc(0,0))
            return
        end if

        call prepare_surface_query(this, "derivative_order_vector", parameters, query_ng, &
            valid, res1, res2, Xt1, Xt2)
        if (.not. valid) then
            allocate(Dgc(0,0))
            return
        end if

        allocate(Dgc(product(query_ng),product(this%nc)), source=0.0_rk)
        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
#if defined(__NVCOMPILER)
        ! NVFORTRAN 26.3 cannot lower the pure local kernel for stdpar GPU execution.
        do point = 1, size(parameters,1)
#else
        do concurrent (point = 1:size(parameters,1))
#endif
            call evaluate_basis_order(this,parameters(point,:),order,Dgc(point,:),weights,side)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one requested mixed derivative at one surface parameter pair.
    !! Without `elem`, the result follows direction-1-fastest global control
    !! order. With `elem`, entries follow the supplied one-based flattened
    !! indices and an out-of-range index returns zero. A finite two-entry `Xt`
    !! is a caller precondition; this overload diagnoses its length and the
    !! derivative orders, but not nonfinite parameter values.
    pure subroutine derivative_order_scalar(this, Xt, order, Dgc, elem, side)
        class(nurbs_surface), target, intent(inout) :: this
        integer, intent(in), optional :: side(:) !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter pair `[u,v]`.
        integer, intent(in) :: order(2)
            !! Nonnegative mixed-derivative orders `[du,dv]`.
        real(rk), allocatable, intent(out) :: Dgc(:)
            !! Dense `[product(nc)]` or selected `[size(elem)]` derivative vector.
        integer, intent(in), contiguous, optional :: elem(:)
            !! Optional ordered one-based flattened control indices.
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call check_derivative_state(this, "derivative_order_scalar", valid, side)
        if (.not. valid) then
            allocate(Dgc(0))
            return
        end if
        if (size(Xt) /= 2 .or. any(order < 0)) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Surface derivative input must contain two parameters and nonnegative orders.",&
                location   = "derivative_order_scalar",&
                suggestion = "Use Xt=[u,v] and order=[du,dv] with nonnegative derivative orders.")
            allocate(Dgc(0))
            return
        end if

        if (present(elem)) then
            allocate(Dgc(size(elem)), source=0.0_rk)
        else
            allocate(Dgc(product(this%nc)), source=0.0_rk)
        end if
        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
        call evaluate_basis_order(this,Xt,order,Dgc,weights,side,elem=elem)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate first parametric derivatives on a tensor-product grid.
    !!
    !! Explicit directional parameter vectors take precedence over resolutions;
    !! omitted directions reuse their caches. With no arguments, cached paired
    !! rows are reused. Flattened grid rows follow
    !! direction-1-fastest `ndgrid` ordering. Query parameters and counts are
    !! local, preserving all cached samples. Missing parameter vectors, invalid
    !! resolutions, and nonfinite explicit parameters produce a diagnostic.
    pure subroutine derivative_vector(this, res1, res2, Xt1, Xt2, dTgc, Tgc, side)
        class(nurbs_surface), target, intent(inout) :: this
        integer, intent(in), optional :: side(:)
        integer, intent(in), optional :: res1, res2
            !! Optional positive uniform-grid sizes in directions 1 and 2.
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:)
            !! Optional directional parameter vectors.
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
            !! Parametric gradients `[ng_total,nc_total,2]`, ordered as \((u,v)\).
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
            !! Optional basis matrix `[ng_total,nc_total]`.
        call query_basis_grid(this, "derivative_vector", res1, res2, Xt1, Xt2, side, Tgc=Tgc, dTgc=dTgc)
    end subroutine derivative_vector
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the parametric basis gradient at one parameter pair.
    !! Optional `elem` selects and orders a local control subset. `Xt` must have
    !! at least two finite entries; later entries are ignored. The spline must
    !! be initialized. This convenience overload does not diagnose those
    !! preconditions.
    pure subroutine derivative_scalar(this, Xt, dTgc, Tgc, elem, side)
        class(nurbs_surface), target, intent(inout) :: this
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter pair \([u,v]\).
        integer, intent(in), optional :: elem(:)
            !! Optional one-based flattened control indices.
        real(rk), allocatable, intent(out) :: dTgc(:,:)
            !! Selected parametric gradients `[nselected,2]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:)
            !! Optional selected basis values `[nselected]`.
        call query_basis_point(this, "derivative_scalar", Xt, elem, side, Tgc=Tgc, dTgc=dTgc)
    end subroutine derivative_scalar
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate complete parametric basis Hessians on a tensor-product grid.
    !!
    !! Let `ncp=product(nc)`. The packed output has shape
    !! `[ng_total,2*ncp,2]`, with
    !! `d2Tgc(g,(a-1)*ncp+A,b)` equal to
    !! \(\partial^2R_A/(\partial\xi_a\partial\xi_b)\). Thus both pure and
    !! mixed derivatives are available and symmetry is represented explicitly.
    !! Parameter selection and validation are the same as for the vector
    !! [[nurbs_surface:derivative]] overload.
    pure subroutine derivative2_vector(this, res1, res2, Xt1, Xt2, d2Tgc, dTgc, Tgc, side)
        class(nurbs_surface), target, intent(inout) :: this
        integer, intent(in), optional :: side(:)
        integer, intent(in), optional :: res1, res2
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:,:)
            !! Packed Hessians `[ng_total,2*nc_total,2]`.
        real(rk), allocatable, intent(out), optional :: dTgc(:,:,:)
            !! Optional first derivatives `[ng_total,nc_total,2]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
            !! Optional basis matrix `[ng_total,nc_total]`.
        call query_basis_grid(this, "derivative2_vector", res1, res2, Xt1, Xt2, side, Tgc=Tgc, dTgc=dTgc, d2Tgc=d2Tgc)
    end subroutine derivative2_vector
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate a complete parametric basis Hessian at one parameter pair.
    !!
    !! For `nselected=product(nc)` or `size(elem)`, the Hessian entry
    !! `d2Tgc((a-1)*nselected+i,b)` stores the selected basis derivative
    !! \(\partial^2R_i/(\partial\xi_a\partial\xi_b)\). Both mixed entries
    !! are populated. Rational normalization uses the complete active support,
    !! including controls omitted by `elem`. `Xt` must have at least
    !! two finite entries; later entries are ignored. Missing spline state or
    !! invalid active parameters produce a diagnostic and empty requested outputs.
    pure subroutine derivative2_scalar(this, Xt, d2Tgc, dTgc, Tgc, elem, side)
        class(nurbs_surface), target, intent(inout) :: this
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter pair \([u,v]\).
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
            !! Packed selected Hessian `[2*nselected,2]`.
        real(rk), allocatable, intent(out), optional :: dTgc(:,:)
            !! Optional selected first derivatives `[nselected,2]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:)
            !! Optional selected basis values `[nselected]`.
        integer, intent(in), optional :: elem(:)
            !! Optional one-based global controls; ordering and repeats are preserved, invalid entries are zero.
        call query_basis_point(this, "derivative2_scalar", Xt, elem, side, Tgc=Tgc, dTgc=dTgc, d2Tgc=d2Tgc)
    end subroutine derivative2_scalar
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the dense surface basis on a tensor-product parameter grid.
    !! Explicit directional vectors take precedence over resolutions; every
    !! in-domain output row sums to one. Omitted directions reuse cached sample
    !! parameters. Query parameters and counts are local, preserving stored
    !! samples. With no arguments, cached paired rows are reused.
    !! Invalid resolutions, missing parameters, and nonfinite explicit
    !! parameters produce a diagnostic.
    pure subroutine basis_vector(this, res1, res2, Xt1, Xt2, Tgc, side)
        class(nurbs_surface), target, intent(inout) :: this
        integer, intent(in), optional :: side(:)
        integer, intent(in), optional :: res1, res2
            !! Optional positive uniform-grid sizes in directions 1 and 2.
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:)
            !! Optional directional parameter vectors.
        real(rk), allocatable, intent(out) :: Tgc(:,:)
            !! Basis matrix `[ng_total,nc_total]`.
        call query_basis_grid(this, "basis_vector", res1, res2, Xt1, Xt2, side, Tgc=Tgc)
    end subroutine basis_vector

    !> Evaluate the requested basis values, gradients, and Hessians on a grid.
    pure subroutine query_basis_grid(this, location, res1, res2, Xt1, Xt2, side, Tgc, dTgc, d2Tgc)
        class(nurbs_surface), target, intent(inout) :: this
        character(*), intent(in) :: location
        integer, intent(in), optional :: res1, res2
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:)
        integer, intent(in), optional :: side(:)
        real(rk), allocatable, intent(out), optional :: Tgc(:,:), dTgc(:,:,:), d2Tgc(:,:,:)
        real(rk), allocatable :: Xt(:,:)
        integer :: query_ng(2)
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call prepare_surface_query(this, location, Xt, query_ng, &
            valid, res1, res2, Xt1, Xt2,side=side)
        if (.not. valid) then
            if (present(d2Tgc)) allocate(d2Tgc(0,0,0))
            if (present(dTgc)) allocate(dTgc(0,0,0))
            if (present(Tgc)) allocate(Tgc(0,0))
            return
        end if

        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
        call compute_basis_data(Xt,this%knot1,this%knot2,this%degree,this%nc, &
            Tgc=Tgc,dTgc=dTgc,d2Tgc=d2Tgc,wrap_parameters=this%wrap_parameters,Wc=weights,side=side)
    end subroutine query_basis_grid
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the surface basis at one parameter pair.
    !! Optional `elem` selects and orders a local control subset. `Xt` must have
    !! at least two finite entries; later entries are ignored. The spline must
    !! be initialized. This convenience overload does not diagnose those
    !! preconditions.
    pure subroutine basis_scalar(this, Xt, Tgc, elem, side)
        class(nurbs_surface), target, intent(inout) :: this
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter pair \([u,v]\).
        integer, intent(in), optional :: elem(:)
            !! Optional one-based flattened control indices.
        real(rk), allocatable, intent(out) :: Tgc(:)
            !! Dense or selected basis vector.

        call query_basis_point(this, "basis_scalar", Xt, elem, side, Tgc=Tgc)
    end subroutine basis_scalar

    !> Evaluate the requested basis values, gradients, and Hessians at one point.
    pure subroutine query_basis_point(this, location, Xt, elem, side, Tgc, dTgc, d2Tgc)
        class(nurbs_surface), target, intent(inout) :: this
        character(*), intent(in) :: location
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: elem(:)
        integer, intent(in), optional :: side(:)
        real(rk), allocatable, intent(out), optional :: Tgc(:), dTgc(:,:), d2Tgc(:,:)
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call check_scalar_query(this,Xt,location,valid,side=side)
        if (.not. valid) then
            if (present(d2Tgc)) allocate(d2Tgc(0,0))
            if (present(dTgc)) allocate(dTgc(0,0))
            if (present(Tgc)) allocate(Tgc(0))
            return
        end if

        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
        call compute_basis_data(Xt,this%knot1,this%knot2,this%degree,this%nc, &
            Tgc=Tgc,dTgc=dTgc,d2Tgc=d2Tgc,elem=elem,wrap_parameters=this%wrap_parameters,Wc=weights,side=side)
    end subroutine query_basis_point
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Insert knots in one tensor direction without changing the surface.
    !!
    !! Requests are applied in array order. The resulting multiplicity may not
    !! exceed `degree(dir)+1`; equality represents a \(C^{-1}\) break in the
    !! selected direction. For a verified periodic direction, the maximum is
    !! `degree(dir)` and insertion also updates the cyclic exterior knots and
    !! repeated control layers. Let \([a,b]\) be its current active interval
    !! and \(\tau\) the value returned by `knot_tolerance`. Writing \(u_i\) for
    !! `Xth(i)`, validation accepts \(u_i\in[a-\tau,b+\tau]\), but inserts the
    !! supplied value exactly: no endpoint snapping or tolerance-based
    !! multiplicity matching is applied.
    !!
    !! Rational data are refined in homogeneous coordinates. Writing
    !! \(\mathbf P\) and \(\mathbf Q\) for the old and new Euclidean control
    !! arrays, after the new weights are fixed `Bs[new_control,old_control]` is
    !! the full flattened tensor-product operator satisfying
    !! \(\mathbf Q_{:,c}=\mathrm{Bs}\,\mathbf P_{:,c}\) for every coordinate
    !! \(c\). `B` is its coordinate-block expansion for control-point-major
    !! vectors. For rational data these are Euclidean fixed-weight maps that
    !! depend on the old and new weights, not homogeneous refinement maps.
    !!
    !! **Algorithm:** W. Boehm, *Computer-Aided Design* 12 (1980), 199--201.
    !! [doi:10.1016/0010-4485(80)90154-2](https://doi.org/10.1016/0010-4485(80)90154-2).
    !! Applied directionally using Piegl and Tiller, Algorithm A5.1.
    pure subroutine insert_knots(this, dir, Xth, r, B, Bs)
        class(nurbs_surface),             intent(inout) :: this
        integer,                          intent(in)    :: dir
            !! Refined parameter direction, 1 or 2.
        real(rk), contiguous,             intent(in)    :: Xth(:)
            !! Finite values in the tolerance-expanded closed active interval; values are not snapped.
        integer,  contiguous,             intent(in)    :: r(:)
            !! Nonnegative counts; final multiplicity is at most `degree(dir)+1`, or `degree(dir)` if periodic.
        real(rk), allocatable,  optional, intent(out)   :: B(:,:)
            !! Optional coordinate-expanded refinement operator.
        real(rk), allocatable,  optional, intent(out)   :: Bs(:,:)
            !! Optional scalar refinement operator.
        real(rk), allocatable :: knots(:), Xc(:,:), Wc(:)

        if (.not.this%err%ok) return
        if (.not.this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_surface",location="insert_knots", &
                message="Refinement requires valid geometry and nonnegative counts.", &
                suggestion="Set complete spline geometry and check the refinement request.")
            return
        end if
        call refine_control_net(this%get_knot(dir),this%get_degree(dir),this%Xc,this%Wc,this%nc,dir, &
            knots,Xc,Wc,this%err,B,Bs,Xth=Xth,r=r)
        if (.not.this%err%ok .or. .not.allocated(Xc)) return
        call set_direction(this,dir,knots,this%get_knot(3-dir),Xc,Wc,[this%degree],[this%wrap_parameters])
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Raise the degree in one tensor direction while preserving active-domain geometry.
    !! Optional `Bs` and `B` are the full scalar and coordinate-expanded
    !! tensor-product mappings from old to elevated Euclidean control values
    !! after fixing the output weights. For rational data they depend on the old
    !! and new weights and are not homogeneous elevation operators.
    !! For `t>0`, a verified periodic direction retains its cyclic knot
    !! extension and repeated control layers. Other non-clamped input in the
    !! selected direction is converted to an equivalent open-clamped
    !! representation on its active interval. Other tensor directions are
    !! unchanged.
    !!
    !! **Algorithm:** H. Prautzsch, *Computer Aided Geometric Design* 1 (1984),
    !! 193--198.
    !! [doi:10.1016/0167-8396(84)90031-1](https://doi.org/10.1016/0167-8396(84)90031-1).
    !! Applied directionally using Piegl and Tiller, Algorithm A5.9.
    pure subroutine elevate_degree(this, dir, t, B, Bs)
        class(nurbs_surface), target,             intent(inout) :: this
        integer,                          intent(in)    :: dir
            !! Elevated parameter direction, 1 or 2.
        integer,                          intent(in)    :: t
            !! Nonnegative degree increment; zero is the identity operation.
        real(rk), allocatable, optional,  intent(out)   :: B(:,:)
            !! Optional coordinate-expanded elevation operator.
        real(rk), allocatable, optional,  intent(out)   :: Bs(:,:)
            !! Optional scalar elevation operator.
        real(rk), allocatable :: knots(:), Xc(:,:), Wc(:)
        integer :: degree(2)

        if (.not.this%err%ok) return
        if (.not.this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_surface",location="elevate_degree", &
                message="Refinement requires valid geometry and nonnegative counts.", &
                suggestion="Set complete spline geometry and check the refinement request.")
            return
        end if
        call refine_control_net(this%get_knot(dir),this%get_degree(dir),this%Xc,this%Wc,this%nc,dir, &
            knots,Xc,Wc,this%err,B,Bs,t=t)
        if (.not.this%err%ok .or. .not.allocated(Xc)) return
        degree=this%degree
        degree(dir)=degree(dir)+t
        call set_direction(this,dir,knots,this%get_knot(3-dir),Xc,Wc,degree,[this%wrap_parameters])
    end subroutine
    !===============================================================================


    !> Reduce degree under one Euclidean error budget for the complete call.
    !! Each candidate is compared with the call-entry geometry. Return the longest
    !! accepted prefix in reduced; over-budget proposals are normal rejections.
    !! Commit once, clearing derived caches only when at least one step succeeds.
    pure subroutine reduce_degree(this, dir, t, tol, reduced, error_bound)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in) :: dir
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
        integer :: nc(2),degree(2)
        real(rk), allocatable :: Xc(:,:),Wc(:)

        if (present(reduced)) reduced=0
        if (present(error_bound)) error_bound=0.0_rk
        if (.not.this%err%ok) return
        if (.not.this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_surface",location="reduce_degree", &
                message="Degree reduction requires valid geometry.",suggestion="Initialize valid knots, controls and weights.")
            return
        end if
        if (dir<1 .or. dir>2) then
            call this%err%set(code=100,severity=1,category="forcad_nurbs_surface",location="reduce_degree", &
                message="Invalid reduction direction.",suggestion="Use a parameter direction from 1 to 2.")
            return
        end if
        knots=this%get_knot(dir)
        call reduce_control_lines(knots,this%degree(dir),this%Xc,this%Wc,this%nc,dir,t,tol, &
            this%err,controls,weight_exponent,accepted,bound)
        if (.not.this%err%ok .or. accepted==0) return
        nc=this%nc
        degree=this%degree
        degree(dir)=degree(dir)-accepted
        nc(dir)=size(controls,1)
        call unpack_refinement_controls(controls,nc,dir,size(this%Xc,2),weight_exponent,Xc,Wc)
        call set_direction(this,dir,knots,this%get_knot(3-dir),Xc,Wc,degree, &
            [this%get_parameter_wrapping(1),this%get_parameter_wrapping(2)])
        if (.not.this%err%ok) return
        if (present(reduced)) reduced=accepted
        if (present(error_bound)) error_bound=bound
    end subroutine reduce_degree

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report whether stored nonuniform weights affect basis evaluation.
    !!
    !! Exactly equal weights select the polynomial basis. Every nonconstant
    !! positive weight list selects rational evaluation; no tolerance discards
    !! weight variation. An active diagnostic or absent weights returns false.
    pure elemental function is_rational(this) result(r)
        class(nurbs_surface), intent(in) :: this
        logical :: r

        r = .false.
        if (.not. this%err%ok) return
        r = this%rational .and. allocated(this%Wc)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Replace control-net visualization connectivity with a copy of `elemConn`.
    pure subroutine set_elem_Xc_vis(this, elemConn)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (.not. this%err%ok) return

        this%elemConn_Xc_vis = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Replace sampled-surface visualization connectivity with a copy of `elemConn`.
    pure subroutine set_elem_Xg_vis(this, elemConn)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (.not. this%err%ok) return

        this%elemConn_Xg_vis = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Replace stored IGA element connectivity with a copy of `elemConn`.
    !!
    !! This low-level setter does not verify consistency with the current knot
    !! vectors or control-net dimensions.
    pure subroutine set_elem(this, elemConn)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (.not. this%err%ok) return

        this%elemConn = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocated copy of control-net visualization connectivity.
    pure function get_elem_Xc_vis(this) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        call get_array(this%elemConn_Xc_vis,elemConn,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocated copy of sampled-surface visualization connectivity.
    pure function get_elem_Xg_vis(this) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        call get_array(this%elemConn_Xg_vis,elemConn,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocated copy of cached IGA element-to-control-point connectivity.
    !! An absent cache or active diagnostic returns shape `[0,0]`.
    !! `cmp_elem()` computes connectivity independently without populating this cache.
    pure function get_elem(this) result(elemConn)
        class(nurbs_surface), intent(in) :: this
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
    pure subroutine remove_knots(this, dir, Xth, r, removed, tol, error_bound)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in) :: dir
            !! Modified parameter direction, 1 through 2.
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
        real(rk), allocatable :: Xc(:,:), Wc(:)
        integer :: nc(2), weight_exponent

        if (present(removed)) allocate(removed(size(r)),source=0)
        if (present(error_bound)) error_bound=0.0_rk
        if (.not.this%err%ok) return
        if (.not.this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_surface",location="remove_knots", &
                message="Knot removal requires valid geometry.",suggestion="Initialize valid controls, knots and weights.")
            return
        end if
        if (dir<1 .or. dir>2) then
            call this%err%set(code=100,severity=1,category="forcad_nurbs_surface",location="remove_knots", &
                message="Invalid knot removal direction.",suggestion="Use dir from 1 through 2.")
            return
        end if
        nc=this%nc
        knots=this%get_knot(dir)
        call coarsen_control_lines(knots,this%degree(dir),this%Xc,this%Wc,nc,dir,Xth,r, &
            this%err,controls,weight_exponent,counts,tol,error_bound)
        if (.not.this%err%ok .or. .not.any(counts>0)) return
        nc(dir)=size(controls,1)
        call unpack_refinement_controls(controls,nc,dir,size(this%Xc,2),weight_exponent,Xc,Wc)
        call set_direction(this,dir,knots,this%get_knot(3-dir),Xc,Wc,[this%degree],[this%wrap_parameters])
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
    !> Construct an axis-aligned rectangular tensor-product surface.
    !!
    !! The control net is uniformly distributed on
    !! `[0,L(1)] x [0,L(2)]` in the plane `z=0`; direction 1 varies fastest.
    !! The generated knot vectors follow the standard `set(nc=...,Xc=...)`
    !! policy.
    !!
    pure subroutine set_tetragon(this, L, nc, Wc)
        class(nurbs_surface), intent(inout) :: this
            !! Surface replaced by the generated rectangle.
        real(rk), intent(in), contiguous :: L(:)
            !! Finite side lengths `[Lx,Ly]`.
        integer, intent(in), contiguous :: nc(:)
            !! Directional control counts `[nu,nv]`, both greater than one.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional finite positive weights `[product(nc)]`.

        if (.not. this%err%ok) return

        call this%set(nc = nc, Xc = tetragon_Xc(L, nc), Wc = Wc)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build IGA connectivity over all nonzero tensor-product knot spans.
    !!
    !! Each row contains \((p+1)(q+1)\) flattened control indices in
    !! direction-1-fastest order. Repeated knots change adjacent overlap
    !! according to their multiplicities. For verified periodic directions,
    !! apply [[nurbs_surface:cmp_dof_map]] before assembling independent cyclic
    !! degrees of freedom.
    pure function cmp_elem(this) result(elemConn)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not.this%err%ok .or. &
            .not.allocated(this%knot1) .or. &
            .not.allocated(this%knot2)) then
            allocate(elemConn(0,0))
        else
            elemConn=spline_connectivity(this%nc,this%degree,this%knot1,knot2=this%knot2)
        end if
    end function cmp_elem
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Map stored controls to independent cyclic IGA degrees of freedom.
    !!
    !! The returned one-based array follows direction-1-fastest storage.
    !! A bounded direction keeps all stored controls; a verified periodic
    !! degree-\(p_d\) direction identifies its final \(p_d\) layers with its
    !! first layers. Apply this map to [[nurbs_surface:cmp_elem]] connectivity.
    pure function cmp_dof_map(this) result(map)
        class(nurbs_surface), intent(in) :: this
        integer, allocatable :: map(:)
        integer :: i, independent_count(2)

        if (.not. this%is_valid()) return

        independent_count = this%nc
        if (this%get_parameter_topology(1) == "periodic") then
            independent_count(1) = this%nc(1) - this%degree(1)
        end if
        if (this%get_parameter_topology(2) == "periodic") then
            independent_count(2) = this%nc(2) - this%degree(2)
        end if

        allocate(map(product(this%nc)))
        do concurrent (i = 0:product(this%nc)-1)
            map(i+1) = modulo(modulo(i,this%nc(1)),independent_count(1)) + 1 + &
                modulo(i/this%nc(1),independent_count(2))*independent_count(1)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Rotate all control points in place using x-y-z Euler angles in degrees.
    !!
    !! The active extrinsic convention is
    !! \(R_z(\theta)R_y(\beta)R_x(\alpha)\). For fewer than three coordinates,
    !! only the leading principal matrix block is applied; coordinates after
    !! the third are unchanged. Cached physical samples are not changed.
    pure subroutine rotate_Xc(this, alpha, beta, theta)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: alpha, beta, theta

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) return

        call rotate_points(this%Xc, this%nc(1)*this%nc(2), alpha, beta, theta)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Rotate all cached physical samples using x-y-z Euler angles in degrees.
    !!
    !! The convention and coordinate-dimension behavior match
    !! [[nurbs_surface:rotate_Xc]]. Control points are not changed.
    pure subroutine rotate_Xg(this, alpha, beta, theta)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in) :: alpha, beta, theta

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xg)) return

        call rotate_points(this%Xg, this%ng(1)*this%ng(2), alpha, beta, theta)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Translate every control point by physical-space vector `vec`.
    !!
    !! `vec` must supply every stored coordinate. Cached samples are not
    !! changed. Extra entries are ignored and finiteness is not checked.
    pure subroutine translate_Xc(this, vec)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: vec(:)
        integer :: i, d

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) return

        do concurrent (i = 1:this%nc(1)*this%nc(2), d = 1:size(this%Xc,2))
            this%Xc(i,d) = this%Xc(i,d) + vec(d)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Translate every cached physical sample by vector `vec`.
    !!
    !! `vec` must supply every stored coordinate. Control points are not
    !! changed. Extra entries are ignored and finiteness is not checked.
    pure subroutine translate_Xg(this, vec)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: vec(:)
        integer :: i, d

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xg)) return

        do concurrent (i = 1:this%ng(1)*this%ng(2), d = 1:size(this%Xg,2))
            this%Xg(i,d) = this%Xg(i,d) + vec(d)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct an exact planar rational annulus.
    !!
    !! Three quadratic rational arcs form each circular boundary and a linear
    !! radial direction connects them. The result has a `7 x 2` control net in
    !! the xy plane at `center(3)`. The two boundary radii are
    !! `abs(radius1)` and `abs(radius2)`. For the named nondegenerate annulus,
    !! use nonzero scales of the same sign and unequal magnitudes (conventionally
    !! `0<radius1<radius2`). Opposite signs make corresponding radial lines
    !! cross the center; these cases are accepted but are not a regular annulus.
    !!
    pure subroutine set_ring(this, center, radius1, radius2)
        class(nurbs_surface), intent(inout) :: this
            !! Surface replaced by the annulus.
        real(rk), intent(in), contiguous :: center(:)
            !! Common center; at least three finite components are required.
        real(rk), intent(in) :: radius1
            !! Finite signed scale of the first circular boundary.
        real(rk), intent(in) :: radius2
            !! Finite signed scale of the second circular boundary.
        real(rk), allocatable :: Xc(:,:), Wc(:), knot1(:), knot2(:)
        integer :: i, d

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
        do concurrent (i = 1:size(Xc, 1), d = 1:size(Xc, 2))
            Xc(i,d) = center(d) + Xc(i,d)
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
    !> Construct an exact planar rational C-shaped annular sector.
    !!
    !! Two quadratic rational arcs span 240 degrees on each circular boundary;
    !! a linear radial direction connects the boundaries. The result has a
    !! `5 x 2` control net in the xy plane at `center(3)`. Boundary radii are
    !! the absolute scale values. Use nonzero scales of the same sign and
    !! unequal magnitudes for a regular C-shaped annular sector; opposite signs
    !! create crossing radial lines.
    !!
    pure subroutine set_C(this, center, radius1, radius2)
        class(nurbs_surface), intent(inout) :: this
            !! Surface replaced by the annular sector.
        real(rk), intent(in), contiguous :: center(:)
            !! Common center; at least three finite components are required.
        real(rk), intent(in) :: radius1
            !! Finite signed scale of the first circular boundary.
        real(rk), intent(in) :: radius2
            !! Finite signed scale of the second circular boundary.
        real(rk), allocatable :: Xc(:,:), Wc(:), knot1(:), knot2(:)
        integer :: i, d

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
        do concurrent (i = 1:size(Xc, 1), d = 1:size(Xc, 2))
            Xc(i,d) = center(d) + Xc(i,d)
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
    !> Construct an exact planar rational half-annulus.
    !!
    !! Two quadratic quarter-arcs form each semicircular boundary and a linear
    !! radial direction connects them. The result has a `5 x 2` control net in
    !! the xy plane at `center(3)`. Use nonzero scales of the same sign and
    !! unequal absolute values for a regular half-annulus. A negative common
    !! sign rotates the template by \(\pi\); opposite signs produce crossing
    !! radial lines.
    !!
    pure subroutine set_half_ring(this, center, radius1, radius2)
        class(nurbs_surface), intent(inout) :: this
            !! Surface replaced by the half-annulus.
        real(rk), intent(in), contiguous :: center(:)
            !! Common center; at least three finite components are required.
        real(rk), intent(in) :: radius1
            !! Finite legacy first-boundary scale; represented radius `abs(radius1)/2`.
        real(rk), intent(in) :: radius2
            !! Finite legacy second-boundary scale; represented radius `abs(radius2)/2`.
        real(rk), allocatable :: Xc(:,:), Wc(:), knot1(:), knot2(:)
        integer :: i, d

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
        do concurrent (i = 1:size(Xc, 1), d = 1:size(Xc, 2))
            Xc(i,d) = center(d) + Xc(i,d)
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
    !> Find the closest cached sample in the shared physical coordinate prefix.
    !! Extra query coordinates are ignored and a longer output is zero-padded.
    !! Ties select the first cached sample. Missing/sticky state, empty data or
    !! nonfinite inputs return zero outputs. Call create() before this discrete search.
    pure subroutine nearest_point(this, point_Xg, nearest_Xg, nearest_Xt, id)
        class(nurbs_surface), intent(in) :: this
        real(rk), intent(in), contiguous :: point_Xg(:)
        real(rk), intent(out), optional :: nearest_Xg(size(point_Xg)), nearest_Xt(2)
        integer, intent(out), optional :: id
        integer :: index, npoint, ncopy

        if (present(nearest_Xg)) nearest_Xg=0.0_rk
        if (present(nearest_Xt)) nearest_Xt=0.0_rk
        if (present(id)) id=0
        if (.not.this%err%ok .or. .not.allocated(this%Xg) .or. .not.allocated(this%Xt)) return
        npoint=min(size(this%Xg,1),size(this%Xt,1))
        index=nearest_sample_index(this%Xg(:npoint,:),point_Xg)
        if (index==0) return
        ncopy=min(size(point_Xg),size(this%Xg,2))
        if (present(nearest_Xg)) nearest_Xg(:ncopy)=this%Xg(index,:ncopy)
        if (present(nearest_Xt)) nearest_Xt=this%Xt(index,1:2)
        if (present(id)) id=index
    end subroutine nearest_point
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Find a bounded local closest-point candidate without changing sample caches.
    !! Span-aware coarse seeds or optional literal initial_Xt initialize a shared
    !! projected Newton solve with Gauss-Newton/descent fallback and backtracking.
    !! Rational geometry jets use local support; no dense basis matrices are formed.
    !! This is a local numerical candidate, without a global-nearest guarantee.
    !! An optional caller err is reset for this call and leaves this%err unchanged;
    !! without it, sticky object diagnostics are respected and failures stored there.
    !! Invalid input returns zero outputs. A failed iteration retains its candidate;
    !! maxit=0 returns the bounded initial candidate with nonconvergence diagnostic109.
    pure subroutine nearest_point2(this, point_Xg, tol, maxit, nearest_Xt, nearest_Xg, initial_Xt, err)
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: point_Xg(:)
            !! Finite query with at least ncoord entries; extra coordinates are ignored.
        real(rk), intent(in) :: tol
            !! Finite nonnegative convergence tolerance.
        integer, intent(in) :: maxit
            !! Nonnegative iteration limit.
        real(rk), intent(out) :: nearest_Xt(2)
        real(rk), intent(out), optional :: nearest_Xg(:)
            !! Optional point output with exactly ncoord entries.
        real(rk), intent(in), optional :: initial_Xt(:)
            !! Optional finite parameter pair inside the literal active box.
        type(debug), intent(out), optional :: err
        type(debug) :: status

        nearest_Xt=0.0_rk
        if (present(nearest_Xg)) nearest_Xg=0.0_rk
        call status%reset()
        if (.not.present(err)) status=this%err
        if (status%ok) then
            if (.not.allocated(this%Xc) .or. .not.allocated(this%knot1) .or. .not.allocated(this%knot2)) then
                call status%set(code=103,severity=1,category="forcad_nurbs_surface",location="nearest_point2", &
                    message="Surface geometry is not initialized.",suggestion="Set complete surface geometry before projection.")
            else
                call nearest_spline_point(this%Xc,this%Wc,this%nc,this%degree,this%knot1,point_Xg,tol,maxit, &
                    nearest_Xt,nearest_Xg,status,knot2=this%knot2,initial_Xt=initial_Xt,wrap=this%wrap_parameters)
            end if
        end if
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
    !! Element and local basis ordering agree with `cmp_elem`. `dA` includes the
    !! tensor Gauss weight, span maps, and either the absolute planar Jacobian or
    !! \(\sqrt{\det G}\) for an embedded surface metric \(G\). For a planar
    !! two-coordinate mapping, gradients are physical Cartesian gradients. For
    !! an embedded surface, they are ambient Cartesian components of the
    !! tangential gradient; optional Hessians are intrinsic covariant surface
    !! Hessians represented as ambient tangential tensors. A singular Jacobian
    !! or rank-deficient metric sets diagnostic code 108. With
    !! `strict=.true.`, a two-coordinate physical map must additionally have a
    !! positive Jacobian determinant. For an embedded surface, strict mode
    !! verifies rank but cannot assign a signed orientation without an external
    !! reference normal. The check is local to the requested quadrature point
    !! and does not prove global injectivity between quadrature points.
    !!
    !! With \(\mathbf J=[\mathbf S_{,u}\ \mathbf S_{,v}]\),
    !! \(\mathbf G=\mathbf J^T\mathbf J\), and
    !! \(\mathbf J^+=\mathbf G^{-1}\mathbf J^T\), the ambient tangential
    !! gradient is
    !!
    !! \[
    !! \nabla_\Gamma R_a=(\mathbf J^+)^T
    !! \begin{bmatrix}R_{a,u}\\R_{a,v}\end{bmatrix},
    !! \qquad
    !! dA=\omega_g\,\Delta u_e\,\Delta v_e\,\sqrt{\det\mathbf G}.
    !! \]
    !!
    !! Let \(\mathbf g^\alpha\) denote the contravariant tangent vectors, the
    !! columns of \(\mathbf J\mathbf G^{-1}\). The optional ambient
    !! representation of the covariant Hessian is
    !!
    !! \[
    !! R_{a\mid\alpha\beta}
    !!   =R_{a,\alpha\beta}
    !!    -\nabla_\Gamma R_a\mathbin{\cdot}\mathbf S_{,\alpha\beta},
    !! \qquad
    !! \operatorname{Hess}_\Gamma R_a
    !!   =\sum_{\alpha,\beta=1}^{2}
    !!      R_{a\mid\alpha\beta}\,
    !!      \mathbf g^\alpha\mathbin{\otimes}\mathbf g^\beta .
    !! \]
    !!
    !! IGA formulation reference: T. J. R. Hughes, J. A. Cottrell, and
    !! Y. Bazilevs, *CMAME* 194 (2005), 4135--4195.
    !! [doi:10.1016/j.cma.2004.10.008](https://doi.org/10.1016/j.cma.2004.10.008).
    !!
    !! Elementwise analysis assumes an invertible geometry map with a smooth
    !! inverse: Y. Bazilevs et al., *MMMAS* 16 (2006), 1031--1090.
    !! [doi:10.1142/S0218202506001455](https://doi.org/10.1142/S0218202506001455).
    pure subroutine ansatz(this, ie, ig, Tgc, dTgc_dXg, dA, ngauss, d2Tgc_dXg2, strict, err, Xg, dXg_dXt, &
        dmeasure_dXc, dmeasure_dlogWc, dgradient_dXc, dgradient_dlogWc, djacobian_dXc, djacobian_dlogWc)
        class(nurbs_surface), intent(inout) :: this
        integer, intent(in) :: ie
            !! One-based direction-1-fastest active element index.
        integer, intent(in) :: ig
            !! One-based flattened tensor quadrature-point index.
        real(rk), intent(out) :: dA
            !! Weighted differential area contribution.
        real(rk), allocatable, intent(out) :: Tgc(:)
            !! Local basis values, size `(degree(1)+1)*(degree(2)+1)`.
        real(rk), allocatable, intent(out) :: dTgc_dXg(:,:)
            !! Local physical or ambient tangential gradients `[nlocal,physical_dim]`.
        integer, intent(in), optional :: ngauss(2)
            !! Optional positive point counts `[2]`; defaults to `degree+1`.
        real(rk), allocatable, target, intent(out), optional :: d2Tgc_dXg2(:,:,:)
            !! Optional planar physical or embedded covariant Hessians `[nlocal,ncoord,ncoord]`.
        logical, intent(in), optional :: strict
            !! Reject orientation reversal in addition to rank deficiency.
        type(debug), intent(out), optional :: err
            !! Caller-owned diagnostic; when present, object diagnostics are neither read nor changed.
        real(rk), allocatable, target, intent(out), optional :: Xg(:), dXg_dXt(:,:)
            !! Physical point `[ncoord]` and raw geometry Jacobian `[ncoord,2]`.
        real(rk), allocatable, target, intent(out), optional :: dmeasure_dXc(:,:), dmeasure_dlogWc(:)
            !! Derivatives of this weighted area w.r.t. local controls and log weights.
        real(rk), allocatable, intent(out), optional :: dgradient_dXc(:,:,:,:), dgradient_dlogWc(:,:,:)
            !! Physical basis gradient derivatives with respect to controls `(A,i,B,j)` and log weights `(A,i,B)`.
        real(rk), allocatable, intent(out), optional :: djacobian_dXc(:,:), djacobian_dlogWc(:,:,:)
            !! Raw Jacobian derivatives for controls (A,alpha) and log weights (i,alpha,A).

        call tensor_ansatz(this%Xc,this%Wc,this%rational,this%degree,this%nc,this%knot1,this%knot2, &
            this%err,err,ie,ig,Tgc,dTgc_dXg,dA,ngauss,strict,Xg,dXg_dXt,d2Tgc_dXg2, &
            dmeasure_dXc,dmeasure_dlogWc,dgradient_dXc,dgradient_dlogWc, &
            djacobian_dXc=djacobian_dXc,djacobian_dlogWc=djacobian_dlogWc)
    end subroutine ansatz
    !===============================================================================

    !===============================================================================
    !> Integrate nonnegative physical area over the active knot-span cells.
    !! One tensor Gauss rule is reused on every cell; degree+1 points per
    !! direction are the default numerical approximation. Local rational
    !! derivatives and the conditioned Jacobian frame match ansatz_element.
    !! Invalid, singular, or strictly inverted points are omitted while valid
    !! contributions remain in area. An active object diagnostic returns zero.
    !! With caller-owned err, object diagnostics are neither read nor changed.
    pure subroutine cmp_area(this, area, ngauss, strict, err)
        class(nurbs_surface), target, intent(inout) :: this
        real(rk), intent(out) :: area
            !! Accumulated nonnegative area, including valid points preceding any failure.
        integer, intent(in), optional :: ngauss(2)
            !! Positive points per direction; default degree+1.
        logical, intent(in), optional :: strict
            !! Reject planar orientation reversal as well as rank deficiency.
        type(debug), target, intent(out), optional :: err
            !! Caller-owned diagnostic, reset at entry and independent of the object diagnostic.
        call tensor_measure(this%Xc,this%degree,this%nc,this%knot1,this%knot2,this%Wc,this%rational, &
            area,this%err,err,ngauss,strict)
    end subroutine cmp_area
    !===============================================================================


    !===============================================================================
    !> Contract support-sized basis values with controls, without a dense design matrix.
    pure function compute_Xg_scalar(Xt, knot1, knot2, degree, nc, Xc, wrap_parameters, Wc, side) result(Xg)
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), Xc(:,:)
        integer, intent(in) :: degree(2), nc(2)
        logical, intent(in), optional :: wrap_parameters(:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk) :: Xg(size(Xc,2)), basis(product(degree+1))
        integer :: first(2), l1, l2, local, global, d

        Xg = 0.0_rk
        call compute_basis_jet_2d(Xt,knot1,knot2,degree,nc,0,first,basis,wrap_parameters=wrap_parameters,Wc=Wc,side=side)
        do l2 = 0, degree(2)
            do l1 = 0, degree(1)
                global = first(1)+l1+nc(1)*(first(2)+l2-1)
                local = l1+(degree(1)+1)*l2+1
                do d = 1, size(Xc,2)
                    Xg(d) = Xg(d)+basis(local)*Xc(global,d)
                end do
            end do
        end do
    end function compute_Xg_scalar
    !===============================================================================


    !===============================================================================
    !> Contract local basis values for paired rows; strided parameter batches are accepted.
    pure function compute_Xg_vector(Xt, knot1, knot2, degree, nc, ng, Xc, wrap_parameters, Wc, side) result(Xg)
        real(rk), intent(in) :: Xt(:,:)
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), Xc(:,:)
        integer, intent(in) :: degree(2), nc(2)
        integer, intent(in), optional :: ng(:)
        logical, intent(in), optional :: wrap_parameters(:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk), allocatable :: Xg(:,:)
        integer :: npoint, point

        npoint=size(Xt,1)
        if (present(ng)) npoint=product(ng)
        allocate(Xg(npoint,size(Xc,2)),source=0.0_rk)
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do point=1,npoint
#else
        do concurrent (point=1:npoint)
#endif
            Xg(point,:)=compute_Xg_scalar(Xt(point,:),knot1,knot2,degree,nc,Xc,wrap_parameters,Wc,side)
        end do
    end function compute_Xg_vector
    !===============================================================================


    !===============================================================================
    !> Evaluate only the requested local derivative order through the shared tensor kernels.
    pure subroutine compute_basis_jet_2d(Xt, knot1, knot2, degree, nc, order, first, basis, gradient, hessian, &
        wrap_parameters, Wc, elem, side)
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2), nc(2), order
        integer, intent(out) :: first(2)
        real(rk), intent(out), contiguous :: basis(:)
        real(rk), intent(out), contiguous, optional :: gradient(:,:), hessian(:,:,:)
        logical, intent(in), optional :: wrap_parameters(:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        integer, intent(in), optional :: elem(:)
        real(rk) :: ders1(0:order,0:degree(1)), ders2(0:order,0:degree(2)), identity(0:order,0:0)
        call tensor_basis_jet(Xt,knot1,knot2,degree=degree,nc=nc,order=order,wrap=wrap_parameters,first=first, &
            basis=basis,gradient=gradient,hessian=hessian,ders1=ders1,ders2=ders2,ders3=identity,Wc=Wc,elem=elem,side=side)
    end subroutine compute_basis_jet_2d
    !===============================================================================


    !===============================================================================
    !> Scatter requested basis data into dense rows, allocating only requested outputs.
    pure subroutine compute_basis_data_vector(Xt, knot1, knot2, degree, nc, Tgc, dTgc, d2Tgc, wrap_parameters, Wc, ng, side)
        real(rk), intent(in) :: Xt(:,:)
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2), nc(2)
        real(rk), allocatable, intent(out), optional :: Tgc(:,:), dTgc(:,:,:), d2Tgc(:,:,:)
        logical, intent(in), optional :: wrap_parameters(:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        integer, intent(in), optional :: ng(2)
        real(rk) :: basis(product(degree+1)), gradient(product(degree+1),2), hessian(product(degree+1),2,2)
        integer :: point, first(2), l1, l2, local, global, ncp, npoint, order

        ncp = product(nc)
        npoint = size(Xt,1)
        if (present(ng)) npoint = product(ng)
        if (npoint<0 .or. npoint>size(Xt,1) .or. size(Xt,2)<2) npoint = 0
        if (present(Tgc)) allocate(Tgc(npoint,ncp),source=0.0_rk)
        if (present(dTgc)) allocate(dTgc(npoint,ncp,2),source=0.0_rk)
        if (present(d2Tgc)) allocate(d2Tgc(npoint,2*ncp,2),source=0.0_rk)
        order = 0
        if (present(dTgc)) order = 1
        if (present(d2Tgc)) order = 2
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do point = 1, npoint
#else
        do concurrent (point = 1:npoint) local(first,l1,l2,local,global,basis,gradient,hessian)
#endif
            call compute_basis_jet_2d(Xt(point,:),knot1,knot2,degree,nc,order,first,basis,gradient,hessian,wrap_parameters,Wc, &
                side=side)
            do l2 = 0, degree(2)
                do l1 = 0, degree(1)
                    global = first(1)+l1+nc(1)*(first(2)+l2-1)
                    local = l1+(degree(1)+1)*l2+1
                    if (present(Tgc)) Tgc(point,global) = basis(local)
                    if (present(dTgc)) dTgc(point,global,:) = gradient(local,:)
                    if (present(d2Tgc)) then
                        d2Tgc(point,global,:) = hessian(local,1,:)
                        d2Tgc(point,ncp+global,:) = hessian(local,2,:)
                    end if
                end do
            end do
        end do
    end subroutine compute_basis_data_vector
    !===============================================================================


    !===============================================================================
    !> Scatter one local calculation in global order or the supplied control selection.
    pure subroutine compute_basis_data_scalar(Xt, knot1, knot2, degree, nc, Tgc, dTgc, d2Tgc, elem, wrap_parameters, Wc, side)
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2), nc(2)
        real(rk), allocatable, intent(out), optional :: Tgc(:), dTgc(:,:), d2Tgc(:,:)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), optional :: wrap_parameters(:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk) :: basis(product(degree+1)), gradient(product(degree+1),2), hessian(product(degree+1),2,2)
        integer :: first(2), l1, l2, local, global, selected, target, nselected, nwork, order

        nselected = product(nc)
        nwork = product(degree+1)
        if (present(elem)) then
            nselected = size(elem)
            nwork = nselected
        end if
        if (present(Tgc)) allocate(Tgc(nselected),source=0.0_rk)
        if (present(dTgc)) allocate(dTgc(nselected,2),source=0.0_rk)
        if (present(d2Tgc)) allocate(d2Tgc(2*nselected,2),source=0.0_rk)
        order = 0
        if (present(dTgc)) order = 1
        if (present(d2Tgc)) order = 2
        call compute_basis_jet_2d(Xt,knot1,knot2,degree,nc,order,first,basis,gradient,hessian,wrap_parameters,Wc,side=side)
        do selected = 1, nwork
            local = selected
            global = first(1)+mod(local-1,degree(1)+1)+nc(1)*(first(2)+(local-1)/(degree(1)+1)-1)
            target = global
            if (present(elem)) then
                global = elem(selected)
                if (global<1 .or. global>product(nc)) cycle
                l1 = mod(global-1,nc(1))+1-first(1)
                l2 = (global-1)/nc(1)+1-first(2)
                if (l1<0 .or. l1>degree(1) .or. l2<0 .or. l2>degree(2)) cycle
                local = l1+(degree(1)+1)*l2+1
                target = selected
            end if
            if (present(Tgc)) Tgc(target) = basis(local)
            if (present(dTgc)) dTgc(target,:) = gradient(local,:)
            if (present(d2Tgc)) then
                d2Tgc(target,:) = hessian(local,1,:)
                d2Tgc(nselected+target,:) = hessian(local,2,:)
            end if
        end do
    end subroutine compute_basis_data_scalar
    !===============================================================================


    !===============================================================================
    !> Compatibility adapter to the shared local surface calculation.
    pure subroutine compute_dTgc_nurbs_2d_vector(Xt, knot1, knot2, degree, nc, ng, Wc, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in), optional :: ng(2)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        call compute_basis_data(Xt,knot1,knot2,degree,nc,Tgc=Tgc,dTgc=dTgc,ng=ng,wrap_parameters=wrap_parameters,Wc=Wc)
    end subroutine compute_dTgc_nurbs_2d_vector
    !===============================================================================


    !===============================================================================
    !> Compatibility adapter to the shared local surface calculation.
    pure subroutine compute_dTgc_nurbs_2d_scalar(Xt, knot1, knot2, degree, nc, Wc, dTgc, Tgc, elem, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        real(rk), intent(in), contiguous :: Wc(:)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        if (present(elem)) then
            if (size(Wc)==size(elem) .and. size(Wc)/=product(nc)) then
                call compute_basis_data(Xt,knot1,knot2,degree,nc,Tgc=Tgc,dTgc=dTgc,elem=elem,wrap_parameters=wrap_parameters)
                call normalize_rational_basis(Tgc,Wc,dTgc)
                return
            end if
        end if
        call compute_basis_data(Xt,knot1,knot2,degree,nc,Tgc=Tgc,dTgc=dTgc,elem=elem,wrap_parameters=wrap_parameters,Wc=Wc)
    end subroutine compute_dTgc_nurbs_2d_scalar
    !===============================================================================


    !===============================================================================
    !> Compatibility adapter to the shared local surface calculation.
    pure subroutine compute_dTgc_bspline_2d_vector(Xt, knot1, knot2, degree, nc, ng, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in), optional :: ng(2)
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        call compute_basis_data(Xt,knot1,knot2,degree,nc,Tgc=Tgc,dTgc=dTgc,ng=ng,wrap_parameters=wrap_parameters)
    end subroutine compute_dTgc_bspline_2d_vector
    !===============================================================================


    !===============================================================================
    !> Compatibility adapter to the shared local surface calculation.
    pure subroutine compute_dTgc_bspline_2d_scalar(Xt, knot1, knot2, degree, nc, dTgc, Tgc, elem, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        call compute_basis_data(Xt,knot1,knot2,degree,nc,Tgc=Tgc,dTgc=dTgc,elem=elem,wrap_parameters=wrap_parameters)
    end subroutine compute_dTgc_bspline_2d_scalar
    !===============================================================================


    !===============================================================================
    !> Compatibility adapter to the shared local surface calculation.
    pure function compute_Tgc_nurbs_2d_vector(Xt, knot1, knot2, degree, nc, ng, Wc, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in), optional :: ng(2)
        real(rk), intent(in), contiguous :: Wc(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Tgc(:,:)
        call compute_basis_data(Xt,knot1,knot2,degree,nc,Tgc=Tgc,ng=ng,wrap_parameters=wrap_parameters,Wc=Wc)
    end function compute_Tgc_nurbs_2d_vector
    !===============================================================================


    !===============================================================================
    !> Compatibility adapter to the shared local surface calculation.
    pure function compute_Tgc_nurbs_2d_scalar(Xt, knot1, knot2, degree, nc, Wc, elem, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        real(rk), intent(in), contiguous :: Wc(:)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Tgc(:)
        integer :: first(2)
        if (present(elem)) then
            if (size(Wc)==size(elem) .and. size(Wc)/=product(nc)) then
                allocate(Tgc(size(elem)))
                call compute_basis_jet_2d(Xt,knot1,knot2,degree,nc,0,first,Tgc, &
                    wrap_parameters=wrap_parameters,Wc=Wc,elem=elem)
                return
            end if
        end if
        call compute_basis_data(Xt,knot1,knot2,degree,nc,Tgc=Tgc,elem=elem,wrap_parameters=wrap_parameters,Wc=Wc)
    end function compute_Tgc_nurbs_2d_scalar
    !===============================================================================


    !===============================================================================
    !> Compatibility adapter to the shared local surface calculation.
    pure function compute_Tgc_bspline_2d_vector(Xt, knot1, knot2, degree, nc, ng, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in), optional :: ng(2)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Tgc(:,:)
        call compute_basis_data(Xt,knot1,knot2,degree,nc,Tgc=Tgc,ng=ng,wrap_parameters=wrap_parameters)
    end function compute_Tgc_bspline_2d_vector
    !===============================================================================


    !===============================================================================
    !> Compatibility adapter to the shared local surface calculation.
    pure function compute_Tgc_bspline_2d_scalar(Xt, knot1, knot2, degree, nc, elem, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        integer, intent(in) :: degree(2)
        integer, intent(in) :: nc(2)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Tgc(:)
        call compute_basis_data(Xt,knot1,knot2,degree,nc,Tgc=Tgc,elem=elem,wrap_parameters=wrap_parameters)
    end function compute_Tgc_bspline_2d_scalar
    !===============================================================================




    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Fit polynomial tensor-product B-spline control points to parameterized data.
    !!
    !! `ndata=[n1,n2]` defines the number of rows used. For a recognized
    !! direction-1-fastest Cartesian parameter grid, Kronecker separability is
    !! exploited through two directional banded normal-equation solves. If a
    !! directional factorization fails, two directional Householder QR solves
    !! preserve the Kronecker structure and avoid a dense tensor design.
    !! The grid test requires exactly `product(ndata)` parameter and data rows.
    !! Repeated coordinates are compared componentwise against 32 times the
    !! larger of their ULPs and `epsilon(rk)` times that parameter direction's
    !! sampled span. No absolute unit floor is used. Otherwise the first
    !! `product(ndata)` rows use one banded tensor-product normal equation, with
    !! a dense pivoted-QR fallback allocated only if Cholesky fails.
    !!
    !! The routine fits a polynomial tensor-product B-spline model. A surface
    !! with nonuniform weights is rejected without changing its controls; use
    !! [[nurbs_surface:lsq_fit_nurbs]] for a rational model. Uniform weights
    !! produce the same polynomial basis. `Xt` and `Xdata` must provide enough finite rows,
    !! `Xt` must have at least two columns, `Xdata` at least one, and each
    !! directional parameter set must give full-rank active-domain coverage.
    !! These shape and finiteness preconditions are not completely diagnosed.
    !! Only `Xc` is replaced; cached `Xg` is not recomputed.
    !!
    !! All paths require numerical full column rank. The QR fallback diagnoses
    !! rank deficiency with code 108 instead of returning nonunique controls.
    pure subroutine lsq_fit_bspline(this, Xt, Xdata, ndata)
        use forcad_utils, only: solve_least_squares_qr, solve_spd_banded
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:,:)
            !! Parameter pairs `[at least product(ndata),2]`.
        real(rk), intent(in), contiguous :: Xdata(:,:)
            !! Data coordinates `[at least product(ndata),ncoord]`.
        integer, intent(in) :: ndata(2)
            !! Directional data counts `[2]`, each at least the corresponding `nc`.
        real(rk), allocatable :: TtT(:,:), TtX(:,:), Xsol(:,:), design(:,:)
        real(rk), allocatable :: G1band(:,:), G2band(:,:), RHS(:,:), rhs12(:,:), rhs21(:,:), sol12(:,:), sol21(:,:)
        real(rk), allocatable :: Phi1(:,:), Phi2(:,:), design1(:,:), design2(:,:)
        real(rk), allocatable :: qr_rhs1(:,:), qr_rhs2(:,:)
        real(rk) :: N1(0:this%degree(1)), N2(0:this%degree(2)), ba, bb
        integer, allocatable :: F1(:), F2(:)
        real(rk) :: grid_tol, grid_scale(2), grid_value, reference_value
        integer :: i, n, dim_, ncp, row, col, i1, i2, bandwidth
        integer :: a1, a2, b1, b2, k, ia, ib, first1, first2
        logical :: tensor_grid, normal_ok

        if (.not. this%err%ok) return

        if (this%is_rational()) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Polynomial fitting cannot retain nonuniform rational weights.",&
                location   = "lsq_fit_bspline",&
                suggestion = "Use lsq_fit_nurbs() or reset the surface with uniform weights.")
            return
        end if

        if (this%nc(1) > ndata(1)) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Invalid number of control points in the first direction.",&
                location   = "lsq_fit_bspline",&
                suggestion = "Ensure that the number of control points does not exceed the number of data points.")
            return
        end if

        if (this%nc(2) > ndata(2)) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Invalid number of control points in the second direction.",&
                location   = "lsq_fit_bspline",&
                suggestion = "Ensure that the number of control points does not exceed the number of data points.")
            return
        end if

        n = ndata(1)*ndata(2)
        dim_ = size(Xdata, 2)
        ncp = this%nc(1)*this%nc(2)

        tensor_grid = size(Xt,1) == n .and. size(Xdata,1) == n
        if (tensor_grid) then
            grid_scale(1) = maxval(Xt(:,1)) - minval(Xt(:,1))
            grid_scale(2) = maxval(Xt(:,2)) - minval(Xt(:,2))
            do i2 = 1, ndata(2)
                do i1 = 1, ndata(1)
                    row = (i2-1)*ndata(1) + i1
                    grid_value = Xt(row,1)
                    reference_value = Xt(i1,1)
                    grid_tol = 32.0_rk*max(&
                        epsilon(1.0_rk)*grid_scale(1),&
                        abs(grid_value-nearest(grid_value,merge(-1.0_rk,1.0_rk,grid_value > 0.0_rk))),&
                        abs(reference_value-nearest(&
                            reference_value,merge(-1.0_rk,1.0_rk,reference_value > 0.0_rk))))
                    if (abs(grid_value - reference_value) > grid_tol) then
                        tensor_grid = .false.
                        exit
                    end if
                    grid_value = Xt(row,2)
                    reference_value = Xt((i2-1)*ndata(1)+1,2)
                    grid_tol = 32.0_rk*max(&
                        epsilon(1.0_rk)*grid_scale(2),&
                        abs(grid_value-nearest(grid_value,merge(-1.0_rk,1.0_rk,grid_value > 0.0_rk))),&
                        abs(reference_value-nearest(&
                            reference_value,merge(-1.0_rk,1.0_rk,reference_value > 0.0_rk))))
                    if (abs(grid_value - reference_value) > grid_tol) then
                        tensor_grid = .false.
                        exit
                    end if
                end do
                if (.not. tensor_grid) exit
            end do
        end if

        if (tensor_grid) then
            allocate(Phi1(0:this%degree(1), ndata(1)), Phi2(0:this%degree(2), ndata(2)))
            allocate(F1(ndata(1)), F2(ndata(2)))
            allocate(G1band(0:this%degree(1), this%nc(1)), G2band(0:this%degree(2), this%nc(2)), &
                RHS(ncp, dim_), source=0.0_rk)

            do i1 = 1, ndata(1)
                call basis_bspline_der(&
                    map_parameter(Xt(i1,1), this%knot1, this%nc(1), this%degree(1), this%wrap_parameters(1)),&
                    this%knot1, this%nc(1), this%degree(1), 0, F1(i1), Phi1(:,i1))
                do a1 = 0, this%degree(1)
                    ia = F1(i1) + a1
                    if (ia < 1 .or. ia > this%nc(1)) cycle
                    ba = Phi1(a1,i1)
                    do b1 = 0, this%degree(1)
                        ib = F1(i1) + b1
                        if (ib < 1 .or. ib > this%nc(1)) cycle
                        if (ia >= ib) G1band(ia-ib,ib) = G1band(ia-ib,ib) + ba*Phi1(b1,i1)
                    end do
                end do
            end do

            do i2 = 1, ndata(2)
                row = (i2-1)*ndata(1) + 1
                call basis_bspline_der(&
                    map_parameter(Xt(row,2), this%knot2, this%nc(2), this%degree(2), this%wrap_parameters(2)),&
                    this%knot2, this%nc(2), this%degree(2), 0, F2(i2), Phi2(:,i2))
                do a2 = 0, this%degree(2)
                    ia = F2(i2) + a2
                    if (ia < 1 .or. ia > this%nc(2)) cycle
                    ba = Phi2(a2,i2)
                    do b2 = 0, this%degree(2)
                        ib = F2(i2) + b2
                        if (ib < 1 .or. ib > this%nc(2)) cycle
                        if (ia >= ib) G2band(ia-ib,ib) = G2band(ia-ib,ib) + ba*Phi2(b2,i2)
                    end do
                end do
            end do

            do i2 = 1, ndata(2)
                do i1 = 1, ndata(1)
                    row = (i2-1)*ndata(1) + i1
                    do a2 = 0, this%degree(2)
                        if (F2(i2)+a2 < 1 .or. F2(i2)+a2 > this%nc(2)) cycle
                        do a1 = 0, this%degree(1)
                            if (F1(i1)+a1 < 1 .or. F1(i1)+a1 > this%nc(1)) cycle
                            ia = (F2(i2)+a2-1)*this%nc(1) + F1(i1) + a1
                            ba = Phi2(a2,i2)*Phi1(a1,i1)
                            do k = 1, dim_
                                RHS(ia,k) = RHS(ia,k) + ba*Xdata(row,k)
                            end do
                        end do
                    end do
                end do
            end do

            allocate(Xsol(ncp,dim_))
            allocate(rhs12(this%nc(1), this%nc(2)), rhs21(this%nc(2), this%nc(1)))
            normal_ok = .true.
            do k = 1, dim_
                do i2 = 1, this%nc(2)
                    do i1 = 1, this%nc(1)
                        rhs12(i1,i2) = RHS((i2-1)*this%nc(1)+i1,k)
                    end do
                end do
                sol12 = solve_spd_banded(G1band, rhs12)
                if (size(sol12,1) /= this%nc(1) .or. size(sol12,2) /= this%nc(2)) then
                    normal_ok = .false.
                    exit
                end if
                do i2 = 1, this%nc(2)
                    do i1 = 1, this%nc(1)
                        rhs21(i2,i1) = sol12(i1,i2)
                    end do
                end do
                sol21 = solve_spd_banded(G2band, rhs21)
                if (size(sol21,1) /= this%nc(2) .or. size(sol21,2) /= this%nc(1)) then
                    normal_ok = .false.
                    exit
                end if
                do i2 = 1, this%nc(2)
                    do i1 = 1, this%nc(1)
                        Xsol((i2-1)*this%nc(1)+i1,k) = sol21(i2,i1)
                    end do
                end do
            end do
            if (normal_ok) then
                this%Xc = Xsol
                return
            end if

            allocate(design1(ndata(1),this%nc(1)), source=0.0_rk)
            allocate(design2(ndata(2),this%nc(2)), source=0.0_rk)
            do i1 = 1, ndata(1)
                do a1 = 0, this%degree(1)
                    ia = F1(i1) + a1
                    if (ia < 1 .or. ia > this%nc(1)) cycle
                    design1(i1,ia) = Phi1(a1,i1)
                end do
            end do
            do i2 = 1, ndata(2)
                do a2 = 0, this%degree(2)
                    ia = F2(i2) + a2
                    if (ia < 1 .or. ia > this%nc(2)) cycle
                    design2(i2,ia) = Phi2(a2,i2)
                end do
            end do

            allocate(qr_rhs1(ndata(1),ndata(2)*dim_))
            do k = 1, dim_
                do i2 = 1, ndata(2)
                    col = (k-1)*ndata(2) + i2
                    do i1 = 1, ndata(1)
                        row = (i2-1)*ndata(1) + i1
                        qr_rhs1(i1,col) = Xdata(row,k)
                    end do
                end do
            end do
            sol12 = solve_least_squares_qr(design1, qr_rhs1)
            if (size(sol12,1) /= this%nc(1) .or. size(sol12,2) /= ndata(2)*dim_) then
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_surface",&
                    message    = "The first directional design is numerically rank deficient.",&
                    location   = "lsq_fit_bspline",&
                    suggestion = "Use distinct first-direction parameters with full spline-space coverage.")
                return
            end if

            allocate(qr_rhs2(ndata(2),this%nc(1)*dim_))
            do k = 1, dim_
                do i1 = 1, this%nc(1)
                    col = (k-1)*this%nc(1) + i1
                    do i2 = 1, ndata(2)
                        qr_rhs2(i2,col) = sol12(i1,(k-1)*ndata(2)+i2)
                    end do
                end do
            end do
            sol21 = solve_least_squares_qr(design2, qr_rhs2)
            if (size(sol21,1) /= this%nc(2) .or. size(sol21,2) /= this%nc(1)*dim_) then
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_surface",&
                    message    = "The second directional design is numerically rank deficient.",&
                    location   = "lsq_fit_bspline",&
                    suggestion = "Use distinct second-direction parameters with full spline-space coverage.")
                return
            end if
            do k = 1, dim_
                do i2 = 1, this%nc(2)
                    do i1 = 1, this%nc(1)
                        Xsol((i2-1)*this%nc(1)+i1,k) = &
                            sol21(i2,(k-1)*this%nc(1)+i1)
                    end do
                end do
            end do
            this%Xc = Xsol
            return
        end if

        bandwidth = this%degree(2)*this%nc(1) + this%degree(1)
        allocate(TtT(0:bandwidth, ncp), TtX(ncp, dim_), source=0.0_rk)
        do i = 1, n
            call basis_bspline_der(&
                map_parameter(Xt(i,1), this%knot1, this%nc(1), this%degree(1), this%wrap_parameters(1)),&
                this%knot1, this%nc(1), this%degree(1), 0, first1, N1)
            call basis_bspline_der(&
                map_parameter(Xt(i,2), this%knot2, this%nc(2), this%degree(2), this%wrap_parameters(2)),&
                this%knot2, this%nc(2), this%degree(2), 0, first2, N2)
            do a2 = 0, this%degree(2)
                if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                do a1 = 0, this%degree(1)
                    if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                    ia = (first2+a2-1)*this%nc(1) + first1 + a1
                    ba = N2(a2)*N1(a1)
                    do k = 1, dim_
                        TtX(ia,k) = TtX(ia,k) + ba*Xdata(i,k)
                    end do
                    do b2 = 0, this%degree(2)
                        if (first2+b2 < 1 .or. first2+b2 > this%nc(2)) cycle
                        do b1 = 0, this%degree(1)
                            if (first1+b1 < 1 .or. first1+b1 > this%nc(1)) cycle
                            ib = (first2+b2-1)*this%nc(1) + first1 + b1
                            bb = N2(b2)*N1(b1)
                            if (ia >= ib) TtT(ia-ib,ib) = TtT(ia-ib,ib) + ba*bb
                        end do
                    end do
                end do
            end do
        end do
        Xsol = solve_spd_banded(TtT, TtX)
        if (size(Xsol,1) /= ncp .or. size(Xsol,2) /= dim_) then
            allocate(design(n,ncp), source=0.0_rk)
            do i = 1, n
                call basis_bspline_der(&
                    map_parameter(Xt(i,1), this%knot1, this%nc(1), this%degree(1), this%wrap_parameters(1)),&
                    this%knot1, this%nc(1), this%degree(1), 0, first1, N1)
                call basis_bspline_der(&
                    map_parameter(Xt(i,2), this%knot2, this%nc(2), this%degree(2), this%wrap_parameters(2)),&
                    this%knot2, this%nc(2), this%degree(2), 0, first2, N2)
                do a2 = 0, this%degree(2)
                    if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                    do a1 = 0, this%degree(1)
                        if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                        ia = (first2+a2-1)*this%nc(1) + first1 + a1
                        design(i,ia) = N2(a2)*N1(a1)
                    end do
                end do
            end do
            Xsol = solve_least_squares_qr(design, Xdata(1:n,:))
        end if
        if (size(Xsol,1) /= ncp .or. size(Xsol,2) /= dim_) then
            call this%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "The tensor-product design is numerically rank deficient.",&
                location   = "lsq_fit_bspline",&
                suggestion = "Use parameters with full spline-space coverage or reduce the control counts.")
            return
        end if
        this%Xc = Xsol
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Alternately fit surface control points and positive rational weights.
    !!
    !! Control points and log-weights are updated from banded normal equations.
    !! Mean log-weight removal fixes projective scale at geometric mean weight
    !! one. Each trial weight vector is used to recompute its least-squares
    !! controls before the residual is tested. A trial is accepted only when
    !! its residual norm does not increase; otherwise controls and weights
    !! remain at the last accepted iterate and the damping is increased.
    !! Accepted trials reduce the damping. Log-weight steps outside the finite
    !! exponential range are rejected before exponentiation.
    !!
    !! With \(n=\mathtt{product(ndata)}\) and residuals
    !! \(r_{ik}=S_k(u_i,v_i)-X_{ik}\), the monitored quantity is
    !!
    !! \[
    !! \rho=\frac{\sqrt{\sum_{i,k}r_{ik}^2}}{n\,n_{coord}}.
    !! \]
    !!
    !! This is not an RMS residual. After an accepted trial, iteration stops
    !! when \(0\le\rho_{old}-\rho\le
    !! \mathtt{tol}\max(1,\rho_{old})\), or when
    !! \(\lVert r\rVert_2^2\le\epsilon_{\mathrm{mach}}
    !! \max(1,\lVert X\rVert_2^2)\). The latter recognizes a fit at the
    !! attainable accuracy of the normal-equation solves without accepting an
    !! increasing trial. The default tolerance is
    !! \(\sqrt{\epsilon_{\mathrm{mach}}}\). Exhausting `maxit` sets diagnostic
    !! code 109 and retains the last accepted controls and weights. Only
    !! controls, weights, and rational-state metadata are updated; cached `Xg`
    !! is not recomputed.
    !! Inputs must provide enough finite rows and at least two parameter columns
    !! and one coordinate column, with `ndata>=nc` componentwise, `maxit>0`,
    !! `tol>=0`, and finite, nonnegative regularization/damping values. The
    !! objective is nonconvex.
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
        class(nurbs_surface), intent(inout) :: this
        real(rk), intent(in), contiguous    :: Xt(:,:)
            !! Parameter pairs `[at least product(ndata),2]`.
        real(rk), intent(in), contiguous    :: Xdata(:,:)
            !! Data coordinates `[at least product(ndata),ncoord]`.
        integer,  intent(in)                :: ndata(2)
            !! Directional data counts `[2]`.
        integer,  intent(in),  optional     :: maxit
            !! Optional positive iteration limit; default 30.
        real(rk), intent(in),  optional     :: tol
            !! Optional nonnegative relative decrease tolerance; default `sqrt(epsilon(rk))`.
        real(rk), intent(in),  optional     :: lambda_xc
            !! Optional nonnegative control-point diagonal regularization.
        real(rk), intent(in),  optional     :: mu0
            !! Optional nonnegative initial log-weight damping.
        real(rk), intent(in),  optional     :: reg_logw
            !! Optional nonnegative additional log-weight regularization.
        real(rk), parameter :: damping_max = sqrt(huge(1.0_rk))
        real(rk), parameter :: log_weight_max = log(huge(1.0_rk)) - 1.0_rk
        real(rk), parameter :: log_weight_min = log(tiny(1.0_rk)) + 1.0_rk
        real(rk), allocatable :: TtT(:,:), TtX(:,:), v(:), Xsol(:,:)
        real(rk), allocatable :: Jband(:,:), Jrhs(:,:), delta_u(:,:)
        real(rk) :: N1(0:this%degree(1)), N2(0:this%degree(2))
        real(rk) :: Cg(size(Xdata, 2))
        real(rk) :: bval, cost_prev, cost_now, cost_sum, cost_reduction, damp, data_scale, denom
        real(rk) :: epss, ha, hb, resid, ta, tb
        real(rk) :: tol_, lamx_, mu, regw
        integer :: dim_, it, maxit_, n, ncp, i, j, k, a1, a2, b1, b2, ia, ib, first1, first2, bandwidth
        logical :: converged, trial_pending

        if (.not. this%err%ok) return

        dim_   = size(Xdata, 2)
        n      = ndata(1) * ndata(2)
        ncp    = this%nc(1) * this%nc(2)
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
                category   = "forcad_nurbs_surface",&
                message    = "Invalid nonlinear fitting controls.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Use maxit > 0 and finite, nonnegative tolerances and damping.")
            return
        end if
        mu = min(mu, damping_max)
        regw = min(regw, damping_max)

        if (this%nc(1) > ndata(1)) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Too few data points in dir-1 for requested control points.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Use nc(1) <= ndata(1).")
            return
        end if
        if (this%nc(2) > ndata(2)) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Too few data points in dir-2 for requested control points.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Use nc(2) <= ndata(2).")
            return
        end if

        if (n <= 0 .or. ncp < 2) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Invalid sizes for LSQ fitting.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Check ndata and nc.")
            return
        end if

        if (allocated(this%Wc)) then
            if (size(this%Wc) /= ncp) deallocate(this%Wc)
        end if
        if (.not. allocated(this%Wc)) allocate(this%Wc(ncp), source=1.0_rk)

        allocate(v(ncp))
        v = log(max(this%Wc, epss))
        v = v - sum(v)/real(ncp, rk)
        this%Wc = exp(v)
        this%rational = nonuniform_weights(this%Wc)

        bandwidth = this%degree(2)*this%nc(1) + this%degree(1)
        allocate(TtT(0:bandwidth, ncp), TtX(ncp, dim_))
        allocate(Jband(0:bandwidth, ncp), Jrhs(ncp,1))

        data_scale = 0.0_rk
        do i = 1, n
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
            do i = 1, n
                call basis_bspline_der(&
                    map_parameter(Xt(i,1), this%knot1, this%nc(1), this%degree(1), this%wrap_parameters(1)),&
                    this%knot1, this%nc(1), this%degree(1), 0, first1, N1)
                call basis_bspline_der(&
                    map_parameter(Xt(i,2), this%knot2, this%nc(2), this%degree(2), this%wrap_parameters(2)),&
                    this%knot2, this%nc(2), this%degree(2), 0, first2, N2)
                denom = 0.0_rk
                do a2 = 0, this%degree(2)
                    if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                    do a1 = 0, this%degree(1)
                        if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                        ia = (first2+a2-1)*this%nc(1) + first1 + a1
                        denom = denom + N2(a2)*N1(a1)*this%Wc(ia)
                    end do
                end do
                if (abs(denom) < epss) denom = sign(epss, denom)
                do a2 = 0, this%degree(2)
                    if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                    do a1 = 0, this%degree(1)
                        if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                        ia = (first2+a2-1)*this%nc(1) + first1 + a1
                        ta = N2(a2)*N1(a1)*this%Wc(ia)/denom
                        do k = 1, dim_
                            TtX(ia,k) = TtX(ia,k) + ta*Xdata(i,k)
                        end do
                        do b2 = 0, this%degree(2)
                            if (first2+b2 < 1 .or. first2+b2 > this%nc(2)) cycle
                            do b1 = 0, this%degree(1)
                                if (first1+b1 < 1 .or. first1+b1 > this%nc(1)) cycle
                                ib = (first2+b2-1)*this%nc(1) + first1 + b1
                                if (ia >= ib) then
                                    tb = N2(b2)*N1(b1)*this%Wc(ib)/denom
                                    TtT(ia-ib,ib) = TtT(ia-ib,ib) + ta*tb
                                end if
                            end do
                        end do
                    end do
                end do
            end do
            if (lamx_ > 0.0_rk) then
                do concurrent (j = 1:ncp)
                    TtT(0,j) = TtT(0,j) + lamx_
                end do
            end if
            Xsol = solve_spd_banded(TtT, TtX)
            if (size(Xsol,1) /= ncp .or. size(Xsol,2) /= dim_ .or. &
                any(.not. ieee_is_finite(Xsol))) then
                if (trial_pending) then
                    do concurrent (j = 1:ncp)
                        this%Wc(j) = exp(v(j))
                    end do
                    mu = min(damping_max, 10.0_rk*max(mu, sqrt(epsilon(1.0_rk))))
                    trial_pending = .false.
                    cycle
                end if
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_surface",&
                    message    = "Banded rational control-point solve failed.",&
                    location   = "lsq_fit_nurbs",&
                    suggestion = "Check parameter coverage or increase lambda_xc.")
                return
            end if

            cost_sum = 0.0_rk
            Jband = 0.0_rk
            Jrhs = 0.0_rk
            do i = 1, n
                call basis_bspline_der(&
                    map_parameter(Xt(i,1), this%knot1, this%nc(1), this%degree(1), this%wrap_parameters(1)),&
                    this%knot1, this%nc(1), this%degree(1), 0, first1, N1)
                call basis_bspline_der(&
                    map_parameter(Xt(i,2), this%knot2, this%nc(2), this%degree(2), this%wrap_parameters(2)),&
                    this%knot2, this%nc(2), this%degree(2), 0, first2, N2)
                denom = 0.0_rk
                do a2 = 0, this%degree(2)
                    if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                    do a1 = 0, this%degree(1)
                        if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                        ia = (first2+a2-1)*this%nc(1) + first1 + a1
                        denom = denom + N2(a2)*N1(a1)*this%Wc(ia)
                    end do
                end do
                if (abs(denom) < epss) denom = sign(epss, denom)
                Cg = 0.0_rk
                do a2 = 0, this%degree(2)
                    if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                    do a1 = 0, this%degree(1)
                        if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                        ia = (first2+a2-1)*this%nc(1) + first1 + a1
                        ta = N2(a2)*N1(a1)*this%Wc(ia)/denom
                        do k = 1, dim_
                            Cg(k) = Cg(k) + ta*Xsol(ia,k)
                        end do
                    end do
                end do
                do k = 1, dim_
                    resid = Cg(k) - Xdata(i,k)
                    cost_sum = cost_sum + resid*resid
                    do a2 = 0, this%degree(2)
                        if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                        do a1 = 0, this%degree(1)
                            if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                            ia = (first2+a2-1)*this%nc(1) + first1 + a1
                            bval = N2(a2)*N1(a1)
                            ha = bval*this%Wc(ia)*(Xsol(ia,k) - Cg(k))/denom
                            Jrhs(ia,1) = Jrhs(ia,1) + ha*resid
                            do b2 = 0, this%degree(2)
                                if (first2+b2 < 1 .or. first2+b2 > this%nc(2)) cycle
                                do b1 = 0, this%degree(1)
                                    if (first1+b1 < 1 .or. first1+b1 > this%nc(1)) cycle
                                    ib = (first2+b2-1)*this%nc(1) + first1 + b1
                                    if (ia >= ib) then
                                        hb = N2(b2)*N1(b1)*this%Wc(ib)*(Xsol(ib,k) - Cg(k))/denom
                                        Jband(ia-ib,ib) = Jband(ia-ib,ib) + ha*hb
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            cost_now = sqrt(cost_sum) / real(n*dim_, rk)

            if (trial_pending) then
                if (.not. ieee_is_finite(cost_now) .or. cost_now > cost_prev) then
                    do concurrent (j = 1:ncp)
                        this%Wc(j) = exp(v(j))
                    end do
                    mu = min(damping_max, 10.0_rk*max(mu, sqrt(epsilon(1.0_rk))))
                    trial_pending = .false.
                    cycle
                end if
                cost_reduction = cost_prev - cost_now
                this%Xc = Xsol
                v = log(this%Wc)
                v = v - sum(v)/real(ncp, rk)
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
                        category   = "forcad_nurbs_surface",&
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
            do concurrent (j = 1:ncp)
                Jband(0,j) = Jband(0,j) + damp
            end do

            delta_u = solve_spd_banded(Jband, Jrhs)
            if (size(delta_u,1) /= ncp .or. size(delta_u,2) /= 1 .or. &
                any(.not. ieee_is_finite(delta_u))) then
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_surface",&
                    message    = "Sparse rational weight solve failed.",&
                    location   = "lsq_fit_nurbs",&
                    suggestion = "Increase mu0 or reg_logw for the rational weight update.")
                return
            end if
            delta_u(:,1) = v - delta_u(:,1)
            delta_u(:,1) = delta_u(:,1) - sum(delta_u(:,1))/real(ncp, rk)
            if (any(.not. ieee_is_finite(delta_u(:,1))) .or. &
                maxval(delta_u(:,1)) > log_weight_max .or. &
                minval(delta_u(:,1)) < log_weight_min) then
                mu = min(damping_max, 10.0_rk*max(mu, sqrt(epsilon(1.0_rk))))
            else
                do concurrent (j = 1:ncp)
                    this%Wc(j) = exp(delta_u(j,1))
                end do
                trial_pending = .true.
            end if
        end do

        if (.not. converged) then
            if (trial_pending) then
                do concurrent (j = 1:ncp)
                    this%Wc(j) = exp(v(j))
                end do
            end if
            this%rational = nonuniform_weights(this%Wc)
            call this%err%set(&
                code       = 109,&
                severity   = 1,&
                category   = "forcad_nurbs_surface",&
                message    = "Nonlinear NURBS fitting did not converge.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Increase maxit or mu0, or relax tol for this data set.")
        end if
    end subroutine
    !===============================================================================
end module forcad_nurbs_surface
