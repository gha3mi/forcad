! author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Trivariate tensor-product B-spline and NURBS volumes.
!!
!! For degrees \(p\), \(q\), and \(r\), the rational volume mapping is
!!
!! \[\mathbf{V}(u,v,w)=
!!   \frac{\sum_{k=1}^{n_3}\sum_{j=1}^{n_2}\sum_{i=1}^{n_1}
!!   N_{i,p}(u)M_{j,q}(v)L_{k,r}(w)w_{ijk}\mathbf{P}_{ijk}}
!!   {\sum_{k=1}^{n_3}\sum_{j=1}^{n_2}\sum_{i=1}^{n_1}
!!   N_{i,p}(u)M_{j,q}(v)L_{k,r}(w)w_{ijk}}.\]
!!
!! Control subscripts in this module's formulas are one-based to match Fortran
!! storage. They are equivalent to conventional zero-based tensor indices after
!! subtracting one in each direction.
!!
!! Let \(B_{ijk}=N_{i,p}M_{j,q}L_{k,r}\) and
!! \(W=\sum_{ijk}B_{ijk}w_{ijk}\). The rational tensor basis
!! \(R_{ijk}=B_{ijk}w_{ijk}/W\) is nonnegative and forms a partition of unity on
!! the active domain because the weights are positive. Every positive-order
!! mixed derivative that exists therefore sums to zero over all basis functions.
!! Polynomial tensor-basis derivatives vanish when an order exceeds the degree
!! in its direction, but rational quotient derivatives generally need not.
!! Uniform positive weight scaling changes only the projective gauge. Rational
!! values scale nonzero weighted basis contributions; derivative recurrences
!! share a weight gauge. Homogeneous refinement selects a common representable
!! gauge for weights and weighted controls and restores the original weights
!! afterward. Geometric derivatives use nested division and never
!! explicitly form \(W^2\).
!!
!! Direction 1 varies fastest in every flattened tensor array:
!! `linear(i,j,k)=i+nc(1)*((j-1)+nc(2)*(k-1))`. Omitting weights, or using
!! equal weights, selects the polynomial tensor-product B-spline form.
!! Every nonconstant positive weight list uses rational normalization.
!!
!! Every direction accepts any finite, nondecreasing knot vector satisfying the
!! spline-space invariants, including uniform or nonuniform, clamped,
!! one-sided clamped, unclamped, and periodic-form vectors. A complete cyclic
!! knot extension together with repeated control and weight layers is verified
!! directionally by [[nurbs_volume:get_parameter_topology]]; malformed
!! periodic-form data remain valid bounded splines. The active interval
!! in direction `d` is `[knot_d(degree(d)+1),knot_d(nc(d)+1)]`. Parameter
!! wrapping is only an evaluation policy and does not create periodic geometry.
!! At an interior knot of multiplicity \(s_d\), the spline space has guaranteed
!! directional continuity \(C^{p_d-s_d}\), although special control data can be
!! smoother. Geometry and basis queries accept one trace selector per direction
!! (−1 left, 0 default, +1 right). The default takes the right trace internally.
!! At a closed endpoint both sides use the inward trace.
!!
!! Directional refinement applies an exact one-dimensional operator to every
!! homogeneous control line, preserving
!! \(\mathbf V_{\mathrm{old}}(u,v,w)=\mathbf V_{\mathrm{new}}(u,v,w)\) to
!! working precision on the active domain. Insertion and degree elevation
!! preserve a verified periodic direction. Other non-clamped directional
!! representations are converted by degree elevation to an equivalent
!! active-domain open-clamped form. Interior knot removal preserves periodic
!! directions by applying the bounded removal test to a guarded cyclic extension.
!! For the square physical Jacobian
!! \(\mathbf J=\partial\mathbf V/\partial(u,v,w)\),
!!
!! \[
!! \mathrm dV=|\det\mathbf J|\,\mathrm du\,\mathrm dv\,\mathrm dw,\qquad
!! \nabla_xR_A=\mathbf J^{-T}\nabla_{\boldsymbol{\xi}}R_A .
!! \]
!!
!! [[nurbs_volume]] supports evaluation, mixed basis derivatives, IGA element
!! and face data, directional refinement, volume quadrature, fitting, discrete
!! and local continuous nearest-point searches, transformations, and legacy VTK
!! export. IGES has no implemented trivariate export path in this API;
!! [[nurbs_volume:export_iges]] reports that condition through
!! [[nurbs_volume:err]].
!!
!! **Primary reference:** L. Piegl and W. Tiller, *The NURBS Book*, 2nd ed.,
!! Springer, 1997.
!! [doi:10.1007/978-3-642-59223-2](https://doi.org/10.1007/978-3-642-59223-2).
!! **Finite-range reference:** IEEE, *IEEE Standard for Floating-Point
!! Arithmetic*, IEEE Std 754-2019, 2019.
!! [IEEE 754-2019](https://standards.ieee.org/ieee/754/6210/).
module forcad_nurbs_volume

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    use forcad_kinds, only: rk
    use forcad_utils, only: valid_trace_side, tensor_ansatz, tensor_ansatz_element, tensor_basis_jet, tensor_measure
    use forcad_utils, only: elemConn_C0, ndgrid, compute_multiplicity, compute_knot_vector, basis_bspline_der, &
        geometry_derivative, tensor_basis_derivative_at, hexahedron_Xc, coarsen_control_lines, reduce_control_lines, &
        refine_control_net, get_array, spline_connectivity, spline_span_plan, spline_element, bezier_element_data, repelem, &
        unique, rotate_points, transform_geometry, export_vtk_legacy, valid_knot_vector, &
        knot_start, knot_end, knot_tolerance, active_knots, active_span_count, active_span, infer_knot_shape, &
        map_parameter, periodic_topology, &
        query_parameter_axis, valid_weights, nonuniform_weights, normalize_rational_basis, &
        unpack_refinement_controls, slice_control_net, isoparametric_grid, restrict_control_net, reparameterize_knots, &
        reorient_patch, map_parameter_values, nearest_sample_index, nearest_spline_point
    use forcad_nurbs_surface, only: nurbs_surface
    use fordebug, only: debug

    implicit none

    private
    public nurbs_volume
#ifndef FORCAD_PUBLIC_DOCS
    public compute_Tgc, compute_dTgc, volume_geometry_derivative
#endif

    !===============================================================================
    ! author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Mutable trivariate B-spline or NURBS volume mapping.
    !!
    !! Initialize the object with [[nurbs_volume:set]]. Directional arrays use
    !! entries 1, 2, and 3 for parameters \(u\), \(v\), and \(w\). Each
    !! `degree(d)` is a polynomial degree; order is `degree(d)+1`. Weights must
    !! be finite and strictly positive. Control coordinates are required to be
    !! finite and have at least one component; coordinate finiteness is a
    !! caller precondition in setter overloads that do not explicitly diagnose
    !! it.
    !!
    !! [[nurbs_volume:create]] samples a Cartesian product grid. Refinement is
    !! applied as one-dimensional homogeneous operators along a selected tensor
    !! direction. Physical IGA gradients, Hessians, and volume integration
    !! require exactly three physical coordinates and a nonsingular Jacobian.
    !!
    !! Successful setters discard sampled points and connectivity while keeping
    !! parameters for create() to refresh. Fixed-shape modify_Xc/modify_Wc
    !! operations discard samples and retain control and IGA connectivity.
    !!
    !! @warning Fitting and the separate
    !! `rotate_Xc` or `translate_Xc` operations do not automatically recompute
    !! every sampled or connectivity cache. [[nurbs_volume:transform]] instead
    !! keeps an existing `Xg` cache affine-consistent with `Xc`. Rebuild
    !! connectivity after a control-count or spline-space change.
    !! @endwarning
    !!
    !! Diagnostic codes are: 100 invalid argument or index; 101 shape mismatch;
    !! 102--105 missing geometry state; 106 underdetermined fitting problem; 107
    !! invalid parameters; 108 singular system or Jacobian; 109 iterative
    !! nonconvergence; and 110 unsupported trivariate IGES export.
    !!
    !! @note Getter generics return allocatable copies of object-owned arrays.
    !! @endnote
    type nurbs_volume
        real(rk), allocatable, private :: Xc(:,:)  !! Flattened control lattice `[product(nc),ncoord]`.
        real(rk), allocatable, private :: Xg(:,:)  !! Cached sampled points `[product(ng),ncoord]`.
        real(rk), allocatable, private :: Wc(:)    !! Positive flattened control weights `[product(nc)]`.
        real(rk), allocatable, private :: Xt1(:)   !! Cached parameters in direction 1 `[ng(1)]`.
        real(rk), allocatable, private :: Xt2(:)   !! Cached parameters in direction 2 `[ng(2)]`.
        real(rk), allocatable, private :: Xt3(:)   !! Cached parameters in direction 3 `[ng(3)]`.
        real(rk), allocatable, private :: Xt(:,:)  !! Flattened parameter grid `[product(ng),3]`.
        real(rk), allocatable, private :: knot1(:) !! Knot vector in direction 1.
        real(rk), allocatable, private :: knot2(:) !! Knot vector in direction 2.
        real(rk), allocatable, private :: knot3(:) !! Knot vector in direction 3.
        integer, private :: degree(3) = 0          !! Polynomial degrees \([p,q,r]\).
        integer, private :: nc(3) = 0              !! Control-point counts by direction.
        integer, private :: ng(3) = 0              !! Retained parameter-grid counts by direction.
        logical, private :: rational = .false.     !! True when nonuniform weights affect the basis.
        logical, private :: wrap_parameters(3) = .false. !! Modulo evaluation policy by direction.
        integer, allocatable, private :: elemConn_Xc_vis(:,:) !! Control-lattice VTK connectivity.
        integer, allocatable, private :: elemConn_Xg_vis(:,:) !! Sampled-volume VTK connectivity.
        integer, allocatable, private :: elemConn(:,:)        !! IGA element-to-control-point map.

        type(debug) :: err !! Recoverable diagnostic state for the most recent failed operation.
    contains
        procedure, private :: set1         !!> Set explicit directional knot vectors, control points, and weights.
        procedure, private :: set2         !!> Build knot vectors from breakpoints, degrees, and continuities.
        procedure, private :: set3         !!> Construct a Bezier or rational Bezier volume.
        procedure, private :: set4         !!> Construct an open uniform volume from degrees and sizes.
        generic :: set => set1, set2, set3, set4 !!> Replace volume geometry or initialize a fitting space.
        procedure :: create                !!> Sample a tensor grid and cache parameters, points, and connectivity.
        procedure, private :: cmp_Xg_scalar, cmp_Xg_vector
        generic :: cmp_Xg => cmp_Xg_scalar, cmp_Xg_vector !!> Evaluate one point or paired parameter rows.
        procedure :: slice                 !!> Extract an exact unsampled surface at one fixed parameter.
        procedure, private :: bezier_element, bezier_all
        generic :: bezier => bezier_element, bezier_all !!> Extract one exact Bezier element or decompose the patch.
        procedure :: restrict              !!> Extract an independent bounded patch on one parameter interval.
        procedure :: split                 !!> Split into two bounded patches at one interior parameter.
        procedure :: reverse               !!> Reverse one parameter direction without changing its active interval.
        procedure, private :: reparameterize_interval, reparameterize_axes
        generic :: reparameterize => reparameterize_interval, reparameterize_axes !!> Rescale or reorder parameter axes.
        procedure, private :: get_Xc_all   !!> Return the flattened control lattice.
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
        generic :: get_knot => get_knoti, get_knot_all !!> Query directional knot-vector data.
        procedure :: get_ng                !!> Return cached sample counts by direction.
        procedure, private :: get_nc_dir   !!> Return the control-point count in one direction.
        procedure, private :: get_nc_all   !!> Return control-point counts in all directions.
        generic :: get_nc => get_nc_all, get_nc_dir !!> Query directional control-point counts.
        procedure :: cmp_degree            !!> Infer all degrees from knots and control-lattice dimensions.
        procedure, private :: get_degree_all !!> Return all polynomial degrees.
        procedure, private :: get_degree_dir !!> Return the degree in one direction.
        generic :: get_degree => get_degree_all, get_degree_dir !!> Query polynomial degrees.
        procedure :: finalize              !!> Release all owned and derived volume data while preserving diagnostics.
        procedure :: cmp_elem_Xc_vis       !!> Build hexahedral connectivity for the control lattice.
        procedure :: cmp_elem_Xg_vis       !!> Build hexahedral connectivity for cached samples.
        procedure :: cmp_elem_Xth          !!> Build hexahedral connectivity for the parameter grid.
        procedure :: cmp_elem              !!> Build nonzero-span tensor-product IGA connectivity.
        procedure :: cmp_dof_map           !!> Map stored periodic controls to independent cyclic degrees of freedom.
        procedure :: get_elem_Xc_vis       !!> Return control-lattice visualization connectivity.
        procedure :: get_elem_Xg_vis       !!> Return sampled-volume visualization connectivity.
        procedure :: get_elem              !!> Return IGA element-to-control-point connectivity.
        procedure :: set_elem_Xc_vis       !!> Replace control-lattice visualization connectivity.
        procedure :: set_elem_Xg_vis       !!> Replace sampled-volume visualization connectivity.
        procedure :: set_elem              !!> Replace stored IGA connectivity.
        procedure :: export_Xc             !!> Export the control lattice in legacy VTK format.
        procedure :: export_Xg             !!> Export cached volume samples in legacy VTK format.
        procedure :: export_Xth            !!> Export the sampled parameter volume in legacy VTK format.
        procedure :: cmp_Xth_in_Xg         !!> Build embedded active isoparametric grid surfaces.
        procedure :: export_Xth_in_Xg      !!> Export parameter-grid cells embedded in physical space.
        procedure :: export_iges           !!> Report that trivariate IGES export is unsupported.
        procedure, private :: modify_Xc_scalar, modify_Xc_all
        generic :: modify_Xc => modify_Xc_scalar, modify_Xc_all !!> Edit one coordinate or the complete control array.
        procedure, private :: modify_Wc_scalar, modify_Wc_all
        generic :: modify_Wc => modify_Wc_scalar, modify_Wc_all !!> Edit one weight or the complete weight vector.
        procedure :: get_multiplicity      !!> Return multiplicities of all distinct stored knots in one direction.
        procedure :: get_continuity        !!> Return \(p_d-s\) metadata for each distinct stored knot.
        procedure :: cmp_nc                !!> Compute control-point counts implied by knots and degrees.
        procedure, private :: basis_vector !!> Evaluate all tensor-product basis functions at many points.
        procedure, private :: basis_scalar !!> Evaluate all tensor-product basis functions at one point.
        generic :: basis => basis_vector, basis_scalar !!> Evaluate the dense rational or polynomial basis.
        procedure, private :: derivative_vector !!> Evaluate parametric basis gradients at many points.
        procedure, private :: derivative_scalar !!> Evaluate the parametric basis gradient at one point.
        generic :: derivative => derivative_vector, derivative_scalar !!> Evaluate first parametric derivatives.
        procedure, private :: derivative2_vector !!> Evaluate packed basis Hessians at many parameter triples.
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
        procedure :: is_valid              !!> Report whether complete, internally consistent volume geometry is stored.
        procedure :: is_rational           !!> Report whether nonuniform weights affect the current basis.
        procedure :: get_parameter_wrapping !!> Report modulo parameter mapping in one direction.
        procedure :: get_parameter_topology !!> Classify one stored parameter topology.
        procedure :: put_to_nurbs          !!> Map structured points to the current parameter domain and cache them.
        procedure :: remove_knots          !!> Remove directional knots that pass an internal geometry test.
        procedure :: transform             !!> Apply one affine map to control points and cached geometry.
        procedure :: rotate_Xc             !!> Rotate control points in three-dimensional physical space.
        procedure :: rotate_Xg             !!> Rotate cached geometry points in three-dimensional physical space.
        procedure :: translate_Xc          !!> Translate the control lattice by a physical-space vector.
        procedure :: translate_Xg          !!> Translate cached geometry points by a physical-space vector.
        procedure :: nearest_point         !!> Select the nearest cached geometry sample.
        procedure :: nearest_point2        !!> Minimize squared distance with bounded Newton iterations and backtracking.
        procedure :: ansatz                !!> Evaluate element-local physical basis data and differential volume.
        procedure :: ansatz_element        !!> Evaluate all quadrature data for one analysis element.
        procedure :: cmp_volume            !!> Integrate mapped volume with elementwise tensor Gauss quadrature.
        procedure :: lsq_fit_bspline       !!> Fit polynomial control points by structured least squares.
        procedure :: lsq_fit_nurbs         !!> Fit control points and positive weights by damped nonlinear least squares.

        ! Faces
        procedure :: cmp_elemFace_Xc_vis   !!> Extract one face from a control-lattice visualization cell.
        procedure :: cmp_elemFace_Xg_vis   !!> Extract one face from a sampled-volume visualization cell.
        procedure :: cmp_elemFace          !!> Extract one face from an IGA element.
        procedure :: cmp_degreeFace        !!> Return tensor degrees retained on one volume face.

        ! Shapes
        procedure :: set_hexahedron        !!> Construct a trilinear hexahedral volume.
        procedure :: set_ring              !!> Construct an exact rational annular volume.
        procedure :: set_half_ring         !!> Construct an exact rational half-annular volume.
        procedure :: set_C                 !!> Construct the predefined C-shaped volume.
    end type
    !===============================================================================

    !> Evaluate polynomial or rational volume points from explicit spline data.
    !! This implementation generic is private; user code should normally call
    !! [[nurbs_volume:cmp_Xg]].
    interface compute_Xg
        module procedure compute_Xg_3d
        module procedure compute_Xg_3d_1point
    end interface

    !> Evaluate dense trivariate tensor-product basis values.
    !! Scalar overloads return `product(nc)` values; vector overloads return one
    !! row per parameter triple. Rational overloads normalize with `Wc`.
    !! For \(B_{ijk}=N_{i,p}(u)M_{j,q}(v)L_{k,r}(w)\),
    !! \(R_{ijk}=B_{ijk}w_{ijk}/\sum_{abc}B_{abc}w_{abc}\).
    interface compute_Tgc
        module procedure compute_Tgc_nurbs_3d_vector
        module procedure compute_Tgc_bspline_3d_vector
        module procedure compute_Tgc_nurbs_3d_scalar
        module procedure compute_Tgc_bspline_3d_scalar
    end interface

    !> Evaluate dense first parametric derivatives of the volume basis.
    !! The last derivative dimension follows parameter directions 1, 2, and 3.
    !! Optional `Tgc` is computed in the same pass.
    !! For \(\alpha\in\{u,v,w\}\), rational derivatives use
    !!
    !! \[
    !! R_{ijk,\alpha}
    !! =\frac{w_{ijk}(B_{ijk,\alpha}W-B_{ijk}W_{,\alpha})}{W^2}.
    !! \]
    interface compute_dTgc
        module procedure compute_dTgc_nurbs_3d_vector
        module procedure compute_dTgc_bspline_3d_vector
        module procedure compute_dTgc_nurbs_3d_scalar
        module procedure compute_dTgc_bspline_3d_scalar
    end interface

    !> Compute requested dense or selected basis data from one local tensor jet.
    interface compute_basis_data
        module procedure compute_basis_data_3d_vector
        module procedure compute_basis_data_3d_scalar
    end interface

contains


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Apply one affine transformation to the volume in place.
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
    !! the transformed control lattice represents the exact affine image of
    !! the original NURBS volume. Knots, weights, parameters, and connectivity
    !! are unchanged. A singular matrix remains an exact map but generally
    !! invalidates the volume as a regular IGA parameterization.
    pure subroutine transform(this, linear, translation, center)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: linear(:,:)
            !! Square linear map `[ncoord,ncoord]`.
        real(rk), intent(in), contiguous, optional :: translation(:)
            !! Translation vector `[ncoord]`; zero when absent.
        real(rk), intent(in), contiguous, optional :: center(:)
            !! Transformation center `[ncoord]`; the origin when absent.

        call transform_geometry(this%Xc,this%Xg,linear,translation,center,this%err,"forcad_nurbs_volume")
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report whether the defining volume geometry is complete and consistent.
    !!
    !! For directional degrees \(p_d\) and control counts \(n_d\), every knot
    !! vector must satisfy
    !! \(\operatorname{size}(U_d)=n_d+p_d+1\). The flattened control lattice
    !! must contain \(n_1n_2n_3\) finite points with at least one physical
    !! coordinate. Explicit weights, when stored, must comprise the same number
    !! of finite positive values and agree with the cached rational
    !! classification. Sample and connectivity caches are derived state and
    !! are not inspected.
    pure function is_valid(this) result(ok)
        class(nurbs_volume), intent(in) :: this
            !! Volume geometry to inspect.
        logical :: ok

        ok = .false.
        if (.not. this%err%ok) return
        if (.not. allocated(this%knot1) .or. .not. allocated(this%knot2) .or. &
            .not. allocated(this%knot3) .or. .not. allocated(this%Xc)) return
        if (.not. valid_knot_vector(this%knot1, this%nc(1), this%degree(1))) return
        if (.not. valid_knot_vector(this%knot2, this%nc(2), this%degree(2))) return
        if (.not. valid_knot_vector(this%knot3, this%nc(3), this%degree(3))) return
        if (size(this%Xc,1) /= product(this%nc) .or. size(this%Xc,2) < 1) return
        if (.not. all(ieee_is_finite(this%Xc))) return

        if (allocated(this%Wc)) then
            if (size(this%Wc) /= product(this%nc) .or. .not. valid_weights(this%Wc)) return
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
    !> Evaluate all quadrature points of one volume analysis element.
    !!
    !! This batched counterpart of [[nurbs_volume:ansatz]] locates the three
    !! active knot spans and constructs the tensor Gauss rule once. The
    !! local-basis index varies fastest in
    !! `Tgc(nlocal,nquadrature)` and
    !! `dTgc_dXg(nlocal,3,nquadrature)`. Differential volumes are returned in
    !! `dV(nquadrature)`. Optional outputs have shapes
    !! `Xg(3,nquadrature)` and
    !! `d2Tgc_dXg2(nlocal,3,3,nquadrature)`.
    !!
    !! Allocatable outputs have `intent(inout)` so repeated element calls reuse
    !! correctly shaped storage. The rational quotient, inverse-Jacobian
    !! gradient, Hessian chain rule, measure, and strict-orientation test are
    !! mathematically identical to [[nurbs_volume:ansatz]]. Any failure zeros
    !! the complete element result before setting `err`.
    pure subroutine ansatz_element(this, ie, Tgc, dTgc_dXg, dV, ngauss, Xg, d2Tgc_dXg2, strict, err, dXg_dXt, &
        dmeasure_dXc, dmeasure_dlogWc, dgradient_dXc, dgradient_dlogWc, djacobian_dXc, djacobian_dlogWc)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: ie
            !! One-based direction-1-fastest active element index.
        real(rk), allocatable, intent(inout) :: Tgc(:,:), dTgc_dXg(:,:,:), dV(:)
            !! Reusable basis `[nlocal,nq]`, physical gradients `[nlocal,3,nq]`, weighted measures `[nq]`.
        integer, intent(in), optional :: ngauss(3)
            !! Positive directional quadrature counts; defaults to `degree+1`.
        real(rk), allocatable, target, intent(inout), optional :: Xg(:,:), d2Tgc_dXg2(:,:,:,:)
            !! Optional physical points `[3,nq]` and physical Hessians `[nlocal,3,3,nq]`.
        logical, intent(in), optional :: strict
            !! Reject orientation reversal in addition to singularity.
        type(debug), intent(out), optional :: err
            !! Caller-owned diagnostic; bypasses and preserves the object diagnostic.
        real(rk), allocatable, target, intent(inout), optional :: dXg_dXt(:,:,:)
            !! Optional raw Jacobians `[3,3,nq]`, coordinate first, parameter direction second.
        real(rk), allocatable, target, intent(inout), optional :: dmeasure_dXc(:,:,:), dmeasure_dlogWc(:,:)
            !! Weighted measure derivatives for the active Euclidean controls and log weights.
        real(rk), allocatable, intent(inout), optional :: dgradient_dXc(:,:,:,:,:), dgradient_dlogWc(:,:,:,:)
            !! Physical basis gradient derivatives with respect to controls `(A,i,B,j,q)` and log weights `(A,i,B,q)`.
        real(rk), allocatable, intent(inout), optional :: djacobian_dXc(:,:,:), djacobian_dlogWc(:,:,:,:)
            !! Raw Jacobian derivatives for controls (A,alpha,q) and log weights (i,alpha,A,q).

        call tensor_ansatz_element(this%Xc,this%Wc,this%rational,this%degree,this%nc,this%knot1,this%knot2, &
            this%err,err,ie,Tgc,dTgc_dXg,dV,ngauss,strict,Xg,dXg_dXt,d2Tgc_dXg2, &
            dmeasure_dXc,dmeasure_dlogWc,dgradient_dXc,dgradient_dlogWc,this%knot3, &
            djacobian_dXc=djacobian_dXc,djacobian_dlogWc=djacobian_dlogWc)
    end subroutine ansatz_element
    !===============================================================================




    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Replace the volume with explicit directional knots and control data.
    !!
    !! `Xc` is direction-1-fastest. If `degree` is absent, an open-clamped
    !! candidate inferred from the initial knot runs takes precedence. If that
    !! candidate is invalid, exhaustive degree/count search must find exactly
    !! one valid triple. Pass explicit degrees whenever this convention is
    !! ambiguous or should not determine the interpretation of non-open data.
    !! `Xc` must contain at least one flattened control row and one finite
    !! coordinate component. Coordinate, knot, and positive-weight invariants
    !! are validated before any core state is replaced.
    !!
    pure subroutine set1(this, knot1, knot2, knot3, Xc, Wc, degree, wrap_parameters)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: knot1(:)
            !! Complete knot vector in direction 1.
        real(rk), intent(in), contiguous :: knot2(:)
            !! Complete knot vector in direction 2.
        real(rk), intent(in), contiguous :: knot3(:)
            !! Complete knot vector in direction 3.
        real(rk), intent(in), contiguous :: Xc(:,:)
            !! Flattened control lattice `[product(nc),ncoord]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional positive flattened weights `[product(nc)]`.
        integer, intent(in), contiguous, optional :: degree(:)
            !! Optional polynomial degrees `[3]`.
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
            !! Optional modulo evaluation flags `[3]`.
        integer :: degree_(3), nc_(3), ncp
        logical :: ok
        logical :: wrap_parameters_(3)

        if (.not. this%err%ok) return

        if (size(Xc,1) < 1 .or. size(Xc,2) < 1) then
            call this%err%set(&
                code       = 101,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Control points need at least one row and one coordinate.",&
                location   = "set1",&
                suggestion = "Provide Xc with shape [product(nc),ncoord], with positive dimensions.")
            return
        end if
        if (.not. all(ieee_is_finite(Xc))) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Control-point coordinates must all be finite.",&
                location   = "set1",&
                suggestion = "Remove NaN and infinite values from Xc before calling set().")
            return
        end if

        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            if (size(wrap_parameters) < 3) then
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Volume parameter wrapping needs one value per parametric direction.",&
                    location   = "set1",&
                    suggestion = "Pass one logical value per parametric direction.")
                return
            end if
            wrap_parameters_ = wrap_parameters(1:3)
        end if

        if (present(degree)) then
            if (size(degree) < 3) then
                call this%err%set(&
                    code       = 103,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Explicit volume degree must provide one entry per parametric direction.",&
                    location   = "set1",&
                    suggestion = "Pass degree=[p1,p2,p3].")
                return
            end if
            degree_ = degree(1:3)
            nc_ = [size(knot1) - degree_(1) - 1, size(knot2) - degree_(2) - 1, size(knot3) - degree_(3) - 1]
            ok = valid_knot_vector(knot1, nc_(1), degree_(1)) .and. &
                valid_knot_vector(knot2, nc_(2), degree_(2)) .and. &
                valid_knot_vector(knot3, nc_(3), degree_(3))
        else
            call infer_knot_shape(knot1, knot2, knot3, size(Xc,1), degree_, nc_, ok)
        end if
        if (.not. ok .or. product(nc_) /= size(Xc,1)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid or ambiguous knot vectors for the supplied control points.",&
                location   = "set1",&
                suggestion = "Use valid knot vectors and pass degree=[p1,p2,p3] for non-open or ambiguous knot vectors.")
            return
        end if
        ncp = product(nc_)

        if (present(Wc)) then
            if (size(Wc) /= ncp) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Weights length mismatch: size(Wc) must equal number of control points.",&
                    location   = "set1",&
                    suggestion = "Provide Wc with size(Wc) == nc(1)*nc(2)*nc(3).")
                return
            else if (.not. valid_weights(Wc)) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Invalid rational weights: all weights must be finite and positive.",&
                    location   = "set1",&
                    suggestion = "Use strictly positive finite Wc values for rational NURBS geometry.")
                return
            end if
        end if

        call invalidate_samples(this,space_changed=.true.)
        this%knot1 = knot1
        this%knot2 = knot2
        this%knot3 = knot3
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
    !> Construct all spline directions from breakpoints and continuities.
    !! Directional continuity arrays are converted to knot multiplicities by
    !! `compute_knot_vector`. Optional controls and weights must match the
    !! resulting tensor-product count. In direction \(d\), breakpoints must be
    !! finite and strictly increasing. Values in \([-1,p_d-1]\) retain a
    !! breakpoint with multiplicity \(p_d-c_i\); `-1` gives full multiplicity
    !! \(p_d+1\), while `c_i=p_d` omits that breakpoint. All generated knot
    !! vectors are validated before object state changes. If `Xc` is absent,
    !! any previous control lattice is released and the result is a
    !! tensor-product spline space for subsequent fitting, not complete volume
    !! geometry. Optional weights are retained only when explicitly supplied.
    pure subroutine set2(this, Xth_dir1, Xth_dir2, Xth_dir3, degree, continuity1, continuity2, continuity3, Xc, Wc)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xth_dir1(:)
            !! Finite, strictly increasing breakpoints in direction 1.
        real(rk), intent(in), contiguous :: Xth_dir2(:)
            !! Finite, strictly increasing breakpoints in direction 2.
        real(rk), intent(in), contiguous :: Xth_dir3(:)
            !! Finite, strictly increasing breakpoints in direction 3.
        integer, intent(in), contiguous :: degree(:)
            !! Polynomial degrees `[3]`.
        integer, intent(in), contiguous :: continuity1(:)
            !! One value per direction-1 breakpoint in `-1:degree(1)`; the upper value omits it.
        integer, intent(in), contiguous :: continuity2(:)
            !! One value per direction-2 breakpoint in `-1:degree(2)`; the upper value omits it.
        integer, intent(in), contiguous :: continuity3(:)
            !! One value per direction-3 breakpoint in `-1:degree(3)`; the upper value omits it.
        real(rk), intent(in), contiguous, optional :: Xc(:,:)
            !! Optional flattened control lattice `[product(computed_nc),ncoord]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional positive weights `[product(computed_nc)]`.
        real(rk), allocatable :: knot1_(:), knot2_(:), knot3_(:)
        integer :: degree_(3), nc_(3), ncp

        if (.not. this%err%ok) return

        if (size(degree) < 3) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Volume degree must provide one entry per parametric direction.",&
                location   = "set2",&
                suggestion = "Pass degree=[p1,p2,p3].")
            return
        end if

        degree_ = degree(1:3)
        knot1_ = compute_knot_vector(Xth_dir1, degree_(1), continuity1)
        knot2_ = compute_knot_vector(Xth_dir2, degree_(2), continuity2)
        knot3_ = compute_knot_vector(Xth_dir3, degree_(3), continuity3)
        nc_ = [size(knot1_) - degree_(1) - 1, size(knot2_) - degree_(2) - 1, size(knot3_) - degree_(3) - 1]
        if (.not. valid_knot_vector(knot1_, nc_(1), degree_(1)) .or. &
            .not. valid_knot_vector(knot2_, nc_(2), degree_(2)) .or. &
            .not. valid_knot_vector(knot3_, nc_(3), degree_(3))) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid knot vector generated from parameter nodes, degree, and continuity.",&
                location   = "set2",&
                suggestion = "Use nondecreasing parameter nodes and knot multiplicity <= degree+1.")
            return
        end if
        ncp = product(nc_)

        if (present(Xc)) then
            call this%set1(knot1_, knot2_, knot3_, Xc, Wc, degree=degree_)
            return
        end if

        if (present(Wc)) then
            if (size(Wc) /= ncp) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume", &
                    message    = "Weights size mismatch in set2",&
                    location   = "set2", &
                    suggestion = "size(Wc) must equal nc(1)*nc(2)*nc(3).")
                return
            else if (.not. valid_weights(Wc)) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Invalid rational weights: all weights must be finite and positive.",&
                    location   = "set2",&
                    suggestion = "Use strictly positive finite Wc values for rational NURBS geometry.")
                return
            end if
        end if

        call invalidate_samples(this,space_changed=.true.)
        this%knot1 = knot1_
        this%knot2 = knot2_
        this%knot3 = knot3_
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
    !> Construct a tensor-product Bezier or rational Bezier volume on \([0,1]^3\).
    !! Directional degrees are `nc(1:3)-1`; endpoint knots are fully clamped.
    pure subroutine set3(this, nc, Xc, Wc)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), contiguous :: nc(:)
            !! Control-point counts `[3]`.
        real(rk), intent(in), contiguous :: Xc(:,:)
            !! Direction-1-fastest control lattice `[product(nc),ncoord]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional positive weights `[product(nc)]`.
        integer :: nc_(3), ncp, i
        real(rk), allocatable :: knot1(:), knot2(:), knot3(:)

        if (.not. this%err%ok) return

        if (size(nc) < 3) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Bezier volume control-point count must provide one entry per parametric direction.",&
                location   = "set3",&
                suggestion = "Pass nc=[nc1,nc2,nc3].")
            return
        end if

        nc_ = nc(1:3)
        ncp = product(nc_)
        if (any(nc_ <= 0) .or. size(Xc,1) /= ncp .or. size(Xc,2) < 1) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid Bezier volume dimensions.",&
                location   = "set3",&
                suggestion = "Require nc > 0 per direction, size(Xc,1) == product(nc), and at least one coordinate dimension.")
            return
        end if

        ! A Bezier direction has p+1 copies of each endpoint and no interior knots.
        knot1=[(0.0_rk,i=1,nc_(1)),(1.0_rk,i=1,nc_(1))]
        knot2=[(0.0_rk,i=1,nc_(2)),(1.0_rk,i=1,nc_(2))]
        knot3=[(0.0_rk,i=1,nc_(3)),(1.0_rk,i=1,nc_(3))]
        call this%set1(knot1,knot2,knot3,Xc,Wc,degree=nc_-1)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct a clamped open-uniform volume on \([0,1]^3\).
    !! Interior knots are uniformly spaced and simple in each direction.
    pure subroutine set4(this, degree, nc, Xc, Wc)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), contiguous :: degree(:)
            !! Polynomial degrees `[3]`, each smaller than `nc`.
        integer, intent(in), contiguous :: nc(:)
            !! Directional control-point counts `[3]`.
        real(rk), intent(in), contiguous :: Xc(:,:)
            !! Direction-1-fastest control lattice `[product(nc),ncoord]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional positive weights `[product(nc)]`.
        integer :: i
        real(rk), allocatable :: knot1(:), knot2(:), knot3(:)

        if (.not. this%err%ok) return

        if (size(degree) < 3 .or. size(nc) < 3) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid volume rank for generated knot-vector construction.",&
                location   = "set4",&
                suggestion = "Provide degree and nc arrays with at least three entries.")
            return
        end if

        if (any(degree(1:3) < 0) .or. any(nc(1:3) <= 0) .or. any(degree(1:3) >= nc(1:3)) .or. &
            size(Xc, 1) /= nc(1)*nc(2)*nc(3) .or. size(Xc, 2) < 1) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid volume dimensions for generated knot-vector construction.",&
                location   = "set4",&
                suggestion = "Require degree >= 0, nc > degree, size(Xc,1) == product(nc), and at least one coordinate.")
            return
        end if

        ! Open-uniform knots have p+1 repeated endpoints and nc-p-1 simple interiors.
        knot1=[(0.0_rk,i=0,degree(1)), &
            (real(i,rk)/real(nc(1)-degree(1),rk),i=1,nc(1)-degree(1)-1), &
            (1.0_rk,i=0,degree(1))]
        knot2=[(0.0_rk,i=0,degree(2)), &
            (real(i,rk)/real(nc(2)-degree(2),rk),i=1,nc(2)-degree(2)-1), &
            (1.0_rk,i=0,degree(2))]
        knot3=[(0.0_rk,i=0,degree(3)), &
            (real(i,rk)/real(nc(3)-degree(3),rk),i=1,nc(3)-degree(3)-1), &
            (1.0_rk,i=0,degree(3))]
        call this%set1(knot1,knot2,knot3,Xc,Wc,degree=degree(1:3))
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Sample and cache volume geometry.
    !!
    !! Paired `Xt` rows take precedence over all directional arguments. Otherwise
    !! explicit directional values override resolutions, with omitted directions
    !! reusing their stored vectors. Calling without arguments reevaluates the
    !! retained parameter grid. Rejected input preserves every existing cache.
    pure subroutine create(this, res1, res2, res3, Xt1, Xt2, Xt3, Xt)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: res1, res2, res3
            !! Positive directional grid sizes, overridden by explicit directional values.
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:), Xt3(:)
            !! Directional grid values; omitted directions reuse the stored vectors.
        real(rk), intent(in), contiguous, optional :: Xt(:,:)
            !! Paired rows with at least three columns; later columns are retained but ignored.
        real(rk), allocatable :: parameters(:,:), points(:,:), axis1(:), axis2(:), axis3(:)
        integer, allocatable :: connectivity(:,:)
        integer :: counts(3)
        logical :: valid, paired

        if (.not. this%err%ok) return
        if (.not. this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_volume", &
                message="Sampling requires valid knots, controls, and positive weights.", &
                location="create",suggestion="Initialize valid volume geometry before sampling.")
            return
        end if
        paired=present(Xt)
        if (present(Xt)) then
            valid=size(Xt,1)>0 .and. size(Xt,2)>=3
            if (valid) valid=all(ieee_is_finite(Xt(:,1:3)))
            if (valid) then
                parameters=Xt
                counts=[size(Xt,1),1,1]
            end if
        else
            call prepare_query_grid(this,"create",parameters,counts,valid,res1,res2,res3,Xt1,Xt2,Xt3, &
                axis1,axis2,axis3)
        end if
        if (.not. valid) then
            call this%err%set(code=107,severity=1,category="forcad_nurbs_volume", &
                message="Sampling parameters must be nonempty and finite in all three directions.", &
                location="create",suggestion="Pass paired Xt rows or valid directional parameters and positive counts.")
            return
        end if
        if (this%rational) then
            points=compute_Xg(parameters,this%knot1,this%knot2,this%knot3,this%degree,this%nc,counts, &
                this%Xc,this%wrap_parameters,Wc=this%Wc)
        else
            points=compute_Xg(parameters,this%knot1,this%knot2,this%knot3,this%degree,this%nc,counts, &
                this%Xc,this%wrap_parameters)
        end if
        if (.not. all(ieee_is_finite(points))) then
            call this%err%set(code=108,severity=1,category="forcad_nurbs_volume", &
                message="Geometry sampling produced a nonfinite point.",location="create", &
                suggestion="Rescale the geometry or choose numerically representable parameters.")
            return
        end if
        connectivity=elemConn_C0(counts(1),counts(2),counts(3),1,1,1)
        if (paired) then
            if (allocated(this%Xt1)) deallocate(this%Xt1)
            if (allocated(this%Xt2)) deallocate(this%Xt2)
            if (allocated(this%Xt3)) deallocate(this%Xt3)
        else if (allocated(axis1)) then
            call move_alloc(axis1,this%Xt1)
            call move_alloc(axis2,this%Xt2)
            call move_alloc(axis3,this%Xt3)
        end if
        call move_alloc(parameters,this%Xt)
        call move_alloc(points,this%Xg)
        call move_alloc(connectivity,this%elemConn_Xg_vis)
        this%ng=counts
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate a point or mixed parametric derivative using the first three parameters.
    !! Later entries are ignored. Unset geometry, active diagnostics or invalid
    !! parameters return an allocated empty result without changing diagnostics.
    pure function cmp_Xg_scalar(this, Xt, order, side) result(Xg)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in), optional :: side(:) !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: order(:) !! Nonnegative parametric derivative order; omitted means position.
        real(rk), allocatable :: Xg(:), rows(:,:)

        rows=cmp_Xg_vector(this,reshape(Xt,[1,size(Xt)]),order,side)
        Xg=reshape(rows,[size(rows)])
    end function cmp_Xg_scalar

    !> Evaluate paired rows, returning `[npoint,ncoord]` without changing caches.
    !! The first three columns must be finite; later columns are ignored.
    !! Strided input is accepted. Valid empty input returns `[0,ncoord]`;
    !! invalid input or unset geometry returns `[0,0]` without changing diagnostics.
    pure function cmp_Xg_vector(this, Xt, order, side) result(Xg)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in), optional :: side(:) !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in) :: Xt(:,:)
        integer, intent(in), optional :: order(:) !! Nonnegative order; omitted means position.
        real(rk), allocatable :: Xg(:,:)

        call volume_geometry_derivative(this,Xt,order,Xg,side)
    end function cmp_Xg_vector

    !> Evaluate patch geometry, optionally retaining a common binary exponent per derivative row.
    pure subroutine volume_geometry_derivative(patch,Xt,order,Xg,side,row_power)
        class(nurbs_volume), intent(in) :: patch
        real(rk), intent(in) :: Xt(:,:)
        integer, optional, intent(in) :: order(:),side(:)
        real(rk), allocatable, intent(out) :: Xg(:,:)
        integer, optional, allocatable, intent(out) :: row_power(:)
        integer :: selected(3)
        logical :: valid

        selected=0
        if (present(order)) then
            selected=-1
            if (size(order)==3) selected=order
        end if
        valid=all(selected>=0) .and. patch%err%ok .and. allocated(patch%Xc) .and. &
            allocated(patch%knot1) .and. allocated(patch%knot2) .and. allocated(patch%knot3)
        if (valid) valid=size(Xt,2) >= 3 .and. valid_trace_side(side,3)
        if (valid) valid=all(ieee_is_finite(Xt(:,1:3)))
        if (.not. valid) then
            allocate(Xg(0,0))
            if (present(row_power)) allocate(row_power(0))
            return
        end if
        if (any(selected>0)) then
            call geometry_derivative(Xt(:,1),patch%nc,patch%degree,patch%knot1,patch%Xc,patch%Wc, &
                selected,patch%wrap_parameters,Xg,knot2=patch%knot2,Xt2=Xt(:,2), &
                knot3=patch%knot3,Xt3=Xt(:,3),side=side,row_power=row_power)
        else if (patch%rational) then
            Xg = compute_Xg(Xt,patch%knot1,patch%knot2,patch%knot3,patch%degree,patch%nc, &
                Xc=patch%Xc,wrap_parameters=patch%wrap_parameters,Wc=patch%Wc,side=side)
        else
            Xg = compute_Xg(Xt,patch%knot1,patch%knot2,patch%knot3,patch%degree,patch%nc, &
                Xc=patch%Xc,wrap_parameters=patch%wrap_parameters,side=side)
        end if
        if (present(row_power).and.all(selected==0)) allocate(row_power(size(Xg,1)),source=0)
    end subroutine

    !> Extract the exact isoparametric surface at fixed `Xt` in direction `dir`.
    !! Retained parameters keep their original order, knots, degrees and wrapping.
    !! Homogeneous control contraction supports interior and unclamped boundary
    !! traces without sampling. Polynomial sources produce polynomial surfaces.
    !! The source is unchanged; failures belong to the returned surface diagnostic.
    pure function slice(this, dir, Xt) result(surface)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
            !! Fixed parameter direction, from one through three.
        real(rk), intent(in) :: Xt
            !! Finite parameter; the source direction's wrapping policy applies.
        type(nurbs_surface) :: surface
        real(rk), allocatable :: Xc(:,:), Wc(:), knot(:)
        integer :: retained(2)

        call check_extraction(this,dir,surface%err,"slice")
        if (.not.surface%err%ok) return
        knot=this%get_knot(dir)
        call slice_control_net(this%Xc,this%Wc,this%nc,dir,knot,this%degree(dir),Xt, &
            this%wrap_parameters(dir),Xc,Wc,surface%err)
        if (.not.surface%err%ok) return
        retained=pack([1,2,3],[1,2,3]/=dir)
        call surface%set(knot1=this%get_knot(retained(1)),knot2=this%get_knot(retained(2)), &
            Xc=Xc,Wc=Wc,degree=this%degree(retained),wrap_parameters=this%wrap_parameters(retained))
    end function slice

    !> Validate geometry-only extraction without changing the source diagnostic.
    pure subroutine check_extraction(this,dir,err,location)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        type(debug), intent(out) :: err
        character(*), intent(in) :: location
        err=this%err
        if (.not.err%ok) return
        if (dir<1 .or. dir>3) then
            call err%set(code=100,severity=1,category="forcad_nurbs_volume", &
                message="The extraction direction must be one, two or three.",location=location, &
                suggestion="Choose dir=1, dir=2 or dir=3.")
        else if (.not.this%is_valid()) then
            call err%set(code=102,severity=1,category="forcad_nurbs_volume", &
                message="Extraction requires complete valid volume geometry.",location=location, &
                suggestion="Initialize knots, controls and any rational weights before extraction.")
        end if
    end subroutine check_extraction

    !> Reconstruct one tensor direction; unchanged knot arrays are supplied independently.
    pure subroutine set_direction(this,dir,knots,other1,other2,Xc,Wc,degree,wrapping)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: dir,degree(3)
        real(rk), intent(in), contiguous :: knots(:),other1(:),other2(:),Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        logical, intent(in) :: wrapping(3)
        select case(dir)
        case(1)
            call this%set(knots,other1,other2,Xc,Wc,degree=degree,wrap_parameters=wrapping)
        case(2)
            call this%set(other1,knots,other2,Xc,Wc,degree=degree,wrap_parameters=wrapping)
        case default
            call this%set(other1,other2,knots,Xc,Wc,degree=degree,wrap_parameters=wrapping)
        end select
    end subroutine set_direction

    !> Reverse a parameter direction, including exterior knots and control layers.
    !! On continuous spans V_rev(u)=V(lo+hi-u). At a discontinuity the new
    !! right trace is the original left trace. Wrapping and explicit weights remain.
    pure function reverse(this,dir) result(patch)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        type(nurbs_volume) :: patch
        integer :: axes(3)

        call check_extraction(this,dir,patch%err,"reverse")
        if (.not.patch%err%ok) return
        axes=[1,2,3]
        axes(dir)=-dir
        call orient(this,axes,patch)
    end function reverse

    !> New axis d takes source axis abs(axes(d)); a negative entry reverses its original interval.
    pure function reparameterize_axes(this,axes) result(patch)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: axes(:)
        type(nurbs_volume) :: patch

        call check_extraction(this,1,patch%err,"reparameterize")
        if (patch%err%ok) call orient(this,axes,patch)
    end function reparameterize_axes

    pure subroutine orient(this,axes,patch)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: axes(:)
        type(nurbs_volume), intent(inout) :: patch
        integer :: source(3)

        call reorient_patch(this%nc,this%degree,axes,this%knot1,this%Xc,this%Wc,patch%knot1,patch%Xc,patch%Wc,patch%err, &
            this%knot2,patch%knot2,this%knot3,patch%knot3)
        if (.not.patch%err%ok) return
        source=abs(axes)
        patch%degree=this%degree(source)
        patch%nc=this%nc(source)
        patch%wrap_parameters=this%wrap_parameters(source)
        patch%rational=this%rational
    end subroutine orient

    !> Rescale the complete selected knot vector to the finite interval `[a,b]`.
    !! The geometry is unchanged under affine parameter correspondence; each
    !! derivative in this direction gains the inverse domain-scale factor.
    pure function reparameterize_interval(this,dir,a,b) result(patch)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), intent(in) :: a,b
        type(nurbs_volume) :: patch
        real(rk), allocatable :: knot(:),knots(:)
        integer :: retained(2)

        call check_extraction(this,dir,patch%err,"reparameterize")
        if (.not.patch%err%ok) return
        knot=this%get_knot(dir)
        call reparameterize_knots(knot,this%nc(dir),this%degree(dir),a,b,knots,patch%err, &
            category="forcad_nurbs_volume")
        if (.not.patch%err%ok) return
        retained=pack([1,2,3],[1,2,3]/=dir)
        call set_direction(patch,dir,knots,this%get_knot(retained(1)),this%get_knot(retained(2)), &
            this%Xc,this%Wc,this%degree,this%wrap_parameters)
    end function reparameterize_interval

    !> Extract one exact element; optional polynomial extraction satisfies N=C*B.
    !! Element and support IDs use direction-1-fastest stored-control ordering.
    !! Bounds have shape [2,3]. The piece retains literal span coordinates,
    !! with clamped knots, right/left endpoint traces and disabled wrapping.
    !! Caller err bypasses the source diagnostic; source geometry and caches are untouched.
    pure subroutine bezier_element(this, ie, piece, C, ids, bounds, err)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: ie
        type(nurbs_volume), intent(out), optional :: piece
        real(rk), allocatable, intent(out), optional :: C(:,:), bounds(:,:)
        integer, allocatable, intent(out), optional :: ids(:)
        type(debug), intent(out), optional :: err

        call bezier_extract(this,ie,piece,C,ids,bounds,err)
    end subroutine bezier_element


    !> Share extraction between direct queries and a decomposition's precomputed span plan.
    pure subroutine bezier_extract(this, ie, piece, C, ids, bounds, err, spans)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: ie
        integer, intent(in), optional :: spans(:)
        type(nurbs_volume), intent(out), optional :: piece
        real(rk), allocatable, intent(out), optional :: C(:,:), bounds(:,:)
        integer, allocatable, intent(out), optional :: ids(:)
        type(debug), intent(out), optional :: err
        real(rk), allocatable :: Xc(:,:), Wc(:)
        integer, allocatable :: nodes(:)
        real(rk) :: limits(2,3)
        integer :: p(3)
        type(debug) :: status

        call bezier_check(this,status,present(err))
        if (status%ok) then
            p=this%degree
            call bezier_element_data(this%knot1,this%nc,p,ie,this%Xc,this%Wc, &
                present(piece),Xc,Wc,nodes,limits,status,C=C,knot2=this%knot2,knot3=this%knot3,rational=this%rational,spans=spans)
            if (status%ok .and. present(piece)) then
                call piece%set(knot1=repelem(limits(:,1),[p(1)+1,p(1)+1]), &
                    knot2=repelem(limits(:,2),[p(2)+1,p(2)+1]), &
                    knot3=repelem(limits(:,3),[p(3)+1,p(3)+1]), &
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
        class(nurbs_volume), intent(in) :: this
        type(nurbs_volume), allocatable, intent(out) :: pieces(:)
        type(debug), intent(out), optional :: err
        type(debug) :: status
        integer, allocatable :: plan(:,:)
        integer :: ie, d, element, ne(3), selected(3)

        call bezier_check(this,status,present(err))
        if (status%ok) then
            plan=spline_span_plan(this%nc,this%degree,this%knot1,knot2=this%knot2,knot3=this%knot3)
            ne=count(plan>0,dim=1)
            allocate(pieces(product(ne)))
            do ie=1,size(pieces)
                element=ie-1
                do d=1,3
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
        class(nurbs_volume), intent(in) :: this
        type(debug), intent(out) :: err
        logical, intent(in) :: bypass

        err=this%err
        if (bypass) call err%reset()
        if (.not.err%ok) return
        if (.not.allocated(this%Xc) .or. &
            .not.allocated(this%knot1) .or. &
            .not.allocated(this%knot2) .or. &
            .not.allocated(this%knot3)) call err%set( &
            code=103,severity=1,category="forcad_nurbs_volume",location="bezier", &
            message="Spline geometry is not initialized.",suggestion="Set knots and controls before Bezier extraction.")
    end subroutine bezier_check


    !> Extract the exact bounded patch on the literal interval `[a,b]`.
    !! Parameters and degrees are unchanged; only the selected direction is
    !! clamped and stops wrapping. At discontinuous cuts the patch keeps the
    !! interior one-sided traces. Geometry is independently owned and unsampled.
    pure function restrict(this,dir,a,b) result(patch)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), intent(in) :: a,b
            !! Finite ordered bounds inside the original active interval; no modulo mapping.
        type(nurbs_volume) :: patch
        real(rk), allocatable :: knot(:),knots(:),Xc(:,:),Wc(:)
        logical :: wrapping(3)
        integer :: retained(2)

        call check_extraction(this,dir,patch%err,"restrict")
        if (.not.patch%err%ok) return
        knot=this%get_knot(dir)
        call restrict_control_net(this%Xc,this%Wc,this%nc,dir,knot,this%degree(dir),a,b,knots,Xc,Wc,patch%err)
        if (.not.patch%err%ok) return
        wrapping=this%wrap_parameters
        wrapping(dir)=.false.
        retained=pack([1,2,3],[1,2,3]/=dir)
        call set_direction(patch,dir,knots,this%get_knot(retained(1)),this%get_knot(retained(2)), &
            Xc,Wc,this%degree,wrapping)
    end function restrict

    !> Split at a strictly interior literal parameter; a failure empties both outputs.
    !! The left patch keeps the trace from below and the right patch from above.
    pure subroutine split(this,dir,Xt,left,right)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), intent(in) :: Xt
        type(nurbs_volume), intent(out) :: left,right
        type(debug) :: err
        real(rk), allocatable :: knot(:)

        call check_extraction(this,dir,err,"split")
        if (err%ok) then
            knot=this%get_knot(dir)
            left=this%restrict(dir,knot_start(knot,this%nc(dir),this%degree(dir)),Xt)
            right=this%restrict(dir,Xt,knot_end(knot,this%nc(dir),this%degree(dir)))
            if (.not.left%err%ok) err=left%err
            if (.not.right%err%ok .and. err%ok) err=right%err
        end if
        if (err%ok) return
        call left%finalize()
        call right%finalize()
        left%err=err
        right%err=err
    end subroutine split
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocated copy of the flattened control lattice.
    !!
    !! The result has shape `[product(nc),ncoord]` and uses
    !! direction-1-fastest ordering. An active diagnostic or unset lattice
    !! produces an allocated `[0,0]` result.
    pure function get_Xc_all(this) result(Xc)
        class(nurbs_volume), intent(in) :: this
        real(rk), allocatable :: Xc(:,:)

        call get_array(this%Xc,Xc,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return one point from the flattened control lattice.
    !!
    !! `n` is the one-based direction-1-fastest control index. Invalid or
    !! unavailable data produces an allocated empty result.
    pure function get_Xci(this, n) result(Xc)
        class(nurbs_volume), intent(in) :: this
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
    !! The elemental query returns zero for an invalid control or coordinate
    !! index and while the object has an active diagnostic.
    pure elemental function get_Xcid(this, n, dir) result(Xc)
        class(nurbs_volume), intent(in) :: this
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
        class(nurbs_volume), intent(in) :: this
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
        class(nurbs_volume), intent(in) :: this
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
        class(nurbs_volume), intent(in) :: this
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
    !! The result has size `product(nc)`. Polynomial volumes, unset
    !! weights, and active diagnostics produce an allocated empty result.
    pure function get_Wc_all(this) result(Wc)
        class(nurbs_volume), intent(in) :: this
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
        class(nurbs_volume), intent(in) :: this
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
    !! Valid directions are 1 through 3. For a Cartesian cache created from
    !! directional vectors, the result contains the corresponding vector.
    !! A row-wise cache created by `put_to_nurbs`, an invalid direction, absent
    !! directional data, or an active diagnostic produces an allocated empty
    !! result.
    pure function get_Xt(this, dir) result(Xt)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), allocatable :: Xt(:)

        select case(dir)
        case(1)
            call get_array(this%Xt1,Xt,this%err%ok)
        case(2)
            call get_array(this%Xt2,Xt,this%err%ok)
        case(3)
            call get_array(this%Xt3,Xt,this%err%ok)
        case default
            allocate(Xt(0))
        end select
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the three directional sample counts.
    !!
    !! A Cartesian cache reports its three directional counts. A row-wise cache
    !! reports `[npoint,1,1]`. The result is zero before sampling and while a
    !! diagnostic is active or geometry samples were invalidated by an edit.
    pure function get_ng(this) result(ng)
        class(nurbs_volume), intent(in) :: this
        integer :: ng(3)

        ng = 0
        if (.not. this%err%ok .or. .not. allocated(this%Xg)) return

        ng = this%ng
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Inspect directional knot sizes without allocating getter copies.
    pure function knot_sizes(this) result(lengths)
        class(nurbs_volume), intent(in) :: this
        integer :: lengths(3)
        lengths=0
        if (allocated(this%knot1)) lengths(1)=size(this%knot1)
        if (allocated(this%knot2)) lengths(2)=size(this%knot2)
        if (allocated(this%knot3)) lengths(3)=size(this%knot3)
    end function knot_sizes

    !> Infer one or all polynomial degrees from stored spline dimensions.
    !!
    !! When a control count is known, direction \(d\) uses
    !! \(p_d=m_d-n_d-1\), where \(m_d\) is the knot-vector size and \(n_d\)
    !! the control count. Otherwise the leading multiplicity is interpreted as
    !! \(p_d+1\). Optional `dir` restricts the update to one direction.
    pure subroutine cmp_degree(this, dir)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: dir
        integer :: d,first,last,inferred(3),lengths(3)

        if (.not.this%err%ok) return
        first=1
        last=3
        if (present(dir)) then
            if (dir<1 .or. dir>3) then
                call this%err%set(code=100,severity=1,category="forcad_nurbs_volume",location="cmp_degree", &
                    message="Invalid direction for computing degree.",suggestion="Choose dir=1, dir=2 or dir=3.")
                return
            end if
            first=dir
            last=dir
        end if
        lengths=knot_sizes(this)
        inferred=this%degree
        do d=first,last
            if (lengths(d)==0) then
                call this%err%set(code=103,severity=1,category="forcad_nurbs_volume",location="cmp_degree", &
                    message="Knot vector is not set.",suggestion="Call set(...) before computing degree.")
                return
            end if
            if (this%nc(d)>0) then
                inferred(d)=lengths(d)-this%nc(d)-1
            else
                select case(d)
                case(1)
                    inferred(d)=compute_multiplicity(this%knot1,this%knot1(1))-1
                case(2)
                    inferred(d)=compute_multiplicity(this%knot2,this%knot2(1))-1
                case default
                    inferred(d)=compute_multiplicity(this%knot3,this%knot3(1))-1
                end select
            end if
        end do
        this%degree=inferred
    end subroutine cmp_degree
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return all polynomial degrees \([p,q,r]\).
    !!
    !! An active diagnostic produces zeros.
    pure function get_degree_all(this) result(degree)
        class(nurbs_volume), intent(in) :: this
        integer :: degree(3)

        degree = 0
        if (.not. this%err%ok) return

        degree(1) = this%degree(1)
        degree(2) = this%degree(2)
        degree(3) = this%degree(3)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the polynomial degree in one parameter direction.
    !!
    !! Valid directions are 1 through 3. Invalid input or an active diagnostic
    !! returns zero.
    pure elemental function get_degree_dir(this, dir) result(degree)
        class(nurbs_volume), intent(in) :: this
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
    !! Valid directions are 1 through 3. Invalid input, unavailable knots, or
    !! an active diagnostic produces an allocated empty result.
    pure function get_knot_all(this, dir) result(knot)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), allocatable :: knot(:)

        select case(dir)
        case(1)
            call get_array(this%knot1,knot,this%err%ok)
        case(2)
            call get_array(this%knot2,knot,this%err%ok)
        case(3)
            call get_array(this%knot3,knot,this%err%ok)
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
        class(nurbs_volume), intent(in) :: this
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
        else if (dir == 3) then
            if (allocated(this%knot3)) then
                if (i >= 1 .and. i <= size(this%knot3)) knot = this%knot3(i)
            end if
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Release all owned volume arrays and reset derived state.
    !!
    !! The diagnostic object is deliberately preserved so callers can inspect
    !! an error after cleanup.
    pure elemental subroutine finalize(this)
        class(nurbs_volume), intent(inout) :: this
        if (allocated(this%Xc)) deallocate(this%Xc)
        if (allocated(this%Xg)) deallocate(this%Xg)
        if (allocated(this%Wc)) deallocate(this%Wc)
        if (allocated(this%Xt1)) deallocate(this%Xt1)
        if (allocated(this%Xt2)) deallocate(this%Xt2)
        if (allocated(this%Xt3)) deallocate(this%Xt3)
        if (allocated(this%Xt)) deallocate(this%Xt)
        if (allocated(this%knot1)) deallocate(this%knot1)
        if (allocated(this%knot2)) deallocate(this%knot2)
        if (allocated(this%knot3)) deallocate(this%knot3)
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
    !> Build structured hexahedral connectivity for the control lattice.
    !!
    !! Optional `p=[p_1,p_2,p_3]` selects visualization cell order;
    !! the default is trilinear. Numbering follows direction-1-fastest storage.
    pure function cmp_elem_Xc_vis(this, p) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), contiguous, optional :: p(:)

        if (.not. this%err%ok) return

        if (present(p)) then
            elemConn = elemConn_C0(this%nc(1), this%nc(2), this%nc(3), p(1), p(2), p(3))
        else
            elemConn = elemConn_C0(this%nc(1), this%nc(2), this%nc(3), 1, 1, 1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build structured hexahedral connectivity for cached geometry samples.
    !!
    !! Optional `p=[p_1,p_2,p_3]` selects visualization cell order;
    !! the default is trilinear. Numbering follows direction-1-fastest storage.
    pure function cmp_elem_Xg_vis(this, p) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), contiguous, optional :: p(:)

        if (.not. this%err%ok .or. .not. allocated(this%Xg)) then
            allocate(elemConn(0,0))
            return
        end if

        if (present(p)) then
            elemConn = elemConn_C0(this%ng(1), this%ng(2), this%ng(3), p(1), p(2), p(3))
        else
            elemConn = elemConn_C0(this%ng(1), this%ng(2), this%ng(3), 1, 1, 1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build hexahedral connectivity for the active parameter-span grid.
    !!
    !! Vertices are triples of distinct active knots. Optional `p` selects
    !! cell order and defaults to `[1,1,1]`.
    pure function cmp_elem_Xth(this, p) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, allocatable :: elemConn(:,:)
        integer, intent(in), contiguous, optional :: p(:)

        if (.not. this%err%ok) return

        if (present(p)) then
            elemConn = elemConn_C0(active_span_count(this%knot1, this%nc(1), this%degree(1)) + 1, &
                active_span_count(this%knot2, this%nc(2), this%degree(2)) + 1, &
                active_span_count(this%knot3, this%nc(3), this%degree(3)) + 1, p(1), p(2), p(3))
        else
            elemConn = elemConn_C0(active_span_count(this%knot1, this%nc(1), this%degree(1)) + 1, &
                active_span_count(this%knot2, this%nc(2), this%degree(2)) + 1, &
                active_span_count(this%knot3, this%nc(3), this%degree(3)) + 1, 1, 1, 1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export the control lattice as legacy VTK hexahedral cells.
    !!
    !! Stored visualization connectivity is used when available; otherwise
    !! trilinear connectivity is generated. Optional `point_data` is
    !! point-associated data with matching `field_names`.
    impure subroutine export_Xc(this, filename, point_data, field_names, encoding)
        class(nurbs_volume), intent(inout) :: this
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
                category   = "forcad_nurbs_volume",&
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

        call export_vtk_legacy(filename=filename, points=this%Xc, elemConn=elemConn, vtkCellType=12, &
                               point_data=point_data, field_names=field_names, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export cached physical samples as legacy VTK hexahedral cells.
    !!
    !! Call [[nurbs_volume:create]] before this operation. Stored visualization
    !! connectivity is used when available; otherwise trilinear connectivity is
    !! generated.
    impure subroutine export_Xg(this, filename, point_data, field_names, encoding)
        class(nurbs_volume), intent(inout) :: this
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
                category   = "forcad_nurbs_volume",&
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

        call export_vtk_legacy(filename=filename, points=this%Xg, elemConn=elemConn, vtkCellType=12, &
                               point_data=point_data, field_names=field_names, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export the active tensor-product parameter grid in parameter space.
    !!
    !! Vertices are triples of distinct active knots and cells span adjacent
    !! values. The written coordinates are \((u,v,w)\), not physical points.
    impure subroutine export_Xth(this, filename, point_data, field_names, encoding)
        class(nurbs_volume), intent(in) :: this
        character(len=*), intent(in) :: filename
        real(rk), intent(in), contiguous, optional :: point_data(:,:)
        character(len=*), intent(in), contiguous, optional :: field_names(:)
        character(len=*), intent(in), optional :: encoding
        integer, allocatable :: elemConn(:,:)
        real(rk), allocatable :: Xth(:,:), Xth1(:), Xth2(:), Xth3(:)

        if (.not. this%err%ok) return

        elemConn = this%cmp_elem_Xth()
        Xth1 = active_knots(this%knot1, this%nc(1), this%degree(1))
        Xth2 = active_knots(this%knot2, this%nc(2), this%degree(2))
        Xth3 = active_knots(this%knot3, this%nc(3), this%degree(3))
        call ndgrid(Xth1, Xth2, Xth3, Xth)

        call export_vtk_legacy(filename=filename, points=Xth, elemConn=elemConn, vtkCellType=12, &
                               point_data=point_data, field_names=field_names, encoding=encoding)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export active isoparametric grid surfaces embedded in physical space.
    !!
    !! Every distinct active knot in one direction defines a surface sampled
    !! over all active spans in the other directions. Optional `res` sets
    !! the minimum directional sampling density and defaults to 10.
    impure subroutine export_Xth_in_Xg(this, filename, res, encoding)
        class(nurbs_volume), intent(inout) :: this
        character(len=*),   intent(in) :: filename
        integer, intent(in), optional  :: res   ! min points per span (>=2)
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
    !> Build active isoparametric grid surfaces embedded in physical space.
    !!
    !! The returned rows of `elemConn` are VTK polyline connectivities into
    !! `Xg`. No file is written.
    impure subroutine cmp_Xth_in_Xg(this, Xg, elemConn, res)
        class(nurbs_volume), intent(inout) :: this
        real(rk), allocatable, intent(out) :: Xg(:,:)
        integer, allocatable, intent(out) :: elemConn(:,:)
        integer, intent(in), optional :: res
        real(rk), allocatable :: Xt_all(:,:)

        if (.not.this%err%ok) return
        if (.not.allocated(this%Xc)) then
            call this%err%set(code=102,severity=1,category="forcad_nurbs_volume", &
                message="Control points are not set.",location="export_Xth_in_Xg", &
                suggestion="Call set(...) before exporting.")
            return
        end if
        if (.not.allocated(this%knot1) .or. .not.allocated(this%knot2) .or. .not.allocated(this%knot3)) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_volume", &
                message="Knot vectors are not set.",location="export_Xth_in_Xg", &
                suggestion="Call set(...) before exporting.")
            return
        end if
        if (size(this%Xc,2)<2 .or. size(this%Xc,2)>3 .or. &
            .not.valid_knot_vector(this%knot1,this%nc(1),this%degree(1)) .or. &
            .not.valid_knot_vector(this%knot2,this%nc(2),this%degree(2)) .or. &
            .not.valid_knot_vector(this%knot3,this%nc(3),this%degree(3))) then
            call this%err%set(code=100,severity=1,category="forcad_nurbs_volume", &
                message="The grid requires valid knots and two or three physical coordinates.", &
                location="export_Xth_in_Xg",suggestion="Check knot vectors and geometry dimension before exporting.")
            return
        end if
        call isoparametric_grid(this%knot1,this%knot2,this%nc,this%degree,Xt_all,elemConn,res,knot3=this%knot3)
        Xg=this%cmp_Xg(Xt_all)
    end subroutine cmp_Xth_in_Xg
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report that direct IGES export of a trivariate NURBS is unsupported.
    !!
    !! IGES entities 126 and 128 represent parametric curves and surfaces; this
    !! API does not approximate a volume or silently replace it by its boundary.
    !! A valid initialized volume and nonempty filename therefore produce
    !! diagnostic code 110 and no file. Use [[nurbs_volume:export_Xg]] for a VTK
    !! representation.
    !!
    !! **Format reference:** U.S. National Bureau of Standards, *Initial
    !! Graphics Exchange Specification (IGES), Version 4.0*, NBSIR 88-3813
    !! (1988), entities 126 and 128.
    !! [NIST publication](https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir88-3813.pdf).
    impure subroutine export_iges(this, filename)
        class(nurbs_volume), intent(inout) :: this
            !! Volume whose diagnostic state receives the recoverable error.
        character(len=*), intent(in) :: filename
            !! Prospective output path; it must be nonempty.

        if (.not. this%err%ok) return

        if (.not. allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Control points are not set.",&
                location   = "export_iges",&
                suggestion = "Call set(...) first before exporting.")
            return
        end if

        if (len_trim(filename) == 0) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Output filename is empty.",&
                location   = "export_iges",&
                suggestion = "Pass a non-empty filename.")
            return
        end if

        call this%err%set(&
            code       = 110,&
            severity   = 1,&
            category   = "forcad_nurbs_volume",&
            message    = "IGES export for trivariate NURBS volumes is not implemented.",&
            location   = "export_iges",&
            suggestion = "Use VTK export for volumes, or export a documented boundary-surface representation once implemented.")
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Edit one finite coordinate without scanning or copying the control lattice.
    !! Successful edits discard physical samples while retaining their parameters
    !! for a subsequent `create()` call.
    pure subroutine modify_Xc_scalar(this,X,num,dir)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: X
        integer, intent(in) :: num,dir
        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) then
            call this%err%set(code=102,severity=1,category="forcad_nurbs_volume", &
                message="Control points are not set.",location="modify_Xc",suggestion="Initialize the control lattice first.")
            return
        end if
        if (num<1 .or. num>size(this%Xc,1) .or. dir<1 .or. dir>size(this%Xc,2) .or. .not. ieee_is_finite(X)) then
            call this%err%set(code=100,severity=1,category="forcad_nurbs_volume", &
                message="Control indices and the replacement coordinate must be valid.",location="modify_Xc", &
                suggestion="Use an existing control and coordinate index with a finite value.")
            return
        end if
        this%Xc(num,dir)=X
        call invalidate_samples(this)
    end subroutine modify_Xc_scalar

    !> Replace all finite control coordinates atomically at the current shape.
    pure subroutine modify_Xc_all(this,Xc)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: Xc(:,:)
        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) then
            call this%err%set(code=102,severity=1,category="forcad_nurbs_volume", &
                message="Control points are not set.",location="modify_Xc",suggestion="Initialize the control lattice first.")
            return
        end if
        if (any(shape(Xc)/=shape(this%Xc)) .or. .not. all(ieee_is_finite(Xc))) then
            call this%err%set(code=101,severity=1,category="forcad_nurbs_volume", &
                message="Replacement controls must have the stored shape and finite coordinates.",location="modify_Xc", &
                suggestion="Keep the control count and coordinate dimension unchanged.")
            return
        end if
        this%Xc=Xc
        call invalidate_samples(this)
    end subroutine modify_Xc_all

    !> Replace one positive weight; polynomial volumes begin with implicit unit weights.
    pure subroutine modify_Wc_scalar(this,W,num)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: W
        integer, intent(in) :: num
        if (.not. this%err%ok) return
        if (any(this%nc<1) .or. num<1 .or. num>product(this%nc)) then
            call this%err%set(code=100,severity=1,category="forcad_nurbs_volume", &
                message="The weight index must belong to an initialized spline space.",location="modify_Wc", &
                suggestion="Initialize the volume and use 1 <= num <= product(nc).")
            return
        end if
        if (.not. ieee_is_finite(W) .or. W<=0.0_rk) then
            call this%err%set(code=105,severity=1,category="forcad_nurbs_volume", &
                message="Replacement weights must be finite and strictly positive.",location="modify_Wc", &
                suggestion="Use a positive finite weight.")
            return
        end if
        if (.not. allocated(this%Wc)) allocate(this%Wc(product(this%nc)),source=1.0_rk)
        this%Wc(num)=W
        this%rational=nonuniform_weights(this%Wc)
        call invalidate_samples(this)
    end subroutine modify_Wc_scalar

    !> Replace all weights atomically without changing the control lattice or knots.
    pure subroutine modify_Wc_all(this,Wc)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: Wc(:)
        if (.not. this%err%ok) return
        if (any(this%nc<1) .or. size(Wc)/=product(this%nc) .or. .not. valid_weights(Wc)) then
            call this%err%set(code=105,severity=1,category="forcad_nurbs_volume", &
                message="Replacement weights must match the control count and be finite and positive.",location="modify_Wc", &
                suggestion="Provide one strictly positive weight for every stored control.")
            return
        end if
        this%Wc=Wc
        this%rational=nonuniform_weights(this%Wc)
        call invalidate_samples(this)
    end subroutine modify_Wc_all

    !> Drop derived geometry, optionally including topology, while keeping parameters.
    pure subroutine invalidate_samples(this,space_changed)
        class(nurbs_volume), intent(inout) :: this
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
    !> Return multiplicities of the distinct knots in one direction.
    !!
    !! Result entries correspond, in nondecreasing order, to the unique values
    !! in the complete selected knot vector. Invalid input or unavailable state
    !! produces an allocated empty result.
    pure function get_multiplicity(this, dir) result(m)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        integer, allocatable :: m(:)

        if (this%err%ok) then
            if (dir==1 .and. allocated(this%knot1)) m=compute_multiplicity(this%knot1)
            if (dir==2 .and. allocated(this%knot2)) m=compute_multiplicity(this%knot2)
            if (dir==3 .and. allocated(this%knot3)) m=compute_multiplicity(this%knot3)
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
    !! continuity; a particular volume mapping can be smoother. Endpoint entries
    !! and entries outside the active interval describe representation
    !! multiplicity, not continuity between two active spans.
    pure function get_continuity(this, dir) result(c)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        integer, allocatable :: c(:)

        c=this%get_multiplicity(dir)
        c=this%get_degree(dir)-c
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Infer one or all directional control counts from knots and degree.
    !!
    !! Each selected direction uses \(n_d=m_d-p_d-1\), where \(m_d\) is the
    !! knot-vector size. Optional `dir` must be in 1 through 3.
    pure subroutine cmp_nc(this, dir)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: dir
        integer :: d,first,last,inferred(3),lengths(3)

        if (.not.this%err%ok) return
        first=1
        last=3
        if (present(dir)) then
            if (dir<1 .or. dir>3) then
                call this%err%set(code=103,severity=1,category="forcad_nurbs_volume",location="cmp_nc", &
                    message="Invalid direction for computing nc.",suggestion="Choose dir=1, dir=2 or dir=3.")
                return
            end if
            first=dir
            last=dir
        end if
        lengths=knot_sizes(this)
        inferred=this%nc
        do d=first,last
            if (lengths(d)==0) then
                call this%err%set(code=103,severity=1,category="forcad_nurbs_volume",location="cmp_nc", &
                    message="Knot vector is not set.",suggestion="Call set(...) before computing nc.")
                return
            end if
            inferred(d)=lengths(d)-this%degree(d)-1
        end do
        this%nc=inferred
    end subroutine cmp_nc
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return control-point counts \([n_1,n_2,n_3]\).
    !!
    !! An active diagnostic produces zeros.
    pure function get_nc_all(this) result(nc)
        class(nurbs_volume), intent(in) :: this
        integer :: nc(3)

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
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        integer :: nc

        nc = 0
        if (.not. this%err%ok) return

        if (dir == 1) then

            if (allocated(this%knot1)) nc = this%nc(1)

        else if (dir == 2) then

            if (allocated(this%knot2)) nc = this%nc(2)

        else if (dir == 3) then

            if (allocated(this%knot3)) nc = this%nc(3)
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Validate the state required by all volume derivative evaluators.
    pure subroutine check_derivative_state(this, location, valid, side)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: side(:)
        character(len=*), intent(in) :: location
        logical, intent(out) :: valid

        valid = .false.
        if (.not. this%err%ok) return
        if (.not. valid_trace_side(side,3)) then
            call this%err%set(code=100,severity=1,category="forcad_nurbs_volume",location=location, &
                message="Invalid knot trace selector.",suggestion="Use -1, 0, or 1 for each parameter direction.")
            return
        end if
        if (.not. allocated(this%knot1) .or. .not. allocated(this%knot2) .or. .not. allocated(this%knot3)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Volume derivatives require three valid knot vectors and control counts.",&
                location   = location,&
                suggestion = "Set the volume before evaluating basis derivatives.")
            return
        end if
        valid = all(this%degree >= 0) .and. all(this%nc > this%degree) .and. &
            size(this%knot1) == this%nc(1) + this%degree(1) + 1 .and. &
            size(this%knot2) == this%nc(2) + this%degree(2) + 1 .and. &
            size(this%knot3) == this%nc(3) + this%degree(3) + 1
        if (valid .and. allocated(this%Wc)) then
            valid = size(this%Wc) == product(this%nc)
        end if
        if (.not. valid) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Volume derivative state is inconsistent.",&
                location   = location,&
                suggestion = "Use matching finite knots, degrees, control counts, and positive weights.")
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> Validate scalar query coordinates after the shared spline-state checks.
    pure subroutine check_query_point(this,Xt,location,valid, side)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: side(:)
        real(rk), intent(in) :: Xt(:)
        character(len=*), intent(in) :: location
        logical, intent(out) :: valid

        call check_derivative_state(this,location,valid,side)
        if (.not. valid) return
        if (size(Xt) < 3) then
            valid = .false.
            call this%err%set(code=100,severity=1,category="forcad_nurbs_volume", &
                message="A volume query requires at least three parameters.",location=location, &
                suggestion="Provide Xt=[u,v,w]; later entries are ignored.")
        else if (.not. all(ieee_is_finite(Xt(1:3)))) then
            valid = .false.
            call this%err%set(code=107,severity=1,category="forcad_nurbs_volume", &
                message="Volume query parameters must be finite.",location=location, &
                suggestion="Provide finite values in Xt(1:3).")
        end if
    end subroutine check_query_point

    !> Select retained paired rows or build one directional tensor grid locally.
    pure subroutine prepare_query_grid(this,location,Xt,ng,valid,res1,res2,res3,Xt1,Xt2,Xt3,axis1,axis2,axis3, side)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: side(:)
        character(len=*), intent(in) :: location
        real(rk), allocatable, intent(out) :: Xt(:,:)
        integer, intent(out) :: ng(3)
        logical, intent(out) :: valid
        integer, intent(in), optional :: res1,res2,res3
        real(rk), intent(in), optional :: Xt1(:),Xt2(:),Xt3(:)
        real(rk), allocatable, intent(out), optional :: axis1(:),axis2(:),axis3(:)
        real(rk), allocatable :: params1(:),params2(:),params3(:)
        logical :: grid_requested

        ng=0
        call check_derivative_state(this,location,valid,side)
        if (.not. valid) return
        grid_requested=present(res1).or.present(res2).or.present(res3).or.present(Xt1).or.present(Xt2).or.present(Xt3)
        if (.not. grid_requested .and. allocated(this%Xt)) then
            Xt=this%Xt
            ng=this%ng
            return
        end if
        call query_parameter_axis(this%knot1,this%nc(1),this%degree(1),this%Xt1,params1,valid,res1,Xt1)
        if (valid) call query_parameter_axis(this%knot2,this%nc(2),this%degree(2),this%Xt2,params2,valid,res2,Xt2)
        if (valid) call query_parameter_axis(this%knot3,this%nc(3),this%degree(3),this%Xt3,params3,valid,res3,Xt3)
        if (.not. valid) then
            call this%err%set(code=107,severity=1,category="forcad_nurbs_volume", &
                message="Each query direction requires finite nonempty parameters or a positive resolution.", &
                location=location,suggestion="Provide directional values/counts or create a parameter grid first.")
            return
        end if
        ng=[size(params1),size(params2),size(params3)]
        call ndgrid(params1,params2,params3,Xt)
        if (present(axis1)) call move_alloc(params1,axis1)
        if (present(axis2)) call move_alloc(params2,axis2)
        if (present(axis3)) call move_alloc(params3,axis3)
    end subroutine prepare_query_grid

    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate one arbitrary mixed volume basis derivative on its local support.
    pure subroutine evaluate_basis_order(this,Xt,order,values,Wc,side,first_active,elem)
        class(nurbs_volume), intent(in) :: this
        real(rk), intent(in) :: Xt(:)
        integer, intent(in) :: order(:)
        real(rk), intent(out) :: values(:)
        real(rk), contiguous, intent(in), optional :: Wc(:)
        integer, intent(in), optional :: side(:)
        integer, intent(out), optional :: first_active(:)
        integer, contiguous, intent(in), optional :: elem(:)

        call tensor_basis_derivative_at(Xt,this%nc,this%degree,this%knot1,order,this%wrap_parameters,values, &
            first_active,elem,Wc,knot2=this%knot2,knot3=this%knot3,side=side)
    end subroutine evaluate_basis_order
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute active mixed derivatives on an explicit volume parameter grid.
    !! Grid points use
    !! \(g=i_1+n_{g,1}[(i_2-1)+n_{g,2}(i_3-1)]\). For local offsets
    !! \(a_d=0,\ldots,p_d\), define
    !!
    !! \[
    !! \ell=1+a_1+(p_1+1)[a_2+(p_2+1)a_3],
    !! \quad
    !! A=f_1+a_1+n_{c,1}[(f_2+a_2-1)+n_{c,2}(f_3+a_3-1)],
    !! \]
    !!
    !! where \(f_d=\mathtt{first\_active(d,g)}\). Then `Dgc(ell,g)` is
    !! \(\partial^{|\mathbf k|}R_A/
    !! (\partial u^{k_1}\partial v^{k_2}\partial w^{k_3})\) for
    !! `order=[k1,k2,k3]`. Wrapping is applied directionally. If an unwrapped
    !! coordinate is outside its active interval, the derivative column is zero
    !! and the corresponding first-active index is 1.
    pure subroutine derivative_order_active_vector(this, Xt1, Xt2, Xt3, order, first_active, Dgc, side)
        class(nurbs_volume), target, intent(inout) :: this
        integer, intent(in), optional :: side(:) !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in), contiguous :: Xt1(:), Xt2(:), Xt3(:)
            !! Finite directional parameter vectors defining a Cartesian grid.
        integer, intent(in) :: order(3)
            !! Nonnegative mixed-derivative orders `[du,dv,dw]`.
        integer, allocatable, intent(out) :: first_active(:,:)
            !! First active control-point indices `[3,npoint]`.
        real(rk), allocatable, intent(out) :: Dgc(:,:)
            !! Local derivatives `[product(degree+1),product(grid sizes)]`.
        integer :: point1, point2, point3
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call check_derivative_state(this, "derivative_order_active_vector", valid, side)
        if (.not. valid) then
            allocate(first_active(0,0))
            allocate(Dgc(0,0))
            return
        end if
        if (any(order < 0) .or. .not. all(ieee_is_finite(Xt1)) .or. .not. all(ieee_is_finite(Xt2)) .or. &
            .not. all(ieee_is_finite(Xt3))) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Active volume derivative inputs must be finite with nonnegative orders.",&
                location   = "derivative_order_active_vector",&
                suggestion = "Use finite parameter grids and order=[du,dv,dw] with nonnegative entries.")
            allocate(first_active(0,0))
            allocate(Dgc(0,0))
            return
        end if

        allocate(Dgc(product(this%degree+1),size(Xt1)*size(Xt2)*size(Xt3)), source=0.0_rk)
        allocate(first_active(3,size(Xt1)*size(Xt2)*size(Xt3)), source=0)
        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
#if defined(__NVCOMPILER)
        do point3 = 1, size(Xt3)
            do point2 = 1, size(Xt2)
                do point1 = 1, size(Xt1)
#else
        do concurrent (point3 = 1:size(Xt3), point2 = 1:size(Xt2), point1 = 1:size(Xt1))
#endif
            call evaluate_basis_order(this,[Xt1(point1),Xt2(point2),Xt3(point3)],order, &
                Dgc(:,point1+size(Xt1)*((point2-1)+size(Xt2)*(point3-1))),weights,side, &
                first_active=first_active(:,point1+size(Xt1)*((point2-1)+size(Xt2)*(point3-1))))
#if defined(__NVCOMPILER)
                end do
            end do
#endif
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one active mixed volume derivative.
    !! Local and global control indices use the same direction-1-fastest
    !! formulas as the vector overload, with one parameter triple.
    pure subroutine derivative_order_active_scalar(this, Xt, order, first_active, Dgc, side)
        class(nurbs_volume), target, intent(inout) :: this
        integer, intent(in), optional :: side(:) !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in), contiguous :: Xt(:)
            !! Finite parameter triple `[u,v,w]`.
        integer, intent(in) :: order(3)
            !! Nonnegative mixed-derivative orders `[du,dv,dw]`.
        integer, intent(out) :: first_active(3)
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
        if (size(Xt) /= 3 .or. any(order < 0)) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Active volume derivative input requires three finite parameters and nonnegative orders.",&
                location   = "derivative_order_active_scalar",&
                suggestion = "Use Xt=[u,v,w] and order=[du,dv,dw] with nonnegative entries.")
            allocate(Dgc(0))
            return
        end if
        if (.not. all(ieee_is_finite(Xt))) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Active volume derivative parameters must be finite.",&
                location   = "derivative_order_active_scalar",&
                suggestion = "Use finite values for Xt=[u,v,w].")
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
    !> Compute one requested mixed derivative on a tensor grid of volume parameters.
    !! For `order=[k1,k2,k3]`, the dense matrix stores
    !!
    !! \[
    !! \mathtt{Dgc(g,A)}=
    !! \frac{\partial^{k_1+k_2+k_3}R_A}
    !!      {\partial u^{k_1}\partial v^{k_2}\partial w^{k_3}},
    !! \quad
    !! g=i_1+n_{g,1}[(i_2-1)+n_{g,2}(i_3-1)].
    !! \]
    !!
    !! Explicit directional vectors take precedence over positive resolutions;
    !! omitted directions reuse sampled parameter vectors. With no directional
    !! arguments, retained paired rows are evaluated directly. Queries never change
    !! sampled parameters, sample counts, or cached geometry. Empty or nonfinite
    !! parameter vectors and unavailable directions produce an empty result
    !! and a recoverable diagnostic.
    pure subroutine derivative_order_vector(this, order, Dgc, res1, res2, res3, Xt1, Xt2, Xt3, side)
        class(nurbs_volume), target, intent(inout) :: this
        integer, intent(in), optional :: side(:) !! Trace per direction: -1 left, 0 default, +1 right.
        integer, intent(in) :: order(3)
            !! Nonnegative mixed-derivative orders `[du,dv,dw]`.
        real(rk), allocatable, intent(out) :: Dgc(:,:)
            !! Dense derivative matrix `[product(ng),product(nc)]`.
        integer, intent(in), optional :: res1, res2, res3
            !! Optional positive generated counts by direction.
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:), Xt3(:)
            !! Optional explicit directional parameter vectors.
        integer :: point, ng(3)
        real(rk), allocatable :: Xt(:,:)
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
                category   = "forcad_nurbs_volume",&
                message    = "Volume derivative orders must be nonnegative.",&
                location   = "derivative_order_vector",&
                suggestion = "Use order=[du,dv,dw] with nonnegative derivative orders.")
            allocate(Dgc(0,0))
            return
        end if

        call prepare_query_grid(this, "derivative_order_vector", Xt, ng, valid, &
            res1, res2, res3, Xt1, Xt2, Xt3)
        if (.not. valid) then
            allocate(Dgc(0,0))
            return
        end if
        allocate(Dgc(product(ng),product(this%nc)), source=0.0_rk)
        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
#if defined(__NVCOMPILER)
        ! NVFORTRAN 26.3 cannot lower the pure local kernel for stdpar GPU execution.
        do point = 1, size(Xt,1)
#else
        do concurrent (point = 1:size(Xt,1))
#endif
            call evaluate_basis_order(this,Xt(point,:),order,Dgc(point,:),weights,side)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one requested mixed derivative at one volume parameter triple.
    !! Without `elem`, the result follows direction-1-fastest global control
    !! order. With `elem`, entries follow the supplied one-based flattened
    !! indices and an out-of-range index returns zero. A finite three-entry `Xt`
    !! is a caller precondition; this overload diagnoses its length and the
    !! derivative orders, but not nonfinite parameter values.
    pure subroutine derivative_order_scalar(this, Xt, order, Dgc, elem, side)
        class(nurbs_volume), target, intent(inout) :: this
        integer, intent(in), optional :: side(:) !! Trace per direction: -1 left, 0 default, +1 right.
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter triple `[u,v,w]`.
        integer, intent(in) :: order(3)
            !! Nonnegative mixed-derivative orders `[du,dv,dw]`.
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
        if (size(Xt) /= 3 .or. any(order < 0)) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Volume derivative input must contain three parameters and nonnegative orders.",&
                location   = "derivative_order_scalar",&
                suggestion = "Use Xt=[u,v,w] and order=[du,dv,dw] with nonnegative derivative orders.")
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
    !! omitted directions reuse their caches. Flattened rows follow
    !! direction-1-fastest `ndgrid` ordering; no arguments reuse retained paired
    !! rows directly. Parameters and counts are local
    !! to the query; sampled metadata and geometry are unchanged. Missing,
    !! empty, or nonfinite directional parameters and nonpositive resolutions
    !! produce an empty result and a recoverable diagnostic.
    pure subroutine derivative_vector(this, res1, res2, res3, Xt1, Xt2, Xt3, dTgc, Tgc, side)
        class(nurbs_volume), target, intent(inout) :: this
        integer, intent(in), optional :: side(:)
        integer, intent(in), optional :: res1, res2, res3
            !! Optional positive uniform-grid sizes in directions 1--3.
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:), Xt3(:)
            !! Optional directional parameter vectors.
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
            !! Parametric gradients `[ng_total,nc_total,3]`, ordered as \((u,v,w)\).
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
            !! Optional basis matrix `[ng_total,nc_total]`.
        call query_basis_grid(this, "derivative_vector", 1, res1, res2, res3, Xt1, Xt2, Xt3, side, Tgc=Tgc, dTgc=dTgc)
    end subroutine derivative_vector
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the parametric basis gradient at one parameter triple.
    !! Optional `elem` selects and orders a local control subset. `Xt` must have
    !! at least three finite entries; later entries are ignored. Invalid state
    !! or parameters return empty requested arrays and a recoverable diagnostic.
    pure subroutine derivative_scalar(this, Xt, dTgc, Tgc, elem, side)
        class(nurbs_volume), target, intent(inout) :: this
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter triple \([u,v,w]\).
        integer, intent(in), optional :: elem(:)
            !! Optional one-based flattened control indices.
        real(rk), allocatable, intent(out) :: dTgc(:,:)
            !! Selected parametric gradients `[nselected,3]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:)
            !! Optional selected basis values `[nselected]`.
        call query_basis_point(this, "derivative_scalar", 1, Xt, elem, side, Tgc=Tgc, dTgc=dTgc)
    end subroutine derivative_scalar
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate complete parametric basis Hessians on a tensor-product grid.
    !!
    !! Let `ncp=product(nc)`. The packed output has shape
    !! `[ng_total,3*ncp,3]`, with
    !! `d2Tgc(g,(a-1)*ncp+A,b)` equal to
    !! \(\partial^2R_A/(\partial\xi_a\partial\xi_b)\). Pure and mixed
    !! derivatives are all available, with symmetry represented explicitly.
    !! Parameter selection and validation are the same as for the vector
    !! [[nurbs_volume:derivative]] overload. Sample caches remain unchanged.
    pure subroutine derivative2_vector(this, res1, res2, res3, Xt1, Xt2, Xt3, d2Tgc, dTgc, Tgc, side)
        class(nurbs_volume), target, intent(inout) :: this
        integer, intent(in), optional :: side(:)
        integer, intent(in), optional :: res1, res2, res3
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:), Xt3(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:,:)
            !! Packed Hessians `[ng_total,3*nc_total,3]`.
        real(rk), allocatable, intent(out), optional :: dTgc(:,:,:)
            !! Optional first derivatives `[ng_total,nc_total,3]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
            !! Optional basis matrix `[ng_total,nc_total]`.
        call query_basis_grid(this, "derivative2_vector", 2, res1, res2, res3, Xt1, Xt2, Xt3, &
            side, Tgc=Tgc, dTgc=dTgc, d2Tgc=d2Tgc)
    end subroutine derivative2_vector
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate a complete parametric basis Hessian at one parameter triple.
    !!
    !! For `ncp=product(nc)`, `d2Tgc((a-1)*ncp+A,b)` stores
    !! \(\partial^2R_A/(\partial\xi_a\partial\xi_b)\). The output shape is
    !! `[3*ncp,3]`; all mixed entries are populated. With `elem`, replace
    !! `ncp` by `size(elem)` and `A` by the position in that list. Rational
    !! normalization always includes the full active support. `Xt` must have
    !! at least three finite entries; later entries are ignored. Invalid state
    !! or parameters return empty requested arrays and a recoverable diagnostic.
    pure subroutine derivative2_scalar(this, Xt, d2Tgc, dTgc, Tgc, elem, side)
        class(nurbs_volume), target, intent(inout) :: this
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter triple \([u,v,w]\).
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
            !! Packed Hessian `[3*nselected,3]`, in selected control order.
        real(rk), allocatable, intent(out), optional :: dTgc(:,:)
            !! Optional selected first derivatives `[nselected,3]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:)
            !! Optional selected basis values `[nselected]`.
        integer, intent(in), optional :: elem(:)
            !! Optional global control IDs; preserves order/repeats and zeros invalid or inactive IDs.
        call query_basis_point(this, "derivative2_scalar", 2, Xt, elem, side, Tgc=Tgc, dTgc=dTgc, d2Tgc=d2Tgc)
    end subroutine derivative2_scalar
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the dense volume basis on a tensor-product parameter grid.
    !! Explicit directional vectors take precedence over resolutions; every
    !! in-domain output row sums to one. No arguments reuse retained paired
    !! rows; omitted directions in a requested grid reuse sampled vectors. Parameters and counts are local to the query; sample
    !! caches remain unchanged. Invalid or unavailable directional parameters
    !! produce an empty result and a recoverable diagnostic.
    pure subroutine basis_vector(this, res1, res2, res3, Xt1, Xt2, Xt3, Tgc, side)
        class(nurbs_volume), target, intent(inout) :: this
        integer, intent(in), optional :: side(:)
        integer, intent(in), optional :: res1, res2, res3
            !! Optional positive uniform-grid sizes in directions 1--3.
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:), Xt3(:)
            !! Optional directional parameter vectors.
        real(rk), allocatable, intent(out) :: Tgc(:,:)
            !! Basis matrix `[ng_total,nc_total]`.
        call query_basis_grid(this, "basis_vector", 0, res1, res2, res3, Xt1, Xt2, Xt3, side, Tgc=Tgc)
    end subroutine basis_vector

    !> Evaluate the requested basis values, gradients, and Hessians on a grid.
    pure subroutine query_basis_grid(this, location, order, res1, res2, res3, Xt1, Xt2, Xt3, &
        side, Tgc, dTgc, d2Tgc)
        class(nurbs_volume), target, intent(inout) :: this
        character(*), intent(in) :: location
        integer, intent(in) :: order
        integer, intent(in), optional :: res1, res2, res3
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:), Xt3(:)
        integer, intent(in), optional :: side(:)
        real(rk), allocatable, intent(out), optional :: Tgc(:,:), dTgc(:,:,:), d2Tgc(:,:,:)
        real(rk), allocatable :: Xt(:,:)
        integer :: ng(3)
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call prepare_query_grid(this, location, Xt, ng, valid, &
            res1, res2, res3, Xt1, Xt2, Xt3,side=side)
        if (.not. valid) then
            if (present(d2Tgc)) allocate(d2Tgc(0,0,0))
            if (present(Tgc)) allocate(Tgc(0,0))
            if (present(dTgc)) allocate(dTgc(0,0,0))
            return
        end if

        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
        call compute_basis_data(Xt,this%knot1,this%knot2,this%knot3,this%degree,this%nc,order, &
            Tgc=Tgc,dTgc=dTgc,d2Tgc=d2Tgc,ng=ng,wrap_parameters=this%wrap_parameters,Wc=weights,side=side)
    end subroutine query_basis_grid
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the volume basis at one parameter triple.
    !! Optional `elem` selects and orders a local control subset. `Xt` must have
    !! at least three finite entries; later entries are ignored. Invalid state
    !! or parameters return empty requested arrays and a recoverable diagnostic.
    pure subroutine basis_scalar(this, Xt, Tgc, elem, side)
        class(nurbs_volume), target, intent(inout) :: this
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter triple \([u,v,w]\).
        integer, intent(in), optional :: elem(:)
            !! Optional one-based flattened control indices.
        real(rk), allocatable, intent(out) :: Tgc(:)
            !! Dense or selected basis vector.
        call query_basis_point(this, "basis_scalar", 0, Xt, elem, side, Tgc=Tgc)
    end subroutine basis_scalar

    !> Evaluate the requested basis values, gradients, and Hessians at one point.
    pure subroutine query_basis_point(this, location, order, Xt, elem, side, Tgc, dTgc, d2Tgc)
        class(nurbs_volume), target, intent(inout) :: this
        character(*), intent(in) :: location
        integer, intent(in) :: order
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: elem(:)
        integer, intent(in), optional :: side(:)
        real(rk), allocatable, intent(out), optional :: Tgc(:), dTgc(:,:), d2Tgc(:,:)
        real(rk), pointer, contiguous :: weights(:)
        logical :: valid

        call check_query_point(this,Xt,location,valid,side=side)
        if (.not. valid) then
            if (present(d2Tgc)) allocate(d2Tgc(0,0))
            if (present(Tgc)) allocate(Tgc(0))
            if (present(dTgc)) allocate(dTgc(0,0))
            return
        end if

        nullify(weights)
        if (this%is_rational()) weights=>this%Wc
        call compute_basis_data(Xt,this%knot1,this%knot2,this%knot3,this%degree,this%nc,order, &
            Tgc=Tgc,dTgc=dTgc,d2Tgc=d2Tgc,elem=elem,wrap_parameters=this%wrap_parameters,Wc=weights,side=side)
    end subroutine query_basis_point
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Insert knots in one tensor direction without changing the volume.
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
        class(nurbs_volume),            intent(inout) :: this
        integer,                        intent(in)    :: dir
            !! Refined parameter direction, 1, 2, or 3.
        real(rk), contiguous,           intent(in)    :: Xth(:)
            !! Finite values in the tolerance-expanded closed active interval; values are not snapped.
        integer,  contiguous,           intent(in)    :: r(:)
            !! Nonnegative counts; final multiplicity is at most `degree(dir)+1`, or `degree(dir)` if periodic.
        real(rk), allocatable, optional,intent(out)   :: B(:,:)
            !! Optional coordinate-expanded refinement operator.
        real(rk), allocatable, optional,intent(out)   :: Bs(:,:)
            !! Optional scalar refinement operator.
        real(rk), allocatable :: knots(:), Xc(:,:), Wc(:)
        integer :: retained(2)

        if (.not.this%err%ok) return
        if (.not.this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_volume",location="insert_knots", &
                message="Refinement requires valid geometry and nonnegative counts.", &
                suggestion="Set complete spline geometry and check the refinement request.")
            return
        end if
        call refine_control_net(this%get_knot(dir),this%get_degree(dir),this%Xc,this%Wc,this%nc,dir, &
            knots,Xc,Wc,this%err,B,Bs,Xth=Xth,r=r)
        if (.not.this%err%ok .or. .not.allocated(Xc)) return
        retained=pack([1,2,3],[1,2,3]/=dir)
        call set_direction(this,dir,knots,this%get_knot(retained(1)),this%get_knot(retained(2)), &
            Xc,Wc,[this%degree],[this%wrap_parameters])
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Raise the degree in one tensor direction while preserving active-domain geometry.
    !! Optional `Bs` and `B` are the full scalar and coordinate-expanded maps
    !! from old to elevated Euclidean flattened control values after fixing the
    !! output weights. For rational data they depend on the old and new weights
    !! and are not homogeneous elevation operators.
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
        class(nurbs_volume), target,             intent(inout) :: this
        integer,                         intent(in)    :: dir
            !! Elevated parameter direction, 1, 2, or 3.
        integer,                         intent(in)    :: t
            !! Nonnegative degree increment; zero is the identity operation.
        real(rk), allocatable, optional, intent(out)   :: B(:,:)
            !! Optional coordinate-expanded elevation operator.
        real(rk), allocatable, optional, intent(out)   :: Bs(:,:)
            !! Optional scalar elevation operator.
        real(rk), allocatable :: knots(:), Xc(:,:), Wc(:)
        integer :: degree(3)
        integer :: retained(2)

        if (.not.this%err%ok) return
        if (.not.this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_volume",location="elevate_degree", &
                message="Refinement requires valid geometry and nonnegative counts.", &
                suggestion="Set complete spline geometry and check the refinement request.")
            return
        end if
        call refine_control_net(this%get_knot(dir),this%get_degree(dir),this%Xc,this%Wc,this%nc,dir, &
            knots,Xc,Wc,this%err,B,Bs,t=t)
        if (.not.this%err%ok .or. .not.allocated(Xc)) return
        degree=this%degree
        degree(dir)=degree(dir)+t
        retained=pack([1,2,3],[1,2,3]/=dir)
        call set_direction(this,dir,knots,this%get_knot(retained(1)),this%get_knot(retained(2)), &
            Xc,Wc,degree,[this%wrap_parameters])
    end subroutine
    !===============================================================================


    !> Reduce degree under one Euclidean error budget for the complete call.
    !! Each candidate is compared with the call-entry geometry. Return the longest
    !! accepted prefix in reduced; over-budget proposals are normal rejections.
    !! Commit once, clearing derived caches only when at least one step succeeds.
    pure subroutine reduce_degree(this, dir, t, tol, reduced, error_bound)
        class(nurbs_volume), intent(inout) :: this
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
        integer :: nc(3),degree(3)
        real(rk), allocatable :: Xc(:,:),Wc(:)
        integer :: retained(2),d

        if (present(reduced)) reduced=0
        if (present(error_bound)) error_bound=0.0_rk
        if (.not.this%err%ok) return
        if (.not.this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_volume",location="reduce_degree", &
                message="Degree reduction requires valid geometry.",suggestion="Initialize valid knots, controls and weights.")
            return
        end if
        if (dir<1 .or. dir>3) then
            call this%err%set(code=100,severity=1,category="forcad_nurbs_volume",location="reduce_degree", &
                message="Invalid reduction direction.",suggestion="Use a parameter direction from 1 to 3.")
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
        retained=pack([1,2,3],[1,2,3]/=dir)
        call set_direction(this,dir,knots,this%get_knot(retained(1)),this%get_knot(retained(2)),Xc,Wc,degree, &
            [(this%get_parameter_wrapping(d),d=1,3)])
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
        class(nurbs_volume), intent(in) :: this
        logical :: r

        r = .false.
        if (.not. this%err%ok) return
        r = this%rational .and. allocated(this%Wc)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report modulo parameter mapping in one volume direction.
    !! Invalid directions or an active diagnostic return `.false.`.
    pure elemental function get_parameter_wrapping(this, dir) result(wrap_parameters)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        logical :: wrap_parameters

        wrap_parameters = .false.
        if (.not. this%err%ok .or. dir < 1 .or. dir > 3) return
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
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        character(len=8) :: topology

        topology = "invalid"
        if (.not. this%is_valid() .or. dir < 1 .or. dir > 3) return

        topology = "bounded"
        if (allocated(this%Wc)) then
            select case (dir)
            case (1)
                if (periodic_topology(this%knot1, this%degree(1), this%Xc, this%Wc, this%nc, 1)) then
                    topology = "periodic"
                end if
            case (2)
                if (periodic_topology(this%knot2, this%degree(2), this%Xc, this%Wc, this%nc, 2)) then
                    topology = "periodic"
                end if
            case (3)
                if (periodic_topology(this%knot3, this%degree(3), this%Xc, this%Wc, this%nc, 3)) then
                    topology = "periodic"
                end if
            case default
                return
            end select
        else
            select case (dir)
            case (1)
                if (periodic_topology(this%knot1, this%degree(1), this%Xc, nc=this%nc, dir=1)) then
                    topology = "periodic"
                end if
            case (2)
                if (periodic_topology(this%knot2, this%degree(2), this%Xc, nc=this%nc, dir=2)) then
                    topology = "periodic"
                end if
            case (3)
                if (periodic_topology(this%knot3, this%degree(3), this%Xc, nc=this%nc, dir=3)) then
                    topology = "periodic"
                end if
            case default
                return
            end select
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Replace control-lattice visualization connectivity with a copy of `elemConn`.
    pure subroutine set_elem_Xc_vis(this, elemConn)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (.not. this%err%ok) return

        this%elemConn_Xc_vis = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Replace sampled-volume visualization connectivity with a copy of `elemConn`.
    pure subroutine set_elem_Xg_vis(this, elemConn)
        class(nurbs_volume), intent(inout) :: this
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
    !! vectors or control-lattice dimensions.
    pure subroutine set_elem(this, elemConn)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), contiguous :: elemConn(:,:)

        if (.not. this%err%ok) return

        this%elemConn = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocated copy of control-lattice visualization connectivity.
    pure function get_elem_Xc_vis(this) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        call get_array(this%elemConn_Xc_vis,elemConn,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocated copy of sampled-volume visualization connectivity.
    pure function get_elem_Xg_vis(this) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        call get_array(this%elemConn_Xg_vis,elemConn,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocated copy of cached IGA element-to-control-point connectivity.
    !! An absent cache or active diagnostic returns shape `[0,0]`.
    !! [[nurbs_volume:cmp_elem]] computes connectivity without storing it;
    !! use [[nurbs_volume:set_elem]] to populate this compatibility cache.
    pure function get_elem(this) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        call get_array(this%elemConn,elemConn,this%err%ok)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct an axis-aligned hexahedral tensor-product volume.
    !!
    !! The control net is uniformly distributed on
    !! `[0,L(1)] x [0,L(2)] x [0,L(3)]`; direction 1 varies fastest. Generated
    !! knot vectors follow the standard `set(nc=...,Xc=...)` policy.
    !!
    pure subroutine set_hexahedron(this, L, nc, Wc)
        class(nurbs_volume), intent(inout) :: this
            !! Volume replaced by the generated hexahedron.
        real(rk), intent(in), contiguous :: L(:)
            !! Finite side lengths `[Lx,Ly,Lz]`.
        integer, intent(in), contiguous :: nc(:)
            !! Directional control counts `[nu,nv,nw]`, all greater than one.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional finite positive weights `[product(nc)]`.

        if (.not. this%err%ok) return

        if (present(Wc)) then
            call this%set(nc, hexahedron_Xc(L, nc), Wc)
        else
            call this%set(nc, hexahedron_Xc(L, nc))
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Parameterize rows by three coordinates and cache the mapped volume.
    !!
    !! Columns 1--3 of `X` are independently mapped from their data ranges
    !! to the three active knot intervals. For direction \(d\),
    !!
    !! \[
    !! \xi_{id}=\xi_{a,d}
    !! +\frac{X_{id}-X_{\min,d}}{X_{\max,d}-X_{\min,d}}
    !!  (\xi_{b,d}-\xi_{a,d}).
    !! \]
    !!
    !! Remaining columns do not affect the parameterization. The routine
    !! atomically replaces the flattened parameter-point pairs
    !! \(((\xi_{i1},\xi_{i2},\xi_{i3}),\mathbf V(\boldsymbol\xi_i))\) and
    !! stores their count as `ng=[size(X,1),1,1]`. Because arbitrary input rows
    !! need not form a Cartesian product, directional parameter vectors are
    !! cleared. When `elemConn` is absent, sampled visualization connectivity is
    !! reset to an empty hexahedral topology.
    pure subroutine put_to_nurbs(this, X, elemConn)
        class(nurbs_volume), intent(inout) :: this
            !! Initialized volume whose sampled geometry cache is replaced.
        real(rk), intent(in), contiguous :: X(:,:)
            !! Ordering data `[npoint,ncolumn]`; columns 1--3 need nonzero ranges.
        integer, intent(in), contiguous, optional :: elemConn(:,:)
            !! Optional one-based visualization connectivity for the cached points.
        real(rk), allocatable :: Xt(:,:), Xg_new(:,:), axis(:),knot(:)
        real(rk) :: lo,hi
        integer :: d

        if (.not. this%err%ok) return

        if (.not.allocated(this%Xc)) then
            call this%err%set(code=102,severity=1,category="forcad_nurbs_volume",location="put_to_nurbs", &
                message="Control points are not set.",suggestion="Call set(...) before mapping points.")
            return
        end if
        if (.not.allocated(this%knot1) .or. .not.allocated(this%knot2) .or. .not.allocated(this%knot3)) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_volume",location="put_to_nurbs", &
                message="Knot vectors are not set.",suggestion="Call set(...) before mapping points.")
            return
        end if
        if (size(X,1)<1 .or. size(X,2)<3) then
            call this%err%set(code=101,severity=1,category="forcad_nurbs_volume",location="put_to_nurbs", &
                message="Input points need at least one row and three coordinates.", &
                suggestion="Provide X with size(X,1)>=1 and size(X,2)>=3.")
            return
        end if

        allocate(Xt(size(X,1),3))
        do d=1,3
            lo=minval(X(:,d))
            hi=maxval(X(:,d))
            if (hi<=lo) then
                call this%err%set(code=100,severity=1,category="forcad_nurbs_volume",location="put_to_nurbs", &
                    message="Input points have zero extent in a parameter direction.", &
                    suggestion="Provide distinct coordinates in every parameter direction.")
                return
            end if
            knot=this%get_knot(d)
            call map_parameter_values(X(:,d),lo,hi,knot_start(knot,this%nc(d),this%degree(d)), &
                knot_end(knot,this%nc(d),this%degree(d)),axis,this%err)
            if (.not.this%err%ok) return
            Xt(:,d)=axis
        end do
        Xg_new=this%cmp_Xg(Xt)
        if (size(Xg_new,1)/=size(X,1) .or. .not.all(ieee_is_finite(Xg_new))) then
            call this%err%set(code=108,severity=1,category="forcad_nurbs_volume",location="put_to_nurbs", &
                message="Mapped geometry is outside the floating-point range.", &
                suggestion="Rescale the control coordinates or input values.")
            return
        end if
        call move_alloc(Xt, this%Xt)
        call move_alloc(Xg_new, this%Xg)
        this%ng = [size(X,1),1,1]
        if (allocated(this%Xt1)) deallocate(this%Xt1)
        if (allocated(this%Xt2)) deallocate(this%Xt2)
        if (allocated(this%Xt3)) deallocate(this%Xt3)

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
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: dir
            !! Modified parameter direction, 1 through 3.
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
        integer :: nc(3), weight_exponent
        integer :: retained(2)

        if (present(removed)) allocate(removed(size(r)),source=0)
        if (present(error_bound)) error_bound=0.0_rk
        if (.not.this%err%ok) return
        if (.not.this%is_valid()) then
            call this%err%set(code=103,severity=1,category="forcad_nurbs_volume",location="remove_knots", &
                message="Knot removal requires valid geometry.",suggestion="Initialize valid controls, knots and weights.")
            return
        end if
        if (dir<1 .or. dir>3) then
            call this%err%set(code=100,severity=1,category="forcad_nurbs_volume",location="remove_knots", &
                message="Invalid knot removal direction.",suggestion="Use dir from 1 through 3.")
            return
        end if
        nc=this%nc
        knots=this%get_knot(dir)
        call coarsen_control_lines(knots,this%degree(dir),this%Xc,this%Wc,nc,dir,Xth,r, &
            this%err,controls,weight_exponent,counts,tol,error_bound)
        if (.not.this%err%ok .or. .not.any(counts>0)) return
        nc(dir)=size(controls,1)
        call unpack_refinement_controls(controls,nc,dir,size(this%Xc,2),weight_exponent,Xc,Wc)
        retained=pack([1,2,3],[1,2,3]/=dir)
        call set_direction(this,dir,knots,this%get_knot(retained(1)),this%get_knot(retained(2)), &
            Xc,Wc,[this%degree],[this%wrap_parameters])
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
    !> Build IGA connectivity over all nonzero tensor-product knot spans.
    !!
    !! Each row contains \((p+1)(q+1)(r+1)\) flattened control indices in
    !! direction-1-fastest order. Knot multiplicities determine overlap
    !! between adjacent rows. For verified periodic directions, apply
    !! [[nurbs_volume:cmp_dof_map]] before assembling independent cyclic
    !! degrees of freedom.
    pure function cmp_elem(this) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not.this%err%ok .or. &
            .not.allocated(this%knot1) .or. &
            .not.allocated(this%knot2) .or. &
            .not.allocated(this%knot3)) then
            allocate(elemConn(0,0))
        else
            elemConn=spline_connectivity(this%nc,this%degree,this%knot1,knot2=this%knot2,knot3=this%knot3)
        end if
    end function cmp_elem
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Map stored controls to independent cyclic IGA degrees of freedom.
    !!
    !! The returned one-based array follows direction-1-fastest storage.
    !! A bounded direction keeps every stored layer; a verified periodic
    !! degree-\(p_d\) direction identifies its final \(p_d\) layers with its
    !! first layers. Apply this map to [[nurbs_volume:cmp_elem]] connectivity.
    pure function cmp_dof_map(this) result(map)
        class(nurbs_volume), intent(in) :: this
        integer, allocatable :: map(:)
        integer :: i, independent_count(3)

        if (.not. this%is_valid()) return

        independent_count = this%nc
        if (this%get_parameter_topology(1) == "periodic") then
            independent_count(1) = this%nc(1) - this%degree(1)
        end if
        if (this%get_parameter_topology(2) == "periodic") then
            independent_count(2) = this%nc(2) - this%degree(2)
        end if
        if (this%get_parameter_topology(3) == "periodic") then
            independent_count(3) = this%nc(3) - this%degree(3)
        end if

        allocate(map(product(this%nc)))
        do concurrent (i = 0:product(this%nc)-1)
            map(i+1) = modulo(modulo(i,this%nc(1)),independent_count(1)) + 1 + &
                modulo(modulo(i/this%nc(1),this%nc(2)),independent_count(2))*independent_count(1) + &
                modulo(i/(this%nc(1)*this%nc(2)),independent_count(3))*&
                independent_count(1)*independent_count(2)
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
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: alpha, beta, theta

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) return

        call rotate_points(this%Xc, this%nc(1)*this%nc(2)*this%nc(3), alpha, beta, theta)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Rotate all cached physical samples using x-y-z Euler angles in degrees.
    !!
    !! The convention and coordinate-dimension behavior match
    !! [[nurbs_volume:rotate_Xc]]. Control points are not changed.
    pure subroutine rotate_Xg(this, alpha, beta, theta)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: alpha, beta, theta

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xg)) return

        call rotate_points(this%Xg, this%ng(1)*this%ng(2)*this%ng(3), alpha, beta, theta)
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
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: vec(:)
        integer :: i, d

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xc)) return

        do concurrent (i = 1:this%nc(1)*this%nc(2)*this%nc(3), d = 1:size(this%Xc,2))
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
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: vec(:)
        integer :: i, d

        if (.not. this%err%ok) return
        if (.not. allocated(this%Xg)) return

        do concurrent (i = 1:this%ng(1)*this%ng(2)*this%ng(3), d = 1:size(this%Xg,2))
            this%Xg(i,d) = this%Xg(i,d) + vec(d)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct an exact rational annular extrusion.
    !!
    !! The first direction contains a three-arc exact circle, the second joins
    !! two circular boundaries linearly, and the third extrudes from
    !! `center(3)` to `center(3)+length`. The control net has shape `7 x 2 x 2`.
    !! Boundary radii are `abs(radius1)` and `abs(radius2)`. For a regular
    !! annular volume, use nonzero radius scales of the same sign and unequal
    !! magnitudes, and a nonzero length. Opposite radius signs make corresponding
    !! radial lines cross the axis. Zero length collapses the third direction.
    !!
    pure subroutine set_ring(this, center, radius1, radius2, length)
        class(nurbs_volume), intent(inout) :: this
            !! Volume replaced by the annular extrusion.
        real(rk), intent(in), contiguous :: center(:)
            !! Base center; at least three finite components are required.
        real(rk), intent(in) :: radius1
            !! Finite signed scale of the first cylindrical boundary.
        real(rk), intent(in) :: radius2
            !! Finite signed scale of the second cylindrical boundary.
        real(rk), intent(in) :: length
            !! Finite signed extrusion length; use a nonzero value for a nondegenerate volume.
        real(rk), allocatable :: Xc(:,:), Wc(:), knot1(:), knot2(:), knot3(:)
        integer :: i, d

        if (.not. this%err%ok) return

        ! Define control points for ring
        allocate(Xc(28, 3))
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

        Xc(15,:)= [ 1.0_rk,  0.0_rk,              length]
        Xc(16,:)= [ 1.0_rk,  sqrt(3.0_rk),        length]
        Xc(17,:)= [-0.5_rk,  sqrt(3.0_rk)/2.0_rk, length]
        Xc(18,:)= [-2.0_rk,  0.0_rk,              length]
        Xc(19,:)= [-0.5_rk, -sqrt(3.0_rk)/2.0_rk, length]
        Xc(20,:)= [ 1.0_rk, -sqrt(3.0_rk),        length]
        Xc(21,:)= [ 1.0_rk,  0.0_rk,              length]

        Xc(15:21,1:2) = Xc(15:21,1:2) * radius1

        Xc(22,:)= [ 1.0_rk,  0.0_rk,              length]
        Xc(23,:)= [ 1.0_rk,  sqrt(3.0_rk),        length]
        Xc(24,:)= [-0.5_rk,  sqrt(3.0_rk)/2.0_rk, length]
        Xc(25,:)= [-2.0_rk,  0.0_rk,              length]
        Xc(26,:)= [-0.5_rk, -sqrt(3.0_rk)/2.0_rk, length]
        Xc(27,:)= [ 1.0_rk, -sqrt(3.0_rk),        length]
        Xc(28,:)= [ 1.0_rk,  0.0_rk,              length]

        Xc(22:28,1:2) = Xc(22:28,1:2) * radius2

        ! Translate the control points
        do concurrent (i = 1:size(Xc, 1), d = 1:size(Xc, 2))
            Xc(i,d) = center(d) + Xc(i,d)
        end do

        ! Define weights for the control points
        Wc = [1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk,&
            1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk,&
            1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk,&
            1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk]

        ! Define knot vector
        knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/3.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot2 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot3 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

        ! Set knot vector, control points, and weights
        call this%set(knot1, knot2, knot3, Xc, Wc)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct an exact rational C-shaped annular extrusion.
    !!
    !! Two quadratic rational arcs span 240 degrees in direction 1, direction 2
    !! joins the circular boundaries, and direction 3 is a linear extrusion.
    !! The control net has shape `5 x 2 x 2`.
    !! For a regular C-shaped annular volume, use nonzero radius scales of the
    !! same sign and unequal magnitudes, together with a nonzero length. A common
    !! negative radius sign rotates the template by \(\pi\); opposite signs make
    !! radial lines cross the axis.
    !!
    pure subroutine set_C(this, center, radius1, radius2, length)
        class(nurbs_volume), intent(inout) :: this
            !! Volume replaced by the C-shaped extrusion.
        real(rk), intent(in), contiguous :: center(:)
            !! Base center; at least three finite components are required.
        real(rk), intent(in) :: radius1
            !! Finite signed scale of the first cylindrical boundary.
        real(rk), intent(in) :: radius2
            !! Finite signed scale of the second cylindrical boundary.
        real(rk), intent(in) :: length
            !! Finite signed extrusion length; use a nonzero value for a nondegenerate volume.
        real(rk), allocatable :: Xc(:,:), Wc(:), knot1(:), knot2(:), knot3(:)
        integer :: i, d

        if (.not. this%err%ok) return

        ! Define control points for C-shape
        allocate(Xc(20, 3))
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

        Xc(11,:)= [ 1.0_rk,  0.0_rk,              length]
        Xc(12,:)= [ 1.0_rk,  sqrt(3.0_rk),        length]
        Xc(13,:)= [-0.5_rk,  sqrt(3.0_rk)/2.0_rk, length]
        Xc(14,:)= [-2.0_rk,  0.0_rk,              length]
        Xc(15,:)= [-0.5_rk, -sqrt(3.0_rk)/2.0_rk, length]

        Xc(11:15,1:2) = Xc(11:15,1:2) * radius1

        Xc(16,:)= [ 1.0_rk,  0.0_rk,              length]
        Xc(17,:)= [ 1.0_rk,  sqrt(3.0_rk),        length]
        Xc(18,:)= [-0.5_rk,  sqrt(3.0_rk)/2.0_rk, length]
        Xc(19,:)= [-2.0_rk,  0.0_rk,              length]
        Xc(20,:)= [-0.5_rk, -sqrt(3.0_rk)/2.0_rk, length]

        Xc(16:20,1:2) = Xc(16:20,1:2) * radius2

        ! Translate the control points
        do concurrent (i = 1:size(Xc, 1), d = 1:size(Xc, 2))
            Xc(i,d) = center(d) + Xc(i,d)
        end do

        ! Define weights for the control points
        Wc = [1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk,&
            1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk,&
            1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk,&
            1.0_rk, 0.5_rk, 1.0_rk, 0.5_rk, 1.0_rk]

        ! Define knot vector
        knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/2.0_rk, 1.0_rk/2.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot2 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot3 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

        ! Set knot vector, control points, and weights
        call this%set(knot1, knot2, knot3, Xc, Wc)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct an exact rational half-annular extrusion.
    !!
    !! Two quadratic quarter-arcs form each semicircular boundary, direction 2
    !! joins the boundaries, and direction 3 is a linear extrusion. The control
    !! net has shape `5 x 2 x 2`.
    !! For a regular half-annular volume, use nonzero radius scales of the same
    !! sign and unequal absolute values, together with a nonzero length. A
    !! negative common radius sign rotates the template by \(\pi\); opposite
    !! signs create crossing radial lines.
    !!
    pure subroutine set_half_ring(this, center, radius1, radius2, length)
        class(nurbs_volume), intent(inout) :: this
            !! Volume replaced by the half-annular extrusion.
        real(rk), intent(in), contiguous :: center(:)
            !! Base center; at least three finite components are required.
        real(rk), intent(in) :: radius1
            !! Finite legacy first-boundary scale; represented radius `abs(radius1)/2`.
        real(rk), intent(in) :: radius2
            !! Finite legacy second-boundary scale; represented radius `abs(radius2)/2`.
        real(rk), intent(in) :: length
            !! Finite signed extrusion length; use a nonzero value for a nondegenerate volume.
        real(rk), allocatable :: Xc(:,:), Wc(:), knot1(:), knot2(:), knot3(:)
        integer :: i, d

        if (.not. this%err%ok) return

        ! Define control points for half ring
        allocate(Xc(20, 3))
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

        Xc(11,:) = [ 0.5_rk, 0.0_rk, length]
        Xc(12,:) = [ 0.5_rk, 0.5_rk, length]
        Xc(13,:) = [ 0.0_rk, 0.5_rk, length]
        Xc(14,:) = [-0.5_rk, 0.5_rk, length]
        Xc(15,:) = [-0.5_rk, 0.0_rk, length]

        Xc(11:15,1:2) = Xc(11:15,1:2) * radius1

        Xc(16,:) = [ 0.5_rk, 0.0_rk, length]
        Xc(17,:) = [ 0.5_rk, 0.5_rk, length]
        Xc(18,:) = [ 0.0_rk, 0.5_rk, length]
        Xc(19,:) = [-0.5_rk, 0.5_rk, length]
        Xc(20,:) = [-0.5_rk, 0.0_rk, length]

        Xc(16:20,1:2) = Xc(16:20,1:2) * radius2

        ! Translate the control points
        do concurrent (i = 1:size(Xc, 1), d = 1:size(Xc, 2))
            Xc(i,d) = center(d) + Xc(i,d)
        end do

        ! Define weights for the control points
        Wc = [1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk,&
            1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk,&
            1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk,&
            1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk, 1.0_rk/sqrt(2.0_rk), 1.0_rk]

        ! Define knot vector
        knot1 = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/2.0_rk, &
            1.0_rk/2.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot2 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        knot3 = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]

        ! Set knot vector, control points, and weights
        call this%set(knot1, knot2, knot3, Xc, Wc)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Find the closest cached sample using a scale-safe Euclidean distance.
    !! Only min(query size,physical dimension) coordinates participate; extra
    !! returned coordinates are zero-padded. Ties retain the first sample.
    !! Missing caches, active diagnostics, empty queries or nonfinite participating
    !! coordinates return zero outputs and id=0. Sampling limits search accuracy.
    pure subroutine nearest_point(this, point_Xg, nearest_Xg, nearest_Xt, id)
        class(nurbs_volume), intent(in) :: this
        real(rk), intent(in), contiguous :: point_Xg(:)
        real(rk), intent(out), optional :: nearest_Xg(size(point_Xg)),nearest_Xt(3)
        integer, intent(out), optional :: id
        integer :: selected,npoint,ncopy

        if (present(nearest_Xg)) nearest_Xg=0.0_rk
        if (present(nearest_Xt)) nearest_Xt=0.0_rk
        if (present(id)) id=0
        if (.not.this%err%ok .or. .not.allocated(this%Xg) .or. .not.allocated(this%Xt)) return
        npoint=min(size(this%Xg,1),size(this%Xt,1))
        selected=nearest_sample_index(this%Xg(:npoint,:),point_Xg)
        if (selected==0) return
        if (present(id)) id=selected
        if (present(nearest_Xg)) then
            ncopy=min(size(point_Xg),size(this%Xg,2))
            nearest_Xg(:ncopy)=this%Xg(selected,:ncopy)
        end if
        if (present(nearest_Xt)) nearest_Xt=this%Xt(selected,1:3)
    end subroutine nearest_point
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute a bounded continuous nearest-point candidate without sampling caches.
    !! Span-aware seeds initialize a safeguarded Newton/Gauss-Newton search;
    !! a finite in-domain initial_Xt bypasses that search. Parameters are scaled
    !! internally to [0,1]. A successful result is a local stationary candidate,
    !! not a certificate of the global nearest point or a zero-distance inverse.
    !! At discontinuities, a left trace uses the representable parameter below
    !! the knot so nearest_Xg agrees with cmp_Xg(nearest_Xt).
    !! Invalid inputs return zero outputs; numerical or convergence failures
    !! retain a valid candidate when available. maxit=0 retains the seed and
    !! reports code109. Caller-owned err bypasses and preserves this%err.
    pure subroutine nearest_point2(this, point_Xg, tol, maxit, nearest_Xt, nearest_Xg, initial_Xt, err)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: point_Xg(:)
            !! Finite query with at least every stored physical coordinate; extras are ignored.
        real(rk), intent(in) :: tol
            !! Finite nonnegative parameter tolerance.
        integer, intent(in) :: maxit
            !! Nonnegative iteration limit.
        real(rk), intent(out) :: nearest_Xt(3)
        real(rk), intent(out), optional :: nearest_Xg(:)
            !! Must match the stored physical dimension when present.
        real(rk), intent(in), optional :: initial_Xt(:)
            !! Optional finite three-entry initial parameter inside the literal active domain.
        type(debug), intent(out), optional :: err
        type(debug) :: status

        nearest_Xt=0.0_rk
        if (present(nearest_Xg)) nearest_Xg=0.0_rk
        if (present(err)) then
            call err%reset()
        else if (.not.this%err%ok) then
            return
        end if
        call status%reset()
        if (.not.allocated(this%Xc) .or. .not.allocated(this%knot1) .or. &
            .not.allocated(this%knot2) .or. .not.allocated(this%knot3)) then
            call status%set(code=103,severity=1,category="forcad_nurbs_volume",location="nearest_point2", &
                message="Nearest-point search requires initialized knots and control points.", &
                suggestion="Set complete volume geometry before searching.")
        else
            call nearest_spline_point(this%Xc,this%Wc,this%nc,this%degree,this%knot1,point_Xg,tol,maxit, &
                nearest_Xt,nearest_Xg,status,knot2=this%knot2,knot3=this%knot3,initial_Xt=initial_Xt, &
                wrap=this%wrap_parameters)
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
    !> Extract global control indices for one face of one IGA element.
    !!
    !! `elem` is a one-based row of current IGA connectivity, computed from
    !! the current knots without requiring sampling or a stored connectivity
    !! cache. Face numbers
    !! are 1 `w_min`, 2 `w_max`, 3 `v_min`, 4 `v_max`,
    !! 5 `u_min`, and 6 `u_max`. The result contains the two
    !! tangential tensor directions with direction 1 varying fastest where it
    !! is tangential. Invalid input or unavailable connectivity returns an
    !! allocated empty vector.
    !> Active controls for one current element, without a global connectivity table.
    pure function element_nodes(this, ie) result(nodes)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: ie
        integer, allocatable :: nodes(:)
        real(rk) :: bounds(2,3)
        type(debug) :: status

        if (.not.allocated(this%knot1) .or. .not.allocated(this%knot2) .or. .not.allocated(this%knot3)) then
            allocate(nodes(0))
            return
        end if
        call spline_element(this%nc,this%degree,this%knot1,ie,nodes,bounds,status, &
            knot2=this%knot2,knot3=this%knot3)
    end function element_nodes

    !> Eight tensor-ordered vertices for one structured trilinear cell.
    pure function cell_nodes(counts, ie) result(nodes)
        integer, intent(in) :: counts(3),ie
        integer, allocatable :: nodes(:)
        integer :: ne(3),e(3),i,j,k
        ne=counts-1
        if (any(ne<1) .or. ie<1 .or. ie>product(ne)) then
            allocate(nodes(0))
            return
        end if
        e=[mod(ie-1,ne(1))+1,mod((ie-1)/ne(1),ne(2))+1,(ie-1)/(ne(1)*ne(2))+1]
        allocate(nodes(8))
        do concurrent(k=0:1,j=0:1,i=0:1)
            nodes(1+i+2*(j+2*k))=e(1)+i+counts(1)*(e(2)+j-1+counts(2)*(e(3)+k-1))
        end do
    end function cell_nodes

    !> Restrict tensor-ordered controls or vertices to one parametric face.
    pure function face_nodes(nodes, n, face) result(selected)
        integer, intent(in) :: nodes(:),n(3),face
        integer, allocatable :: selected(:)
        integer :: k
        if (size(nodes)/=product(n) .or. any(n<1)) then
            allocate(selected(0))
            return
        end if
        select case(face)
        case(1)
            selected=nodes(1:n(1)*n(2))
        case(2)
            selected=nodes(size(nodes)-n(1)*n(2)+1:)
        case(3,4)
            allocate(selected(n(1)*n(3)))
            do k=0,n(3)-1
                if(face==3) then
                    selected(k*n(1)+1:(k+1)*n(1))=nodes(k*n(1)*n(2)+1:k*n(1)*n(2)+n(1))
                else
                    selected(k*n(1)+1:(k+1)*n(1))=nodes((k+1)*n(1)*n(2)-n(1)+1:(k+1)*n(1)*n(2))
                end if
            end do
        case(5)
            selected=nodes(1::n(1))
        case(6)
            selected=nodes(n(1)::n(1))
        case default
            allocate(selected(0))
        end select
    end function face_nodes

    pure function cmp_elemFace(this, elem, face) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: elem, face
        integer, allocatable :: elemConn(:)
        if (.not. this%err%ok) then
            allocate(elemConn(0))
            return
        end if
        elemConn=face_nodes(element_nodes(this,elem),this%degree+1,face)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Extract one quadrilateral face from a control-lattice VTK cell.
    !!
    !! `elem` selects a cell in the current default trilinear visualization
    !! connectivity, generated independently of stored connectivity overrides.
    !! Face numbering is
    !! `[w_min,w_max,v_min,v_max,u_min,u_max]`. The result has four
    !! global point indices, or is empty for invalid or unavailable state.
    pure function cmp_elemFace_Xc_vis(this, elem, face) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: elem, face
        integer, allocatable :: elemConn(:)
        if (.not. this%err%ok) then
            allocate(elemConn(0))
            return
        end if
        elemConn=face_nodes(cell_nodes(this%nc,elem),[2,2,2],face)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Extract one quadrilateral face from a sampled-volume VTK cell.
    !!
    !! `elem` selects a cell in the current default trilinear visualization
    !! connectivity, generated independently of stored connectivity overrides.
    !! Face numbering is
    !! `[w_min,w_max,v_min,v_max,u_min,u_max]`. The result has four
    !! global point indices, or is empty for invalid or unavailable state.
    pure function cmp_elemFace_Xg_vis(this, elem, face) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: elem, face
        integer, allocatable :: elemConn(:)
        if (.not. this%err%ok .or. .not. allocated(this%Xg)) then
            allocate(elemConn(0))
            return
        end if
        elemConn=face_nodes(cell_nodes(this%ng,elem),[2,2,2],face)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the tensor degree retained on one face.
    !!
    !! Face numbering is `[w_min,w_max,v_min,v_max,u_min,u_max]`. The
    !! normal-direction entry is zero; for example a `w`-face returns
    !! \([p,q,0]\). Invalid input or an active diagnostic returns zeros.
    pure function cmp_degreeFace(this, face) result(degree)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: face
        integer :: degree(3)

        degree = 0
        if (.not. this%err%ok) return

        select case (face)
          case(1)
            degree = [this%degree(1), this%degree(2), 0]
          case(2)
            degree = [this%degree(1), this%degree(2), 0]
          case(3)
            degree = [this%degree(1), 0, this%degree(3)]
          case(4)
            degree = [this%degree(1), 0, this%degree(3)]
          case(5)
            degree = [0, this%degree(2), this%degree(3)]
          case(6)
            degree = [0, this%degree(2), this%degree(3)]
          case default
            degree = 0
        end select
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute element-local IGA basis data at one quadrature point.
    !!
    !! The volume must have exactly three physical coordinates. Element and
    !! local basis order agree with `cmp_elem`. `dV` includes tensor quadrature
    !! weight, all span maps, and the absolute physical Jacobian determinant.
    !! Gradients and optional Hessians are transformed to Cartesian coordinates;
    !! a singular Jacobian sets diagnostic code 108. With `strict=.true.`, the
    !! Jacobian determinant must also be positive, so an orientation-reversing
    !! integration point is rejected instead of contributing positive measure.
    !! The check is local to the requested quadrature point and does not prove
    !! global injectivity between quadrature points.
    !!
    !! With \(\mathbf J=\partial\mathbf V/\partial(u,v,w)\),
    !!
    !! \[
    !! \nabla_xR_a=\mathbf J^{-T}\nabla_{\boldsymbol{\xi}}R_a,\qquad
    !! dV=\omega_g\,\Delta u_e\,\Delta v_e\,\Delta w_e
    !!     |\det\mathbf J|.
    !! \]
    !!
    !! For Cartesian indices \(d,e\), the optional physical Hessian applies the
    !! inverse-map chain rule,
    !!
    !! \[
    !! \frac{\partial^2R_a}{\partial x_d\,\partial x_e}
    !!   =\sum_{\alpha,\beta=1}^{3}
    !!      \left(R_{a,\alpha\beta}
    !!       -\nabla_xR_a\mathbin{\cdot}\mathbf V_{,\alpha\beta}\right)
    !!      (\mathbf J^{-1})_{\alpha d}
    !!      (\mathbf J^{-1})_{\beta e}.
    !! \]
    !!
    !! IGA formulation reference: T. J. R. Hughes, J. A. Cottrell, and
    !! Y. Bazilevs, *CMAME* 194 (2005), 4135--4195.
    !! [doi:10.1016/j.cma.2004.10.008](https://doi.org/10.1016/j.cma.2004.10.008).
    !!
    !! Elementwise analysis assumes an invertible geometry map with a smooth
    !! inverse: Y. Bazilevs et al., *MMMAS* 16 (2006), 1031--1090.
    !! [doi:10.1142/S0218202506001455](https://doi.org/10.1142/S0218202506001455).
    pure subroutine ansatz(this, ie, ig, Tgc, dTgc_dXg, dV, ngauss, d2Tgc_dXg2, strict, err, Xg, dXg_dXt, &
        dmeasure_dXc, dmeasure_dlogWc, dgradient_dXc, dgradient_dlogWc, djacobian_dXc, djacobian_dlogWc)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: ie, ig
            !! Active element and flattened tensor quadrature-point indices, both one-based.
        real(rk), allocatable, target, intent(out) :: Tgc(:), dTgc_dXg(:,:)
            !! Local basis values `[nlocal]` and physical gradients `[nlocal,3]`.
        real(rk), intent(out) :: dV
            !! Weighted differential volume contribution.
        integer, intent(in), optional :: ngauss(3)
            !! Positive directional quadrature counts; defaults to `degree+1`.
        real(rk), allocatable, target, intent(out), optional :: d2Tgc_dXg2(:,:,:)
            !! Optional physical Hessians `[nlocal,3,3]`.
        logical, intent(in), optional :: strict
            !! Reject orientation reversal in addition to singularity.
        type(debug), intent(out), optional :: err
            !! Caller-owned diagnostic; bypasses and preserves the object diagnostic.
        real(rk), allocatable, target, intent(out), optional :: Xg(:), dXg_dXt(:,:)
            !! Optional physical point `[3]` and raw Jacobian `[3,3]`, coordinate first.
        real(rk), allocatable, target, intent(out), optional :: dmeasure_dXc(:,:), dmeasure_dlogWc(:)
            !! Weighted measure derivatives for the active Euclidean controls and log weights.
        real(rk), allocatable, intent(out), optional :: dgradient_dXc(:,:,:,:), dgradient_dlogWc(:,:,:)
            !! Physical basis gradient derivatives with respect to controls `(A,i,B,j)` and log weights `(A,i,B)`.
        real(rk), allocatable, intent(out), optional :: djacobian_dXc(:,:), djacobian_dlogWc(:,:,:)
            !! Raw Jacobian derivatives for controls (A,alpha) and log weights (i,alpha,A).

        call tensor_ansatz(this%Xc,this%Wc,this%rational,this%degree,this%nc,this%knot1,this%knot2, &
            this%err,err,ie,ig,Tgc,dTgc_dXg,dV,ngauss,strict,Xg,dXg_dXt,d2Tgc_dXg2, &
            dmeasure_dXc,dmeasure_dlogWc,dgradient_dXc,dgradient_dlogWc,this%knot3, &
            djacobian_dXc=djacobian_dXc,djacobian_dlogWc=djacobian_dlogWc)
    end subroutine ansatz


    !===============================================================================

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Integrate physical volume over all nonzero knot-span cells.
    !! Each cell uses the same tensor Gauss rule, defaulting to `degree+1`
    !! points per direction. This is an approximation for nonpolynomial maps.
    !! The local rational basis and scaled Jacobian frame are shared with IGA.
    !! Invalid quadrature points or Jacobians are omitted from the partial
    !! result and reported through the diagnostic. Reversed maps contribute
    !! positively unless `strict` is true. This does not certify injectivity.
    !! Caller-owned `err` is reset and leaves the object diagnostic untouched;
    !! without it an active object diagnostic returns zero. Caches are unchanged.
    pure subroutine cmp_volume(this, volume, ngauss, strict, err)
        class(nurbs_volume), target, intent(inout) :: this
        real(rk), intent(out) :: volume
            !! Accumulated nonnegative physical volume; partial on quadrature failure.
        integer, intent(in), optional :: ngauss(3)
            !! Positive directional quadrature counts; defaults to `degree+1`.
        logical, intent(in), optional :: strict
            !! Omit and diagnose orientation-reversed quadrature points.
        type(debug), target, intent(out), optional :: err
            !! Caller-owned diagnostic for shared-geometry integration.
        call tensor_measure(this%Xc,this%degree,this%nc,this%knot1,this%knot2,this%Wc,this%rational, &
            volume,this%err,err,ngauss,strict,this%knot3)
    end subroutine cmp_volume
    !===============================================================================


    !===============================================================================
    !===============================================================================

    !> Evaluate requested dense basis data with support-sized working arrays.
    !! Packed Hessian component `(a,b)` is at `(point,(a-1)*ncp+global,b)`.
    !! Dense output rows are scattered directly, avoiding strided helper temporaries.
    pure subroutine compute_basis_data_3d_vector(Xt,knot1,knot2,knot3,degree,nc,order, &
        Tgc,dTgc,d2Tgc,ng,wrap_parameters,Wc, side)
        real(rk), intent(in) :: Xt(:,:)
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3), nc(3), order
        real(rk), allocatable, intent(out), optional :: Tgc(:,:), dTgc(:,:,:), d2Tgc(:,:,:)
        integer, intent(in), optional :: ng(3)
        logical, intent(in), optional :: wrap_parameters(:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk) :: basis(0:product(degree+1)-1)
        real(rk) :: ders1(0:order,0:degree(1)), ders2(0:order,0:degree(2)), ders3(0:order,0:degree(3))
        real(rk) :: gradient(0:merge(product(degree+1),0,order>0)-1,3)
        real(rk) :: hessian(0:merge(product(degree+1),0,order>1)-1,3,3)
        logical :: wrap(3)
        integer :: point, npoint, ncp, first(3), local, global, l1, l2, l3, plane, a, b

        npoint=size(Xt,1)
        if (present(ng)) npoint=product(ng)
        ncp=product(nc)
        if (present(Tgc)) allocate(Tgc(npoint,ncp),source=0.0_rk)
        if (present(dTgc)) allocate(dTgc(npoint,ncp,3),source=0.0_rk)
        if (present(d2Tgc)) allocate(d2Tgc(npoint,3*ncp,3),source=0.0_rk)
        wrap=.false.
        if (present(wrap_parameters)) wrap(1:min(3,size(wrap_parameters)))= &
            wrap_parameters(1:min(3,size(wrap_parameters)))
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do point=1,npoint
#else
        do concurrent (point=1:npoint) local(basis,gradient,hessian,ders1,ders2,ders3,first,local,global,l1,l2,l3,plane,a,b)
#endif
            call tensor_basis_jet(Xt(point,:),knot1,knot2,knot3,degree,nc,order,wrap, &
                first,basis,gradient,hessian,ders1,ders2,ders3,Wc,side=side)
            do local=0,product(degree+1)-1
                l1=mod(local,degree(1)+1)
                plane=local/(degree(1)+1)
                l2=mod(plane,degree(2)+1)
                l3=plane/(degree(2)+1)
                global=first(1)+l1+nc(1)*((first(2)+l2-1)+nc(2)*(first(3)+l3-1))
                if (global<1 .or. global>ncp) cycle
                if (present(Tgc)) Tgc(point,global)=basis(local)
                if (present(dTgc)) dTgc(point,global,:)=gradient(local,:)
                if (present(d2Tgc)) then
                    do b=1,3
                        do a=1,b
                            d2Tgc(point,(a-1)*ncp+global,b)=hessian(local,a,b)
                            d2Tgc(point,(b-1)*ncp+global,a)=hessian(local,a,b)
                        end do
                    end do
                end if
            end do
        end do
    end subroutine compute_basis_data_3d_vector
    !===============================================================================

    !> Evaluate requested dense or selected basis data from a single local jet.
    !! Selected IDs preserve ordering/repeats; inactive and invalid entries are zero.
    !! Packed component `(a,b)` is at `((a-1)*nselected+i,b)`.
    pure subroutine compute_basis_data_3d_scalar(Xt,knot1,knot2,knot3,degree,nc,order, &
        Tgc,dTgc,d2Tgc,elem,wrap_parameters,Wc, side)
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3), nc(3), order
        real(rk), allocatable, intent(out), optional :: Tgc(:), dTgc(:,:), d2Tgc(:,:)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), optional :: wrap_parameters(:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk) :: basis(0:product(degree+1)-1)
        real(rk) :: ders1(0:order,0:degree(1)), ders2(0:order,0:degree(2)), ders3(0:order,0:degree(3))
        real(rk) :: gradient(0:merge(product(degree+1),0,order>0)-1,3)
        real(rk) :: hessian(0:merge(product(degree+1),0,order>1)-1,3,3)
        logical :: wrap(3)
        integer :: ncp, nselected, nwork, first(3), i, local, global, target, l1, l2, l3, plane, a, b

        ncp=product(nc)
        nselected=ncp
        nwork=product(degree+1)
        if (present(elem)) then
            nselected=size(elem)
            nwork=nselected
        end if
        if (present(Tgc)) allocate(Tgc(nselected),source=0.0_rk)
        if (present(dTgc)) allocate(dTgc(nselected,3),source=0.0_rk)
        if (present(d2Tgc)) allocate(d2Tgc(3*nselected,3),source=0.0_rk)
        wrap=.false.
        if (present(wrap_parameters)) wrap(1:min(3,size(wrap_parameters)))= &
            wrap_parameters(1:min(3,size(wrap_parameters)))
        call tensor_basis_jet(Xt,knot1,knot2,knot3,degree,nc,order,wrap,first,basis,gradient,hessian,ders1,ders2,ders3,Wc, &
                side=side)
        do i=1,nwork
            if (present(elem)) then
                global=elem(i)
                if (global<1 .or. global>ncp) cycle
                l1=mod(global-1,nc(1))+1-first(1)
                l2=mod((global-1)/nc(1),nc(2))+1-first(2)
                l3=(global-1)/(nc(1)*nc(2))+1-first(3)
                if (l1<0 .or. l1>degree(1) .or. l2<0 .or. l2>degree(2) .or. l3<0 .or. l3>degree(3)) cycle
                local=l1+(degree(1)+1)*(l2+(degree(2)+1)*l3)
                target=i
            else
                local=i-1
                l1=mod(local,degree(1)+1)
                plane=local/(degree(1)+1)
                l2=mod(plane,degree(2)+1)
                l3=plane/(degree(2)+1)
                target=first(1)+l1+nc(1)*((first(2)+l2-1)+nc(2)*(first(3)+l3-1))
                if (target<1 .or. target>ncp) cycle
            end if
            if (present(Tgc)) Tgc(target)=basis(local)
            if (present(dTgc)) dTgc(target,:)=gradient(local,:)
            if (present(d2Tgc)) then
                do b=1,3
                    do a=1,b
                        d2Tgc((a-1)*nselected+target,b)=hessian(local,a,b)
                        d2Tgc((b-1)*nselected+target,a)=hessian(local,a,b)
                    end do
                end do
            end if
        end do
    end subroutine compute_basis_data_3d_scalar
    !===============================================================================

    !> Contract the local basis directly with controls for paired parameter rows.
    pure function compute_Xg_3d(Xt,knot1,knot2,knot3,degree,nc,ng,Xc,wrap_parameters,Wc, side) result(Xg)
        real(rk), intent(in) :: Xt(:,:)
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:), Xc(:,:)
        integer, intent(in) :: degree(3), nc(3)
        integer, intent(in), optional :: ng(3)
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
            Xg(point,:)=compute_Xg_3d_1point(Xt(point,:),knot1,knot2,knot3,degree,nc,Xc,wrap_parameters,Wc,side)
        end do
    end function compute_Xg_3d
    !===============================================================================

    !> Contract one local basis with the control coordinates.
    pure function compute_Xg_3d_1point(Xt,knot1,knot2,knot3,degree,nc,Xc,wrap_parameters,Wc, side) result(Xg)
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:), Xc(:,:)
        integer, intent(in) :: degree(3), nc(3)
        logical, intent(in), optional :: wrap_parameters(:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk) :: Xg(size(Xc,2)), basis(0:product(degree+1)-1), unused_gradient(0,3), unused_hessian(0,3,3)
        real(rk) :: ders1(0:0,0:degree(1)), ders2(0:0,0:degree(2)), ders3(0:0,0:degree(3))
        logical :: wrap(3)
        integer :: first(3), local, global, l1, l2, l3, plane, d

        Xg=0.0_rk
        wrap=.false.
        if (present(wrap_parameters)) wrap(1:min(3,size(wrap_parameters)))= &
            wrap_parameters(1:min(3,size(wrap_parameters)))
        call tensor_basis_jet(Xt,knot1,knot2,knot3,degree,nc,0,wrap, &
            first,basis,unused_gradient,unused_hessian,ders1,ders2,ders3,Wc,side=side)
        do local=0,product(degree+1)-1
            l1=mod(local,degree(1)+1)
            plane=local/(degree(1)+1)
            l2=mod(plane,degree(2)+1)
            l3=plane/(degree(2)+1)
            global=first(1)+l1+nc(1)*((first(2)+l2-1)+nc(2)*(first(3)+l3-1))
            if (global<1 .or. global>product(nc)) cycle
            do d=1,size(Xc,2)
                Xg(d)=Xg(d)+basis(local)*Xc(global,d)
            end do
        end do
    end function compute_Xg_3d_1point
    !===============================================================================

    !> Positional adapter for the exported rational basis overload.
    pure function compute_Tgc_nurbs_3d_vector(Xt, knot1, knot2, knot3, degree, nc, ng, Wc, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: ng(3)
        real(rk), intent(in), contiguous :: Wc(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Tgc(:,:)

        call compute_basis_data(Xt,knot1,knot2,knot3,degree,nc,0, &
            Tgc=Tgc,ng=ng,wrap_parameters=wrap_parameters,Wc=Wc)
    end function compute_Tgc_nurbs_3d_vector
    !===============================================================================

    !> Positional adapter for the exported rational basis overload.
    !! Local weights normalize the supplied list; global weights normalize full support.
    pure function compute_Tgc_nurbs_3d_scalar(Xt, knot1, knot2, knot3, degree, nc, Wc, elem, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        real(rk), intent(in), contiguous :: Wc(:)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Tgc(:)
        real(rk) :: ders1(0:0,0:degree(1)), ders2(0:0,0:degree(2)), ders3(0:0,0:degree(3)), unused_gradient(0:-1,3)
        integer :: first(3)
        logical :: wrap(3)

        if (present(elem)) then
            if (size(Wc)==size(elem) .and. size(Wc)/=product(nc)) then
                allocate(Tgc(size(elem)))
                wrap=.false.
                if (present(wrap_parameters)) wrap(1:min(3,size(wrap_parameters)))= &
                    wrap_parameters(1:min(3,size(wrap_parameters)))
                call tensor_basis_jet(Xt,knot1,knot2,knot3,degree,nc,0,wrap,first,Tgc,unused_gradient, &
                    ders1=ders1,ders2=ders2,ders3=ders3,Wc=Wc,elem=elem)
                return
            end if
        end if
        call compute_basis_data(Xt,knot1,knot2,knot3,degree,nc,0, &
            Tgc=Tgc,elem=elem,wrap_parameters=wrap_parameters,Wc=Wc)
    end function compute_Tgc_nurbs_3d_scalar
    !===============================================================================

    !> Positional adapter for the exported polynomial basis overload.
    pure function compute_Tgc_bspline_3d_vector(Xt, knot1, knot2, knot3, degree, nc, ng, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: ng(3)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Tgc(:,:)

        call compute_basis_data(Xt,knot1,knot2,knot3,degree,nc,0, &
            Tgc=Tgc,ng=ng,wrap_parameters=wrap_parameters)
    end function compute_Tgc_bspline_3d_vector
    !===============================================================================

    !> Positional adapter for the exported polynomial basis overload.
    pure function compute_Tgc_bspline_3d_scalar(Xt, knot1, knot2, knot3, degree, nc, elem, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Tgc(:)

        call compute_basis_data(Xt,knot1,knot2,knot3,degree,nc,0, &
            Tgc=Tgc,elem=elem,wrap_parameters=wrap_parameters)
    end function compute_Tgc_bspline_3d_scalar
    !===============================================================================

    !> Positional adapter for the exported rational gradient overload.
    pure subroutine compute_dTgc_nurbs_3d_vector(Xt, knot1, knot2, knot3, degree, nc, ng, Wc, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: ng(3)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)

        call compute_basis_data(Xt,knot1,knot2,knot3,degree,nc,1, &
            Tgc=Tgc,dTgc=dTgc,ng=ng,wrap_parameters=wrap_parameters,Wc=Wc)
    end subroutine compute_dTgc_nurbs_3d_vector
    !===============================================================================

    !> Positional adapter for the exported rational gradient overload.
    !! Local weights normalize the supplied list; global weights normalize full support.
    pure subroutine compute_dTgc_nurbs_3d_scalar(Xt, knot1, knot2, knot3, degree, nc, Wc, dTgc, Tgc, elem, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        real(rk), intent(in), contiguous :: Wc(:)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:)

        if (present(elem)) then
            if (size(Wc)==size(elem) .and. size(Wc)/=product(nc)) then
                call compute_basis_data(Xt,knot1,knot2,knot3,degree,nc,1, &
                    Tgc=Tgc,dTgc=dTgc,elem=elem,wrap_parameters=wrap_parameters)
                call normalize_rational_basis(Tgc,Wc,dTgc)
                return
            end if
        end if
        call compute_basis_data(Xt,knot1,knot2,knot3,degree,nc,1, &
            Tgc=Tgc,dTgc=dTgc,elem=elem,wrap_parameters=wrap_parameters,Wc=Wc)
    end subroutine compute_dTgc_nurbs_3d_scalar
    !===============================================================================

    !> Positional adapter for the exported polynomial gradient overload.
    pure subroutine compute_dTgc_bspline_3d_vector(Xt, knot1, knot2, knot3, degree, nc, ng, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: ng(3)
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)

        call compute_basis_data(Xt,knot1,knot2,knot3,degree,nc,1, &
            Tgc=Tgc,dTgc=dTgc,ng=ng,wrap_parameters=wrap_parameters)
    end subroutine compute_dTgc_bspline_3d_vector
    !===============================================================================

    !> Positional adapter for the exported polynomial gradient overload.
    pure subroutine compute_dTgc_bspline_3d_scalar(Xt, knot1, knot2, knot3, degree, nc, dTgc, Tgc, elem, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:)

        call compute_basis_data(Xt,knot1,knot2,knot3,degree,nc,1, &
            Tgc=Tgc,dTgc=dTgc,elem=elem,wrap_parameters=wrap_parameters)
    end subroutine compute_dTgc_bspline_3d_scalar
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Fit polynomial tensor-product B-spline controls to parameterized volume data.
    !!
    !! `ndata=[n1,n2,n3]` defines the rows used. A recognized Cartesian grid is
    !! solved by three directional banded normal-equation factorizations. If a
    !! directional factorization fails, three directional Householder QR solves
    !! preserve Kronecker structure without forming the dense tensor design.
    !! The grid test requires exactly `product(ndata)` parameter and data rows.
    !! Repeated coordinates are compared componentwise against 32 times the
    !! larger of their ULPs and `epsilon(rk)` times that direction's sampled
    !! span. No absolute unit floor is used. Otherwise the first
    !! `product(ndata)` rows use one banded tensor-product normal system, with a
    !! dense pivoted-QR fallback allocated only after Cholesky failure.
    !!
    !! The routine fits a polynomial tensor-product B-spline model. A volume
    !! with nonuniform weights is rejected without changing its controls; use
    !! [[nurbs_volume:lsq_fit_nurbs]] for a rational model. Uniform weights
    !! produce the same polynomial basis. `Xt` and `Xdata` must provide enough finite rows,
    !! `Xt` must have at least three columns, `Xdata` at least one, and each
    !! directional parameter set must give full-rank active-domain coverage.
    !! These shape and finiteness preconditions are not completely diagnosed.
    !! Only `Xc` is replaced; cached `Xg` is not recomputed.
    !!
    !! All paths require numerical full column rank. The QR fallback diagnoses
    !! rank deficiency with code 108 instead of returning nonunique controls.
    pure subroutine lsq_fit_bspline(this, Xt, Xdata, ndata)
        use forcad_utils, only: solve_least_squares_qr, solve_spd_banded
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:,:)
            !! Parameter triples `[at least product(ndata),3]`.
        real(rk), intent(in), contiguous :: Xdata(:,:)
            !! Data coordinates `[at least product(ndata),ncoord]`.
        integer, intent(in) :: ndata(3)
            !! Directional data counts `[3]`, each at least the corresponding `nc`.
        real(rk), allocatable :: TtT(:,:), TtX(:,:), Xsol(:,:), design(:,:)
        real(rk), allocatable :: G1band(:,:), G2band(:,:), G3band(:,:), RHS(:,:)
        real(rk), allocatable :: rhs1(:,:), rhs2(:,:), rhs3(:,:), sol1(:,:), sol2(:,:), sol3(:,:)
        real(rk), allocatable :: Phi1(:,:), Phi2(:,:), Phi3(:,:)
        real(rk), allocatable :: design1(:,:), design2(:,:), design3(:,:)
        real(rk), allocatable :: qr_rhs1(:,:), qr_rhs2(:,:), qr_rhs3(:,:)
        real(rk) :: N1(0:this%degree(1)), N2(0:this%degree(2)), N3(0:this%degree(3))
        real(rk) :: ba, bb, n32a, n32b
        integer, allocatable :: F1(:), F2(:), F3(:)
        real(rk) :: grid_tol, grid_scale(3), grid_value, reference_value
        integer :: i, n, ncp, dim_, row, col, i1, i2, i3, bandwidth
        integer :: a1, a2, a3, b1, b2, b3, k, ia, ib, first1, first2, first3
        logical :: tensor_grid, normal_ok

        if (.not. this%err%ok) return

        if (this%is_rational()) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Polynomial fitting cannot retain nonuniform rational weights.",&
                location   = "lsq_fit_bspline",&
                suggestion = "Use lsq_fit_nurbs() or reset the volume with uniform weights.")
            return
        end if

        if (this%nc(1) > ndata(1)) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid number of control points in the first direction.",&
                location   = "lsq_fit_bspline",&
                suggestion = "Ensure that the number of control points does not exceed the number of data points.")
            return
        end if
        if (this%nc(2) > ndata(2)) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid number of control points in the second direction.",&
                location   = "lsq_fit_bspline",&
                suggestion = "Ensure that the number of control points does not exceed the number of data points.")
            return
        end if
        if (this%nc(3) > ndata(3)) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid number of control points in the third direction.",&
                location   = "lsq_fit_bspline",&
                suggestion = "Ensure that the number of control points does not exceed the number of data points.")
            return
        end if

        n = ndata(1)*ndata(2)*ndata(3)
        ncp = this%nc(1)*this%nc(2)*this%nc(3)
        dim_ = size(Xdata, 2)

        tensor_grid = size(Xt,1) == n .and. size(Xdata,1) == n
        if (tensor_grid) then
            grid_scale(1) = maxval(Xt(:,1)) - minval(Xt(:,1))
            grid_scale(2) = maxval(Xt(:,2)) - minval(Xt(:,2))
            grid_scale(3) = maxval(Xt(:,3)) - minval(Xt(:,3))
            do i3 = 1, ndata(3)
                do i2 = 1, ndata(2)
                    do i1 = 1, ndata(1)
                        row = ((i3-1)*ndata(2) + i2-1)*ndata(1) + i1
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
                        grid_value = Xt(row,3)
                        reference_value = Xt(((i3-1)*ndata(2))*ndata(1)+1,3)
                        grid_tol = 32.0_rk*max(&
                            epsilon(1.0_rk)*grid_scale(3),&
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
                if (.not. tensor_grid) exit
            end do
        end if

        if (tensor_grid) then
            allocate(Phi1(0:this%degree(1), ndata(1)), Phi2(0:this%degree(2), ndata(2)), Phi3(0:this%degree(3), ndata(3)))
            allocate(F1(ndata(1)), F2(ndata(2)), F3(ndata(3)))
            allocate(G1band(0:this%degree(1), this%nc(1)), G2band(0:this%degree(2), this%nc(2)), &
                G3band(0:this%degree(3), this%nc(3)), RHS(ncp, dim_), source=0.0_rk)

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

            do i3 = 1, ndata(3)
                row = ((i3-1)*ndata(2))*ndata(1) + 1
                call basis_bspline_der(&
                    map_parameter(Xt(row,3), this%knot3, this%nc(3), this%degree(3), this%wrap_parameters(3)),&
                    this%knot3, this%nc(3), this%degree(3), 0, F3(i3), Phi3(:,i3))
                do a3 = 0, this%degree(3)
                    ia = F3(i3) + a3
                    if (ia < 1 .or. ia > this%nc(3)) cycle
                    ba = Phi3(a3,i3)
                    do b3 = 0, this%degree(3)
                        ib = F3(i3) + b3
                        if (ib < 1 .or. ib > this%nc(3)) cycle
                        if (ia >= ib) G3band(ia-ib,ib) = G3band(ia-ib,ib) + ba*Phi3(b3,i3)
                    end do
                end do
            end do

            do i3 = 1, ndata(3)
                do i2 = 1, ndata(2)
                    do i1 = 1, ndata(1)
                        row = ((i3-1)*ndata(2) + i2-1)*ndata(1) + i1
                        do a3 = 0, this%degree(3)
                            if (F3(i3)+a3 < 1 .or. F3(i3)+a3 > this%nc(3)) cycle
                            do a2 = 0, this%degree(2)
                                if (F2(i2)+a2 < 1 .or. F2(i2)+a2 > this%nc(2)) cycle
                                n32a = Phi3(a3,i3)*Phi2(a2,i2)
                                do a1 = 0, this%degree(1)
                                    if (F1(i1)+a1 < 1 .or. F1(i1)+a1 > this%nc(1)) cycle
                                    ia = ((F3(i3)+a3-1)*this%nc(2) + F2(i2)+a2-1)*this%nc(1) + F1(i1) + a1
                                    ba = n32a*Phi1(a1,i1)
                                    do k = 1, dim_
                                        RHS(ia,k) = RHS(ia,k) + ba*Xdata(row,k)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do

            allocate(Xsol(ncp,dim_))
            allocate(rhs1(this%nc(1), this%nc(2)*this%nc(3)))
            allocate(rhs2(this%nc(2), this%nc(1)*this%nc(3)))
            allocate(rhs3(this%nc(3), this%nc(1)*this%nc(2)))

            normal_ok = .true.
            do k = 1, dim_
                do i3 = 1, this%nc(3)
                    do i2 = 1, this%nc(2)
                        col = (i3-1)*this%nc(2) + i2
                        do i1 = 1, this%nc(1)
                            rhs1(i1,col) = RHS(((i3-1)*this%nc(2)+i2-1)*this%nc(1)+i1,k)
                        end do
                    end do
                end do

                sol1 = solve_spd_banded(G1band, rhs1)
                if (size(sol1,1) /= this%nc(1) .or. size(sol1,2) /= this%nc(2)*this%nc(3)) then
                    normal_ok = .false.
                    exit
                end if
                do i3 = 1, this%nc(3)
                    do i1 = 1, this%nc(1)
                        col = (i3-1)*this%nc(1) + i1
                        do i2 = 1, this%nc(2)
                            rhs2(i2,col) = sol1(i1,(i3-1)*this%nc(2)+i2)
                        end do
                    end do
                end do

                sol2 = solve_spd_banded(G2band, rhs2)
                if (size(sol2,1) /= this%nc(2) .or. size(sol2,2) /= this%nc(1)*this%nc(3)) then
                    normal_ok = .false.
                    exit
                end if
                do i2 = 1, this%nc(2)
                    do i1 = 1, this%nc(1)
                        col = (i2-1)*this%nc(1) + i1
                        do i3 = 1, this%nc(3)
                            rhs3(i3,col) = sol2(i2,(i3-1)*this%nc(1)+i1)
                        end do
                    end do
                end do

                sol3 = solve_spd_banded(G3band, rhs3)
                if (size(sol3,1) /= this%nc(3) .or. size(sol3,2) /= this%nc(1)*this%nc(2)) then
                    normal_ok = .false.
                    exit
                end if
                do i3 = 1, this%nc(3)
                    do i2 = 1, this%nc(2)
                        do i1 = 1, this%nc(1)
                            Xsol(((i3-1)*this%nc(2)+i2-1)*this%nc(1)+i1,k) = &
                                sol3(i3,(i2-1)*this%nc(1)+i1)
                        end do
                    end do
                end do
            end do
            if (normal_ok) then
                this%Xc = Xsol
                return
            end if

            allocate(design1(ndata(1),this%nc(1)), source=0.0_rk)
            allocate(design2(ndata(2),this%nc(2)), source=0.0_rk)
            allocate(design3(ndata(3),this%nc(3)), source=0.0_rk)
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
            do i3 = 1, ndata(3)
                do a3 = 0, this%degree(3)
                    ia = F3(i3) + a3
                    if (ia < 1 .or. ia > this%nc(3)) cycle
                    design3(i3,ia) = Phi3(a3,i3)
                end do
            end do

            allocate(qr_rhs1(ndata(1),ndata(2)*ndata(3)*dim_))
            do k = 1, dim_
                do i3 = 1, ndata(3)
                    do i2 = 1, ndata(2)
                        col = ((k-1)*ndata(3)+i3-1)*ndata(2) + i2
                        do i1 = 1, ndata(1)
                            row = ((i3-1)*ndata(2)+i2-1)*ndata(1) + i1
                            qr_rhs1(i1,col) = Xdata(row,k)
                        end do
                    end do
                end do
            end do
            sol1 = solve_least_squares_qr(design1, qr_rhs1)
            if (size(sol1,1) /= this%nc(1) .or. size(sol1,2) /= ndata(2)*ndata(3)*dim_) then
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "The first directional design is numerically rank deficient.",&
                    location   = "lsq_fit_bspline",&
                    suggestion = "Use distinct first-direction parameters with full spline-space coverage.")
                return
            end if

            allocate(qr_rhs2(ndata(2),this%nc(1)*ndata(3)*dim_))
            do k = 1, dim_
                do i3 = 1, ndata(3)
                    do i1 = 1, this%nc(1)
                        col = ((k-1)*ndata(3)+i3-1)*this%nc(1) + i1
                        do i2 = 1, ndata(2)
                            qr_rhs2(i2,col) = sol1(i1,&
                                ((k-1)*ndata(3)+i3-1)*ndata(2)+i2)
                        end do
                    end do
                end do
            end do
            sol2 = solve_least_squares_qr(design2, qr_rhs2)
            if (size(sol2,1) /= this%nc(2) .or. size(sol2,2) /= this%nc(1)*ndata(3)*dim_) then
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "The second directional design is numerically rank deficient.",&
                    location   = "lsq_fit_bspline",&
                    suggestion = "Use distinct second-direction parameters with full spline-space coverage.")
                return
            end if

            allocate(qr_rhs3(ndata(3),this%nc(1)*this%nc(2)*dim_))
            do k = 1, dim_
                do i2 = 1, this%nc(2)
                    do i1 = 1, this%nc(1)
                        col = ((k-1)*this%nc(2)+i2-1)*this%nc(1) + i1
                        do i3 = 1, ndata(3)
                            qr_rhs3(i3,col) = sol2(i2,&
                                ((k-1)*ndata(3)+i3-1)*this%nc(1)+i1)
                        end do
                    end do
                end do
            end do
            sol3 = solve_least_squares_qr(design3, qr_rhs3)
            if (size(sol3,1) /= this%nc(3) .or. size(sol3,2) /= this%nc(1)*this%nc(2)*dim_) then
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "The third directional design is numerically rank deficient.",&
                    location   = "lsq_fit_bspline",&
                    suggestion = "Use distinct third-direction parameters with full spline-space coverage.")
                return
            end if
            do k = 1, dim_
                do i3 = 1, this%nc(3)
                    do i2 = 1, this%nc(2)
                        do i1 = 1, this%nc(1)
                            Xsol(((i3-1)*this%nc(2)+i2-1)*this%nc(1)+i1,k) = &
                                sol3(i3,((k-1)*this%nc(2)+i2-1)*this%nc(1)+i1)
                        end do
                    end do
                end do
            end do
            this%Xc = Xsol
            return
        end if

        bandwidth = this%degree(3)*this%nc(1)*this%nc(2) + this%degree(2)*this%nc(1) + this%degree(1)
        allocate(TtT(0:bandwidth, ncp), TtX(ncp, dim_), source=0.0_rk)
        do i = 1, n
            call basis_bspline_der(&
                map_parameter(Xt(i,1), this%knot1, this%nc(1), this%degree(1), this%wrap_parameters(1)),&
                this%knot1, this%nc(1), this%degree(1), 0, first1, N1)
            call basis_bspline_der(&
                map_parameter(Xt(i,2), this%knot2, this%nc(2), this%degree(2), this%wrap_parameters(2)),&
                this%knot2, this%nc(2), this%degree(2), 0, first2, N2)
            call basis_bspline_der(&
                map_parameter(Xt(i,3), this%knot3, this%nc(3), this%degree(3), this%wrap_parameters(3)),&
                this%knot3, this%nc(3), this%degree(3), 0, first3, N3)
            do a3 = 0, this%degree(3)
                if (first3+a3 < 1 .or. first3+a3 > this%nc(3)) cycle
                do a2 = 0, this%degree(2)
                    if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                    n32a = N3(a3)*N2(a2)
                    do a1 = 0, this%degree(1)
                        if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                        ia = ((first3+a3-1)*this%nc(2) + first2+a2-1)*this%nc(1) + first1 + a1
                        ba = n32a*N1(a1)
                        do k = 1, dim_
                            TtX(ia,k) = TtX(ia,k) + ba*Xdata(i,k)
                        end do
                        do b3 = 0, this%degree(3)
                            if (first3+b3 < 1 .or. first3+b3 > this%nc(3)) cycle
                            do b2 = 0, this%degree(2)
                                if (first2+b2 < 1 .or. first2+b2 > this%nc(2)) cycle
                                n32b = N3(b3)*N2(b2)
                                do b1 = 0, this%degree(1)
                                    if (first1+b1 < 1 .or. first1+b1 > this%nc(1)) cycle
                                    ib = ((first3+b3-1)*this%nc(2) + first2+b2-1)*this%nc(1) + first1 + b1
                                    bb = n32b*N1(b1)
                                    if (ia >= ib) TtT(ia-ib,ib) = TtT(ia-ib,ib) + ba*bb
                                end do
                            end do
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
                call basis_bspline_der(&
                    map_parameter(Xt(i,3), this%knot3, this%nc(3), this%degree(3), this%wrap_parameters(3)),&
                    this%knot3, this%nc(3), this%degree(3), 0, first3, N3)
                do a3 = 0, this%degree(3)
                    if (first3+a3 < 1 .or. first3+a3 > this%nc(3)) cycle
                    do a2 = 0, this%degree(2)
                        if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                        n32a = N3(a3)*N2(a2)
                        do a1 = 0, this%degree(1)
                            if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                            ia = ((first3+a3-1)*this%nc(2)+first2+a2-1)*this%nc(1) + first1 + a1
                            design(i,ia) = n32a*N1(a1)
                        end do
                    end do
                end do
            end do
            Xsol = solve_least_squares_qr(design, Xdata(1:n,:))
        end if
        if (size(Xsol,1) /= ncp .or. size(Xsol,2) /= dim_) then
            call this%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
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
    !> Alternately fit volume control points and positive rational weights.
    !!
    !! Banded normal equations update controls and log-weights. Mean log-weight
    !! removal fixes projective scale at geometric mean weight one. Each trial
    !! weight vector is used to recompute its least-squares controls before the
    !! residual is tested. A trial is accepted only when its residual norm does
    !! not increase; otherwise controls and weights remain at the last accepted
    !! iterate and the damping is increased. Accepted trials reduce the
    !! damping. Log-weight steps outside the finite exponential range are
    !! rejected before exponentiation.
    !!
    !! With \(n=\mathtt{product(ndata)}\) and residuals
    !! \(r_{ik}=V_k(u_i,v_i,w_i)-X_{ik}\), the monitored quantity is
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
    !! Inputs must provide enough finite rows and at least three parameter
    !! columns and one coordinate column, with `ndata>=nc` componentwise,
    !! `maxit>0`, `tol>=0`, and finite, nonnegative
    !! regularization/damping values. The objective is nonconvex.
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
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous  :: Xt(:,:)
            !! Parameter triples `[at least product(ndata),3]`.
        real(rk), intent(in), contiguous  :: Xdata(:,:)
            !! Data coordinates `[at least product(ndata),ncoord]`.
        integer,  intent(in)              :: ndata(3)
            !! Directional data counts `[3]`.
        integer,  intent(in),  optional   :: maxit
            !! Optional positive iteration limit; default 30.
        real(rk), intent(in),  optional   :: tol
            !! Optional nonnegative relative decrease tolerance; default `sqrt(epsilon(rk))`.
        real(rk), intent(in),  optional   :: lambda_xc
            !! Optional nonnegative control-point diagonal regularization.
        real(rk), intent(in),  optional   :: mu0
            !! Optional nonnegative initial log-weight damping.
        real(rk), intent(in),  optional   :: reg_logw
            !! Optional nonnegative additional log-weight regularization.
        real(rk), parameter :: damping_max = sqrt(huge(1.0_rk))
        real(rk), parameter :: log_weight_max = log(huge(1.0_rk)) - 1.0_rk
        real(rk), parameter :: log_weight_min = log(tiny(1.0_rk)) + 1.0_rk
        real(rk), allocatable :: TtT(:,:), TtX(:,:), v(:), Xsol(:,:)
        real(rk), allocatable :: Jband(:,:), Jrhs(:,:), delta_u(:,:)
        real(rk) :: N1(0:this%degree(1)), N2(0:this%degree(2)), N3(0:this%degree(3))
        real(rk) :: Cg(size(Xdata, 2))
        real(rk) :: bval, cost_prev, cost_now, cost_sum, cost_reduction, damp, data_scale, denom
        real(rk) :: epss, ha, hb, n32a, n32b
        real(rk) :: resid, ta, tb, tol_, lamx_, mu, regw
        integer :: dim_, it, maxit_, n, ncp, i, j, k, a1, a2, a3, b1, b2, b3
        integer :: ia, ib, first1, first2, first3, bandwidth, n12
        logical :: converged, trial_pending

        if (.not. this%err%ok) return

        dim_   = size(Xdata, 2)
        n      = ndata(1) * ndata(2) * ndata(3)
        n12    = this%nc(1) * this%nc(2)
        ncp    = n12 * this%nc(3)
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
                category   = "forcad_nurbs_volume",&
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
                category   = "forcad_nurbs_volume",&
                message    = "Too few samples in dir-1.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Use nc(1) <= ndata(1).")
            return
        end if
        if (this%nc(2) > ndata(2)) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Too few samples in dir-2.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Use nc(2) <= ndata(2).")
            return
        end if
        if (this%nc(3) > ndata(3)) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Too few samples in dir-3.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Use nc(3) <= ndata(3).")
            return
        end if
        if (n <= 0 .or. ncp < 2) then
            call this%err%set(&
                code       = 106,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
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

        bandwidth = this%degree(3)*n12 + this%degree(2)*this%nc(1) + this%degree(1)
        allocate(TtT(0:bandwidth,ncp), TtX(ncp,dim_))
        allocate(Jband(0:bandwidth,ncp), Jrhs(ncp,1))

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
                call basis_bspline_der(&
                    map_parameter(Xt(i,3), this%knot3, this%nc(3), this%degree(3), this%wrap_parameters(3)),&
                    this%knot3, this%nc(3), this%degree(3), 0, first3, N3)
                denom = 0.0_rk
                do a3 = 0, this%degree(3)
                    if (first3+a3 < 1 .or. first3+a3 > this%nc(3)) cycle
                    do a2 = 0, this%degree(2)
                        if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                        n32a = N3(a3)*N2(a2)
                        do a1 = 0, this%degree(1)
                            if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                            ia = ((first3+a3-1)*this%nc(2) + first2+a2-1)*this%nc(1) + first1 + a1
                            denom = denom + n32a*N1(a1)*this%Wc(ia)
                        end do
                    end do
                end do
                if (abs(denom) < epss) denom = sign(epss, denom)
                do a3 = 0, this%degree(3)
                    if (first3+a3 < 1 .or. first3+a3 > this%nc(3)) cycle
                    do a2 = 0, this%degree(2)
                        if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                        n32a = N3(a3)*N2(a2)
                        do a1 = 0, this%degree(1)
                            if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                            ia = ((first3+a3-1)*this%nc(2) + first2+a2-1)*this%nc(1) + first1 + a1
                            ta = n32a*N1(a1)*this%Wc(ia)/denom
                            do k = 1, dim_
                                TtX(ia,k) = TtX(ia,k) + ta*Xdata(i,k)
                            end do
                            do b3 = 0, this%degree(3)
                                if (first3+b3 < 1 .or. first3+b3 > this%nc(3)) cycle
                                do b2 = 0, this%degree(2)
                                    if (first2+b2 < 1 .or. first2+b2 > this%nc(2)) cycle
                                    n32b = N3(b3)*N2(b2)
                                    do b1 = 0, this%degree(1)
                                        if (first1+b1 < 1 .or. first1+b1 > this%nc(1)) cycle
                                        ib = ((first3+b3-1)*this%nc(2) + first2+b2-1)*this%nc(1) + first1 + b1
                                        if (ia >= ib) then
                                            tb = n32b*N1(b1)*this%Wc(ib)/denom
                                            TtT(ia-ib,ib) = TtT(ia-ib,ib) + ta*tb
                                        end if
                                    end do
                                end do
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
                    category   = "forcad_nurbs_volume",&
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
                call basis_bspline_der(&
                    map_parameter(Xt(i,3), this%knot3, this%nc(3), this%degree(3), this%wrap_parameters(3)),&
                    this%knot3, this%nc(3), this%degree(3), 0, first3, N3)
                denom = 0.0_rk
                do a3 = 0, this%degree(3)
                    if (first3+a3 < 1 .or. first3+a3 > this%nc(3)) cycle
                    do a2 = 0, this%degree(2)
                        if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                        n32a = N3(a3)*N2(a2)
                        do a1 = 0, this%degree(1)
                            if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                            ia = ((first3+a3-1)*this%nc(2) + first2+a2-1)*this%nc(1) + first1 + a1
                            denom = denom + n32a*N1(a1)*this%Wc(ia)
                        end do
                    end do
                end do
                if (abs(denom) < epss) denom = sign(epss, denom)
                Cg = 0.0_rk
                do a3 = 0, this%degree(3)
                    if (first3+a3 < 1 .or. first3+a3 > this%nc(3)) cycle
                    do a2 = 0, this%degree(2)
                        if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                        n32a = N3(a3)*N2(a2)
                        do a1 = 0, this%degree(1)
                            if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                            ia = ((first3+a3-1)*this%nc(2) + first2+a2-1)*this%nc(1) + first1 + a1
                            ta = n32a*N1(a1)*this%Wc(ia)/denom
                            do k = 1, dim_
                                Cg(k) = Cg(k) + ta*Xsol(ia,k)
                            end do
                        end do
                    end do
                end do
                do k = 1, dim_
                    resid = Cg(k) - Xdata(i,k)
                    cost_sum = cost_sum + resid*resid
                    do a3 = 0, this%degree(3)
                        if (first3+a3 < 1 .or. first3+a3 > this%nc(3)) cycle
                        do a2 = 0, this%degree(2)
                            if (first2+a2 < 1 .or. first2+a2 > this%nc(2)) cycle
                            n32a = N3(a3)*N2(a2)
                            do a1 = 0, this%degree(1)
                                if (first1+a1 < 1 .or. first1+a1 > this%nc(1)) cycle
                                ia = ((first3+a3-1)*this%nc(2) + first2+a2-1)*this%nc(1) + first1 + a1
                                bval = n32a*N1(a1)
                                ha = bval*this%Wc(ia)*(Xsol(ia,k) - Cg(k))/denom
                                Jrhs(ia,1) = Jrhs(ia,1) + ha*resid
                                do b3 = 0, this%degree(3)
                                    if (first3+b3 < 1 .or. first3+b3 > this%nc(3)) cycle
                                    do b2 = 0, this%degree(2)
                                        if (first2+b2 < 1 .or. first2+b2 > this%nc(2)) cycle
                                        n32b = N3(b3)*N2(b2)
                                        do b1 = 0, this%degree(1)
                                            if (first1+b1 < 1 .or. first1+b1 > this%nc(1)) cycle
                                            ib = ((first3+b3-1)*this%nc(2) + first2+b2-1)*this%nc(1) + first1 + b1
                                            if (ia >= ib) then
                                                hb = n32b*N1(b1)*this%Wc(ib)*(Xsol(ib,k) - Cg(k))/denom
                                                Jband(ia-ib,ib) = Jband(ia-ib,ib) + ha*hb
                                            end if
                                        end do
                                    end do
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
                        category   = "forcad_nurbs_volume",&
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
                    category   = "forcad_nurbs_volume",&
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
                category   = "forcad_nurbs_volume",&
                message    = "Nonlinear NURBS fitting did not converge.",&
                location   = "lsq_fit_nurbs",&
                suggestion = "Increase maxit or mu0, or relax tol for this data set.")
        end if
    end subroutine
    !===============================================================================
end module forcad_nurbs_volume
