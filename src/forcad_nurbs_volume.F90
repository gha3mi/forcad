!> author: Seyed Ali Ghasemi
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
!! sums multiply active weights by one exact radix power selected from their
!! largest model exponent. Directional homogeneous refinement applies this
!! normalization to the complete weight vector and restores its original gauge
!! after dehomogenization. Geometric derivatives use nested division and never
!! explicitly form \(W^2\).
!!
!! Direction 1 varies fastest in every flattened tensor array:
!! `linear(i,j,k)=i+nc(1)*((j-1)+nc(2)*(k-1))`. Omitting weights, or using
!! equal weights, selects the polynomial tensor-product B-spline form.
!! More precisely, supplied weights are classified as uniform when their
!! maximum deviation from the first weight is at most
!! `32*epsilon(rk)*maxval(Wc)`; that path treats them as exactly equal and
!! omits rational normalization.
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
!! smoother. At a discontinuous derivative, evaluation follows the selected
!! one-sided span convention.
!!
!! Directional refinement applies an exact one-dimensional operator to every
!! homogeneous control line, preserving
!! \(\mathbf V_{\mathrm{old}}(u,v,w)=\mathbf V_{\mathrm{new}}(u,v,w)\) to
!! working precision on the active domain. Insertion and degree elevation
!! preserve a verified periodic direction. Other non-clamped directional
!! representations are converted by degree elevation to an equivalent
!! active-domain open-clamped form. Periodic directional knot removal is
!! rejected because the bounded removal kernel cannot preserve cyclic
!! identities.
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
    use forcad_utils, only: basis_bspline, elemConn_C0, kron, ndgrid, compute_multiplicity, compute_knot_vector, &
        basis_bspline_der, insert_knot_A_5_1, insert_knot_periodic_A_5_1, findspan, elevate_degree_A_5_9, &
        hexahedron_Xc, remove_knots_A_5_8, &
        elemConn_Cn, unique, rotate_points, gauss_leg, export_vtk_legacy, basis_bspline_2der, fill_uniform, &
        valid_knot_vector, knot_start, knot_end, knot_tolerance, active_knots, active_span_count, active_knot_multiplicity, &
        infer_knot_shape, sparse_left_matmul, show_pyvista_singlepatch, basis_bspline_der_all_active, &
        tensor_basis_derivative_local, tensor_basis_derivatives2_local, map_parameter, periodic_topology
    use fordebug, only: debug

    implicit none

    private
    public nurbs_volume, compute_Tgc, compute_dTgc

    !===============================================================================
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
    !! @warning Core-data setters, fitting, control-lattice transformations,
    !! and weight changes do not automatically clear or recompute every sampled
    !! or connectivity cache. Call `create` after a geometry change before
    !! using cached `Xg`, and rebuild connectivity after a control-count or
    !! spline-space change.
    !! @endwarning
    !!
    !! Diagnostic codes are: 100 invalid argument or index; 101 shape mismatch;
    !! 102--105 missing geometry state; 106 underdetermined fitting problem; 107
    !! invalid parameters; 108 singular system or Jacobian; 109 iterative
    !! nonconvergence; and 110 unsupported trivariate IGES export.
    !!
    !! @note Getter generics return allocatable copies of object-owned arrays.
    !! @endnote
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
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
        integer, private :: ng(3) = 0              !! Cached sample counts by direction.
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
        procedure :: cmp_Xg                !!> Evaluate physical points at one or more parameter triples.
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
        procedure :: export_Xth_in_Xg      !!> Export parameter-grid cells embedded in physical space.
        procedure :: export_iges           !!> Report that trivariate IGES export is unsupported.
        procedure :: modify_Xc             !!> Replace one stored control coordinate.
        procedure :: modify_Wc             !!> Replace one positive weight and update rational state.
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
        procedure :: is_valid              !!> Report whether complete, internally consistent volume geometry is stored.
        procedure :: is_rational           !!> Report whether nonuniform weights affect the current basis.
        procedure :: get_parameter_wrapping !!> Report modulo parameter mapping in one direction.
        procedure :: get_parameter_topology !!> Classify one stored parameter topology.
        procedure :: put_to_nurbs          !!> Map structured points to the current parameter domain and cache them.
        procedure :: remove_knots          !!> Remove directional knots that pass an internal geometry test.
        procedure :: rotate_Xc             !!> Rotate control points in three-dimensional physical space.
        procedure :: rotate_Xg             !!> Rotate cached geometry points in three-dimensional physical space.
        procedure :: translate_Xc          !!> Translate the control lattice by a physical-space vector.
        procedure :: translate_Xg          !!> Translate cached geometry points by a physical-space vector.
        procedure :: show                  !!> Display previously exported representations with PyVista.
        procedure :: nearest_point         !!> Select the nearest cached geometry sample.
        procedure :: nearest_point2        !!> Minimize squared distance with bounded Newton iterations and backtracking.
        procedure :: ansatz                !!> Evaluate element-local physical basis data and differential volume.
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
        module procedure compute_Xg_nurbs_3d
        module procedure compute_Xg_bspline_3d
        module procedure compute_Xg_nurbs_3d_1point
        module procedure compute_Xg_bspline_3d_1point
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

    !> Evaluate dense parametric Hessians of the volume basis.
    !!
    !! For `ncp=product(nc)`, scalar output uses shape `[3*ncp,3]` and stores
    !! \(\partial^2R_A/(\partial\xi_a\partial\xi_b)\) at
    !! `d2Tgc((a-1)*ncp+A,b)`. Vector output prepends the parameter-point
    !! dimension. All mixed entries are represented explicitly. This
    !! implementation generic is private; use [[nurbs_volume:derivative_order]]
    !! when only one derivative multi-index is required.
    interface compute_d2Tgc
        module procedure compute_d2Tgc_nurbs_3d_vector
        module procedure compute_d2Tgc_bspline_3d_vector
        module procedure compute_d2Tgc_nurbs_3d_scalar
        module procedure compute_d2Tgc_bspline_3d_scalar
    end interface

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Test whether every supplied volume weight is finite and strictly positive.
    !! Public setters separately validate the required control count.
    pure logical function valid_weights(Wc) result(ok)
        real(rk), intent(in), contiguous :: Wc(:)
            !! Candidate flattened control weights.

        ok = all(ieee_is_finite(Wc) .and. Wc > 0.0_rk)
    end function
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
    !> Pack volume control lines for one-dimensional refinement.
    !!
    !! The selected direction becomes the first dimension of `Xdir`; both
    !! orthogonal control directions and all physical coordinates are packed
    !! into columns. Rational input stores \(\widehat w\mathbf P\) and
    !! \(\widehat w\), with one exact radix-power normalization shared by the
    !! complete control lattice. `weight_exponent` records that gauge for
    !! restoration by [[unpack_refine_volume]]. `dir` must be in `1:3` and `Wc`
    !! must be present when `rational` is true.
    pure subroutine pack_refine_volume(Xc, Wc, nc, dir, rational, Xdir, weight_exponent)
        real(rk), intent(in), contiguous :: Xc(:,:)
            !! Direction-1-fastest control lattice `[product(nc),ncoord]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Flattened weights `[product(nc)]` for rational input.
        integer, intent(in), contiguous :: nc(:)
            !! Directional control counts `[n1,n2,n3]`.
        integer, intent(in) :: dir
            !! Refinement direction in `1:3`.
        logical, intent(in) :: rational
            !! Pack homogeneous controls when true.
        real(rk), allocatable, intent(out) :: Xdir(:,:)
            !! Packed directional control lines.
        integer, intent(out), optional :: weight_exponent
            !! Radix exponent removed from rational weights.
        integer :: n1, n2, n3, ncoord, nval, i1, i2, i3, c, weight_exponent_

        n1 = nc(1)
        n2 = nc(2)
        n3 = nc(3)
        ncoord = size(Xc,2)
        nval = ncoord + merge(1, 0, rational)
        weight_exponent_ = 0
        if (rational) weight_exponent_ = exponent(maxval(Wc))
        if (present(weight_exponent)) weight_exponent = weight_exponent_

        select case (dir)
        case (1)
            allocate(Xdir(n1, n2*n3*nval))
            if (rational) then
                do concurrent (c = 1:ncoord, i3 = 1:n3, i2 = 1:n2, i1 = 1:n1)
                    Xdir(i1, i2 + (i3-1)*n2 + (c-1)*n2*n3) = &
                        Xc((i3-1)*n1*n2 + (i2-1)*n1 + i1, c)*&
                        scale(Wc((i3-1)*n1*n2 + (i2-1)*n1 + i1), -weight_exponent_)
                end do
                do concurrent (i3 = 1:n3, i2 = 1:n2, i1 = 1:n1)
                    Xdir(i1, i2 + (i3-1)*n2 + ncoord*n2*n3) = &
                        scale(Wc((i3-1)*n1*n2 + (i2-1)*n1 + i1), -weight_exponent_)
                end do
            else
                do concurrent (c = 1:ncoord, i3 = 1:n3, i2 = 1:n2, i1 = 1:n1)
                    Xdir(i1, i2 + (i3-1)*n2 + (c-1)*n2*n3) = Xc((i3-1)*n1*n2 + (i2-1)*n1 + i1, c)
                end do
            end if
        case (2)
            allocate(Xdir(n2, n1*n3*nval))
            if (rational) then
                do concurrent (c = 1:ncoord, i3 = 1:n3, i1 = 1:n1, i2 = 1:n2)
                    Xdir(i2, i1 + (i3-1)*n1 + (c-1)*n1*n3) = &
                        Xc((i3-1)*n1*n2 + (i2-1)*n1 + i1, c)*&
                        scale(Wc((i3-1)*n1*n2 + (i2-1)*n1 + i1), -weight_exponent_)
                end do
                do concurrent (i3 = 1:n3, i1 = 1:n1, i2 = 1:n2)
                    Xdir(i2, i1 + (i3-1)*n1 + ncoord*n1*n3) = &
                        scale(Wc((i3-1)*n1*n2 + (i2-1)*n1 + i1), -weight_exponent_)
                end do
            else
                do concurrent (c = 1:ncoord, i3 = 1:n3, i1 = 1:n1, i2 = 1:n2)
                    Xdir(i2, i1 + (i3-1)*n1 + (c-1)*n1*n3) = Xc((i3-1)*n1*n2 + (i2-1)*n1 + i1, c)
                end do
            end if
        case (3)
            allocate(Xdir(n3, n1*n2*nval))
            if (rational) then
                do concurrent (c = 1:ncoord, i2 = 1:n2, i1 = 1:n1, i3 = 1:n3)
                    Xdir(i3, i1 + (i2-1)*n1 + (c-1)*n1*n2) = &
                        Xc((i3-1)*n1*n2 + (i2-1)*n1 + i1, c)*&
                        scale(Wc((i3-1)*n1*n2 + (i2-1)*n1 + i1), -weight_exponent_)
                end do
                do concurrent (i2 = 1:n2, i1 = 1:n1, i3 = 1:n3)
                    Xdir(i3, i1 + (i2-1)*n1 + ncoord*n1*n2) = &
                        scale(Wc((i3-1)*n1*n2 + (i2-1)*n1 + i1), -weight_exponent_)
                end do
            else
                do concurrent (c = 1:ncoord, i2 = 1:n2, i1 = 1:n1, i3 = 1:n3)
                    Xdir(i3, i1 + (i2-1)*n1 + (c-1)*n1*n2) = Xc((i3-1)*n1*n2 + (i2-1)*n1 + i1, c)
                end do
            end if
        case default
            return
        end select
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Unpack refined volume control lines into the flattened control lattice.
    !!
    !! This reverses [[pack_refine_volume]]. Rational rows are dehomogenized by
    !! their normalized stored weights, then `weight_exponent` restores the
    !! original projective gauge. The caller must supply `dir` in `1:3` and
    !! nonzero homogeneous weights.
    pure subroutine unpack_refine_volume(Xdir, nc, dir, ncoord, rational, Xc, Wc, weight_exponent)
        real(rk), intent(in), contiguous :: Xdir(:,:)
            !! Refined directional control lines.
        integer, intent(in), contiguous :: nc(:)
            !! Refined directional control counts `[n1,n2,n3]`.
        integer, intent(in) :: dir
            !! Direction represented by the first dimension of `Xdir`.
        integer, intent(in) :: ncoord
            !! Number of physical coordinates per control point.
        logical, intent(in) :: rational
            !! Dehomogenize the final packed component when true.
        real(rk), allocatable, intent(out) :: Xc(:,:)
            !! Direction-1-fastest control lattice `[product(nc),ncoord]`.
        real(rk), allocatable, intent(out), optional :: Wc(:)
            !! Optional recovered flattened weights.
        integer, intent(in), optional :: weight_exponent
            !! Radix exponent restored to rational weights.
        integer :: n1, n2, n3, i1, i2, i3, c, weight_exponent_

        n1 = nc(1)
        n2 = nc(2)
        n3 = nc(3)
        weight_exponent_ = 0
        if (present(weight_exponent)) weight_exponent_ = weight_exponent
        allocate(Xc(n1*n2*n3, ncoord))
        if (rational .and. present(Wc)) allocate(Wc(n1*n2*n3))

        select case (dir)
        case (1)
            if (rational) then
                do concurrent (c = 1:ncoord, i3 = 1:n3, i2 = 1:n2, i1 = 1:n1)
                    Xc((i3-1)*n1*n2 + (i2-1)*n1 + i1, c) = &
                        Xdir(i1, i2 + (i3-1)*n2 + (c-1)*n2*n3)/Xdir(i1, i2 + (i3-1)*n2 + ncoord*n2*n3)
                end do
                if (present(Wc)) then
                    do concurrent (i3 = 1:n3, i2 = 1:n2, i1 = 1:n1)
                        Wc((i3-1)*n1*n2 + (i2-1)*n1 + i1) = &
                            scale(Xdir(i1, i2 + (i3-1)*n2 + ncoord*n2*n3), weight_exponent_)
                    end do
                end if
            else
                do concurrent (c = 1:ncoord, i3 = 1:n3, i2 = 1:n2, i1 = 1:n1)
                    Xc((i3-1)*n1*n2 + (i2-1)*n1 + i1, c) = Xdir(i1, i2 + (i3-1)*n2 + (c-1)*n2*n3)
                end do
            end if
        case (2)
            if (rational) then
                do concurrent (c = 1:ncoord, i3 = 1:n3, i1 = 1:n1, i2 = 1:n2)
                    Xc((i3-1)*n1*n2 + (i2-1)*n1 + i1, c) = &
                        Xdir(i2, i1 + (i3-1)*n1 + (c-1)*n1*n3)/Xdir(i2, i1 + (i3-1)*n1 + ncoord*n1*n3)
                end do
                if (present(Wc)) then
                    do concurrent (i3 = 1:n3, i1 = 1:n1, i2 = 1:n2)
                        Wc((i3-1)*n1*n2 + (i2-1)*n1 + i1) = &
                            scale(Xdir(i2, i1 + (i3-1)*n1 + ncoord*n1*n3), weight_exponent_)
                    end do
                end if
            else
                do concurrent (c = 1:ncoord, i3 = 1:n3, i1 = 1:n1, i2 = 1:n2)
                    Xc((i3-1)*n1*n2 + (i2-1)*n1 + i1, c) = Xdir(i2, i1 + (i3-1)*n1 + (c-1)*n1*n3)
                end do
            end if
        case (3)
            if (rational) then
                do concurrent (c = 1:ncoord, i2 = 1:n2, i1 = 1:n1, i3 = 1:n3)
                    Xc((i3-1)*n1*n2 + (i2-1)*n1 + i1, c) = &
                        Xdir(i3, i1 + (i2-1)*n1 + (c-1)*n1*n2)/Xdir(i3, i1 + (i2-1)*n1 + ncoord*n1*n2)
                end do
                if (present(Wc)) then
                    do concurrent (i2 = 1:n2, i1 = 1:n1, i3 = 1:n3)
                        Wc((i3-1)*n1*n2 + (i2-1)*n1 + i1) = &
                            scale(Xdir(i3, i1 + (i2-1)*n1 + ncoord*n1*n2), weight_exponent_)
                    end do
                end if
            else
                do concurrent (c = 1:ncoord, i2 = 1:n2, i1 = 1:n1, i3 = 1:n3)
                    Xc((i3-1)*n1*n2 + (i2-1)*n1 + i1, c) = Xdir(i3, i1 + (i2-1)*n1 + (c-1)*n1*n2)
                end do
            end if
        case default
            return
        end select
    end subroutine
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
    !! At least one finite coordinate per control point is a caller
    !! precondition; this overload validates the flattened row count but does
    !! not diagnose a zero coordinate dimension or nonfinite coordinates. Knot
    !! and positive-weight invariants are validated before state is replaced.
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

        this%knot1 = knot1
        this%knot2 = knot2
        this%knot3 = knot3
        this%degree = degree_
        this%nc = nc_
        this%Xc = Xc
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
    !> Construct all spline directions from breakpoints and continuities.
    !! Directional continuity arrays are converted to knot multiplicities by
    !! [[compute_knot_vector]]. Optional controls and weights must match the
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
            if (size(Xc,1) /= ncp) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume", &
                    message    = "Control points size mismatch in set2",&
                    location   = "set2", &
                    suggestion = "size(Xc,1) must equal nc(1)*nc(2)*nc(3).")
                return
            end if
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

        this%knot1 = knot1_
        this%knot2 = knot2_
        this%knot3 = knot3_
        this%degree = degree_
        this%nc = nc_
        this%wrap_parameters = .false.
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
        integer :: nc_(3), ncp

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

        if (present(Wc)) then
            if (size(Wc) /= ncp) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Weights length mismatch: size(Wc) must equal number of control points.",&
                    location   = "set3",&
                    suggestion = "Provide Wc with size(Wc) == nc(1)*nc(2)*nc(3).")
                return
            else if (.not. valid_weights(Wc)) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Invalid rational weights: all weights must be finite and positive.",&
                    location   = "set3",&
                    suggestion = "Use strictly positive finite Wc values for rational NURBS geometry.")
                return
            end if
        end if

        this%Xc = Xc
        this%nc = nc_

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
                deallocate(this%knot2)
                allocate(this%knot2(2*this%nc(2)))
            end if
        else
            allocate(this%knot2(2*this%nc(2)))
        end if
        this%knot2(1:this%nc(2)) = 0.0_rk
        this%knot2(this%nc(2)+1:2*this%nc(2)) = 1.0_rk

        if (allocated(this%knot3)) then
            if (size(this%knot3) /= 2*this%nc(3)) then
                deallocate(this%knot3)
                allocate(this%knot3(2*this%nc(3)))
            end if
        else
            allocate(this%knot3(2*this%nc(3)))
        end if
        this%knot3(1:this%nc(3)) = 0.0_rk
        this%knot3(this%nc(3)+1:2*this%nc(3)) = 1.0_rk

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
        integer :: m(3), i

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
                suggestion = "Require degree >= 0, nc > degree per direction, size(Xc,1) == product(nc), and at least one coordinate dimension.")
            return
        end if

        if (present(Wc)) then
            if (size(Wc) /= nc(1)*nc(2)*nc(3)) then
                call this%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Weights length mismatch: size(Wc) must equal number of control points.",&
                    location   = "set4",&
                    suggestion = "Provide Wc with size(Wc) == nc(1)*nc(2)*nc(3).")
                return
            else if (.not. valid_weights(Wc)) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Invalid rational weights: all weights must be finite and positive.",&
                    location   = "set4",&
                    suggestion = "Use strictly positive finite Wc values for rational NURBS geometry.")
                return
            end if
        end if

        if (allocated(this%Xc)) then
            if (size(this%Xc,1) /= size(Xc,1) .or. size(this%Xc,2) /= size(Xc,2)) deallocate(this%Xc)
        end if

        this%Xc = Xc
        this%nc = nc(1:3)
        this%degree = degree(1:3)

        ! Size of knot vectors
        m = nc(1:3) + degree(1:3) + 1

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
                deallocate(this%knot2)
                allocate(this%knot2(m(2)))
            end if
        else
            allocate(this%knot2(m(2)))
        end if
        this%knot2(1:degree(2)+1) = 0.0_rk
        this%knot2(degree(2)+2:m(2)-degree(2)-1) = [(real(i, rk)/(m(2)-2*degree(2)-1), i=1, m(2)-2*degree(2)-2)]
        this%knot2(m(2)-degree(2):m(2)) = 1.0_rk

        if (allocated(this%knot3)) then
            if (size(this%knot3) /= m(3)) then
                deallocate(this%knot3)
                allocate(this%knot3(m(3)))
            end if
        else
            allocate(this%knot3(m(3)))
        end if
        this%knot3(1:degree(3)+1) = 0.0_rk
        this%knot3(degree(3)+2:m(3)-degree(3)-1) = [(real(i, rk)/(m(3)-2*degree(3)-1), i=1, m(3)-2*degree(3)-2)]
        this%knot3(m(3)-degree(3):m(3)) = 1.0_rk
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
    !> Sample and cache volume geometry.
    !!
    !! Explicit `Xt1:3` replace directional parameter vectors; corresponding
    !! `res1:3` generate uniform vectors when explicit values are absent.
    !! Existing vectors are reused when neither is provided. Unless arbitrary
    !! rows `Xt(:,1:3)` are supplied, `ndgrid` forms a direction-1-fastest
    !! Cartesian product.
    pure subroutine create(this, res1, res2, res3, Xt1, Xt2, Xt3, Xt)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: res1
            !! Optional positive generated count in direction 1.
        integer, intent(in), optional :: res2
            !! Optional positive generated count in direction 2.
        integer, intent(in), optional :: res3
            !! Optional positive generated count in direction 3.
        real(rk), intent(in), contiguous, optional :: Xt1(:)
            !! Optional explicit direction-1 parameters.
        real(rk), intent(in), contiguous, optional :: Xt2(:)
            !! Optional explicit direction-2 parameters.
        real(rk), intent(in), contiguous, optional :: Xt3(:)
            !! Optional explicit direction-3 parameters.
        real(rk), intent(in), contiguous, optional :: Xt(:,:)
            !! Optional arbitrary parameter points `[npoint,at least 3]`.
        real(rk), allocatable :: Xg_new(:,:)

        if (.not. this%err%ok) return

        ! check
        if (.not.allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Control points are not set.",&
                location   = "create",&
                suggestion = "Call set(...) first before create().")
            return
        end if

        if (.not.allocated(this%knot1) .or. .not.allocated(this%knot2) .or. .not.allocated(this%knot3)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Knot vector is not set.",&
                location   = "create",&
                suggestion = "Call set(...) first before create().")
            return
        end if

        ! Set parameter values
        if (present(Xt1)) then
            if (size(Xt1) < 1) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Parameter values in direction 1 are empty.",&
                    location   = "create",&
                    suggestion = "Pass nonempty Xt1 or use res1 >= 1.")
                return
            end if
            if (allocated(this%Xt1)) then
                if (size(this%Xt1) /= size(Xt1)) deallocate(this%Xt1)
            end if
            this%Xt1 = Xt1
        else if (present(res1)) then
            if (res1 < 1) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Resolution in direction 1 must be at least one.",&
                    location   = "create",&
                    suggestion = "Use res1 >= 1 or pass explicit nonempty Xt1 values.")
                return
            end if
            if (allocated(this%Xt1)) then
                if (size(this%Xt1) /= res1) then
                    deallocate(this%Xt1)
                    allocate(this%Xt1(res1))
                end if
            else
                allocate(this%Xt1(res1))
            end if
            call fill_uniform(knot_start(this%knot1, this%nc(1), this%degree(1)), &
                knot_end(this%knot1, this%nc(1), this%degree(1)), this%Xt1)
        end if

        ! Set parameter values
        if (present(Xt2)) then
            if (size(Xt2) < 1) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Parameter values in direction 2 are empty.",&
                    location   = "create",&
                    suggestion = "Pass nonempty Xt2 or use res2 >= 1.")
                return
            end if
            if (allocated(this%Xt2)) then
                if (size(this%Xt2) /= size(Xt2)) deallocate(this%Xt2)
            end if
            this%Xt2 = Xt2
        else if (present(res2)) then
            if (res2 < 1) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Resolution in direction 2 must be at least one.",&
                    location   = "create",&
                    suggestion = "Use res2 >= 1 or pass explicit nonempty Xt2 values.")
                return
            end if
            if (allocated(this%Xt2)) then
                if (size(this%Xt2) /= res2) then
                    deallocate(this%Xt2)
                    allocate(this%Xt2(res2))
                end if
            else
                allocate(this%Xt2(res2))
            end if
            call fill_uniform(knot_start(this%knot2, this%nc(2), this%degree(2)), &
                knot_end(this%knot2, this%nc(2), this%degree(2)), this%Xt2)
        end if

        ! Set parameter values
        if (present(Xt3)) then
            if (size(Xt3) < 1) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Parameter values in direction 3 are empty.",&
                    location   = "create",&
                    suggestion = "Pass nonempty Xt3 or use res3 >= 1.")
                return
            end if
            if (allocated(this%Xt3)) then
                if (size(this%Xt3) /= size(Xt3)) deallocate(this%Xt3)
            end if
            this%Xt3 = Xt3
        else if (present(res3)) then
            if (res3 < 1) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Resolution in direction 3 must be at least one.",&
                    location   = "create",&
                    suggestion = "Use res3 >= 1 or pass explicit nonempty Xt3 values.")
                return
            end if
            if (allocated(this%Xt3)) then
                if (size(this%Xt3) /= res3) then
                    deallocate(this%Xt3)
                    allocate(this%Xt3(res3))
                end if
            else
                allocate(this%Xt3(res3))
            end if
            call fill_uniform(knot_start(this%knot3, this%nc(3), this%degree(3)), &
                knot_end(this%knot3, this%nc(3), this%degree(3)), this%Xt3)
        end if

        if (present(Xt)) then
            if (size(Xt,1) < 1 .or. size(Xt,2) < 3) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Parameter point array has invalid shape.",&
                    location   = "create",&
                    suggestion = "Pass Xt with at least one row and three parameter columns.")
                return
            end if
            this%Xt = Xt
            if (allocated(this%Xt1) .and. allocated(this%Xt2) .and. allocated(this%Xt3)) then
                if (size(this%Xt,1) == size(this%Xt1)*size(this%Xt2)*size(this%Xt3)) then
                    this%ng(1) = size(this%Xt1,1)
                    this%ng(2) = size(this%Xt2,1)
                    this%ng(3) = size(this%Xt3,1)
                else
                    this%ng = [size(this%Xt,1), 1, 1]
                end if
            else
                this%ng = [size(this%Xt,1), 1, 1]
            end if
        else

            if (.not. allocated(this%Xt1) .or. .not. allocated(this%Xt2) .or. .not. allocated(this%Xt3)) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Parameter values are not set.",&
                    location   = "create",&
                    suggestion = "Pass Xt, or pass Xt1/Xt2/Xt3 or res1/res2/res3.")
                return
            end if
            ! Set number of geometry points
            this%ng(1) = size(this%Xt1,1)
            this%ng(2) = size(this%Xt2,1)
            this%ng(3) = size(this%Xt3,1)

            call ndgrid(this%Xt1, this%Xt2, this%Xt3, this%Xt)
        end if
        if (this%is_rational()) then ! NURBS
            Xg_new = compute_Xg(this%Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                this%Xc, this%Wc, this%wrap_parameters)
        else ! B-Spline
            Xg_new = compute_Xg(this%Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                this%Xc, this%wrap_parameters)
        end if
        call move_alloc(Xg_new, this%Xg)
        this%elemConn_Xg_vis = this%cmp_elem_Xg_vis()
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate one physical volume point.
    pure function cmp_Xg(this, Xt) result(Xg)
        class(nurbs_volume), intent(in) :: this
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter triple; at least three entries are required.
        real(rk), allocatable :: Xg(:)

        if (.not. this%err%ok .or. .not. allocated(this%Xc) .or. &
            .not. allocated(this%knot1) .or. .not. allocated(this%knot2) .or. .not. allocated(this%knot3)) then
            allocate(Xg(0))
            return
        end if

        if (this%is_rational()) then ! NURBS
            Xg = compute_Xg(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                this%Xc, this%Wc, this%wrap_parameters)
        else ! B-Spline
            Xg = compute_Xg(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%Xc, this%wrap_parameters)
        end if
    end function
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

        if (.not. this%err%ok) then
            allocate(Xt(0))
            return
        end if

        if (dir == 1) then
            if (allocated(this%Xt1)) then
                Xt = this%Xt1
            else
                allocate(Xt(0))
            end if
        else if (dir == 2) then
            if (allocated(this%Xt2)) then
                Xt = this%Xt2
            else
                allocate(Xt(0))
            end if
        else if (dir == 3) then
            if (allocated(this%Xt3)) then
                Xt = this%Xt3
            else
                allocate(Xt(0))
            end if
        else
            allocate(Xt(0))
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the three directional sample counts.
    !!
    !! A Cartesian cache reports its three directional counts. A row-wise cache
    !! reports `[npoint,1,1]`. The result is zero before sampling and while a
    !! diagnostic is active.
    pure function get_ng(this) result(ng)
        class(nurbs_volume), intent(in) :: this
        integer :: ng(3)

        ng = 0
        if (.not. this%err%ok) return

        ng = this%ng
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Infer one or all polynomial degrees from stored spline dimensions.
    !!
    !! When a control count is known, direction \(d\) uses
    !! \(p_d=m_d-n_d-1\), where \(m_d\) is the knot-vector size and \(n_d\)
    !! the control count. Otherwise the leading multiplicity is interpreted as
    !! \(p_d+1\). Optional `dir` restricts the update to one direction.
    pure subroutine cmp_degree(this, dir)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: dir
        integer, allocatable :: m1(:), m2(:), m3(:)

        if (.not. this%err%ok) return

        if (present(dir)) then
            if (dir == 1) then
                if (this%nc(1) > 0) then
                    this%degree(1) = size(this%knot1) - this%nc(1) - 1
                else
                    m1 = this%get_multiplicity(1)
                    this%degree(1) = m1(1) - 1
                end if
            else if (dir == 2) then
                if (this%nc(2) > 0) then
                    this%degree(2) = size(this%knot2) - this%nc(2) - 1
                else
                    m2 = this%get_multiplicity(2)
                    this%degree(2) = m2(1) - 1
                end if
            else if (dir == 3) then
                if (this%nc(3) > 0) then
                    this%degree(3) = size(this%knot3) - this%nc(3) - 1
                else
                    m3 = this%get_multiplicity(3)
                    this%degree(3) = m3(1) - 1
                end if
            else
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Invalid direction for degree.",&
                    location   = "cmp_degree",&
                    suggestion = "Check the direction argument.")
                return
            end if
        else
            if (this%nc(1) > 0) then
                this%degree(1) = size(this%knot1) - this%nc(1) - 1
            else
                m1 = this%get_multiplicity(1)
                this%degree(1) = m1(1) - 1
            end if

            if (this%nc(2) > 0) then
                this%degree(2) = size(this%knot2) - this%nc(2) - 1
            else
                m2 = this%get_multiplicity(2)
                this%degree(2) = m2(1) - 1
            end if

            if (this%nc(3) > 0) then
                this%degree(3) = size(this%knot3) - this%nc(3) - 1
            else
                m3 = this%get_multiplicity(3)
                this%degree(3) = m3(1) - 1
            end if
        end if

    end subroutine
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

        degree = 0
        if (.not. this%err%ok) return

        if (dir == 1) then
            degree = this%degree(1)
        else if (dir == 2) then
            degree = this%degree(2)
        else if (dir == 3) then
            degree = this%degree(3)
        end if
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

        if (.not. this%err%ok) then
            allocate(knot(0))
            return
        end if

        if (dir == 1) then
            if (allocated(this%knot1)) then
                knot = this%knot1
            else
                allocate(knot(0))
            end if
        else if (dir == 2) then
            if (allocated(this%knot2)) then
                knot = this%knot2
            else
                allocate(knot(0))
            end if
        else if (dir == 3) then
            if (allocated(this%knot3)) then
                knot = this%knot3
            else
                allocate(knot(0))
            end if
        else
            allocate(knot(0))
        end if

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

        if (.not. this%err%ok) return

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

        integer :: a, b, t, s, o, r, res_min, ncoord, i, j, k, ne_total, np, idx, m, N1sp, N2sp, N3sp, L, N, res1, res2, res3
        integer :: ng_eval(3)

        real(rk), allocatable :: U1(:), U2(:), U3(:)        ! unique knots
        real(rk), allocatable :: U1r(:), U2r(:), U3r(:)     ! refined per dir (length N)
        real(rk), allocatable :: Xt_all(:,:), Xg_all(:,:)   ! batched params & geometry
        integer,  allocatable :: elemConn(:,:)              ! [ne_total, N]

        if (.not. this%err%ok) return

        if (.not. allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Control points are not set.",&
                location   = "export_Xth_in_Xg",&
                suggestion = "Call set(...) first before exporting.")
            return
        end if

        if (.not. allocated(this%knot1) .or. .not. allocated(this%knot2) .or. .not. allocated(this%knot3)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Knot vector is not set.",&
                location   = "export_Xth_in_Xg",&
                suggestion = "Call set(...) first before exporting.")
            return
        end if

        res_min = 10
        if (present(res)) res_min = max(2, res)

        U1 = active_knots(this%knot1, this%nc(1), this%degree(1))
        if (size(U1) < 2) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "knot1 needs >= 2 unique values.",&
                location   = "export_Xth_in_Xg",&
                suggestion = "Check the knot vector for sufficient unique values.")
            return
        end if
        U2 = active_knots(this%knot2, this%nc(2), this%degree(2))
        if (size(U2) < 2) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "knot2 needs >= 2 unique values.",&
                location   = "export_Xth_in_Xg",&
                suggestion = "Check the knot vector for sufficient unique values.")
            return
        end if
        U3 = active_knots(this%knot3, this%nc(3), this%degree(3))
        if (size(U3) < 2) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "knot3 needs >= 2 unique values.",&
                location   = "export_Xth_in_Xg",&
                suggestion = "Check the knot vector for sufficient unique values.")
            return
        end if

        ! spans per direction
        N1sp = size(U1)-1
        N2sp = size(U2)-1
        N3sp = size(U3)-1

        L = N1sp
        a = L
        b = N2sp
        do
            t = mod(a,b)
            if (t==0) exit
            a = b
            b = t
        end do
        L = (L / b) * N2sp

        a = L
        b = N3sp
        do
            t = mod(a,b)
            if (t==0) exit
            a = b
            b = t
        end do
        L = (L / b) * N3sp

        L = L * max(1, res_min-1)

        N = L + 1
        res1 = L / N1sp + 1
        res2 = L / N2sp + 1
        res3 = L / N3sp + 1

        ncoord = size(this%Xc,2)
        if (ncoord < 2 .or. ncoord > 3) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid geometry dimension.",&
                location   = "export_Xth_in_Xg",&
                suggestion = "Check the geometry dimension before exporting the NURBS volume.")
            return
        end if

        ! Allocate refined knot vectors
        allocate(U1r( (size(U1)-1)*(res1-1) + 1 ))
        allocate(U2r( (size(U2)-1)*(res2-1) + 1 ))
        allocate(U3r( (size(U3)-1)*(res3-1) + 1 ))

        do concurrent (s = 1:size(U1)-1, r = 1:res1-1) local(o)
            o = (s-1)*(res1-1)
            U1r(o+r) = U1(s) + (U1(s+1)-U1(s)) * real(r-1, rk) / real(res1-1, rk)
        end do
        U1r(size(U1r)) = U1(size(U1))
        do concurrent (s = 1:size(U2)-1, r = 1:res2-1) local(o)
            o = (s-1)*(res2-1)
            U2r(o+r) = U2(s) + (U2(s+1)-U2(s)) * real(r-1, rk) / real(res2-1, rk)
        end do
        U2r(size(U2r)) = U2(size(U2))
        do concurrent (s = 1:size(U3)-1, r = 1:res3-1) local(o)
            o = (s-1)*(res3-1)
            U3r(o+r) = U3(s) + (U3(s+1)-U3(s)) * real(r-1, rk) / real(res3-1, rk)
        end do
        U3r(size(U3r)) = U3(size(U3))

        if (size(U1r) /= N .or. size(U2r) /= N .or. size(U3r) /= N) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Refinement size mismatch.",&
                location   = "export_Xth_in_Xg",&
                suggestion = "Check the refinement process for consistency.")
            return
        end if

        ! total element count and node count
        ne_total = size(U2)*size(U3) + size(U1)*size(U3) + size(U1)*size(U2)
        np = ne_total * N

        ! Allocate global arrays
        allocate(Xt_all(np,3), elemConn(ne_total, N))

        ! dir-1: u varies (v=U2(j), w=U3(k))
        do concurrent (k = 1:size(U3), j = 1:size(U2), m = 1:N) local(idx)
            idx = ((k-1)*size(U2) + j - 1)*N + m
            Xt_all(idx,1) = U1r(m)
            Xt_all(idx,2) = U2(j)
            Xt_all(idx,3) = U3(k)
        end do

        ! dir-2: v varies (u=U1(i), w=U3(k))
        do concurrent (k = 1:size(U3), i = 1:size(U1), m = 1:N) local(idx)
            idx = (size(U2)*size(U3) + (k-1)*size(U1) + i - 1)*N + m
            Xt_all(idx,1) = U1(i)
            Xt_all(idx,2) = U2r(m)
            Xt_all(idx,3) = U3(k)
        end do

        ! dir-3: w varies (u=U1(i), v=U2(j))
        do concurrent (j = 1:size(U2), i = 1:size(U1), m = 1:N) local(idx)
            idx = (size(U2)*size(U3) + size(U1)*size(U3) + (j-1)*size(U1) + i - 1)*N + m
            Xt_all(idx,1) = U1(i)
            Xt_all(idx,2) = U2(j)
            Xt_all(idx,3) = U3r(m)
        end do
        ! compute global points
        ng_eval(1) = np
        ng_eval(2) = 1
        ng_eval(3) = 1
        if (this%is_rational()) then
            Xg_all = compute_Xg(Xt_all, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                ng_eval, this%Xc, this%Wc, this%wrap_parameters)
        else
            Xg_all = compute_Xg(Xt_all, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                ng_eval, this%Xc, this%wrap_parameters)
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
    !> Replace one coordinate of one flattened control point.
    !!
    !! `num` is the one-based direction-1-fastest control index and
    !! `dir` selects the physical coordinate. `X` is expected to be finite but
    !! is not checked. Cached geometry is not
    !! recomputed; call [[nurbs_volume:create]] when samples must follow the
    !! modified control lattice.
    pure subroutine modify_Xc(this,X,num,dir)
        class(nurbs_volume), intent(inout) :: this
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
                    category   = "forcad_nurbs_volume",&
                    message    = "Invalid control point index or direction.",&
                    location   = "modify_Xc",&
                    suggestion = "Use 1 <= num <= product(nc) and 1 <= dir <= spatial dimension.")
                return
            end if
            this%Xc(num,dir) = X
        else
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
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
    !> Replace one finite, strictly positive flattened control weight.
    !!
    !! The rational-state flag is recomputed from the stored weights. Cached
    !! geometry is unchanged; call [[nurbs_volume:create]] to resample it.
    pure subroutine modify_Wc(this,W,num)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: W
        integer, intent(in) :: num

        if (.not. this%err%ok) return

        if (allocated(this%Wc)) then
            if (num < lbound(this%Wc,1) .or. num > ubound(this%Wc,1)) then
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Invalid weight index.",&
                    location   = "modify_Wc",&
                    suggestion = "Use 1 <= num <= number of control points.")
                return
            end if
            if (.not. ieee_is_finite(W) .or. W <= 0.0_rk) then
                call this%err%set(&
                    code       = 105,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
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
                category   = "forcad_nurbs_volume",&
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
    !> Return multiplicities of the distinct knots in one direction.
    !!
    !! Result entries correspond, in nondecreasing order, to the unique values
    !! in the complete selected knot vector. Invalid input or unavailable state
    !! produces an allocated empty result.
    pure function get_multiplicity(this, dir) result(m)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        integer, allocatable :: m(:)

        if (.not. this%err%ok) then
            allocate(m(0))
            return
        end if

        if (dir == 1) then

            if (.not.allocated(this%knot1)) then
                allocate(m(0))
            else
                m = compute_multiplicity(this%knot1)
            end if

        else if (dir == 2) then

            if (.not.allocated(this%knot2)) then
                allocate(m(0))
            else
                m = compute_multiplicity(this%knot2)
            end if

        else if (dir == 3) then

            if (.not.allocated(this%knot3)) then
                allocate(m(0))
            else
                m = compute_multiplicity(this%knot3)
            end if

        else
            allocate(m(0))
        end if

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

        if (.not. this%err%ok) then
            allocate(c(0))
            return
        end if

        if (dir == 1) then

            if (.not.allocated(this%knot1)) then
                allocate(c(0))
            else
                c = this%degree(1) - compute_multiplicity(this%knot1)
            end if

        else if (dir == 2) then

            if (.not.allocated(this%knot2)) then
                allocate(c(0))
            else
                c = this%degree(2) - compute_multiplicity(this%knot2)
            end if

        else if (dir == 3) then

            if (.not.allocated(this%knot3)) then
                allocate(c(0))
            else
                c = this%degree(3) - compute_multiplicity(this%knot3)
            end if

        else
            allocate(c(0))
        end if

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

        if (.not. this%err%ok) return

        if (present(dir)) then

            if (dir == 1) then

                ! check
                if (.not.allocated(this%knot1)) then
                    call this%err%set(&
                        code       = 103,&
                        severity   = 1,&
                        category   = "forcad_nurbs_volume",&
                        message    = "Knot vector is not set.",&
                        location   = "cmp_nc",&
                        suggestion = "Call set(...) first before computing nc." )
                    return
                else
                    this%nc(1) = size(this%knot1) - this%degree(1) - 1
                end if

            else if (dir == 2) then

                ! check
                if (.not.allocated(this%knot2)) then
                    call this%err%set(&
                        code       = 103,&
                        severity   = 1,&
                        category   = "forcad_nurbs_volume",&
                        message    = "Knot vector is not set.",&
                        location   = "cmp_nc",&
                        suggestion = "Call set(...) first before computing nc." )
                    return
                else
                    this%nc(2) = size(this%knot2) - this%degree(2) - 1
                end if

            else if (dir == 3) then

                ! check
                if (.not.allocated(this%knot3)) then
                    call this%err%set(&
                        code       = 103,&
                        severity   = 1,&
                        category   = "forcad_nurbs_volume",&
                        message    = "Knot vector is not set.",&
                        location   = "cmp_nc",&
                        suggestion = "Call set(...) first before computing nc." )
                    return
                else
                    this%nc(3) = size(this%knot3) - this%degree(3) - 1
                end if

            else
                call this%err%set(&
                    code       = 103,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Invalid direction for computing number of control points.",&
                    location   = "cmp_nc",&
                    suggestion = "Use dir=1 or dir=2 to specify the direction.")
                return
            end if

        else

            ! check
            if (.not.allocated(this%knot1)) then
                call this%err%set(&
                    code       = 103,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Knot vector is not set.",&
                    location   = "cmp_nc",&
                    suggestion = "Call set(...) first before computing nc.")
                return
            else
                this%nc(1) = size(this%knot1) - this%degree(1) - 1
            end if

            ! check
            if (.not.allocated(this%knot2)) then
                call this%err%set(&
                    code       = 103,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Knot vector is not set.",&
                    location   = "cmp_nc",&
                    suggestion = "Call set(...) first before computing nc.")
                return
            else
                this%nc(2) = size(this%knot2) - this%degree(2) - 1
            end if

            ! check
            if (.not.allocated(this%knot3)) then
                call this%err%set(&
                    code       = 103,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Knot vector is not set.",&
                    location   = "cmp_nc",&
                    suggestion = "Call set(...) first before computing nc.")
                return
            else
                this%nc(3) = size(this%knot3) - this%degree(3) - 1
            end if

        end if

    end subroutine
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
    pure subroutine check_derivative_state(this, location, valid)
        class(nurbs_volume), intent(inout) :: this
        character(len=*), intent(in) :: location
        logical, intent(out) :: valid

        valid = .false.
        if (.not. this%err%ok) return
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
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate one arbitrary mixed volume basis derivative on its local support.
    pure subroutine compute_derivative_order_3d_local(&
        Xt1,&
        Xt2,&
        Xt3,&
        knot1,&
        knot2,&
        knot3,&
        degree,&
        nc,&
        order,&
        wrap_parameters,&
        values,&
        first_active,&
        elem,&
        Wc)
        real(rk), intent(in) :: Xt1, Xt2, Xt3
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3), nc(3), order(3)
        logical, intent(in) :: wrap_parameters(3)
        real(rk), intent(out) :: values(:)
        integer, intent(out), contiguous, optional :: first_active(:)
        integer, intent(in), contiguous, optional :: elem(:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk) :: ders1(0:max(0,order(1)),0:degree(1))
        real(rk) :: ders2(0:max(0,order(2)),0:degree(2))
        real(rk) :: ders3(0:max(0,order(3)),0:degree(3))
        real(rk) :: local_values(0:(degree(1)+1)*(degree(2)+1)*(degree(3)+1)-1)
        integer :: first1, first2, first3, i, i1, i2, i3, l1, l2, l3, local, global, plane_index

        values = 0.0_rk
        call basis_bspline_der_all_active(&
            Xt     = map_parameter(Xt1, knot1, nc(1), degree(1), wrap_parameters(1)),&
            knot   = knot1,&
            nc     = nc(1),&
            degree = degree(1),&
            nder   = order(1),&
            first  = first1,&
            ders   = ders1)
        call basis_bspline_der_all_active(&
            Xt     = map_parameter(Xt2, knot2, nc(2), degree(2), wrap_parameters(2)),&
            knot   = knot2,&
            nc     = nc(2),&
            degree = degree(2),&
            nder   = order(2),&
            first  = first2,&
            ders   = ders2)
        call basis_bspline_der_all_active(&
            Xt     = map_parameter(Xt3, knot3, nc(3), degree(3), wrap_parameters(3)),&
            knot   = knot3,&
            nc     = nc(3),&
            degree = degree(3),&
            nder   = order(3),&
            first  = first3,&
            ders   = ders3)
        if (present(Wc)) then
            call tensor_basis_derivative_local(&
                first1 = first1,&
                nc1    = nc(1),&
                ders1  = ders1,&
                first2 = first2,&
                nc2    = nc(2),&
                ders2  = ders2,&
                first3 = first3,&
                nc3    = nc(3),&
                ders3  = ders3,&
                values = local_values,&
                Wc     = Wc)
        else
            call tensor_basis_derivative_local(&
                first1 = first1,&
                nc1    = nc(1),&
                ders1  = ders1,&
                first2 = first2,&
                nc2    = nc(2),&
                ders2  = ders2,&
                first3 = first3,&
                nc3    = nc(3),&
                ders3  = ders3,&
                values = local_values)
        end if

        if (present(first_active)) then
            if (size(first_active) >= 3) first_active(1:3) = [first1,first2,first3]
            do i = 0, min(size(local_values),size(values))-1
                values(i+1) = local_values(i)
            end do
            return
        end if

        if (present(elem)) then
            do i = 1, min(size(values),size(elem))
                global = elem(i)
                if (global < 1 .or. global > product(nc)) cycle
                i3 = (global-1)/(nc(1)*nc(2)) + 1
                plane_index = global - (i3-1)*nc(1)*nc(2)
                i2 = (plane_index-1)/nc(1) + 1
                i1 = plane_index - (i2-1)*nc(1)
                l1 = i1 - first1
                l2 = i2 - first2
                l3 = i3 - first3
                if (l1 >= 0 .and. l1 <= degree(1) .and. l2 >= 0 .and. l2 <= degree(2) .and. &
                    l3 >= 0 .and. l3 <= degree(3)) then
                    local = l1 + (degree(1)+1)*(l2 + (degree(2)+1)*l3)
                    values(i) = local_values(local)
                end if
            end do
        else
            do l3 = 0, degree(3)
                do l2 = 0, degree(2)
                    do l1 = 0, degree(1)
                        i1 = first1 + l1
                        i2 = first2 + l2
                        i3 = first3 + l3
                        if (i1 < 1 .or. i1 > nc(1) .or. i2 < 1 .or. i2 > nc(2) .or. &
                            i3 < 1 .or. i3 > nc(3)) cycle
                        local = l1 + (degree(1)+1)*(l2 + (degree(2)+1)*l3)
                        global = i1 + nc(1)*((i2-1) + nc(2)*(i3-1))
                        if (global <= size(values)) values(global) = local_values(local)
                    end do
                end do
            end do
        end if
    end subroutine
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
    pure subroutine derivative_order_active_vector(this, Xt1, Xt2, Xt3, order, first_active, Dgc)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt1(:), Xt2(:), Xt3(:)
            !! Finite directional parameter vectors defining a Cartesian grid.
        integer, intent(in) :: order(3)
            !! Nonnegative mixed-derivative orders `[du,dv,dw]`.
        integer, allocatable, intent(out) :: first_active(:,:)
            !! First active control-point indices `[3,npoint]`.
        real(rk), allocatable, intent(out) :: Dgc(:,:)
            !! Local derivatives `[product(degree+1),product(grid sizes)]`.
        integer :: point1, point2, point3
        logical :: valid

        call check_derivative_state(this, "derivative_order_active_vector", valid)
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
        if (this%is_rational()) then
#if defined(__NVCOMPILER)
            do point3 = 1, size(Xt3)
                do point2 = 1, size(Xt2)
                    do point1 = 1, size(Xt1)
#else
            do concurrent (point3 = 1:size(Xt3), point2 = 1:size(Xt2), point1 = 1:size(Xt1))
#endif
                call compute_derivative_order_3d_local(&
                    Xt1             = Xt1(point1),&
                    Xt2             = Xt2(point2),&
                    Xt3             = Xt3(point3),&
                    knot1           = this%knot1,&
                    knot2           = this%knot2,&
                    knot3           = this%knot3,&
                    degree          = this%degree,&
                    nc              = this%nc,&
                    order           = order,&
                    wrap_parameters = this%wrap_parameters,&
                    values          = Dgc(:,point1+size(Xt1)*((point2-1)+size(Xt2)*(point3-1))),&
                    first_active    = first_active(:,point1+size(Xt1)*((point2-1)+size(Xt2)*(point3-1))),&
                    Wc              = this%Wc)
#if defined(__NVCOMPILER)
                    end do
                end do
#endif
            end do
        else
#if defined(__NVCOMPILER)
            do point3 = 1, size(Xt3)
                do point2 = 1, size(Xt2)
                    do point1 = 1, size(Xt1)
#else
            do concurrent (point3 = 1:size(Xt3), point2 = 1:size(Xt2), point1 = 1:size(Xt1))
#endif
                call compute_derivative_order_3d_local(&
                    Xt1             = Xt1(point1),&
                    Xt2             = Xt2(point2),&
                    Xt3             = Xt3(point3),&
                    knot1           = this%knot1,&
                    knot2           = this%knot2,&
                    knot3           = this%knot3,&
                    degree          = this%degree,&
                    nc              = this%nc,&
                    order           = order,&
                    wrap_parameters = this%wrap_parameters,&
                    values          = Dgc(:,point1+size(Xt1)*((point2-1)+size(Xt2)*(point3-1))),&
                    first_active    = first_active(:,point1+size(Xt1)*((point2-1)+size(Xt2)*(point3-1))))
#if defined(__NVCOMPILER)
                    end do
                end do
#endif
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one active mixed volume derivative.
    !! Local and global control indices use the same direction-1-fastest
    !! formulas as the vector overload, with one parameter triple.
    pure subroutine derivative_order_active_scalar(this, Xt, order, first_active, Dgc)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:)
            !! Finite parameter triple `[u,v,w]`.
        integer, intent(in) :: order(3)
            !! Nonnegative mixed-derivative orders `[du,dv,dw]`.
        integer, intent(out) :: first_active(3)
            !! First active control-point index in each direction.
        real(rk), allocatable, intent(out) :: Dgc(:)
            !! Local derivatives `[product(degree+1)]`.
        logical :: valid

        first_active = 0
        call check_derivative_state(this, "derivative_order_active_scalar", valid)
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
        if (this%is_rational()) then
            call compute_derivative_order_3d_local(&
                Xt1             = Xt(1),&
                Xt2             = Xt(2),&
                Xt3             = Xt(3),&
                knot1           = this%knot1,&
                knot2           = this%knot2,&
                knot3           = this%knot3,&
                degree          = this%degree,&
                nc              = this%nc,&
                order           = order,&
                wrap_parameters = this%wrap_parameters,&
                values          = Dgc,&
                first_active    = first_active,&
                Wc              = this%Wc)
        else
            call compute_derivative_order_3d_local(&
                Xt1             = Xt(1),&
                Xt2             = Xt(2),&
                Xt3             = Xt(3),&
                knot1           = this%knot1,&
                knot2           = this%knot2,&
                knot3           = this%knot3,&
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
    !! Explicit directional vectors replace their caches; otherwise a positive
    !! resolution generates an inclusive active-interval grid, and an omitted
    !! direction reuses its cache. The method updates `ng` and directional
    !! parameter caches but not `Xt` or cached geometry `Xg`. Nonfinite explicit
    !! parameters are not diagnosed here and yield zero rows when no span is
    !! selected.
    pure subroutine derivative_order_vector(this, order, Dgc, res1, res2, res3, Xt1, Xt2, Xt3)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: order(3)
            !! Nonnegative mixed-derivative orders `[du,dv,dw]`.
        real(rk), allocatable, intent(out) :: Dgc(:,:)
            !! Dense derivative matrix `[product(ng),product(nc)]`.
        integer, intent(in), optional :: res1, res2, res3
            !! Optional positive generated counts by direction.
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:), Xt3(:)
            !! Optional explicit directional parameter vectors.
        integer :: point1, point2, point3
        logical :: valid

        call check_derivative_state(this, "derivative_order_vector", valid)
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

        if (present(Xt1)) then
            this%Xt1 = Xt1
        else if (present(res1)) then
            if (res1 < 1) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "First volume derivative resolution must be positive.",&
                    location   = "derivative_order_vector",&
                    suggestion = "Use res1 >= 1 or supply a nonempty Xt1 array.")
                allocate(Dgc(0,0))
                return
            end if
            if (.not. allocated(this%Xt1)) allocate(this%Xt1(res1))
            if (size(this%Xt1) /= res1) then
                deallocate(this%Xt1)
                allocate(this%Xt1(res1))
            end if
            call fill_uniform(&
                a = knot_start(this%knot1, this%nc(1), this%degree(1)),&
                b = knot_end(this%knot1, this%nc(1), this%degree(1)),&
                x = this%Xt1)
        else if (.not. allocated(this%Xt1)) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "No first-direction volume parameters are available.",&
                location   = "derivative_order_vector",&
                suggestion = "Supply Xt1 or res1 before requesting vector derivatives.")
            allocate(Dgc(0,0))
            return
        end if

        if (present(Xt2)) then
            this%Xt2 = Xt2
        else if (present(res2)) then
            if (res2 < 1) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Second volume derivative resolution must be positive.",&
                    location   = "derivative_order_vector",&
                    suggestion = "Use res2 >= 1 or supply a nonempty Xt2 array.")
                allocate(Dgc(0,0))
                return
            end if
            if (.not. allocated(this%Xt2)) allocate(this%Xt2(res2))
            if (size(this%Xt2) /= res2) then
                deallocate(this%Xt2)
                allocate(this%Xt2(res2))
            end if
            call fill_uniform(&
                a = knot_start(this%knot2, this%nc(2), this%degree(2)),&
                b = knot_end(this%knot2, this%nc(2), this%degree(2)),&
                x = this%Xt2)
        else if (.not. allocated(this%Xt2)) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "No second-direction volume parameters are available.",&
                location   = "derivative_order_vector",&
                suggestion = "Supply Xt2 or res2 before requesting vector derivatives.")
            allocate(Dgc(0,0))
            return
        end if

        if (present(Xt3)) then
            this%Xt3 = Xt3
        else if (present(res3)) then
            if (res3 < 1) then
                call this%err%set(&
                    code       = 107,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Third volume derivative resolution must be positive.",&
                    location   = "derivative_order_vector",&
                    suggestion = "Use res3 >= 1 or supply a nonempty Xt3 array.")
                allocate(Dgc(0,0))
                return
            end if
            if (.not. allocated(this%Xt3)) allocate(this%Xt3(res3))
            if (size(this%Xt3) /= res3) then
                deallocate(this%Xt3)
                allocate(this%Xt3(res3))
            end if
            call fill_uniform(&
                a = knot_start(this%knot3, this%nc(3), this%degree(3)),&
                b = knot_end(this%knot3, this%nc(3), this%degree(3)),&
                x = this%Xt3)
        else if (.not. allocated(this%Xt3)) then
            call this%err%set(&
                code       = 107,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "No third-direction volume parameters are available.",&
                location   = "derivative_order_vector",&
                suggestion = "Supply Xt3 or res3 before requesting vector derivatives.")
            allocate(Dgc(0,0))
            return
        end if

        this%ng = [size(this%Xt1),size(this%Xt2),size(this%Xt3)]
        allocate(Dgc(product(this%ng),product(this%nc)), source=0.0_rk)
        if (this%is_rational()) then
#if defined(__NVCOMPILER)
            ! NVFORTRAN 26.3 cannot lower the pure local kernel for stdpar GPU execution.
            do point3 = 1, this%ng(3)
                do point2 = 1, this%ng(2)
                    do point1 = 1, this%ng(1)
#else
            do concurrent (point3 = 1:this%ng(3), point2 = 1:this%ng(2), point1 = 1:this%ng(1))
#endif
                call compute_derivative_order_3d_local(&
                    Xt1    = this%Xt1(point1),&
                    Xt2    = this%Xt2(point2),&
                    Xt3    = this%Xt3(point3),&
                    knot1  = this%knot1,&
                    knot2  = this%knot2,&
                    knot3  = this%knot3,&
                    degree = this%degree,&
                    nc     = this%nc,&
                    order  = order,&
                    wrap_parameters = this%wrap_parameters,&
                    values = Dgc(point1+this%ng(1)*((point2-1)+this%ng(2)*(point3-1)),:),&
                    Wc     = this%Wc)
#if defined(__NVCOMPILER)
                    end do
                end do
#endif
            end do
        else
#if defined(__NVCOMPILER)
            do point3 = 1, this%ng(3)
                do point2 = 1, this%ng(2)
                    do point1 = 1, this%ng(1)
#else
            do concurrent (point3 = 1:this%ng(3), point2 = 1:this%ng(2), point1 = 1:this%ng(1))
#endif
                call compute_derivative_order_3d_local(&
                    Xt1    = this%Xt1(point1),&
                    Xt2    = this%Xt2(point2),&
                    Xt3    = this%Xt3(point3),&
                    knot1  = this%knot1,&
                    knot2  = this%knot2,&
                    knot3  = this%knot3,&
                    degree = this%degree,&
                    nc     = this%nc,&
                    order  = order,&
                    wrap_parameters = this%wrap_parameters,&
                    values = Dgc(point1+this%ng(1)*((point2-1)+this%ng(2)*(point3-1)),:))
#if defined(__NVCOMPILER)
                    end do
                end do
#endif
            end do
        end if
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
    pure subroutine derivative_order_scalar(this, Xt, order, Dgc, elem)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter triple `[u,v,w]`.
        integer, intent(in) :: order(3)
            !! Nonnegative mixed-derivative orders `[du,dv,dw]`.
        real(rk), allocatable, intent(out) :: Dgc(:)
            !! Dense `[product(nc)]` or selected `[size(elem)]` derivative vector.
        integer, intent(in), contiguous, optional :: elem(:)
            !! Optional ordered one-based flattened control indices.
        logical :: valid

        call check_derivative_state(this, "derivative_order_scalar", valid)
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
        if (this%is_rational()) then
            call compute_derivative_order_3d_local(&
                Xt1    = Xt(1),&
                Xt2    = Xt(2),&
                Xt3    = Xt(3),&
                knot1  = this%knot1,&
                knot2  = this%knot2,&
                knot3  = this%knot3,&
                degree = this%degree,&
                nc     = this%nc,&
                order  = order,&
                wrap_parameters = this%wrap_parameters,&
                values = Dgc,&
                elem   = elem,&
                Wc     = this%Wc)
        else
            call compute_derivative_order_3d_local(&
                Xt1    = Xt(1),&
                Xt2    = Xt(2),&
                Xt3    = Xt(3),&
                knot1  = this%knot1,&
                knot2  = this%knot2,&
                knot3  = this%knot3,&
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
    !> Evaluate first parametric derivatives on a tensor-product grid.
    !!
    !! Explicit directional parameter vectors take precedence over resolutions;
    !! omitted directions reuse their caches. Flattened rows follow
    !! direction-1-fastest `ndgrid` ordering. This convenience overload updates
    !! directional caches and `ng`, but not the flattened cache `Xt` or geometry
    !! cache `Xg`. It assumes a valid spline state and an explicit, generated,
    !! or allocated parameter vector in every direction; those preconditions
    !! are not diagnosed here.
    pure subroutine derivative_vector(this, res1, res2, res3, Xt1, Xt2, Xt3, dTgc, Tgc)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: res1, res2, res3
            !! Optional positive uniform-grid sizes in directions 1--3.
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:), Xt3(:)
            !! Optional directional parameter vectors.
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
            !! Parametric gradients `[ng_total,nc_total,3]`, ordered as \((u,v,w)\).
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
            !! Optional basis matrix `[ng_total,nc_total]`.
        real(rk), allocatable :: Tgc_(:,:)
        real(rk), allocatable :: Xt(:,:)

        if (.not. this%err%ok) return

        ! Set parameter values
        if (present(Xt1)) then
            if (allocated(this%Xt1)) then
                if (size(Xt1,1) /= size(this%Xt1,1)) deallocate(this%Xt1)
            end if
            this%Xt1 = Xt1
        else if (present(res1)) then
            if (allocated(this%Xt1)) then
                if (size(this%Xt1,1) /= res1) then
                    deallocate(this%Xt1)
                    allocate(this%Xt1(res1))
                end if
            else
                allocate(this%Xt1(res1))
            end if
            call fill_uniform(knot_start(this%knot1, this%nc(1), this%degree(1)), &
                knot_end(this%knot1, this%nc(1), this%degree(1)), this%Xt1)
        end if

        ! Set parameter values
        if (present(Xt2)) then
            if (allocated(this%Xt2)) then
                if (size(Xt2,1) /= size(this%Xt2,1)) deallocate(this%Xt2)
            end if
            this%Xt2 = Xt2
        else if (present(res2)) then
            if (allocated(this%Xt2)) then
                if (size(this%Xt2,1) /= res2) then
                    deallocate(this%Xt2)
                    allocate(this%Xt2(res2))
                end if
            else
                allocate(this%Xt2(res2))
            end if
            call fill_uniform(knot_start(this%knot2, this%nc(2), this%degree(2)), &
                knot_end(this%knot2, this%nc(2), this%degree(2)), this%Xt2)
        end if

        ! Set parameter values
        if (present(Xt3)) then
            if (allocated(this%Xt3)) then
                if (size(Xt3,1) /= size(this%Xt3,1)) deallocate(this%Xt3)
            end if
            this%Xt3 = Xt3
        else if (present(res3)) then
            if (allocated(this%Xt3)) then
                if (size(this%Xt3,1) /= res3) then
                    deallocate(this%Xt3)
                    allocate(this%Xt3(res3))
                end if
            else
                allocate(this%Xt3(res3))
            end if
            call fill_uniform(knot_start(this%knot3, this%nc(3), this%degree(3)), &
                knot_end(this%knot3, this%nc(3), this%degree(3)), this%Xt3)
        end if

        ! Set number of geometry points
        this%ng(1) = size(this%Xt1,1)
        this%ng(2) = size(this%Xt2,1)
        this%ng(3) = size(this%Xt3,1)

        call ndgrid(this%Xt1, this%Xt2, this%Xt3, Xt)

        if (this%is_rational()) then ! NURBS
            if (present(Tgc)) then
                call compute_dTgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                    this%Wc, dTgc, Tgc, this%wrap_parameters)
            else
                call compute_dTgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                    this%Wc, dTgc, Tgc_, this%wrap_parameters)
            end if
        else
            if (present(Tgc)) then
                call compute_dTgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                    dTgc, Tgc, this%wrap_parameters)
            else
                call compute_dTgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                    dTgc, Tgc_, this%wrap_parameters)
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the parametric basis gradient at one parameter triple.
    !! Optional `elem` selects and orders a local control subset. `Xt` must have
    !! at least three finite entries; later entries are ignored. The spline must
    !! be initialized. This convenience overload does not diagnose those
    !! preconditions.
    pure subroutine derivative_scalar(this, Xt, dTgc, Tgc, elem)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter triple \([u,v,w]\).
        integer, intent(in), optional :: elem(:)
            !! Optional one-based flattened control indices.
        real(rk), allocatable, intent(out) :: dTgc(:,:)
            !! Selected parametric gradients `[nselected,3]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:)
            !! Optional selected basis values `[nselected]`.
        real(rk), allocatable :: Tgc_(:)

        if (.not. this%err%ok) return

        if (this%is_rational()) then ! NURBS
            if (present(elem)) then
                if (present(Tgc)) then
                    call compute_dTgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                        this%Wc, dTgc, Tgc, elem, this%wrap_parameters)
                else
                    call compute_dTgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                        this%Wc, dTgc, Tgc_, elem, this%wrap_parameters)
                end if
            else
                if (present(Tgc)) then
                    call compute_dTgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                        this%Wc, dTgc, Tgc, wrap_parameters=this%wrap_parameters)
                else
                    call compute_dTgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                        this%Wc, dTgc, Tgc_, wrap_parameters=this%wrap_parameters)
                end if
            end if
        else
            if (present(Tgc)) then
                call compute_dTgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                    dTgc, Tgc, elem, this%wrap_parameters)
            else
                call compute_dTgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                    dTgc, Tgc_, elem, this%wrap_parameters)
            end if
        end if
    end subroutine
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
    !! Parameter selection, cache side effects, and unchecked convenience-
    !! overload preconditions are the same as for the vector
    !! [[nurbs_volume:derivative]] overload.
    pure subroutine derivative2_vector(this, res1, res2, res3, Xt1, Xt2, Xt3, d2Tgc, dTgc, Tgc)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: res1, res2, res3
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:), Xt3(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:,:)
            !! Packed Hessians `[ng_total,3*nc_total,3]`.
        real(rk), allocatable, intent(out), optional :: dTgc(:,:,:)
            !! Optional first derivatives `[ng_total,nc_total,3]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:,:)
            !! Optional basis matrix `[ng_total,nc_total]`.
        real(rk), allocatable :: dTgc_(:,:,:), Tgc_(:,:)
        real(rk), allocatable :: Xt(:,:)

        if (.not. this%err%ok) return

        ! Set parameter values
        if (present(Xt1)) then
            if (allocated(this%Xt1)) then
                if (size(Xt1,1) /= size(this%Xt1,1)) deallocate(this%Xt1)
            end if
            this%Xt1 = Xt1
        else if (present(res1)) then
            if (allocated(this%Xt1)) then
                if (size(this%Xt1,1) /= res1) then
                    deallocate(this%Xt1)
                    allocate(this%Xt1(res1))
                end if
            else
                allocate(this%Xt1(res1))
            end if
            call fill_uniform(knot_start(this%knot1, this%nc(1), this%degree(1)), &
                knot_end(this%knot1, this%nc(1), this%degree(1)), this%Xt1)
        end if

        ! Set parameter values
        if (present(Xt2)) then
            if (allocated(this%Xt2)) then
                if (size(Xt2,1) /= size(this%Xt2,1)) deallocate(this%Xt2)
            end if
            this%Xt2 = Xt2
        else if (present(res2)) then
            if (allocated(this%Xt2)) then
                if (size(this%Xt2,1) /= res2) then
                    deallocate(this%Xt2)
                    allocate(this%Xt2(res2))
                end if
            else
                allocate(this%Xt2(res2))
            end if
            call fill_uniform(knot_start(this%knot2, this%nc(2), this%degree(2)), &
                knot_end(this%knot2, this%nc(2), this%degree(2)), this%Xt2)
        end if

        ! Set parameter values
        if (present(Xt3)) then
            if (allocated(this%Xt3)) then
                if (size(Xt3,1) /= size(this%Xt3,1)) deallocate(this%Xt3)
            end if
            this%Xt3 = Xt3
        else if (present(res3)) then
            if (allocated(this%Xt3)) then
                if (size(this%Xt3,1) /= res3) then
                    deallocate(this%Xt3)
                    allocate(this%Xt3(res3))
                end if
            else
                allocate(this%Xt3(res3))
            end if
            call fill_uniform(knot_start(this%knot3, this%nc(3), this%degree(3)), &
                knot_end(this%knot3, this%nc(3), this%degree(3)), this%Xt3)
        end if

        ! Set number of geometry points
        this%ng(1) = size(this%Xt1,1)
        this%ng(2) = size(this%Xt2,1)
        this%ng(3) = size(this%Xt3,1)

        call ndgrid(this%Xt1, this%Xt2, this%Xt3, Xt)

        if (this%is_rational()) then ! NURBS
            if (present(dTgc)) then
                if (present(Tgc)) then
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                        this%Wc, d2Tgc, dTgc, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                        this%Wc, d2Tgc, dTgc, Tgc_, this%wrap_parameters)
                end if
            else
                if (present(Tgc)) then
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                        this%Wc, d2Tgc, dTgc_, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                        this%Wc, d2Tgc, dTgc_, Tgc_, this%wrap_parameters)
                end if
            end if
        else
            if (present(dTgc)) then
                if (present(Tgc)) then
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                        d2Tgc, dTgc, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                        d2Tgc, dTgc, Tgc_, this%wrap_parameters)
                end if
            else
                if (present(Tgc)) then
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                        d2Tgc, dTgc_, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, &
                        d2Tgc, dTgc_, Tgc_, this%wrap_parameters)
                end if
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate a complete parametric basis Hessian at one parameter triple.
    !!
    !! For `ncp=product(nc)`, `d2Tgc((a-1)*ncp+A,b)` stores
    !! \(\partial^2R_A/(\partial\xi_a\partial\xi_b)\). The output shape is
    !! `[3*ncp,3]`; all mixed entries are populated. `Xt` must have at least
    !! three finite entries; later entries are ignored. The spline must be
    !! initialized. This convenience overload does not diagnose those
    !! preconditions.
    pure subroutine derivative2_scalar(this, Xt, d2Tgc, dTgc, Tgc)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter triple \([u,v,w]\).
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
            !! Packed dense Hessian `[3*nc_total,3]`.
        real(rk), allocatable, intent(out), optional :: dTgc(:,:)
            !! Optional dense first derivatives `[nc_total,3]`.
        real(rk), allocatable, intent(out), optional :: Tgc(:)
            !! Optional dense basis values `[nc_total]`.
        real(rk), allocatable :: dTgc_(:,:), Tgc_(:)

        if (.not. this%err%ok) return

        if (this%is_rational()) then ! NURBS
            if (present(dTgc)) then
                if (present(Tgc)) then
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                        this%Wc, d2Tgc, dTgc, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                        this%Wc, d2Tgc, dTgc, Tgc_, this%wrap_parameters)
                end if
            else
                if (present(Tgc)) then
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                        this%Wc, d2Tgc, dTgc_, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                        this%Wc, d2Tgc, dTgc_, Tgc_, this%wrap_parameters)
                end if
            end if
        else
            if (present(dTgc)) then
                if (present(Tgc)) then
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                        d2Tgc, dTgc, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                        d2Tgc, dTgc, Tgc_, this%wrap_parameters)
                end if
            else
                if (present(Tgc)) then
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                        d2Tgc, dTgc_, Tgc, this%wrap_parameters)
                else
                    call compute_d2Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                        d2Tgc, dTgc_, Tgc_, this%wrap_parameters)
                end if
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the dense volume basis on a tensor-product parameter grid.
    !! Explicit directional vectors take precedence over resolutions; every
    !! in-domain output row sums to one. This convenience overload updates
    !! directional caches and `ng`, but not `Xt` or `Xg`, and assumes an
    !! initialized spline plus available parameters in all three directions;
    !! those preconditions are not diagnosed here.
    pure subroutine basis_vector(this, res1, res2, res3, Xt1, Xt2, Xt3, Tgc)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: res1, res2, res3
            !! Optional positive uniform-grid sizes in directions 1--3.
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:), Xt3(:)
            !! Optional directional parameter vectors.
        real(rk), allocatable, intent(out) :: Tgc(:,:)
            !! Basis matrix `[ng_total,nc_total]`.
        real(rk), allocatable :: Xt(:,:)

        if (.not. this%err%ok) return

        ! Set parameter values
        if (present(Xt1)) then
            if (allocated(this%Xt1)) then
                if (size(Xt1,1) /= size(this%Xt1,1)) deallocate(this%Xt1)
            end if
            this%Xt1 = Xt1
        else if (present(res1)) then
            if (allocated(this%Xt1)) then
                if (size(this%Xt1,1) /= res1) then
                    deallocate(this%Xt1)
                    allocate(this%Xt1(res1))
                end if
            else
                allocate(this%Xt1(res1))
            end if
            call fill_uniform(knot_start(this%knot1, this%nc(1), this%degree(1)), &
                knot_end(this%knot1, this%nc(1), this%degree(1)), this%Xt1)
        end if

        ! Set parameter values
        if (present(Xt2)) then
            if (allocated(this%Xt2)) then
                if (size(Xt2,1) /= size(this%Xt2,1)) deallocate(this%Xt2)
            end if
            this%Xt2 = Xt2
        else if (present(res2)) then
            if (allocated(this%Xt2)) then
                if (size(this%Xt2,1) /= res2) then
                    deallocate(this%Xt2)
                    allocate(this%Xt2(res2))
                end if
            else
                allocate(this%Xt2(res2))
            end if
            call fill_uniform(knot_start(this%knot2, this%nc(2), this%degree(2)), &
                knot_end(this%knot2, this%nc(2), this%degree(2)), this%Xt2)
        end if

        ! Set parameter values
        if (present(Xt3)) then
            if (allocated(this%Xt3)) then
                if (size(Xt3,1) /= size(this%Xt3,1)) deallocate(this%Xt3)
            end if
            this%Xt3 = Xt3
        else if (present(res3)) then
            if (allocated(this%Xt3)) then
                if (size(this%Xt3,1) /= res3) then
                    deallocate(this%Xt3)
                    allocate(this%Xt3(res3))
                end if
            else
                allocate(this%Xt3(res3))
            end if
            call fill_uniform(knot_start(this%knot3, this%nc(3), this%degree(3)), &
                knot_end(this%knot3, this%nc(3), this%degree(3)), this%Xt3)
        end if

        ! Set number of geometry points
        this%ng(1) = size(this%Xt1,1)
        this%ng(2) = size(this%Xt2,1)
        this%ng(3) = size(this%Xt3,1)

        call ndgrid(this%Xt1, this%Xt2, this%Xt3, Xt)

        if (this%is_rational()) then ! NURBS
            Tgc = compute_Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                this%ng, this%Wc, this%wrap_parameters)
        else
            Tgc = compute_Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, this%ng, this%wrap_parameters)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the volume basis at one parameter triple.
    !! Optional `elem` selects and orders a local control subset. `Xt` must have
    !! at least three finite entries; later entries are ignored. The spline must
    !! be initialized. This convenience overload does not diagnose those
    !! preconditions.
    pure subroutine basis_scalar(this, Xt, Tgc, elem)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:)
            !! Parameter triple \([u,v,w]\).
        integer, intent(in), optional :: elem(:)
            !! Optional one-based flattened control indices.
        real(rk), allocatable, intent(out) :: Tgc(:)
            !! Dense or selected basis vector.

        if (.not. this%err%ok) return

        if (this%is_rational()) then ! NURBS
            if (present(elem)) then
                Tgc = compute_Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                    this%Wc, elem, this%wrap_parameters)
            else
                Tgc = compute_Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                    this%Wc, wrap_parameters=this%wrap_parameters)
            end if
        else
            if (present(elem)) then
                Tgc = compute_Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, elem, this%wrap_parameters)
            else
                Tgc = compute_Tgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                    wrap_parameters=this%wrap_parameters)
            end if
        end if
    end subroutine
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

        integer :: k, i, s, j, n_new, n1_old, weight_exponent
        real(rk) :: lo, hi, ktol
        real(rk), allocatable :: Xdir(:,:), Xdir_new(:,:), Xc_new(:,:), Wc_new(:)
        real(rk), allocatable :: knot_work(:), knot_new(:)
        integer :: nc_old(3), nc_work(3), ncoord, ncp_old
        real(rk), allocatable :: Wc_old(:)
        real(rk), allocatable :: A1(:,:), A_re_loc(:,:), A_work(:,:)
        logical :: rational, need_map, periodic

        if (.not. this%err%ok) return
        if (size(Xth) /= size(r)) then
            call this%err%set(&
                code       = 101,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Knot insertion input sizes do not match.",&
                location   = "insert_knots",&
                suggestion = "Pass Xth and r arrays with the same size.")
            return
        end if
        if (any(.not. ieee_is_finite(Xth)) .or. any(r < 0)) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid knot insertion request.",&
                location   = "insert_knots",&
                suggestion = "Use finite knot values and nonnegative insertion counts.")
            return
        end if

        ncoord     = size(this%Xc,2)
        rational = this%is_rational()
        need_map = present(B) .or. present(Bs)
        nc_old = this%nc
        nc_work = this%nc

        select case(dir)
        case(1)
            n1_old = nc_old(1)
            knot_work = this%knot1
        case(2)
            n1_old = nc_old(2)
            knot_work = this%knot2
        case(3)
            n1_old = nc_old(3)
            knot_work = this%knot3
        case default
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid direction for inserting knots.",&
                location   = "insert_knots",&
                suggestion = "Use dir=1 or dir=2 or dir=3 to specify the direction.")
            return
        end select

        if (need_map) then
            ncp_old = size(this%Xc,1)
            if (rational) then
                allocate(Wc_old(ncp_old))
                Wc_old = this%Wc
            end if
            allocate(A1(n1_old,n1_old), source=0.0_rk)
            do concurrent (j = 1:n1_old)
                A1(j,j) = 1.0_rk
            end do
        end if

        if (rational) then
            call pack_refine_volume(this%Xc, this%Wc, nc_work, dir, .true., Xdir, weight_exponent)
        else
            call pack_refine_volume(this%Xc, nc=nc_work, dir=dir, rational=.false., Xdir=Xdir)
        end if
        periodic = periodic_topology(knot_work, this%degree(dir), Xdir)
        do i = 1, size(Xth)
            if (r(i) == 0) cycle
            lo = knot_start(knot_work, nc_work(dir), this%degree(dir))
            hi = knot_end(knot_work, nc_work(dir), this%degree(dir))
            ktol = knot_tolerance(knot_work, this%degree(dir)+1, nc_work(dir)+1)
            if (Xth(i) < lo - ktol .or. Xth(i) > hi + ktol) then
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Knot insertion value is outside the active parameter interval.",&
                    location   = "insert_knots",&
                    suggestion = "Insert knots only inside [knot_start, knot_end].")
                return
            end if
            s = compute_multiplicity(knot_work, Xth(i))
            if (s + r(i) > this%degree(dir) + merge(0,1,periodic)) then
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Knot insertion would exceed the supported knot multiplicity.",&
                    location   = "insert_knots",&
                    suggestion = "Keep the final multiplicity within the current topology's limit.")
                return
            end if
            k = findspan(nc_work(dir)-1, this%degree(dir), Xth(i), knot_work)

            if (periodic) then
                if (need_map) then
                    call insert_knot_periodic_A_5_1(this%degree(dir), knot_work, Xdir, Xth(i), k, s, r(i), &
                        n_new, knot_new, Xdir_new, A_re_loc)
                else
                    call insert_knot_periodic_A_5_1(this%degree(dir), knot_work, Xdir, Xth(i), k, s, r(i), &
                        n_new, knot_new, Xdir_new)
                end if
            else
                if (need_map) then
                    call insert_knot_A_5_1(this%degree(dir), knot_work, Xdir, Xth(i), k, s, r(i), &
                        n_new, knot_new, Xdir_new, A_re_loc)
                else
                    call insert_knot_A_5_1(this%degree(dir), knot_work, Xdir, Xth(i), k, s, r(i), &
                        n_new, knot_new, Xdir_new)
                end if
            end if
            if (n_new < 0) then
                call this%err%set(&
                    code       = 108,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Knot insertion failed to preserve the spline topology.",&
                    location   = "insert_knots",&
                    suggestion = "Check the knot value, multiplicity, and periodic representation.")
                return
            end if
            if (need_map) then
                call sparse_left_matmul(A_re_loc, A1, A_work)
                call move_alloc(A_work, A1)
            end if

            call move_alloc(Xdir_new, Xdir)
            call move_alloc(knot_new, knot_work)
            nc_work(dir) = n_new + 1
        end do

        if (rational) then
            call unpack_refine_volume(Xdir, nc_work, dir, ncoord, .true., Xc_new, Wc_new, weight_exponent)
            select case(dir)
            case(1)
                call this%set(knot1=knot_work, knot2=this%get_knot(2), knot3=this%get_knot(3), Xc=Xc_new, Wc=Wc_new, &
                    degree=this%degree, wrap_parameters=this%wrap_parameters)
            case(2)
                call this%set(knot1=this%get_knot(1), knot2=knot_work, knot3=this%get_knot(3), Xc=Xc_new, Wc=Wc_new, &
                    degree=this%degree, wrap_parameters=this%wrap_parameters)
            case(3)
                call this%set(knot1=this%get_knot(1), knot2=this%get_knot(2), knot3=knot_work, Xc=Xc_new, Wc=Wc_new, &
                    degree=this%degree, wrap_parameters=this%wrap_parameters)
            case default
                return
            end select
        else
            call unpack_refine_volume(Xdir, nc_work, dir, ncoord, .false., Xc_new)
            select case(dir)
            case(1)
                call this%set(knot1=knot_work, knot2=this%get_knot(2), knot3=this%get_knot(3), &
                    Xc=Xc_new, degree=this%degree, wrap_parameters=this%wrap_parameters)
            case(2)
                call this%set(knot1=this%get_knot(1), knot2=knot_work, knot3=this%get_knot(3), &
                    Xc=Xc_new, degree=this%degree, wrap_parameters=this%wrap_parameters)
            case(3)
                call this%set(knot1=this%get_knot(1), knot2=this%get_knot(2), knot3=knot_work, &
                    Xc=Xc_new, degree=this%degree, wrap_parameters=this%wrap_parameters)
            case default
                return
            end select
        end if
        if (need_map) then
            block
                real(rk), allocatable :: S_loc(:,:)
                integer :: nc1, nc2, nc3, n1_new
                integer :: i1, j1, i2, i3, i2_old, i2_new, ii, i3_old, i3_new
                integer :: mS, nS, c

                nc1 = this%nc(1)
                nc2 = this%nc(2)
                nc3 = this%nc(3)

                select case(dir)
                case (1)
                    n1_new = this%nc(1)
                    mS = n1_new*nc2*nc3
                    nS = nc_old(1)*nc2*nc3
                    if (present(Bs)) then
                        allocate(S_loc(mS,nS), source=0.0_rk)
                    else
                        allocate(B(mS*ncoord, nS*ncoord), source=0.0_rk)
                    end if
                    if (rational) then
                        if (present(Bs)) then
                            do concurrent (i3=0:nc3-1, i2=0:nc2-1, j1=1:nc_old(1), i1=1:n1_new, &
                                A1(i1,j1) < 0.0_rk .or. A1(i1,j1) > 0.0_rk .or. ieee_is_nan(A1(i1,j1)))
                                S_loc((i3*nc2 + i2)*n1_new + i1,(i3*nc2 + i2)*nc_old(1) + j1) = &
                                    A1(i1,j1) * Wc_old((i3*nc2 + i2)*nc_old(1) + j1) / &
                                    this%Wc((i3*nc2 + i2)*n1_new + i1)
                            end do
                        else
                            do concurrent (i3=0:nc3-1, i2=0:nc2-1, j1=1:nc_old(1), i1=1:n1_new, c=1:ncoord, &
                                A1(i1,j1) < 0.0_rk .or. A1(i1,j1) > 0.0_rk .or. ieee_is_nan(A1(i1,j1)))
                                B(((i3*nc2 + i2)*n1_new + i1 - 1)*ncoord + c, &
                                    ((i3*nc2 + i2)*nc_old(1) + j1 - 1)*ncoord + c) = &
                                    A1(i1,j1) * Wc_old((i3*nc2 + i2)*nc_old(1) + j1) / &
                                    this%Wc((i3*nc2 + i2)*n1_new + i1)
                            end do
                        end if
                    else
                        if (present(Bs)) then
                            do concurrent (i3=0:nc3-1, i2=0:nc2-1, j1=1:nc_old(1), i1=1:n1_new, &
                                A1(i1,j1) < 0.0_rk .or. A1(i1,j1) > 0.0_rk .or. ieee_is_nan(A1(i1,j1)))
                                S_loc((i3*nc2 + i2)*n1_new + i1,(i3*nc2 + i2)*nc_old(1) + j1) = A1(i1,j1)
                            end do
                        else
                            do concurrent (i3=0:nc3-1, i2=0:nc2-1, j1=1:nc_old(1), i1=1:n1_new, c=1:ncoord, &
                                A1(i1,j1) < 0.0_rk .or. A1(i1,j1) > 0.0_rk .or. ieee_is_nan(A1(i1,j1)))
                                B(((i3*nc2 + i2)*n1_new + i1 - 1)*ncoord + c, &
                                    ((i3*nc2 + i2)*nc_old(1) + j1 - 1)*ncoord + c) = A1(i1,j1)
                            end do
                        end if
                    end if

                case (2)
                    n1_new = this%nc(2)
                    mS = n1_new*nc1*nc3
                    nS = nc_old(2)*nc1*nc3
                    if (present(Bs)) then
                        allocate(S_loc(mS,nS), source=0.0_rk)
                    else
                        allocate(B(mS*ncoord, nS*ncoord), source=0.0_rk)
                    end if
                    if (rational) then
                        if (present(Bs)) then
                            do concurrent (i3=0:nc3-1, i2_old=1:nc_old(2), i2_new=1:n1_new, ii=0:nc1-1, &
                                A1(i2_new,i2_old) < 0.0_rk .or. A1(i2_new,i2_old) > 0.0_rk .or. &
                                ieee_is_nan(A1(i2_new,i2_old)))
                                S_loc(ii + 1 + (i2_new-1)*nc1 + i3*nc1*n1_new, &
                                    ii + 1 + (i2_old-1)*nc1 + i3*nc1*nc_old(2)) = &
                                    A1(i2_new,i2_old) * Wc_old(ii + 1 + (i2_old-1)*nc1 + i3*nc1*nc_old(2)) / &
                                    this%Wc(ii + 1 + (i2_new-1)*nc1 + i3*nc1*n1_new)
                            end do
                        else
                            do concurrent (i3=0:nc3-1, i2_old=1:nc_old(2), i2_new=1:n1_new, ii=0:nc1-1, c=1:ncoord, &
                                A1(i2_new,i2_old) < 0.0_rk .or. A1(i2_new,i2_old) > 0.0_rk .or. &
                                ieee_is_nan(A1(i2_new,i2_old)))
                                B((ii + (i2_new-1)*nc1 + i3*nc1*n1_new)*ncoord + c, &
                                    (ii + (i2_old-1)*nc1 + i3*nc1*nc_old(2))*ncoord + c) = &
                                    A1(i2_new,i2_old) * Wc_old(ii + 1 + (i2_old-1)*nc1 + i3*nc1*nc_old(2)) / &
                                    this%Wc(ii + 1 + (i2_new-1)*nc1 + i3*nc1*n1_new)
                            end do
                        end if
                    else
                        if (present(Bs)) then
                            do concurrent (i3=0:nc3-1, i2_old=1:nc_old(2), i2_new=1:n1_new, ii=0:nc1-1, &
                                A1(i2_new,i2_old) < 0.0_rk .or. A1(i2_new,i2_old) > 0.0_rk .or. &
                                ieee_is_nan(A1(i2_new,i2_old)))
                                S_loc(ii + 1 + (i2_new-1)*nc1 + i3*nc1*n1_new, &
                                    ii + 1 + (i2_old-1)*nc1 + i3*nc1*nc_old(2)) = A1(i2_new,i2_old)
                            end do
                        else
                            do concurrent (i3=0:nc3-1, i2_old=1:nc_old(2), i2_new=1:n1_new, ii=0:nc1-1, c=1:ncoord, &
                                A1(i2_new,i2_old) < 0.0_rk .or. A1(i2_new,i2_old) > 0.0_rk .or. &
                                ieee_is_nan(A1(i2_new,i2_old)))
                                B((ii + (i2_new-1)*nc1 + i3*nc1*n1_new)*ncoord + c, &
                                    (ii + (i2_old-1)*nc1 + i3*nc1*nc_old(2))*ncoord + c) = A1(i2_new,i2_old)
                            end do
                        end if
                    end if

                case (3)
                    n1_new = this%nc(3)
                    mS = n1_new*nc1*nc2
                    nS = nc_old(3)*nc1*nc2
                    if (present(Bs)) then
                        allocate(S_loc(mS,nS), source=0.0_rk)
                    else
                        allocate(B(mS*ncoord, nS*ncoord), source=0.0_rk)
                    end if
                    if (rational) then
                        if (present(Bs)) then
                            do concurrent (i3_old=1:nc_old(3), i3_new=1:n1_new, i2=0:nc2-1, ii=0:nc1-1, &
                                A1(i3_new,i3_old) < 0.0_rk .or. A1(i3_new,i3_old) > 0.0_rk .or. &
                                ieee_is_nan(A1(i3_new,i3_old)))
                                S_loc(ii + 1 + i2*nc1 + (i3_new-1)*nc1*nc2, &
                                    ii + 1 + i2*nc1 + (i3_old-1)*nc1*nc2) = &
                                    A1(i3_new,i3_old) * Wc_old(ii + 1 + i2*nc1 + (i3_old-1)*nc1*nc2) / &
                                    this%Wc(ii + 1 + i2*nc1 + (i3_new-1)*nc1*nc2)
                            end do
                        else
                            do concurrent (i3_old=1:nc_old(3), i3_new=1:n1_new, i2=0:nc2-1, ii=0:nc1-1, c=1:ncoord, &
                                A1(i3_new,i3_old) < 0.0_rk .or. A1(i3_new,i3_old) > 0.0_rk .or. &
                                ieee_is_nan(A1(i3_new,i3_old)))
                                B((ii + i2*nc1 + (i3_new-1)*nc1*nc2)*ncoord + c, &
                                    (ii + i2*nc1 + (i3_old-1)*nc1*nc2)*ncoord + c) = &
                                    A1(i3_new,i3_old) * Wc_old(ii + 1 + i2*nc1 + (i3_old-1)*nc1*nc2) / &
                                    this%Wc(ii + 1 + i2*nc1 + (i3_new-1)*nc1*nc2)
                            end do
                        end if
                    else
                        if (present(Bs)) then
                            do concurrent (i3_old=1:nc_old(3), i3_new=1:n1_new, i2=0:nc2-1, ii=0:nc1-1, &
                                A1(i3_new,i3_old) < 0.0_rk .or. A1(i3_new,i3_old) > 0.0_rk .or. &
                                ieee_is_nan(A1(i3_new,i3_old)))
                                S_loc(ii + 1 + i2*nc1 + (i3_new-1)*nc1*nc2, &
                                    ii + 1 + i2*nc1 + (i3_old-1)*nc1*nc2) = A1(i3_new,i3_old)
                            end do
                        else
                            do concurrent (i3_old=1:nc_old(3), i3_new=1:n1_new, i2=0:nc2-1, ii=0:nc1-1, c=1:ncoord, &
                                A1(i3_new,i3_old) < 0.0_rk .or. A1(i3_new,i3_old) > 0.0_rk .or. &
                                ieee_is_nan(A1(i3_new,i3_old)))
                                B((ii + i2*nc1 + (i3_new-1)*nc1*nc2)*ncoord + c, &
                                    (ii + i2*nc1 + (i3_old-1)*nc1*nc2)*ncoord + c) = A1(i3_new,i3_old)
                            end do
                        end if
                    end if
                case default
                    return
                end select

                if (present(B) .and. present(Bs)) then
                    allocate(B(mS*ncoord, nS*ncoord), source=0.0_rk)
                    do c = 1, ncoord
                        B(c:mS*ncoord:ncoord, c:nS*ncoord:ncoord) = S_loc
                    end do
                end if

                if (present(Bs)) then
                    call move_alloc(S_loc, Bs)
                end if
            end block
        end if

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
        class(nurbs_volume),             intent(inout) :: this
        integer,                         intent(in)    :: dir
            !! Elevated parameter direction, 1, 2, or 3.
        integer,                         intent(in)    :: t
            !! Nonnegative degree increment; zero is the identity operation.
        real(rk), allocatable, optional, intent(out)   :: B(:,:)
            !! Optional coordinate-expanded elevation operator.
        real(rk), allocatable, optional, intent(out)   :: Bs(:,:)
            !! Optional scalar elevation operator.

        integer :: n1_new, weight_exponent
        real(rk), allocatable :: knot_new(:), Xdir(:,:), Xdir_new(:,:), Xc_new(:,:), Wc_new(:)
        real(rk), allocatable :: Tdir(:,:)
        integer :: nc_old(3), d, mS, nS, c
        real(rk), allocatable :: Wc_old(:), S_loc(:,:)
        integer :: i1, j1, i2, i3, i2_old, i2_new, ii, i3_old, i3_new
        integer :: degree_new(3), nc_new(3)
        logical :: rational, need_map, success

        if (.not. this%err%ok) return
        if (t < 0) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid degree elevation request.",&
                location   = "elevate_degree",&
                suggestion = "Use a nonnegative degree elevation count.")
            return
        end if

        select case(dir)
        case (1, 2, 3)
        case default
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid direction for elevating degree.",&
                location   = "elevate_degree",&
                suggestion = "Use dir=1 or dir=2 or dir=3 to specify the direction.")
            return
        end select

        d = size(this%Xc,2)
        rational = this%is_rational()
        need_map = present(B) .or. present(Bs)
        nc_old = this%nc

        if (need_map .and. rational) then
            allocate(Wc_old(size(this%Wc)))
            Wc_old = this%Wc
        end if

        if (rational) then
            call pack_refine_volume(this%Xc, this%Wc, this%nc, dir, .true., Xdir, weight_exponent)
        else
            call pack_refine_volume(this%Xc, nc=this%nc, dir=dir, rational=.false., Xdir=Xdir)
        end if

        if (need_map) then
            select case(dir)
            case(1)
                call elevate_degree_A_5_9(&
                    t        = t,&
                    knot     = this%knot1,&
                    degree   = this%degree(1),&
                    Xcw      = Xdir,&
                    nc_new   = n1_new,&
                    knot_new = knot_new,&
                    Xcw_new  = Xdir_new,&
                    Tmap     = Tdir,&
                    success  = success)
            case(2)
                call elevate_degree_A_5_9(&
                    t        = t,&
                    knot     = this%knot2,&
                    degree   = this%degree(2),&
                    Xcw      = Xdir,&
                    nc_new   = n1_new,&
                    knot_new = knot_new,&
                    Xcw_new  = Xdir_new,&
                    Tmap     = Tdir,&
                    success  = success)
            case(3)
                call elevate_degree_A_5_9(&
                    t        = t,&
                    knot     = this%knot3,&
                    degree   = this%degree(3),&
                    Xcw      = Xdir,&
                    nc_new   = n1_new,&
                    knot_new = knot_new,&
                    Xcw_new  = Xdir_new,&
                    Tmap     = Tdir,&
                    success  = success)
            case default
                return
            end select
        else
            select case(dir)
            case(1)
                call elevate_degree_A_5_9(&
                    t        = t,&
                    knot     = this%knot1,&
                    degree   = this%degree(1),&
                    Xcw      = Xdir,&
                    nc_new   = n1_new,&
                    knot_new = knot_new,&
                    Xcw_new  = Xdir_new,&
                    success  = success)
            case(2)
                call elevate_degree_A_5_9(&
                    t        = t,&
                    knot     = this%knot2,&
                    degree   = this%degree(2),&
                    Xcw      = Xdir,&
                    nc_new   = n1_new,&
                    knot_new = knot_new,&
                    Xcw_new  = Xdir_new,&
                    success  = success)
            case(3)
                call elevate_degree_A_5_9(&
                    t        = t,&
                    knot     = this%knot3,&
                    degree   = this%degree(3),&
                    Xcw      = Xdir,&
                    nc_new   = n1_new,&
                    knot_new = knot_new,&
                    Xcw_new  = Xdir_new,&
                    success  = success)
            case default
                return
            end select
        end if

        if (.not. success) then
            call this%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Degree elevation failed to construct a valid directional transformation.",&
                location   = "elevate_degree",&
                suggestion = "Check the knot vector, degree, and active spline domain.")
            return
        end if

        nc_new = this%nc
        nc_new(dir) = n1_new
        degree_new = this%degree
        degree_new(dir) = degree_new(dir) + t

        if (rational) then
            call unpack_refine_volume(Xdir_new, nc_new, dir, d, .true., Xc_new, Wc_new, weight_exponent)
            select case(dir)
            case(1)
                call this%set(knot1=knot_new, knot2=this%get_knot(2), knot3=this%get_knot(3), &
                    Xc=Xc_new, Wc=Wc_new, degree=degree_new, wrap_parameters=this%wrap_parameters)
            case(2)
                call this%set(knot1=this%get_knot(1), knot2=knot_new, knot3=this%get_knot(3), &
                    Xc=Xc_new, Wc=Wc_new, degree=degree_new, wrap_parameters=this%wrap_parameters)
            case(3)
                call this%set(knot1=this%get_knot(1), knot2=this%get_knot(2), knot3=knot_new, &
                    Xc=Xc_new, Wc=Wc_new, degree=degree_new, wrap_parameters=this%wrap_parameters)
            case default
                return
            end select
        else
            call unpack_refine_volume(Xdir_new, nc_new, dir, d, .false., Xc_new)
            select case(dir)
            case(1)
                call this%set(knot1=knot_new, knot2=this%get_knot(2), knot3=this%get_knot(3), &
                    Xc=Xc_new, degree=degree_new, wrap_parameters=this%wrap_parameters)
            case(2)
                call this%set(knot1=this%get_knot(1), knot2=knot_new, knot3=this%get_knot(3), &
                    Xc=Xc_new, degree=degree_new, wrap_parameters=this%wrap_parameters)
            case(3)
                call this%set(knot1=this%get_knot(1), knot2=this%get_knot(2), knot3=knot_new, &
                    Xc=Xc_new, degree=degree_new, wrap_parameters=this%wrap_parameters)
            case default
                return
            end select
        end if
        if (need_map) then
            select case(dir)
            case(1)
                mS = this%nc(1)*this%nc(2)*this%nc(3)
                nS = nc_old(1) *this%nc(2)*this%nc(3)
                if (present(Bs)) then
                    allocate(S_loc(mS,nS), source=0.0_rk)
                else
                    allocate(B(mS*d, nS*d), source=0.0_rk)
                end if

                if (rational) then
                    if (present(Bs)) then
                        do concurrent (i3=0:this%nc(3)-1, i2=0:this%nc(2)-1, j1=1:nc_old(1), i1=1:this%nc(1), &
                            Tdir(i1,j1) < 0.0_rk .or. Tdir(i1,j1) > 0.0_rk .or. ieee_is_nan(Tdir(i1,j1)))
                            S_loc((i3*this%nc(2) + i2)*this%nc(1) + i1, &
                                (i3*this%nc(2) + i2)*nc_old(1) + j1) = &
                                Tdir(i1,j1) * Wc_old((i3*this%nc(2) + i2)*nc_old(1) + j1) / &
                                this%Wc((i3*this%nc(2) + i2)*this%nc(1) + i1)
                        end do
                    else
                        do concurrent (i3=0:this%nc(3)-1, i2=0:this%nc(2)-1, j1=1:nc_old(1), i1=1:this%nc(1), c=1:d, &
                            Tdir(i1,j1) < 0.0_rk .or. Tdir(i1,j1) > 0.0_rk .or. ieee_is_nan(Tdir(i1,j1)))
                            B(((i3*this%nc(2) + i2)*this%nc(1) + i1 - 1)*d + c, &
                                ((i3*this%nc(2) + i2)*nc_old(1) + j1 - 1)*d + c) = &
                                Tdir(i1,j1) * Wc_old((i3*this%nc(2) + i2)*nc_old(1) + j1) / &
                                this%Wc((i3*this%nc(2) + i2)*this%nc(1) + i1)
                        end do
                    end if
                else
                    if (present(Bs)) then
                        do concurrent (i3=0:this%nc(3)-1, i2=0:this%nc(2)-1, j1=1:nc_old(1), i1=1:this%nc(1), &
                            Tdir(i1,j1) < 0.0_rk .or. Tdir(i1,j1) > 0.0_rk .or. ieee_is_nan(Tdir(i1,j1)))
                            S_loc((i3*this%nc(2) + i2)*this%nc(1) + i1, &
                                (i3*this%nc(2) + i2)*nc_old(1) + j1) = Tdir(i1,j1)
                        end do
                    else
                        do concurrent (i3=0:this%nc(3)-1, i2=0:this%nc(2)-1, j1=1:nc_old(1), i1=1:this%nc(1), c=1:d, &
                            Tdir(i1,j1) < 0.0_rk .or. Tdir(i1,j1) > 0.0_rk .or. ieee_is_nan(Tdir(i1,j1)))
                            B(((i3*this%nc(2) + i2)*this%nc(1) + i1 - 1)*d + c, &
                                ((i3*this%nc(2) + i2)*nc_old(1) + j1 - 1)*d + c) = Tdir(i1,j1)
                        end do
                    end if
                end if

            case(2)
                mS = this%nc(2)*this%nc(1)*this%nc(3)
                nS = nc_old(2) *this%nc(1)*this%nc(3)
                if (present(Bs)) then
                    allocate(S_loc(mS,nS), source=0.0_rk)
                else
                    allocate(B(mS*d, nS*d), source=0.0_rk)
                end if

                if (rational) then
                    if (present(Bs)) then
                        do concurrent (i3=0:this%nc(3)-1, i2_old=1:nc_old(2), i2_new=1:this%nc(2), ii=0:this%nc(1)-1, &
                            Tdir(i2_new,i2_old) < 0.0_rk .or. Tdir(i2_new,i2_old) > 0.0_rk .or. &
                            ieee_is_nan(Tdir(i2_new,i2_old)))
                            S_loc(ii + 1 + (i2_new-1)*this%nc(1) + i3*this%nc(1)*this%nc(2), &
                                ii + 1 + (i2_old-1)*this%nc(1) + i3*this%nc(1)*nc_old(2)) = &
                                Tdir(i2_new,i2_old) * Wc_old(ii + 1 + (i2_old-1)*this%nc(1) + i3*this%nc(1)*nc_old(2)) / &
                                this%Wc(ii + 1 + (i2_new-1)*this%nc(1) + i3*this%nc(1)*this%nc(2))
                        end do
                    else
                        do concurrent (i3=0:this%nc(3)-1, i2_old=1:nc_old(2), i2_new=1:this%nc(2), ii=0:this%nc(1)-1, c=1:d, &
                            Tdir(i2_new,i2_old) < 0.0_rk .or. Tdir(i2_new,i2_old) > 0.0_rk .or. &
                            ieee_is_nan(Tdir(i2_new,i2_old)))
                            B((ii + (i2_new-1)*this%nc(1) + i3*this%nc(1)*this%nc(2))*d + c, &
                                (ii + (i2_old-1)*this%nc(1) + i3*this%nc(1)*nc_old(2))*d + c) = &
                                Tdir(i2_new,i2_old) * Wc_old(ii + 1 + (i2_old-1)*this%nc(1) + i3*this%nc(1)*nc_old(2)) / &
                                this%Wc(ii + 1 + (i2_new-1)*this%nc(1) + i3*this%nc(1)*this%nc(2))
                        end do
                    end if
                else
                    if (present(Bs)) then
                        do concurrent (i3=0:this%nc(3)-1, i2_old=1:nc_old(2), i2_new=1:this%nc(2), ii=0:this%nc(1)-1, &
                            Tdir(i2_new,i2_old) < 0.0_rk .or. Tdir(i2_new,i2_old) > 0.0_rk .or. &
                            ieee_is_nan(Tdir(i2_new,i2_old)))
                            S_loc(ii + 1 + (i2_new-1)*this%nc(1) + i3*this%nc(1)*this%nc(2), &
                                ii + 1 + (i2_old-1)*this%nc(1) + i3*this%nc(1)*nc_old(2)) = Tdir(i2_new,i2_old)
                        end do
                    else
                        do concurrent (i3=0:this%nc(3)-1, i2_old=1:nc_old(2), i2_new=1:this%nc(2), ii=0:this%nc(1)-1, c=1:d, &
                            Tdir(i2_new,i2_old) < 0.0_rk .or. Tdir(i2_new,i2_old) > 0.0_rk .or. &
                            ieee_is_nan(Tdir(i2_new,i2_old)))
                            B((ii + (i2_new-1)*this%nc(1) + i3*this%nc(1)*this%nc(2))*d + c, &
                                (ii + (i2_old-1)*this%nc(1) + i3*this%nc(1)*nc_old(2))*d + c) = Tdir(i2_new,i2_old)
                        end do
                    end if
                end if

            case(3)
                mS = this%nc(3)*this%nc(1)*this%nc(2)
                nS = nc_old(3) *this%nc(1)*this%nc(2)
                if (present(Bs)) then
                    allocate(S_loc(mS,nS), source=0.0_rk)
                else
                    allocate(B(mS*d, nS*d), source=0.0_rk)
                end if

                if (rational) then
                    if (present(Bs)) then
                        do concurrent (i3_old=1:nc_old(3), i3_new=1:this%nc(3), i2=0:this%nc(2)-1, ii=0:this%nc(1)-1, &
                            Tdir(i3_new,i3_old) < 0.0_rk .or. Tdir(i3_new,i3_old) > 0.0_rk .or. &
                            ieee_is_nan(Tdir(i3_new,i3_old)))
                            S_loc(ii + 1 + i2*this%nc(1) + (i3_new-1)*this%nc(1)*this%nc(2), &
                                ii + 1 + i2*this%nc(1) + (i3_old-1)*this%nc(1)*this%nc(2)) = &
                                Tdir(i3_new,i3_old) * Wc_old(ii + 1 + i2*this%nc(1) + (i3_old-1)*this%nc(1)*this%nc(2)) / &
                                this%Wc(ii + 1 + i2*this%nc(1) + (i3_new-1)*this%nc(1)*this%nc(2))
                        end do
                    else
                        do concurrent (i3_old=1:nc_old(3), i3_new=1:this%nc(3), i2=0:this%nc(2)-1, ii=0:this%nc(1)-1, c=1:d, &
                            Tdir(i3_new,i3_old) < 0.0_rk .or. Tdir(i3_new,i3_old) > 0.0_rk .or. &
                            ieee_is_nan(Tdir(i3_new,i3_old)))
                            B((ii + i2*this%nc(1) + (i3_new-1)*this%nc(1)*this%nc(2))*d + c, &
                                (ii + i2*this%nc(1) + (i3_old-1)*this%nc(1)*this%nc(2))*d + c) = &
                                Tdir(i3_new,i3_old) * Wc_old(ii + 1 + i2*this%nc(1) + (i3_old-1)*this%nc(1)*this%nc(2)) / &
                                this%Wc(ii + 1 + i2*this%nc(1) + (i3_new-1)*this%nc(1)*this%nc(2))
                        end do
                    end if
                else
                    if (present(Bs)) then
                        do concurrent (i3_old=1:nc_old(3), i3_new=1:this%nc(3), i2=0:this%nc(2)-1, ii=0:this%nc(1)-1, &
                            Tdir(i3_new,i3_old) < 0.0_rk .or. Tdir(i3_new,i3_old) > 0.0_rk .or. &
                            ieee_is_nan(Tdir(i3_new,i3_old)))
                            S_loc(ii + 1 + i2*this%nc(1) + (i3_new-1)*this%nc(1)*this%nc(2), &
                                ii + 1 + i2*this%nc(1) + (i3_old-1)*this%nc(1)*this%nc(2)) = Tdir(i3_new,i3_old)
                        end do
                    else
                        do concurrent (i3_old=1:nc_old(3), i3_new=1:this%nc(3), i2=0:this%nc(2)-1, ii=0:this%nc(1)-1, c=1:d, &
                            Tdir(i3_new,i3_old) < 0.0_rk .or. Tdir(i3_new,i3_old) > 0.0_rk .or. &
                            ieee_is_nan(Tdir(i3_new,i3_old)))
                            B((ii + i2*this%nc(1) + (i3_new-1)*this%nc(1)*this%nc(2))*d + c, &
                                (ii + i2*this%nc(1) + (i3_old-1)*this%nc(1)*this%nc(2))*d + c) = Tdir(i3_new,i3_old)
                        end do
                    end if
                end if
            case default
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Invalid direction for elevating degree.",&
                    location   = "elevate_degree",&
                    suggestion = "Use dir=1 or dir=2 or dir=3 to specify the direction.")
                return
            end select

            if (present(Bs)) call move_alloc(S_loc, Bs)

            if (present(B) .and. present(Bs)) then
                d = size(this%Xc,2)
                allocate(B(mS*d, nS*d), source=0.0_rk)
                do c = 1, d
                    B(c:mS*d:d, c:nS*d:d) = Bs
                end do
            end if

            if (allocated(Wc_old)) deallocate(Wc_old)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report whether stored nonuniform weights affect basis evaluation.
    !!
    !! With positive stored weights, the result is true precisely when
    !!
    !! \[
    !! \max_A|w_A-w_1|>
    !! 32\,\epsilon_{rk}\max_Aw_A.
    !! \]
    !!
    !! Weights within that threshold are evaluated through the polynomial path,
    !! not merely reported as approximately uniform. An active diagnostic or
    !! absent explicit weights returns false.
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
    !! cyclic identities in [[periodic_topology]]. Every other valid direction
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
    !> Replace sampled-volume visualization connectivity with a copy of `elemConn`.
    pure subroutine set_elem_Xg_vis(this, elemConn)
        class(nurbs_volume), intent(inout) :: this
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
    !> Replace stored IGA element connectivity with a copy of `elemConn`.
    !!
    !! This low-level setter does not verify consistency with the current knot
    !! vectors or control-lattice dimensions.
    pure subroutine set_elem(this, elemConn)
        class(nurbs_volume), intent(inout) :: this
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
    !> Return an allocated copy of control-lattice visualization connectivity.
    pure function get_elem_Xc_vis(this) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        elemConn = this%elemConn_Xc_vis
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocated copy of sampled-volume visualization connectivity.
    pure function get_elem_Xg_vis(this) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        elemConn = this%elemConn_Xg_vis
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an allocated copy of IGA element-to-control-point connectivity.
    pure function get_elem(this) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, allocatable :: elemConn(:,:)

        if (.not. this%err%ok) return

        elemConn = this%elemConn
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
        real(rk), allocatable :: Xt(:,:), Xg_new(:,:)
        real(rk) :: min_X1, max_X1, min_X2, max_X2, min_X3, max_X3
        real(rk) :: start1, end1, start2, end2, start3, end3, scale1, scale2, scale3
        integer :: i, ng_eval(3)

        if (.not. this%err%ok) return

        if (.not. allocated(this%Xc)) then
            call this%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Control points are not set.",&
                location   = "put_to_nurbs",&
                suggestion = "Call set(...) first before mapping points to the NURBS volume.")
            return
        end if
        if (.not. allocated(this%knot1) .or. .not. allocated(this%knot2) .or. .not. allocated(this%knot3)) then
            call this%err%set(&
                code       = 103,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Knot vector is not set.",&
                location   = "put_to_nurbs",&
                suggestion = "Call set(...) first before mapping points to the NURBS volume.")
            return
        end if
        if (size(X,1) < 1 .or. size(X,2) < 3) then
            call this%err%set(&
                code       = 101,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Input points must have at least one point and three coordinates.",&
                location   = "put_to_nurbs",&
                suggestion = "Provide X(:,:) with size(X,1) >= 1 and size(X,2) >= 3.")
            return
        end if

        min_X1 = X(1,1)
        max_X1 = X(1,1)
        min_X2 = X(1,2)
        max_X2 = X(1,2)
        min_X3 = X(1,3)
        max_X3 = X(1,3)
        do i = 2, size(X,1)
            min_X1 = min(min_X1, X(i,1))
            max_X1 = max(max_X1, X(i,1))
            min_X2 = min(min_X2, X(i,2))
            max_X2 = max(max_X2, X(i,2))
            min_X3 = min(min_X3, X(i,3))
            max_X3 = max(max_X3, X(i,3))
        end do
        if (max_X1 <= min_X1 .or. max_X2 <= min_X2 .or. max_X3 <= min_X3) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Input points have zero extent in at least one volume parameter direction.",&
                location   = "put_to_nurbs",&
                suggestion = "Provide points with distinct first, second, and third coordinates.")
            return
        end if

        start1 = knot_start(this%knot1, this%nc(1), this%degree(1))
        end1 = knot_end(this%knot1, this%nc(1), this%degree(1))
        start2 = knot_start(this%knot2, this%nc(2), this%degree(2))
        end2 = knot_end(this%knot2, this%nc(2), this%degree(2))
        start3 = knot_start(this%knot3, this%nc(3), this%degree(3))
        end3 = knot_end(this%knot3, this%nc(3), this%degree(3))

        allocate(Xt(size(X,1),3))
        scale1 = (end1 - start1)/(max_X1 - min_X1)
        scale2 = (end2 - start2)/(max_X2 - min_X2)
        scale3 = (end3 - start3)/(max_X3 - min_X3)
        do concurrent (i = 1:size(X,1))
            Xt(i,1) = (X(i,1) - min_X1)*scale1 + start1
            Xt(i,2) = (X(i,2) - min_X2)*scale2 + start2
            Xt(i,3) = (X(i,3) - min_X3)*scale3 + start3
        end do
        ng_eval(1) = size(X,1)
        ng_eval(2) = 1
        ng_eval(3) = 1
        if (this%is_rational()) then
            Xg_new = compute_Xg(Xt=Xt, knot1=this%knot1, knot2=this%knot2, knot3=this%knot3, degree=this%degree, &
                nc=this%nc, ng=ng_eval, Xc=this%Xc, Wc=this%Wc, wrap_parameters=this%wrap_parameters)
        else
            Xg_new = compute_Xg(Xt=Xt, knot1=this%knot1, knot2=this%knot2, knot3=this%knot3, degree=this%degree, &
                nc=this%nc, ng=ng_eval, Xc=this%Xc, wrap_parameters=this%wrap_parameters)
        end if
        call move_alloc(Xt, this%Xt)
        call move_alloc(Xg_new, this%Xg)
        this%ng = ng_eval
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
    !> Attempt geometry-preserving removal of interior directional knots.
    !! Up to `r(i)` copies are attempted in array order. `Xth(i)` must equal an
    !! existing stored knot exactly and lie more than the current
    !! `knot_tolerance` from both active endpoints. The longest removable prefix
    !! is accepted; the first failed componentwise, roundoff-scaled control test
    !! ends that request. Rational data are tested in homogeneous coordinates,
    !! so the criterion is not a certified Euclidean geometry-error bound.
    !! A verified periodic direction is rejected when any removal is requested,
    !! because bounded Algorithm A5.8 does not update the cyclic knot extension
    !! and repeated control layers.
    !!
    !! **Algorithm:** W. Tiller, *Computer-Aided Design* 24 (1992), 445--453.
    !! [doi:10.1016/0010-4485(92)90012-Y](https://doi.org/10.1016/0010-4485(92)90012-Y).
    !! Applied directionally using Piegl and Tiller, Algorithm A5.8.
    pure subroutine remove_knots(this, dir ,Xth,r)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: dir
            !! Modified parameter direction, 1, 2, or 3.
        real(rk), intent(in), contiguous :: Xth(:)
            !! Exact stored values of interior knots outside the endpoint-tolerance bands.
        integer, intent(in), contiguous :: r(:)
            !! Nonnegative maximum removal counts not exceeding current multiplicity.
        integer :: k, i, s, d, t, weight_exponent
        real(rk) :: lo, hi, ktol
        real(rk), allocatable :: Xdir(:,:), Xdir_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:), knot_work(:)
        integer :: nc_work(3)
        logical :: rational, changed

        if (.not. this%err%ok) return
        if (size(Xth) /= size(r)) then
            call this%err%set(&
                code       = 101,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Knot removal input sizes do not match.",&
                location   = "remove_knots",&
                suggestion = "Pass Xth and r arrays with the same size.")
            return
        end if
        if (any(.not. ieee_is_finite(Xth)) .or. any(r < 0)) then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid knot removal request.",&
                location   = "remove_knots",&
                suggestion = "Use finite knot values and nonnegative removal counts.")
            return
        end if

        select case(dir)
        case(1)
            knot_work = this%knot1
        case(2)
            knot_work = this%knot2
        case(3)
            knot_work = this%knot3
        case default
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid direction for removing knots.",&
                location   = "remove_knots",&
                suggestion = "Use dir=1 or dir=2 or dir=3 to specify the direction.")
            return
        end select
        if (any(r > 0) .and. this%get_parameter_topology(dir) == "periodic") then
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Periodic knot removal is not supported without preserving cyclic topology.",&
                location   = "remove_knots",&
                suggestion = "Retain the periodic knots or first construct an explicit bounded representation.")
            return
        end if

        d = size(this%Xc,2)
        rational = this%is_rational()
        nc_work = this%nc
        changed = .false.

        if (rational) then
            call pack_refine_volume(this%Xc, this%Wc, nc_work, dir, .true., Xdir, weight_exponent)
        else
            call pack_refine_volume(this%Xc, nc=nc_work, dir=dir, rational=.false., Xdir=Xdir)
        end if

        do i = 1, size(Xth)
            if (r(i) == 0) cycle
            lo = knot_start(knot_work, nc_work(dir), this%degree(dir))
            hi = knot_end(knot_work, nc_work(dir), this%degree(dir))
            ktol = knot_tolerance(knot_work, this%degree(dir)+1, nc_work(dir)+1)
            if (Xth(i) < lo - ktol .or. Xth(i) > hi + ktol) then
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Knot removal value is outside the active parameter interval.",&
                    location   = "remove_knots",&
                    suggestion = "Remove knots only inside [knot_start, knot_end].")
                return
            end if
            if (Xth(i) <= lo + ktol .or. Xth(i) >= hi - ktol) then
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Active-boundary knots cannot be removed by geometry-preserving knot removal.",&
                    location   = "remove_knots",&
                    suggestion = "Remove only knots strictly inside the active parameter interval.")
                return
            end if
            s = compute_multiplicity(knot_work, Xth(i))
            if (s == 0 .or. r(i) > s) then
                call this%err%set(&
                    code       = 100,&
                    severity   = 1,&
                    category   = "forcad_nurbs_volume",&
                    message    = "Requested knot removal multiplicity is not available.",&
                    location   = "remove_knots",&
                    suggestion = "Remove only knots present in the knot vector and use r <= current multiplicity.")
                return
            end if
            k = findspan(nc_work(dir)-1, this%degree(dir), Xth(i), knot_work)
            k = k + 1

            call remove_knots_A_5_8(&
                p        = this%degree(dir),&
                knot     = knot_work,&
                Pw       = Xdir,&
                u        = Xth(i),&
                r        = k,&
                s        = s,&
                num      = r(i),&
                t        = t,&
                knot_new = knot_new,&
                Pw_new   = Xdir_new)

            if (t > 0) then
                call move_alloc(Xdir_new, Xdir)
                call move_alloc(knot_new, knot_work)
                nc_work(dir) = size(Xdir,1)
                changed = .true.
            end if
        end do

        if (.not. changed) return

        if (rational) then
            call unpack_refine_volume(Xdir, nc_work, dir, d, .true., Xc_new, Wc_new, weight_exponent)
            select case(dir)
            case(1)
                call this%set(knot1=knot_work, knot2=this%get_knot(2), knot3=this%get_knot(3), &
                    Xc=Xc_new, Wc=Wc_new, degree=this%degree, wrap_parameters=this%wrap_parameters)
            case(2)
                call this%set(knot1=this%get_knot(1), knot2=knot_work, knot3=this%get_knot(3), &
                    Xc=Xc_new, Wc=Wc_new, degree=this%degree, wrap_parameters=this%wrap_parameters)
            case(3)
                call this%set(knot1=this%get_knot(1), knot2=this%get_knot(2), knot3=knot_work, &
                    Xc=Xc_new, Wc=Wc_new, degree=this%degree, wrap_parameters=this%wrap_parameters)
            case default
                return
            end select
        else
            call unpack_refine_volume(Xdir, nc_work, dir, d, .false., Xc_new)
            select case(dir)
            case(1)
                call this%set(knot1=knot_work, knot2=this%get_knot(2), knot3=this%get_knot(3), &
                    Xc=Xc_new, degree=this%degree, wrap_parameters=this%wrap_parameters)
            case(2)
                call this%set(knot1=this%get_knot(1), knot2=knot_work, knot3=this%get_knot(3), &
                    Xc=Xc_new, degree=this%degree, wrap_parameters=this%wrap_parameters)
            case(3)
                call this%set(knot1=this%get_knot(1), knot2=this%get_knot(2), knot3=knot_work, &
                    Xc=Xc_new, degree=this%degree, wrap_parameters=this%wrap_parameters)
            case default
                return
            end select
        end if
    end subroutine
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

        if (.not. this%err%ok) return

        call elemConn_Cn(this%nc(1), this%nc(2), this%nc(3),&
            this%degree(1),this%degree(2),this%degree(3),&
            active_knots(this%knot1, this%nc(1), this%degree(1)),&
            active_knots(this%knot2, this%nc(2), this%degree(2)),&
            active_knots(this%knot3, this%nc(3), this%degree(3)),&
            active_knot_multiplicity(this%knot1, this%nc(1), this%degree(1)),&
            active_knot_multiplicity(this%knot2, this%nc(2), this%degree(2)),&
            active_knot_multiplicity(this%knot3, this%nc(3), this%degree(3)),&
            elemConn)
    end function
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
    !> Display previously exported volume representations with PyVista.
    !!
    !! The filenames are consumed by the optional Python viewer; this procedure
    !! does not export them. Scalar `Xg` point fields and available colormaps
    !! are selected in the viewer. Defining `NOSHOW_PYVISTA` makes the
    !! operation a no-op.
    impure subroutine show(this, vtkfile_Xc, vtkfile_Xg, vtkfile_Xth_in_Xg)
        class(nurbs_volume), intent(inout) :: this
        character(len=*), intent(in) :: vtkfile_Xc
            !! Existing control-geometry VTK path.
        character(len=*), intent(in) :: vtkfile_Xg
            !! Existing sampled-geometry VTK path containing optional point fields.
        character(len=*), intent(in), optional :: vtkfile_Xth_in_Xg
            !! Optional existing parameter-grid VTK path.

        if (.not. this%err%ok) return

#ifndef NOSHOW_PYVISTA
        call show_pyvista_singlepatch(&
            vtkfile_Xc        = vtkfile_Xc,&
            vtkfile_Xg        = vtkfile_Xg,&
            vtkfile_Xth_in_Xg = vtkfile_Xth_in_Xg,&
            rank_name         = "volume")
#endif
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
    !> Find the closest point among cached volume samples.
    !! This discrete search requires `create`; accuracy is limited by the cached
    !! parameter grid. Distance uses the first
    !! `min(size(point_Xg),physical_dimension)` coordinates. Extra query
    !! coordinates are ignored, missing physical coordinates are omitted, and
    !! a longer returned vector is zero-padded. Equal-distance ties select the
    !! earliest direction-1-fastest cached sample. A valid allocated cache and a
    !! nonempty query are unchecked preconditions. If `err` is already active,
    !! optional outputs are not defined.
    pure subroutine nearest_point(this, point_Xg, nearest_Xg, nearest_Xt, id)
        class(nurbs_volume), intent(in) :: this
        real(rk), intent(in), contiguous :: point_Xg(:)
            !! Query point.
        real(rk), intent(out), optional :: nearest_Xg(size(point_Xg))
            !! Optional closest cached physical point.
        real(rk), intent(out), optional :: nearest_Xt(3)
            !! Optional parameter triple of that sample.
        integer, intent(out), optional :: id
            !! Optional one-based flattened sample identifier.
        integer :: id_, i, d, dim_, ncopy, npoint
        real(rk) :: dist, best_dist, dx

        if (.not. this%err%ok) return

        id_ = 1
        npoint = min(size(this%Xg,1), size(this%Xt,1))
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
        if (present(nearest_Xt)) then
            do concurrent (d = 1:3)
                nearest_Xt(d) = this%Xt(id_,d)
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute a continuous nearest-point candidate on the active parameter box.
    !!
    !! A streamed \(10^3\) seed search initializes projected Newton iteration
    !! for
    !! \(f(\boldsymbol{\xi})=
    !! \tfrac12\|\mathbf{V}(\boldsymbol{\xi})-\mathbf{x}\|^2\) on the active
    !! parameter box \(\Omega\). Convergence is tested with
    !!
    !! \[
    !! \mathbf g_P(\boldsymbol{\xi})=
    !! \boldsymbol{\xi}-
    !! \Pi_\Omega\!\left(\boldsymbol{\xi}-\nabla f(\boldsymbol{\xi})\right).
    !! \]
    !!
    !! Vanishing \(\mathbf g_P\) is equivalent to the box-constrained
    !! first-order KKT conditions: a gradient component is zero in a free
    !! direction, nonnegative at its lower bound, or nonpositive at its upper
    !! bound. The Newton proposal is projected onto \(\Omega\); a non-descent
    !! proposal is replaced by the projected negative-gradient displacement.
    !! A trial changes the iterate only after satisfying the Armijo
    !! sufficient-decrease inequality. Failure to accept one of 50 reductions
    !! retains the current iterate and sets code 109. A rejected Hessian sets
    !! code 108, and exhausting `maxit` sets code 109. On successful return, the
    !! result is a KKT-stationary local candidate and need not be the global
    !! nearest point.
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

        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: point_Xg(:)
            !! Finite query with at least the volume's physical dimension; extra components are ignored.
        real(rk), intent(in) :: tol
            !! Nonnegative absolute gradient-mapping norm tolerance; not internally validated.
        integer, intent(in) :: maxit
            !! Positive maximum iteration count; not internally validated.
        real(rk), intent(out) :: nearest_Xt(3)
            !! Returned parameter triple.
        real(rk), intent(out), optional :: nearest_Xg(size(this%Xc,2))
            !! Optional physical point at `nearest_Xt`.

        real(rk) :: obj, obj_trial, grad(3), hess(3,3), dk(3), detH
        real(rk) :: alphak, tau, beta
        real(rk) :: lower_bounds(3), upper_bounds(3), xt(3)
        real(rk), allocatable :: dTgc(:,:), d2Tgc(:,:)
        real(rk) :: Xg(size(this%Xc,2)), residual(size(this%Xc,2)), xk(3)
        real(rk) :: dXg(size(this%Xc,2),3), d2Xg(size(this%Xc,2),3,3)
        integer :: k, l, i, i1, i2, i3, nseed, d, dim_, ncp, a, b
        logical :: convergenz, line_search_accepted
        real(rk) :: dist, best_dist, grad_norm, grad_step

        nearest_Xt = 0.0_rk
        if (present(nearest_Xg)) nearest_Xg = 0.0_rk
        if (.not. this%err%ok) return

        alphak = 0.0_rk
        dk     = 0.0_rk
        k      = 0
        dim_   = size(this%Xc,2)

        ! bounds
        lower_bounds = [knot_start(this%knot1, this%nc(1), this%degree(1)), &
            knot_start(this%knot2, this%nc(2), this%degree(2)), &
            knot_start(this%knot3, this%nc(3), this%degree(3))]
        upper_bounds = [knot_end(this%knot1, this%nc(1), this%degree(1)), &
            knot_end(this%knot2, this%nc(2), this%degree(2)), &
            knot_end(this%knot3, this%nc(3), this%degree(3))]

        ! initial guess (streamed coarse search)
        nseed = 10
        xk = lower_bounds
        best_dist = huge(1.0_rk)
        do i3 = 1, nseed
            if (nseed > 1) then
                xt(3) = lower_bounds(3) + (upper_bounds(3) - lower_bounds(3))*real(i3-1, rk)/real(nseed-1, rk)
            else
                xt(3) = lower_bounds(3)
            end if
            do i2 = 1, nseed
                if (nseed > 1) then
                    xt(2) = lower_bounds(2) + (upper_bounds(2) - lower_bounds(2))*real(i2-1, rk)/real(nseed-1, rk)
                else
                    xt(2) = lower_bounds(2)
                end if
                do i1 = 1, nseed
                    if (nseed > 1) then
                        xt(1) = lower_bounds(1) + (upper_bounds(1) - lower_bounds(1))*real(i1-1, rk)/real(nseed-1, rk)
                    else
                        xt(1) = lower_bounds(1)
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
            end do
        end do

        ! clamp initial guess to bounds
        do d = 1, 3
            xk(d) = max(min(xk(d), upper_bounds(d)), lower_bounds(d))
        end do
        nearest_Xt = xk
        if (present(nearest_Xg)) nearest_Xg = this%cmp_Xg(xk)

        convergenz = .false.

        do while (.not. convergenz .and. k < maxit)

            ! objective, gradient, hessian
            Xg = this%cmp_Xg(xk)
            call this%derivative2(Xt=xk, d2Tgc=d2Tgc, dTgc=dTgc)

            obj = 0.0_rk
            do d = 1, dim_
                residual(d) = Xg(d) - point_Xg(d)
                obj = obj + residual(d)*residual(d)
            end do
            obj = 0.5_rk*obj

            ncp = this%nc(1)*this%nc(2)*this%nc(3)
            dXg = 0.0_rk
            d2Xg = 0.0_rk
            do i = 1, ncp
                do d = 1, dim_
                    dXg(d,1) = dXg(d,1) + dTgc(i,1)*this%Xc(i,d)
                    dXg(d,2) = dXg(d,2) + dTgc(i,2)*this%Xc(i,d)
                    dXg(d,3) = dXg(d,3) + dTgc(i,3)*this%Xc(i,d)
                    d2Xg(d,1,1) = d2Xg(d,1,1) + d2Tgc(i,1)*this%Xc(i,d)
                    d2Xg(d,2,1) = d2Xg(d,2,1) + d2Tgc(ncp+i,1)*this%Xc(i,d)
                    d2Xg(d,3,1) = d2Xg(d,3,1) + d2Tgc(2*ncp+i,1)*this%Xc(i,d)
                    d2Xg(d,1,2) = d2Xg(d,1,2) + d2Tgc(i,2)*this%Xc(i,d)
                    d2Xg(d,2,2) = d2Xg(d,2,2) + d2Tgc(ncp+i,2)*this%Xc(i,d)
                    d2Xg(d,3,2) = d2Xg(d,3,2) + d2Tgc(2*ncp+i,2)*this%Xc(i,d)
                    d2Xg(d,1,3) = d2Xg(d,1,3) + d2Tgc(i,3)*this%Xc(i,d)
                    d2Xg(d,2,3) = d2Xg(d,2,3) + d2Tgc(ncp+i,3)*this%Xc(i,d)
                    d2Xg(d,3,3) = d2Xg(d,3,3) + d2Tgc(2*ncp+i,3)*this%Xc(i,d)
                end do
            end do

            grad = 0.0_rk
            hess = 0.0_rk
            do d = 1, dim_
                grad(1) = grad(1) + residual(d)*dXg(d,1)
                grad(2) = grad(2) + residual(d)*dXg(d,2)
                grad(3) = grad(3) + residual(d)*dXg(d,3)
                do b = 1, 3
                    do a = 1, 3
                        hess(a,b) = hess(a,b) + dXg(d,a)*dXg(d,b) + residual(d)*d2Xg(d,a,b)
                    end do
                end do
            end do
            grad_norm = 0.0_rk
            do i = 1, 3
                grad_norm = grad_norm + &
                    (xk(i)-max(min(xk(i)-grad(i),upper_bounds(i)),lower_bounds(i)))**2
            end do
            grad_norm = sqrt(grad_norm)

            if (grad_norm <= tol) then
                convergenz  = .true.
                nearest_Xt  = xk
                if (present(nearest_Xg)) nearest_Xg = Xg
            else
                ! Newton step
                detH = hess(1,1)*(hess(2,2)*hess(3,3) - hess(2,3)*hess(3,2)) &
                     - hess(1,2)*(hess(2,1)*hess(3,3) - hess(2,3)*hess(3,1)) &
                     + hess(1,3)*(hess(2,1)*hess(3,2) - hess(2,2)*hess(3,1))
                if (abs(detH) <= 32.0_rk*epsilon(1.0_rk)*max(1.0_rk, &
                    abs(hess(1,1)), abs(hess(1,2)), abs(hess(1,3)), &
                    abs(hess(2,1)), abs(hess(2,2)), abs(hess(2,3)), &
                    abs(hess(3,1)), abs(hess(3,2)), abs(hess(3,3)))) then
                    call this%err%set(&
                        code       = 108,&
                        severity   = 1,&
                        category   = "forcad_nurbs_volume",&
                        message    = "Singular Hessian in nearest-point iteration.",&
                        location   = "nearest_point2",&
                        suggestion = "Check the volume geometry, tolerance, or initial sampling grid.")
                    nearest_Xt = xk
                    if (present(nearest_Xg)) nearest_Xg = Xg
                    return
                end if
                dk(1) = -((hess(2,2)*hess(3,3) - hess(2,3)*hess(3,2))*grad(1) &
                        + (hess(1,3)*hess(3,2) - hess(1,2)*hess(3,3))*grad(2) &
                        + (hess(1,2)*hess(2,3) - hess(1,3)*hess(2,2))*grad(3))/detH
                dk(2) = -((hess(2,3)*hess(3,1) - hess(2,1)*hess(3,3))*grad(1) &
                        + (hess(1,1)*hess(3,3) - hess(1,3)*hess(3,1))*grad(2) &
                        + (hess(1,3)*hess(2,1) - hess(1,1)*hess(2,3))*grad(3))/detH
                dk(3) = -((hess(2,1)*hess(3,2) - hess(2,2)*hess(3,1))*grad(1) &
                        + (hess(1,2)*hess(3,1) - hess(1,1)*hess(3,2))*grad(2) &
                        + (hess(1,1)*hess(2,2) - hess(1,2)*hess(2,1))*grad(3))/detH
                do i = 1, 3
                    dk(i) = max(min(xk(i)+dk(i),upper_bounds(i)),lower_bounds(i)) - xk(i)
                end do
                grad_step = grad(1)*dk(1) + grad(2)*dk(2) + grad(3)*dk(3)
                if (grad_step >= 0.0_rk) then
                    do i = 1, 3
                        dk(i) = max(min(xk(i)-grad(i),upper_bounds(i)),lower_bounds(i)) - xk(i)
                    end do
                    grad_step = grad(1)*dk(1) + grad(2)*dk(2) + grad(3)*dk(3)
                end if
                if (.not. ieee_is_finite(grad_step) .or. grad_step >= 0.0_rk) then
                    nearest_Xt = xk
                    if (present(nearest_Xg)) nearest_Xg = Xg
                    call this%err%set(&
                        code       = 109,&
                        severity   = 1,&
                        category   = "forcad_nurbs_volume",&
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
                    xt(1) = xk(1) + alphak*dk(1)
                    xt(2) = xk(2) + alphak*dk(2)
                    xt(3) = xk(3) + alphak*dk(3)
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
                        category   = "forcad_nurbs_volume",&
                        message    = "Armijo line search failed in nearest-point iteration.",&
                        location   = "nearest_point2",&
                        suggestion = "Relax tolerance or check geometry smoothness and conditioning.")
                    return
                end if

                do d = 1, 3
                    xk(d) = max(min(xt(d), upper_bounds(d)), lower_bounds(d))
                end do

                k = k + 1
            end if
        end do

        if (.not. convergenz) then
            nearest_Xt = xk
            if (present(nearest_Xg)) nearest_Xg = this%cmp_Xg(xk)
            call this%err%set(&
                code       = 109,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Nearest-point iteration did not converge.",&
                location   = "nearest_point2",&
                suggestion = "Increase maxit, relax tolerance, or check geometry conditioning.")
        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Extract global control indices for one face of one IGA element.
    !!
    !! `elem` is a one-based row of stored IGA connectivity. Face numbers
    !! are 1 `w_min`, 2 `w_max`, 3 `v_min`, 4 `v_max`,
    !! 5 `u_min`, and 6 `u_max`. The result contains the two
    !! tangential tensor directions with direction 1 varying fastest where it
    !! is tangential. Invalid input or unavailable connectivity returns an
    !! allocated empty vector.
    pure function cmp_elemFace(this, elem, face) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: elem
        integer, intent(in) :: face
        integer, allocatable :: elemConn(:)
        integer :: n(3), ii, jj, k

        if (.not. this%err%ok .or. face < 1 .or. face > 6) then
            allocate(elemConn(0))
            return
        end if

        if (.not. allocated(this%elemConn)) then
            allocate(elemConn(0))
            return
        end if

        if (elem < 1 .or. elem > size(this%elemConn,1)) then
            allocate(elemConn(0))
            return
        end if

        !> number of nodes in each direction
        n = this%degree + 1

        select case(face)
          case(1)
            allocate(elemConn(n(1)*n(2)))
            elemConn = this%elemConn(elem, 1 : n(1)*n(2))
          case(2)
            allocate(elemConn(n(1)*n(2)))
            elemConn = this%elemConn(elem, n(1)*n(2)*n(3)-n(1)*n(2)+1 : n(1)*n(2)*n(3))
          case(3)
            allocate(elemConn(n(1)*n(3)))
            k = 0
            do jj = 1, n(3)
                do ii = 1, n(1)
                    k = k+1
                    elemConn(k) = this%elemConn(elem, n(1)*n(2)*jj + ii - n(1)*n(2))
                end do
            end do
          case(4)
            allocate(elemConn(n(1)*n(3)))
            k = 0
            do jj = 1, n(3)
                do ii = 1, n(1)
                    k = k+1
                    elemConn(k) = this%elemConn(elem, n(1)*n(2)*jj + ii - n(1))
                end do
            end do
          case(5)
            allocate(elemConn(n(2)*n(3)))
            do ii = 1, n(2)*n(3)
                elemConn(ii) = this%elemConn(elem, n(1)*ii - n(1) + 1)
            end do
          case(6)
            allocate(elemConn(n(2)*n(3)))
            do ii = 1, n(2)*n(3)
                elemConn(ii) = this%elemConn(elem, n(1)*ii)
            end do
          case default
            allocate(elemConn(0))
        end select
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Extract one quadrilateral face from a control-lattice VTK cell.
    !!
    !! `elem` selects a stored visualization cell. Face numbering is
    !! `[w_min,w_max,v_min,v_max,u_min,u_max]`. The result has four
    !! global point indices, or is empty for invalid or unavailable state.
    pure function cmp_elemFace_Xc_vis(this, elem, face) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: elem
        integer, intent(in) :: face
        integer, allocatable :: elemConn(:)
        integer :: n(3), ii, jj, k

        if (.not. this%err%ok .or. face < 1 .or. face > 6) then
            allocate(elemConn(0))
            return
        end if

        if (.not. allocated(this%elemConn_Xc_vis)) then
            allocate(elemConn(0))
            return
        end if

        if (elem < 1 .or. elem > size(this%elemConn_Xc_vis,1)) then
            allocate(elemConn(0))
            return
        end if

        !> number of nodes in each direction
        n = [2,2,2]

        select case(face)
          case(1)
            allocate(elemConn(n(1)*n(2)))
            elemConn = this%elemConn_Xc_vis(elem, 1 : n(1)*n(2))
          case(2)
            allocate(elemConn(n(1)*n(2)))
            elemConn = this%elemConn_Xc_vis(elem, n(1)*n(2)*n(3)-n(1)*n(2)+1 : n(1)*n(2)*n(3))
          case(3)
            allocate(elemConn(n(1)*n(3)))
            k = 0
            do jj = 1, n(3)
                do ii = 1, n(1)
                    k = k+1
                    elemConn(k) = this%elemConn_Xc_vis(elem, n(1)*n(2)*jj + ii - n(1)*n(2))
                end do
            end do
          case(4)
            allocate(elemConn(n(1)*n(3)))
            k = 0
            do jj = 1, n(3)
                do ii = 1, n(1)
                    k = k+1
                    elemConn(k) = this%elemConn_Xc_vis(elem, n(1)*n(2)*jj + ii - n(1))
                end do
            end do
          case(5)
            allocate(elemConn(n(2)*n(3)))
            do ii = 1, n(2)*n(3)
                elemConn(ii) = this%elemConn_Xc_vis(elem, n(1)*ii - n(1) + 1)
            end do
          case(6)
            allocate(elemConn(n(2)*n(3)))
            do ii = 1, n(2)*n(3)
                elemConn(ii) = this%elemConn_Xc_vis(elem, n(1)*ii)
            end do
          case default
            allocate(elemConn(0))
        end select
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Extract one quadrilateral face from a sampled-volume VTK cell.
    !!
    !! `elem` selects a stored visualization cell. Face numbering is
    !! `[w_min,w_max,v_min,v_max,u_min,u_max]`. The result has four
    !! global point indices, or is empty for invalid or unavailable state.
    pure function cmp_elemFace_Xg_vis(this, elem, face) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: elem
        integer, intent(in) :: face
        integer, allocatable :: elemConn(:)
        integer :: n(3), ii, jj, k

        if (.not. this%err%ok .or. face < 1 .or. face > 6) then
            allocate(elemConn(0))
            return
        end if

        if (.not. allocated(this%elemConn_Xg_vis)) then
            allocate(elemConn(0))
            return
        end if

        if (elem < 1 .or. elem > size(this%elemConn_Xg_vis,1)) then
            allocate(elemConn(0))
            return
        end if

        !> number of nodes in each direction
        n = [2,2,2]

        select case(face)
          case(1)
            allocate(elemConn(n(1)*n(2)))
            elemConn = this%elemConn_Xg_vis(elem, 1 : n(1)*n(2))
          case(2)
            allocate(elemConn(n(1)*n(2)))
            elemConn = this%elemConn_Xg_vis(elem, n(1)*n(2)*n(3)-n(1)*n(2)+1 : n(1)*n(2)*n(3))
          case(3)
            allocate(elemConn(n(1)*n(3)))
            k = 0
            do jj = 1, n(3)
                do ii = 1, n(1)
                    k = k+1
                    elemConn(k) = this%elemConn_Xg_vis(elem, n(1)*n(2)*jj + ii - n(1)*n(2))
                end do
            end do
          case(4)
            allocate(elemConn(n(1)*n(3)))
            k = 0
            do jj = 1, n(3)
                do ii = 1, n(1)
                    k = k+1
                    elemConn(k) = this%elemConn_Xg_vis(elem, n(1)*n(2)*jj + ii - n(1))
                end do
            end do
          case(5)
            allocate(elemConn(n(2)*n(3)))
            do ii = 1, n(2)*n(3)
                elemConn(ii) = this%elemConn_Xg_vis(elem, n(1)*ii - n(1) + 1)
            end do
          case(6)
            allocate(elemConn(n(2)*n(3)))
            do ii = 1, n(2)*n(3)
                elemConn(ii) = this%elemConn_Xg_vis(elem, n(1)*ii)
            end do
          case default
            allocate(elemConn(0))
        end select
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
    pure subroutine ansatz(this, ie, ig, Tgc, dTgc_dXg, dV, ngauss, d2Tgc_dXg2, strict)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: ie
            !! One-based direction-1-fastest active element index.
        integer, intent(in) :: ig
            !! One-based flattened tensor quadrature-point index.
        real(rk), intent(out) :: dV
            !! Weighted differential volume contribution.
        real(rk), allocatable, intent(out) :: Tgc(:)
            !! Local basis values, size `product(degree+1)`.
        real(rk), allocatable, intent(out) :: dTgc_dXg(:,:)
            !! Local physical gradients `[nlocal,3]`.
        integer, intent(in), optional :: ngauss(3)
            !! Optional positive point counts `[3]`; defaults to `degree+1`.
        real(rk), allocatable, intent(out), optional :: d2Tgc_dXg2(:,:,:)
            !! Optional physical Hessians `[nlocal,3,3]`.
        logical, intent(in), optional :: strict
            !! Reject orientation reversal in addition to singularity.
        real(rk), allocatable :: Xth1(:), Xth2(:), Xth3(:), Xksi(:,:), Wksi(:)
        real(rk) :: Xth1_e(2), Xth2_e(2), Xth3_e(2)
        integer, allocatable :: elem_c(:,:)
        integer, allocatable :: mul1(:), mul2(:), mul3(:)
        real(rk), allocatable :: dTgc_dXt(:,:), hessian(:,:,:)
        real(rk) :: Xt(3), dXt_dXksi(3)
        real(rk) :: dXg_dXt(size(this%Xc, 2), 3), d2Xg_dXt2(size(this%Xc,2),3,3), Jinv(3,3)
        real(rk) :: corrected(3,3)
        real(rk) :: ders1(0:2,0:this%degree(1)), ders2(0:2,0:this%degree(2)), ders3(0:2,0:this%degree(3))
        real(rk) :: detA, g11, g22, g33, tol
        integer :: degree_(3)
        integer :: a, d, e, alpha, beta, e1, e2, e3, i1, i2, i3, nelem1, nelem2, nelem3
        integer :: dim_, nloc, pref1, pref2, pref3, tmp, first1, first2, first3
        logical :: strict_

        if (.not. this%err%ok) return

        strict_ = .false.
        if (present(strict)) strict_ = strict

        if (present(ngauss)) then
            degree_ = ngauss - 1
        else
            degree_ = this%degree
        end if
        if (any(degree_ < 0)) then
            allocate(Tgc(0), dTgc_dXg(0,0))
            if (present(d2Tgc_dXg2)) allocate(d2Tgc_dXg2(0,0,0))
            dV = 0.0_rk
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Every quadrature count must be positive.",&
                location   = "ansatz",&
                suggestion = "Set every ngauss component to at least one.")
            return
        end if
        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], degree_, Xksi, Wksi)

        Xth1 = active_knots(this%knot1, this%nc(1), this%degree(1))
        Xth2 = active_knots(this%knot2, this%nc(2), this%degree(2))
        Xth3 = active_knots(this%knot3, this%nc(3), this%degree(3))
        nelem1 = size(Xth1) - 1
        nelem2 = size(Xth2) - 1
        nelem3 = size(Xth3) - 1
        if (ie < 1 .or. ie > nelem1*nelem2*nelem3 .or. ig < 1 .or. ig > size(Wksi)) then
            allocate(Tgc(0), dTgc_dXg(0,0))
            if (present(d2Tgc_dXg2)) allocate(d2Tgc_dXg2(0,0,0))
            dV = 0.0_rk
            call this%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Invalid ansatz element or quadrature index.",&
                location   = "ansatz",&
                suggestion = "Use an element index in the active element range and a valid quadrature index.")
            return
        end if

        e1 = mod(ie - 1, nelem1) + 1
        tmp = (ie - 1)/nelem1
        e2 = mod(tmp, nelem2) + 1
        e3 = tmp/nelem2 + 1
        Xth1_e(1) = Xth1(e1)
        Xth1_e(2) = Xth1(e1+1)
        Xth2_e(1) = Xth2(e2)
        Xth2_e(2) = Xth2(e2+1)
        Xth3_e(1) = Xth3(e3)
        Xth3_e(2) = Xth3(e3+1)
        mul1 = active_knot_multiplicity(this%knot1, this%nc(1), this%degree(1))
        mul2 = active_knot_multiplicity(this%knot2, this%nc(2), this%degree(2))
        mul3 = active_knot_multiplicity(this%knot3, this%nc(3), this%degree(3))
        pref1 = this%degree(1) + 1
        do i1 = 2, e1
            pref1 = pref1 + mul1(i1)
        end do
        pref2 = this%degree(2) + 1
        do i2 = 2, e2
            pref2 = pref2 + mul2(i2)
        end do
        pref3 = this%degree(3) + 1
        do i3 = 2, e3
            pref3 = pref3 + mul3(i3)
        end do

        nloc = (this%degree(1) + 1)*(this%degree(2) + 1)*(this%degree(3) + 1)
        allocate(elem_c(1,nloc))
        do concurrent (i3 = 0:this%degree(3), i2 = 0:this%degree(2), i1 = 0:this%degree(1))
            elem_c(1, i3*(this%degree(2)+1)*(this%degree(1)+1) + i2*(this%degree(1)+1) + i1 + 1) = &
                ((pref3 - this%degree(3) + i3) - 1)*this%nc(1)*this%nc(2) + &
                ((pref2 - this%degree(2) + i2) - 1)*this%nc(1) + pref1 - this%degree(1) + i1
        end do

        dXt_dXksi(1) = Xth1_e(2) - Xth1_e(1)
        dXt_dXksi(2) = Xth2_e(2) - Xth2_e(1)
        dXt_dXksi(3) = Xth3_e(2) - Xth3_e(1)
        Xt(1) = Xth1_e(1) + Xksi(ig,1)*dXt_dXksi(1)
        Xt(2) = Xth2_e(1) + Xksi(ig,2)*dXt_dXksi(2)
        Xt(3) = Xth3_e(1) + Xksi(ig,3)*dXt_dXksi(3)

        if (present(d2Tgc_dXg2)) then
            allocate(Tgc(nloc), dTgc_dXt(nloc,3), hessian(nloc,3,3))
            call basis_bspline_der_all_active(&
                Xt     = Xt(1),&
                knot   = this%knot1,&
                nc     = this%nc(1),&
                degree = this%degree(1),&
                nder   = 2,&
                first  = first1,&
                ders   = ders1)
            call basis_bspline_der_all_active(&
                Xt     = Xt(2),&
                knot   = this%knot2,&
                nc     = this%nc(2),&
                degree = this%degree(2),&
                nder   = 2,&
                first  = first2,&
                ders   = ders2)
            call basis_bspline_der_all_active(&
                Xt     = Xt(3),&
                knot   = this%knot3,&
                nc     = this%nc(3),&
                degree = this%degree(3),&
                nder   = 2,&
                first  = first3,&
                ders   = ders3)
            if (this%is_rational()) then
                call tensor_basis_derivatives2_local(&
                    first1   = first1,&
                    nc1      = this%nc(1),&
                    ders1    = ders1,&
                    first2   = first2,&
                    nc2      = this%nc(2),&
                    ders2    = ders2,&
                    first3   = first3,&
                    nc3      = this%nc(3),&
                    ders3    = ders3,&
                    basis    = Tgc,&
                    gradient = dTgc_dXt,&
                    hessian  = hessian,&
                    Wc       = this%Wc)
            else
                call tensor_basis_derivatives2_local(&
                    first1   = first1,&
                    nc1      = this%nc(1),&
                    ders1    = ders1,&
                    first2   = first2,&
                    nc2      = this%nc(2),&
                    ders2    = ders2,&
                    first3   = first3,&
                    nc3      = this%nc(3),&
                    ders3    = ders3,&
                    basis    = Tgc,&
                    gradient = dTgc_dXt,&
                    hessian  = hessian)
            end if
        else if (this%is_rational()) then
            call compute_dTgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                this%Wc, dTgc_dXt, Tgc, elem_c(1,:), this%wrap_parameters)
        else
            call compute_dTgc(Xt, this%knot1, this%knot2, this%knot3, this%degree, this%nc, &
                dTgc_dXt, Tgc, elem_c(1,:), this%wrap_parameters)
        end if

        nloc = size(Tgc)
        dim_ = size(this%Xc, 2)
        if (dim_ /= 3) then
            call this%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Volume ansatz requires three-dimensional physical coordinates.",&
                location   = "ansatz",&
                suggestion = "Use volume control points with exactly three coordinates.")
            dV = 0.0_rk
            return
        end if
        allocate(dTgc_dXg(nloc, 3), source=0.0_rk)
        if (present(d2Tgc_dXg2)) allocate(d2Tgc_dXg2(nloc,3,3), source=0.0_rk)
        dXg_dXt = 0.0_rk
        d2Xg_dXt2 = 0.0_rk
        do a = 1, nloc
            do d = 1, dim_
                dXg_dXt(d,1) = dXg_dXt(d,1) + this%Xc(elem_c(1,a),d)*dTgc_dXt(a,1)
                dXg_dXt(d,2) = dXg_dXt(d,2) + this%Xc(elem_c(1,a),d)*dTgc_dXt(a,2)
                dXg_dXt(d,3) = dXg_dXt(d,3) + this%Xc(elem_c(1,a),d)*dTgc_dXt(a,3)
                if (present(d2Tgc_dXg2)) then
                    do beta = 1, 3
                        do alpha = 1, 3
                            d2Xg_dXt2(d,alpha,beta) = d2Xg_dXt2(d,alpha,beta) + &
                                this%Xc(elem_c(1,a),d)*hessian(a,alpha,beta)
                        end do
                    end do
                end if
            end do
        end do

        detA = &
            + dXg_dXt(1,1)*( dXg_dXt(2,2)*dXg_dXt(3,3) - dXg_dXt(2,3)*dXg_dXt(3,2) ) &
            - dXg_dXt(1,2)*( dXg_dXt(2,1)*dXg_dXt(3,3) - dXg_dXt(2,3)*dXg_dXt(3,1) ) &
            + dXg_dXt(1,3)*( dXg_dXt(2,1)*dXg_dXt(3,2) - dXg_dXt(2,2)*dXg_dXt(3,1) )
        g11 = dot_product(dXg_dXt(:,1), dXg_dXt(:,1))
        g22 = dot_product(dXg_dXt(:,2), dXg_dXt(:,2))
        g33 = dot_product(dXg_dXt(:,3), dXg_dXt(:,3))
        tol = 64.0_rk*epsilon(1.0_rk)*sqrt(g11)*sqrt(g22)*sqrt(g33)
        if (g11 <= 0.0_rk .or. g22 <= 0.0_rk .or. g33 <= 0.0_rk .or. abs(detA) <= tol) then
            call this%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Singular volume Jacobian.",&
                location   = "ansatz",&
                suggestion = "Check for degenerate control points or invalid element geometry.")
            dV = 0.0_rk
            return
        end if
        if (strict_ .and. detA < 0.0_rk) then
            call this%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Inverted volume Jacobian in strict analysis mode.",&
                location   = "ansatz",&
                suggestion = "Orient the volume parameter directions consistently.")
            dV = 0.0_rk
            return
        end if

        Jinv(1,1) = dXg_dXt(2,2)*dXg_dXt(3,3) - dXg_dXt(2,3)*dXg_dXt(3,2)
        Jinv(1,2) = dXg_dXt(1,3)*dXg_dXt(3,2) - dXg_dXt(1,2)*dXg_dXt(3,3)
        Jinv(1,3) = dXg_dXt(1,2)*dXg_dXt(2,3) - dXg_dXt(1,3)*dXg_dXt(2,2)
        Jinv(2,1) = dXg_dXt(2,3)*dXg_dXt(3,1) - dXg_dXt(2,1)*dXg_dXt(3,3)
        Jinv(2,2) = dXg_dXt(1,1)*dXg_dXt(3,3) - dXg_dXt(1,3)*dXg_dXt(3,1)
        Jinv(2,3) = dXg_dXt(1,3)*dXg_dXt(2,1) - dXg_dXt(1,1)*dXg_dXt(2,3)
        Jinv(3,1) = dXg_dXt(2,1)*dXg_dXt(3,2) - dXg_dXt(2,2)*dXg_dXt(3,1)
        Jinv(3,2) = dXg_dXt(1,2)*dXg_dXt(3,1) - dXg_dXt(1,1)*dXg_dXt(3,2)
        Jinv(3,3) = dXg_dXt(1,1)*dXg_dXt(2,2) - dXg_dXt(1,2)*dXg_dXt(2,1)
        Jinv = Jinv/detA

        do d = 1, 3
            do a = 1, nloc
                dTgc_dXg(a,d) = dTgc_dXt(a,1)*Jinv(1,d) &
                    + dTgc_dXt(a,2)*Jinv(2,d) &
                    + dTgc_dXt(a,3)*Jinv(3,d)
            end do
        end do
        if (present(d2Tgc_dXg2)) then
            do a = 1, nloc
                do beta = 1, 3
                    do alpha = 1, 3
                        corrected(alpha,beta) = hessian(a,alpha,beta)
                        do d = 1, 3
                            corrected(alpha,beta) = corrected(alpha,beta) - &
                                dTgc_dXg(a,d)*d2Xg_dXt2(d,alpha,beta)
                        end do
                    end do
                end do
                do e = 1, 3
                    do d = 1, 3
                        do beta = 1, 3
                            do alpha = 1, 3
                                d2Tgc_dXg2(a,d,e) = d2Tgc_dXg2(a,d,e) + &
                                    Jinv(alpha,d)*corrected(alpha,beta)*Jinv(beta,e)
                            end do
                        end do
                    end do
                end do
            end do
        end if
        dV = abs(detA*dXt_dXksi(1)*dXt_dXksi(2)*dXt_dXksi(3))*Wksi(ig)
    end subroutine
    !===============================================================================

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Integrate physical volume over all nonzero knot-span cells.
    !!
    !! The mapping must have exactly three physical coordinates. The computed
    !! approximation is
    !!
    !! \[V_h=\sum_e\sum_g\omega_g\,
    !!   |\det J(u_{eg},v_{eg},w_{eg})|\,|\det J_{e,\xi}|.\]
    !!
    !! Every nonzero active knot-span cell receives an independent tensor
    !! Gauss-Legendre rule. Because the absolute Jacobian determinant is
    !! generally nonpolynomial, especially for rational mappings, the default
    !! `degree+1` rule is an approximation and not an exact-volume guarantee.
    !! A quadrature point is treated as singular when the determinant test falls
    !! below a `64*epsilon(rk)` relative scale; diagnostic code 108 is set and
    !! that point's contribution is omitted from the partial result. Taking the
    !! absolute determinant measures folded or orientation-reversed cells
    !! positively; this routine does not certify global injectivity.
    !!
    !! A valid initialized volume and positive entries of `ngauss` are caller
    !! preconditions; quadrature counts are not diagnosed here. If an existing
    !! diagnostic is active, the routine returns `volume=0`. A mapping with a
    !! physical dimension other than three also returns zero and sets code 108.
    !! Rational Jacobian columns are evaluated in a local exact radix-power
    !! weight gauge using \((S_{,a}-\mathbf V W_{,a})/W\). Avoiding explicit
    !! \(W^2\) preserves the integrated measure under common finite weight
    !! rescaling.
    pure subroutine cmp_volume(this, volume, ngauss)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(out) :: volume
            !! Accumulated nonnegative physical volume.
        integer, intent(in), optional :: ngauss(3)
            !! Optional positive points per direction; defaults to `degree+1`.
        real(rk), allocatable :: Xth1(:), Xth2(:), Xth3(:), Xksi(:,:), Wksi(:)
        real(rk) :: N1(0:this%degree(1)), N2(0:this%degree(2)), N3(0:this%degree(3))
        real(rk) :: dN1(0:this%degree(1)), dN2(0:this%degree(2)), dN3(0:this%degree(3))
        real(rk) :: Xt(3), dXt_dXksi(3)
        real(rk) :: bval, db1v, db2v, db3v, b32, denom, ddenom1, ddenom2, ddenom3, wi
        real(rk) :: dV, S_d, dS1_d, dS2_d, dS3_d, dx1_d, dx2_d, dx3_d
        real(rk) :: dx11, dx12, dx13, dx21, dx22, dx23, dx31, dx32, dx33
        real(rk) :: detA, g11, g22, g33, tol
        integer :: ie, ig, e1, e2, e3, nelem1, nelem2, tmp, dim_, n12
        integer :: first1, first2, first3, j1, j2, j3, idx1, idx2, idx3, idx, d
        integer :: j10, j11, j20, j21, j30, j31
        integer :: ngauss_(3), nelem, nquad, singular_state, weight_exponent
        logical :: rational

        volume = 0.0_rk
        if (.not. this%err%ok) return

        if (present(ngauss)) then
            ngauss_ = ngauss
        else
            ngauss_ = this%degree + 1
        end if

        call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], ngauss_-1, Xksi, Wksi)
        Xth1 = active_knots(this%knot1, this%nc(1), this%degree(1))
        Xth2 = active_knots(this%knot2, this%nc(2), this%degree(2))
        Xth3 = active_knots(this%knot3, this%nc(3), this%degree(3))
        nelem1 = size(Xth1) - 1
        nelem2 = size(Xth2) - 1
        nelem = nelem1*nelem2*(size(Xth3) - 1)
        nquad = product(ngauss_)
        n12 = this%nc(1)*this%nc(2)
        rational = this%is_rational()
        dim_ = size(this%Xc, 2)
        if (dim_ /= 3) then
            call this%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Volume ansatz requires three-dimensional physical coordinates.",&
                location   = "cmp_volume",&
                suggestion = "Use volume control points with exactly three coordinates.")
            return
        end if
        singular_state = 0
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do ie = 1, nelem
#else
        do concurrent (ie = 1: nelem) &
            local(N1, N2, N3, dN1, dN2, dN3, Xt, dXt_dXksi, bval, db1v, db2v, db3v, &
                b32, denom, ddenom1, ddenom2, ddenom3, wi, dV, S_d, dS1_d, &
                dS2_d, dS3_d, dx1_d, dx2_d, dx3_d, dx11, dx12, dx13, dx21, dx22, dx23, &
                dx31, dx32, dx33, detA, g11, g22, g33, tol, ig, e1, e2, e3, tmp, first1, first2, first3, &
                j1, j2, j3, idx1, idx2, idx3, idx, d, j10, j11, j20, j21, j30, j31, weight_exponent) &
            reduce(+:volume) &
            reduce(max:singular_state)
#endif
            e1 = mod(ie - 1, nelem1) + 1
            tmp = (ie - 1)/nelem1
            e2 = mod(tmp, nelem2) + 1
            e3 = tmp/nelem2 + 1
            dXt_dXksi(1) = Xth1(e1+1) - Xth1(e1)
            dXt_dXksi(2) = Xth2(e2+1) - Xth2(e2)
            dXt_dXksi(3) = Xth3(e3+1) - Xth3(e3)
            dV = 0.0_rk
            do ig = 1, nquad
                Xt(1) = Xth1(e1) + Xksi(ig,1)*dXt_dXksi(1)
                Xt(2) = Xth2(e2) + Xksi(ig,2)*dXt_dXksi(2)
                Xt(3) = Xth3(e3) + Xksi(ig,3)*dXt_dXksi(3)

                call basis_bspline_der(Xt(1), this%knot1, this%nc(1), this%degree(1), 1, first1, N1, dN1)
                call basis_bspline_der(Xt(2), this%knot2, this%nc(2), this%degree(2), 1, first2, N2, dN2)
                call basis_bspline_der(Xt(3), this%knot3, this%nc(3), this%degree(3), 1, first3, N3, dN3)

                j10 = max(0, 1 - first1)
                j11 = min(this%degree(1), this%nc(1) - first1)
                j20 = max(0, 1 - first2)
                j21 = min(this%degree(2), this%nc(2) - first2)
                j30 = max(0, 1 - first3)
                j31 = min(this%degree(3), this%nc(3) - first3)
                dx11 = 0.0_rk
                dx12 = 0.0_rk
                dx13 = 0.0_rk
                dx21 = 0.0_rk
                dx22 = 0.0_rk
                dx23 = 0.0_rk
                dx31 = 0.0_rk
                dx32 = 0.0_rk
                dx33 = 0.0_rk
                g11 = 0.0_rk
                g22 = 0.0_rk
                g33 = 0.0_rk

                if (rational) then
                    idx = ((first3+j30-1)*this%nc(2) + first2+j20-1)*this%nc(1) + first1+j10
                    weight_exponent = exponent(this%Wc(idx))
                    do j3 = j30, j31
                        idx3 = first3 + j3
                        do j2 = j20, j21
                            idx2 = first2 + j2
                            do j1 = j10, j11
                                idx1 = first1 + j1
                                idx = (idx3 - 1)*n12 + (idx2 - 1)*this%nc(1) + idx1
                                weight_exponent = max(weight_exponent, exponent(this%Wc(idx)))
                            end do
                        end do
                    end do
                    denom = 0.0_rk
                    ddenom1 = 0.0_rk
                    ddenom2 = 0.0_rk
                    ddenom3 = 0.0_rk
                    do j3 = j30, j31
                        idx3 = first3 + j3
                        do j2 = j20, j21
                            idx2 = first2 + j2
                            b32 = N3(j3)*N2(j2)
                            do j1 = j10, j11
                                idx1 = first1 + j1
                                idx = (idx3 - 1)*n12 + (idx2 - 1)*this%nc(1) + idx1
                                wi = scale(this%Wc(idx), -weight_exponent)
                                bval = b32*N1(j1)
                                db1v = b32*dN1(j1)
                                db2v = N3(j3)*dN2(j2)*N1(j1)
                                db3v = dN3(j3)*N2(j2)*N1(j1)
                                denom = denom + bval*wi
                                ddenom1 = ddenom1 + db1v*wi
                                ddenom2 = ddenom2 + db2v*wi
                                ddenom3 = ddenom3 + db3v*wi
                            end do
                        end do
                    end do

                    if (denom > 0.0_rk) then
                        do d = 1, 3
                            S_d = 0.0_rk
                            dS1_d = 0.0_rk
                            dS2_d = 0.0_rk
                            dS3_d = 0.0_rk
                            do j3 = j30, j31
                                idx3 = first3 + j3
                                do j2 = j20, j21
                                    idx2 = first2 + j2
                                    b32 = N3(j3)*N2(j2)
                                    do j1 = j10, j11
                                        idx1 = first1 + j1
                                        idx = (idx3 - 1)*n12 + (idx2 - 1)*this%nc(1) + idx1
                                        wi = scale(this%Wc(idx), -weight_exponent)
                                        bval = b32*N1(j1)
                                        db1v = b32*dN1(j1)
                                        db2v = N3(j3)*dN2(j2)*N1(j1)
                                        db3v = dN3(j3)*N2(j2)*N1(j1)
                                        S_d = S_d + bval*wi*this%Xc(idx,d)
                                        dS1_d = dS1_d + db1v*wi*this%Xc(idx,d)
                                        dS2_d = dS2_d + db2v*wi*this%Xc(idx,d)
                                        dS3_d = dS3_d + db3v*wi*this%Xc(idx,d)
                                    end do
                                end do
                            end do
                            S_d = S_d/denom
                            dx1_d = (dS1_d - S_d*ddenom1)/denom
                            dx2_d = (dS2_d - S_d*ddenom2)/denom
                            dx3_d = (dS3_d - S_d*ddenom3)/denom
                            g11 = g11 + dx1_d*dx1_d
                            g22 = g22 + dx2_d*dx2_d
                            g33 = g33 + dx3_d*dx3_d
                            if (d == 1) then
                                dx11 = dx1_d
                                dx12 = dx2_d
                                dx13 = dx3_d
                            else if (d == 2) then
                                dx21 = dx1_d
                                dx22 = dx2_d
                                dx23 = dx3_d
                            else
                                dx31 = dx1_d
                                dx32 = dx2_d
                                dx33 = dx3_d
                            end if
                        end do
                    end if
                else
                    do d = 1, 3
                        dx1_d = 0.0_rk
                        dx2_d = 0.0_rk
                        dx3_d = 0.0_rk
                        do j3 = j30, j31
                            idx3 = first3 + j3
                            do j2 = j20, j21
                                idx2 = first2 + j2
                                b32 = N3(j3)*N2(j2)
                                do j1 = j10, j11
                                    idx1 = first1 + j1
                                    idx = (idx3 - 1)*n12 + (idx2 - 1)*this%nc(1) + idx1
                                    db1v = b32*dN1(j1)
                                    db2v = N3(j3)*dN2(j2)*N1(j1)
                                    db3v = dN3(j3)*N2(j2)*N1(j1)
                                    dx1_d = dx1_d + db1v*this%Xc(idx,d)
                                    dx2_d = dx2_d + db2v*this%Xc(idx,d)
                                    dx3_d = dx3_d + db3v*this%Xc(idx,d)
                                end do
                            end do
                        end do
                        g11 = g11 + dx1_d*dx1_d
                        g22 = g22 + dx2_d*dx2_d
                        g33 = g33 + dx3_d*dx3_d
                        if (d == 1) then
                            dx11 = dx1_d
                            dx12 = dx2_d
                            dx13 = dx3_d
                        else if (d == 2) then
                            dx21 = dx1_d
                            dx22 = dx2_d
                            dx23 = dx3_d
                        else
                            dx31 = dx1_d
                            dx32 = dx2_d
                            dx33 = dx3_d
                        end if
                    end do
                end if

                detA = &
                    + dx11*(dx22*dx33 - dx23*dx32) &
                    - dx12*(dx21*dx33 - dx23*dx31) &
                    + dx13*(dx21*dx32 - dx22*dx31)
                tol = 64.0_rk*epsilon(1.0_rk)*sqrt(g11)*sqrt(g22)*sqrt(g33)
                if (g11 <= 0.0_rk .or. g22 <= 0.0_rk .or. g33 <= 0.0_rk .or. abs(detA) <= tol) then
                    singular_state = max(singular_state, 1)
                else
                    dV = dV + abs(detA*dXt_dXksi(1)*dXt_dXksi(2)*dXt_dXksi(3))*Wksi(ig)
                end if
            end do
            volume = volume + dV
        end do

        if (singular_state > 0) then
            call this%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
                message    = "Singular volume Jacobian.",&
                location   = "cmp_volume",&
                suggestion = "Check for degenerate control points or invalid element geometry.")
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate a rational volume at a list of parameter triples.
    !!
    !! `Xt` contains one `[u,v,w]` triple per row. The result has shape
    !! `[ng_,size(Xc,2)]`, where `ng_=product(ng)` when supplied and otherwise
    !! `size(Xt,1)`. Direction 1 varies fastest in the flattened control lattice;
    !! only local tensor-product support contributes to each row.
    pure function compute_Xg_nurbs_3d(Xt, knot1, knot2, knot3, degree, nc, ng, Xc, Wc, wrap_parameters) result(Xg)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: ng(3)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous :: Wc(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Xg(:,:)
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3))
        real(rk) :: bw, wi, wsum, n32, wtol
        integer :: i, j1, j2, j3, d, idx, first1, first2, first3, ng_, dim_
        integer :: j10, j11, j20, j21, j30, j31, weight_exponent
        logical :: wrap_parameters_(3)

        if (present(ng)) then
            ng_ = product(ng)
        else
            ng_ = size(Xt, 1)
        end if

        dim_ = size(Xc,2)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        allocate(Xg(ng_, dim_))
        wtol = 0.0_rk
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N1, N2, N3, bw, wi, wsum, n32, j1, j2, j3, d, idx, &
            first1, first2, first3, j10, j11, j20, j21, j30, j31, weight_exponent)
#endif
            do d = 1, dim_
                Xg(i,d) = 0.0_rk
            end do
            call basis_bspline_der(map_parameter(Xt(i,1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
                knot1, nc(1), degree(1), 0, first1, N1)
            call basis_bspline_der(map_parameter(Xt(i,2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
                knot2, nc(2), degree(2), 0, first2, N2)
            call basis_bspline_der(map_parameter(Xt(i,3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
                knot3, nc(3), degree(3), 0, first3, N3)
            j10 = max(0, 1 - first1)
            j11 = min(degree(1), nc(1) - first1)
            j20 = max(0, 1 - first2)
            j21 = min(degree(2), nc(2) - first2)
            j30 = max(0, 1 - first3)
            j31 = min(degree(3), nc(3) - first3)
            idx = ((first3+j30-1)*nc(2) + first2+j20-1)*nc(1) + first1+j10
            weight_exponent = exponent(Wc(idx))
            do j3 = j30, j31
                do j2 = j20, j21
                    do j1 = j10, j11
                        idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1+j1
                        weight_exponent = max(weight_exponent, exponent(Wc(idx)))
                    end do
                end do
            end do
            wsum = 0.0_rk
            do j3 = j30, j31
                do j2 = j20, j21
                    n32 = N3(j3)*N2(j2)
                    do j1 = j10, j11
                        idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1+j1
                        wi = scale(Wc(idx), -weight_exponent)
                        bw = n32*N1(j1)*wi
                        wsum = wsum + bw
                        do d = 1, dim_
                            Xg(i,d) = Xg(i,d) + bw*Xc(idx,d)
                        end do
                    end do
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
    !> Evaluate a rational volume at one parameter triple.
    !!
    !! The fixed-size result contains one value per physical coordinate.
    !! Directional wrapping is optional; a zero rational denominator leaves the
    !! initialized result zero.
    pure function compute_Xg_nurbs_3d_1point(Xt, knot1, knot2, knot3, degree, nc, Xc, Wc, wrap_parameters) result(Xg)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        real(rk), intent(in), contiguous :: Xc(:,:)
        real(rk), intent(in), contiguous :: Wc(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk) :: Xg(size(Xc,2))
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3))
        real(rk) :: bw, wi, wsum, n32, wtol
        integer :: j1, j2, j3, d, idx, first1, first2, first3, dim_
        integer :: j10, j11, j20, j21, j30, j31, weight_exponent
        logical :: wrap_parameters_(3)

        Xg = 0.0_rk
        dim_ = size(Xc,2)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        wtol = 0.0_rk
        call basis_bspline_der(map_parameter(Xt(1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
            knot1, nc(1), degree(1), 0, first1, N1)
        call basis_bspline_der(map_parameter(Xt(2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
            knot2, nc(2), degree(2), 0, first2, N2)
        call basis_bspline_der(map_parameter(Xt(3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
            knot3, nc(3), degree(3), 0, first3, N3)
        j10 = max(0, 1 - first1)
        j11 = min(degree(1), nc(1) - first1)
        j20 = max(0, 1 - first2)
        j21 = min(degree(2), nc(2) - first2)
        j30 = max(0, 1 - first3)
        j31 = min(degree(3), nc(3) - first3)
        idx = ((first3+j30-1)*nc(2) + first2+j20-1)*nc(1) + first1+j10
        weight_exponent = exponent(Wc(idx))
        do j3 = j30, j31
            do j2 = j20, j21
                do j1 = j10, j11
                    idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1+j1
                    weight_exponent = max(weight_exponent, exponent(Wc(idx)))
                end do
            end do
        end do
        wsum = 0.0_rk
        do j3 = j30, j31
            do j2 = j20, j21
                n32 = N3(j3)*N2(j2)
                do j1 = j10, j11
                    idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1+j1
                    wi = scale(Wc(idx), -weight_exponent)
                    bw = n32*N1(j1)*wi
                    wsum = wsum + bw
                    do d = 1, dim_
                        Xg(d) = Xg(d) + bw*Xc(idx,d)
                    end do
                end do
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
    !> Evaluate a polynomial B-spline volume at a list of parameter triples.
    !!
    !! Evaluation accumulates only
    !! `(degree(1)+1)*(degree(2)+1)*(degree(3)+1)` local controls per point and
    !! avoids constructing a dense basis matrix.
    pure function compute_Xg_bspline_3d(Xt, knot1, knot2, knot3, degree, nc, ng, Xc, wrap_parameters) result(Xg)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: ng(3)
        real(rk), intent(in), contiguous :: Xc(:,:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Xg(:,:)
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3)), n32
        integer :: i, j1, j2, j3, d, idx, first1, first2, first3, ng_, dim_
        logical :: wrap_parameters_(3)

        if (present(ng)) then
            ng_ = product(ng)
        else
            ng_ = size(Xt, 1)
        end if

        dim_ = size(Xc,2)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        allocate(Xg(ng_, dim_))
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N1, N2, N3, n32, j1, j2, j3, d, idx, first1, first2, first3)
#endif
            do d = 1, dim_
                Xg(i,d) = 0.0_rk
            end do
            call basis_bspline_der(map_parameter(Xt(i,1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
                knot1, nc(1), degree(1), 0, first1, N1)
            call basis_bspline_der(map_parameter(Xt(i,2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
                knot2, nc(2), degree(2), 0, first2, N2)
            call basis_bspline_der(map_parameter(Xt(i,3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
                knot3, nc(3), degree(3), 0, first3, N3)
            do j3 = max(0, 1 - first3), min(degree(3), nc(3) - first3)
                do j2 = max(0, 1 - first2), min(degree(2), nc(2) - first2)
                    n32 = N3(j3)*N2(j2)
                    do j1 = max(0, 1 - first1), min(degree(1), nc(1) - first1)
                        idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1+j1
                        do d = 1, dim_
                            Xg(i,d) = Xg(i,d) + n32*N1(j1)*Xc(idx,d)
                        end do
                    end do
                end do
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate a polynomial B-spline volume at one parameter triple.
    !!
    !! `Xt(1:3)` contains `[u,v,w]`; the result has length `size(Xc,2)` and uses
    !! direction-1-fastest control-lattice indexing.
    pure function compute_Xg_bspline_3d_1point(Xt, knot1, knot2, knot3, degree, nc, Xc, wrap_parameters) result(Xg)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        real(rk), intent(in), contiguous :: Xc(:,:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk) :: Xg(size(Xc,2))
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3)), n32
        integer :: j1, j2, j3, d, idx, first1, first2, first3, dim_
        logical :: wrap_parameters_(3)

        Xg = 0.0_rk
        dim_ = size(Xc,2)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        call basis_bspline_der(map_parameter(Xt(1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
            knot1, nc(1), degree(1), 0, first1, N1)
        call basis_bspline_der(map_parameter(Xt(2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
            knot2, nc(2), degree(2), 0, first2, N2)
        call basis_bspline_der(map_parameter(Xt(3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
            knot3, nc(3), degree(3), 0, first3, N3)
        do j3 = max(0, 1 - first3), min(degree(3), nc(3) - first3)
            do j2 = max(0, 1 - first2), min(degree(2), nc(2) - first2)
                n32 = N3(j3)*N2(j2)
                do j1 = max(0, 1 - first1), min(degree(1), nc(1) - first1)
                    idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1+j1
                    do d = 1, dim_
                        Xg(d) = Xg(d) + n32*N1(j1)*Xc(idx,d)
                    end do
                end do
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate rational volume basis values and gradients at many parameter triples.
    !!
    !! `Tgc` has shape `[ng_,ncp]` and `dTgc` shape `[ng_,ncp,3]`, where
    !! `ncp=product(nc)`. The final derivative index follows `[u,v,w]` and all
    !! quotient-rule terms share the same local polynomial basis pass.
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
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3))
        real(rk) :: dN1(0:degree(1)), dN2(0:degree(2)), dN3(0:degree(3))
        real(rk) :: denom, dden1, dden2, dden3, b, db1, db2, db3, wi, tval, wtol
        integer :: i, j1, j2, j3, idx1, idx2, idx3, idx, first1, first2, first3, n21, ncp, ng_
        integer :: j10, j11, j20, j21, j30, j31, weight_exponent
        logical :: wrap_parameters_(3)

        if (present(ng)) then
            ng_ = product(ng)
        else
            ng_ = size(Xt, 1)
        end if

        n21 = nc(1)*nc(2)
        ncp = n21*nc(3)
        allocate(dTgc(ng_, ncp, 3), Tgc(ng_, ncp), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        wtol = 0.0_rk
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N1, N2, N3, dN1, dN2, dN3, denom, dden1, dden2, dden3, &
            b, db1, db2, db3, wi, tval, j1, j2, j3, idx1, idx2, idx3, idx, first1, first2, first3, &
            j10, j11, j20, j21, j30, j31, weight_exponent)
#endif
            call basis_bspline_der(map_parameter(Xt(i,1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
                knot1, nc(1), degree(1), 1, first1, N1, dN1)
            call basis_bspline_der(map_parameter(Xt(i,2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
                knot2, nc(2), degree(2), 1, first2, N2, dN2)
            call basis_bspline_der(map_parameter(Xt(i,3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
                knot3, nc(3), degree(3), 1, first3, N3, dN3)
            j10 = max(0, 1 - first1)
            j11 = min(degree(1), nc(1) - first1)
            j20 = max(0, 1 - first2)
            j21 = min(degree(2), nc(2) - first2)
            j30 = max(0, 1 - first3)
            j31 = min(degree(3), nc(3) - first3)
            idx = ((first3+j30-1)*nc(2) + first2+j20-1)*nc(1) + first1+j10
            weight_exponent = exponent(Wc(idx))
            do j3 = j30, j31
                do j2 = j20, j21
                    do j1 = j10, j11
                        idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1+j1
                        weight_exponent = max(weight_exponent, exponent(Wc(idx)))
                    end do
                end do
            end do
            denom = 0.0_rk
            dden1 = 0.0_rk
            dden2 = 0.0_rk
            dden3 = 0.0_rk
            do j3 = j30, j31
                idx3 = first3 + j3
                do j2 = j20, j21
                    idx2 = first2 + j2
                    do j1 = j10, j11
                        idx1 = first1 + j1
                        idx = (idx3 - 1)*n21 + (idx2 - 1)*nc(1) + idx1
                        wi = scale(Wc(idx), -weight_exponent)
                        b = N3(j3)*N2(j2)*N1(j1)
                        db1 = N3(j3)*N2(j2)*dN1(j1)
                        db2 = N3(j3)*dN2(j2)*N1(j1)
                        db3 = dN3(j3)*N2(j2)*N1(j1)
                        denom = denom + b*wi
                        dden1 = dden1 + db1*wi
                        dden2 = dden2 + db2*wi
                        dden3 = dden3 + db3*wi
                    end do
                end do
            end do
            if (abs(denom) > wtol) then
                do j3 = j30, j31
                    idx3 = first3 + j3
                    do j2 = j20, j21
                        idx2 = first2 + j2
                        do j1 = j10, j11
                            idx1 = first1 + j1
                            idx = (idx3 - 1)*n21 + (idx2 - 1)*nc(1) + idx1
                            wi = scale(Wc(idx), -weight_exponent)
                            b = N3(j3)*N2(j2)*N1(j1)
                            db1 = N3(j3)*N2(j2)*dN1(j1)
                            db2 = N3(j3)*dN2(j2)*N1(j1)
                            db3 = dN3(j3)*N2(j2)*N1(j1)
                            tval = b*wi/denom
                            Tgc(i,idx) = tval
                            dTgc(i,idx,1) = (db1*wi - tval*dden1)/denom
                            dTgc(i,idx,2) = (db2*wi - tval*dden2)/denom
                            dTgc(i,idx,3) = (db3*wi - tval*dden3)/denom
                        end do
                    end do
                end do
            end if
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate rational volume basis values and gradients at one parameter triple.
    !!
    !! Dense output uses `ncp=product(nc)` rows. With `elem`, output order
    !! follows that one-based control-index list and `Wc` may be either the full
    !! flattened vector or an element-local vector of length `size(elem)`.
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
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3))
        real(rk) :: dN1(0:degree(1)), dN2(0:degree(2)), dN3(0:degree(3))
        real(rk) :: bval, b32, db1v, db2v, db3v, denom, ddenom1, ddenom2, ddenom3, wi, wtol
        integer :: a, i1, i2, i3, idx, idx1, idx2, idx3, flat, ncp, nloc, n12, first1, first2, first3, l1, l2, l3
        integer :: l10, l11, l20, l21, l30, l31, weight_exponent
        logical :: local_weights
        logical :: wrap_parameters_(3)

        n12 = nc(1)*nc(2)
        ncp = n12*nc(3)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        wtol = 0.0_rk
        if (.not. present(elem)) then
            allocate(dTgc(ncp, 3), Tgc(ncp), source=0.0_rk)
            call basis_bspline_der(map_parameter(Xt(1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
                knot1, nc(1), degree(1), 1, first1, N1, dN1)
            call basis_bspline_der(map_parameter(Xt(2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
                knot2, nc(2), degree(2), 1, first2, N2, dN2)
            call basis_bspline_der(map_parameter(Xt(3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
                knot3, nc(3), degree(3), 1, first3, N3, dN3)
            l10 = max(0, 1 - first1)
            l11 = min(degree(1), nc(1) - first1)
            l20 = max(0, 1 - first2)
            l21 = min(degree(2), nc(2) - first2)
            l30 = max(0, 1 - first3)
            l31 = min(degree(3), nc(3) - first3)
            idx = ((first3+l30-1)*nc(2) + first2+l20-1)*nc(1) + first1+l10
            weight_exponent = exponent(Wc(idx))
            do l3 = l30, l31
                do l2 = l20, l21
                    do l1 = l10, l11
                        idx = ((first3+l3-1)*nc(2) + first2+l2-1)*nc(1) + first1+l1
                        weight_exponent = max(weight_exponent, exponent(Wc(idx)))
                    end do
                end do
            end do
            denom = 0.0_rk
            ddenom1 = 0.0_rk
            ddenom2 = 0.0_rk
            ddenom3 = 0.0_rk
            do l3 = l30, l31
                idx3 = first3 + l3
                do l2 = l20, l21
                    idx2 = first2 + l2
                    b32 = N3(l3)*N2(l2)
                    do l1 = l10, l11
                        idx1 = first1 + l1
                        idx = (idx3 - 1)*n12 + (idx2 - 1)*nc(1) + idx1
                        bval = b32*N1(l1)
                        db1v = b32*dN1(l1)
                        db2v = N3(l3)*dN2(l2)*N1(l1)
                        db3v = dN3(l3)*N2(l2)*N1(l1)
                        wi = scale(Wc(idx), -weight_exponent)
                        denom = denom + bval*wi
                        ddenom1 = ddenom1 + db1v*wi
                        ddenom2 = ddenom2 + db2v*wi
                        ddenom3 = ddenom3 + db3v*wi
                    end do
                end do
            end do
            if (abs(denom) > wtol) then
                do l3 = l30, l31
                    idx3 = first3 + l3
                    do l2 = l20, l21
                        idx2 = first2 + l2
                        b32 = N3(l3)*N2(l2)
                        do l1 = l10, l11
                            idx1 = first1 + l1
                            idx = (idx3 - 1)*n12 + (idx2 - 1)*nc(1) + idx1
                            bval = b32*N1(l1)
                            db1v = b32*dN1(l1)
                            db2v = N3(l3)*dN2(l2)*N1(l1)
                            db3v = dN3(l3)*N2(l2)*N1(l1)
                            wi = scale(Wc(idx), -weight_exponent)
                            Tgc(idx) = bval*wi/denom
                            dTgc(idx,1) = (db1v*wi - Tgc(idx)*ddenom1)/denom
                            dTgc(idx,2) = (db2v*wi - Tgc(idx)*ddenom2)/denom
                            dTgc(idx,3) = (db3v*wi - Tgc(idx)*ddenom3)/denom
                        end do
                    end do
                end do
            end if
        else
            nloc = size(elem)
            allocate(dTgc(nloc, 3), source=0.0_rk)
            allocate(Tgc(nloc), source=0.0_rk)
            call basis_bspline_der(map_parameter(Xt(1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
                knot1, nc(1), degree(1), 1, first1, N1, dN1)
            call basis_bspline_der(map_parameter(Xt(2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
                knot2, nc(2), degree(2), 1, first2, N2, dN2)
            call basis_bspline_der(map_parameter(Xt(3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
                knot3, nc(3), degree(3), 1, first3, N3, dN3)
            denom = 0.0_rk
            ddenom1 = 0.0_rk
            ddenom2 = 0.0_rk
            ddenom3 = 0.0_rk
            local_weights = size(Wc) == nloc .and. size(Wc) /= ncp
            if (local_weights) then
                weight_exponent = exponent(maxval(Wc))
                do a = 1, nloc
                    if (elem(a) < 1 .or. elem(a) > ncp) cycle
                    flat = elem(a) - 1
                    i1 = mod(flat, nc(1)) + 1
                    i2 = mod(flat/nc(1), nc(2)) + 1
                    i3 = flat/n12 + 1
                    l1 = i1 - first1
                    l2 = i2 - first2
                    l3 = i3 - first3
                    if (l1 < 0 .or. l1 > degree(1) .or. l2 < 0 .or. l2 > degree(2) .or. &
                        l3 < 0 .or. l3 > degree(3)) cycle
                    b32 = N3(l3)*N2(l2)
                    bval = b32*N1(l1)
                    db1v = b32*dN1(l1)
                    db2v = N3(l3)*dN2(l2)*N1(l1)
                    db3v = dN3(l3)*N2(l2)*N1(l1)
                    wi = scale(Wc(a), -weight_exponent)
                    denom = denom + bval*wi
                    ddenom1 = ddenom1 + db1v*wi
                    ddenom2 = ddenom2 + db2v*wi
                    ddenom3 = ddenom3 + db3v*wi
                end do
            else
                l10 = max(0, 1 - first1)
                l11 = min(degree(1), nc(1) - first1)
                l20 = max(0, 1 - first2)
                l21 = min(degree(2), nc(2) - first2)
                l30 = max(0, 1 - first3)
                l31 = min(degree(3), nc(3) - first3)
                idx = ((first3+l30-1)*nc(2) + first2+l20-1)*nc(1) + first1+l10
                weight_exponent = exponent(Wc(idx))
                do l3 = l30, l31
                    do l2 = l20, l21
                        do l1 = l10, l11
                            idx = ((first3+l3-1)*nc(2) + first2+l2-1)*nc(1) + first1+l1
                            weight_exponent = max(weight_exponent, exponent(Wc(idx)))
                        end do
                    end do
                end do
                do l3 = l30, l31
                    idx3 = first3 + l3
                    do l2 = l20, l21
                        idx2 = first2 + l2
                        b32 = N3(l3)*N2(l2)
                        do l1 = l10, l11
                            idx1 = first1 + l1
                            idx = (idx3 - 1)*n12 + (idx2 - 1)*nc(1) + idx1
                            bval = b32*N1(l1)
                            db1v = b32*dN1(l1)
                            db2v = N3(l3)*dN2(l2)*N1(l1)
                            db3v = dN3(l3)*N2(l2)*N1(l1)
                            wi = scale(Wc(idx), -weight_exponent)
                            denom = denom + bval*wi
                            ddenom1 = ddenom1 + db1v*wi
                            ddenom2 = ddenom2 + db2v*wi
                            ddenom3 = ddenom3 + db3v*wi
                        end do
                    end do
                end do
            end if
            if (abs(denom) > wtol) then
                do a = 1, nloc
                    if (elem(a) < 1 .or. elem(a) > ncp) cycle
                    flat = elem(a) - 1
                    i1 = mod(flat, nc(1)) + 1
                    i2 = mod(flat/nc(1), nc(2)) + 1
                    i3 = flat/n12 + 1
                    l1 = i1 - first1
                    l2 = i2 - first2
                    l3 = i3 - first3
                    if (local_weights) then
                        wi = scale(Wc(a), -weight_exponent)
                    else
                        wi = scale(Wc(elem(a)), -weight_exponent)
                    end if
                    if (l1 >= 0 .and. l1 <= degree(1) .and. l2 >= 0 .and. l2 <= degree(2) .and. &
                        l3 >= 0 .and. l3 <= degree(3)) then
                        b32 = N3(l3)*N2(l2)
                        bval = b32*N1(l1)
                        db1v = b32*dN1(l1)
                        db2v = N3(l3)*dN2(l2)*N1(l1)
                        db3v = dN3(l3)*N2(l2)*N1(l1)
                        Tgc(a) = bval*wi/denom
                        dTgc(a,1) = (db1v*wi - Tgc(a)*ddenom1)/denom
                        dTgc(a,2) = (db2v*wi - Tgc(a)*ddenom2)/denom
                        dTgc(a,3) = (db3v*wi - Tgc(a)*ddenom3)/denom
                    end if
                end do
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate polynomial volume basis values and gradients at many parameter triples.
    !!
    !! Outputs have shapes `[ng_,ncp]` and `[ng_,ncp,3]`. Entries outside each
    !! point's local tensor-product support remain zero.
    pure subroutine compute_dTgc_bspline_3d_vector(Xt, knot1, knot2, knot3, degree, nc, ng, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: ng(3)
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3))
        real(rk) :: dN1(0:degree(1)), dN2(0:degree(2)), dN3(0:degree(3))
        real(rk) :: b32
        integer :: i, j1, j2, j3, idx1, idx2, idx3, idx, first1, first2, first3, n21, ncp, ng_
        logical :: wrap_parameters_(3)

        if (present(ng)) then
            ng_ = product(ng)
        else
            ng_ = size(Xt, 1)
        end if

        n21 = nc(1)*nc(2)
        ncp = n21*nc(3)
        allocate(dTgc(ng_, ncp, 3), Tgc(ng_, ncp), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if

#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N1, N2, N3, dN1, dN2, dN3, b32, &
            j1, j2, j3, idx1, idx2, idx3, idx, first1, first2, first3)
#endif
            call basis_bspline_der(map_parameter(Xt(i,1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
                knot1, nc(1), degree(1), 1, first1, N1, dN1)
            call basis_bspline_der(map_parameter(Xt(i,2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
                knot2, nc(2), degree(2), 1, first2, N2, dN2)
            call basis_bspline_der(map_parameter(Xt(i,3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
                knot3, nc(3), degree(3), 1, first3, N3, dN3)
            do j3 = max(0, 1 - first3), min(degree(3), nc(3) - first3)
                idx3 = first3 + j3
                do j2 = max(0, 1 - first2), min(degree(2), nc(2) - first2)
                    idx2 = first2 + j2
                    b32 = N3(j3)*N2(j2)
                    do j1 = max(0, 1 - first1), min(degree(1), nc(1) - first1)
                        idx1 = first1 + j1
                        idx = (idx3 - 1)*n21 + (idx2 - 1)*nc(1) + idx1
                        Tgc(i,idx) = b32*N1(j1)
                        dTgc(i,idx,1) = b32*dN1(j1)
                        dTgc(i,idx,2) = N3(j3)*dN2(j2)*N1(j1)
                        dTgc(i,idx,3) = dN3(j3)*N2(j2)*N1(j1)
                    end do
                end do
            end do
        end do

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate polynomial volume basis values and gradients at one parameter triple.
    !!
    !! Without `elem`, outputs use flattened global control order. With `elem`,
    !! they follow the supplied one-based element order; inactive or invalid
    !! listed controls remain zero. Gradient columns correspond to `[u,v,w]`.
    pure subroutine compute_dTgc_bspline_3d_scalar(Xt, knot1, knot2, knot3, degree, nc, dTgc, Tgc, elem, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3))
        real(rk) :: dN1(0:degree(1)), dN2(0:degree(2)), dN3(0:degree(3))
        real(rk) :: b32
        integer :: a, i1, i2, i3, idx, idx1, idx2, idx3, flat, nloc, n12, first1, first2, first3, l1, l2, l3
        logical :: wrap_parameters_(3)

        n12 = nc(1)*nc(2)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        if (.not. present(elem)) then
            allocate(dTgc(nc(1)*nc(2)*nc(3), 3), Tgc(nc(1)*nc(2)*nc(3)), source=0.0_rk)
            call basis_bspline_der(map_parameter(Xt(1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
                knot1, nc(1), degree(1), 1, first1, N1, dN1)
            call basis_bspline_der(map_parameter(Xt(2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
                knot2, nc(2), degree(2), 1, first2, N2, dN2)
            call basis_bspline_der(map_parameter(Xt(3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
                knot3, nc(3), degree(3), 1, first3, N3, dN3)

            do l3 = max(0, 1 - first3), min(degree(3), nc(3) - first3)
                idx3 = first3 + l3
                do l2 = max(0, 1 - first2), min(degree(2), nc(2) - first2)
                    idx2 = first2 + l2
                    b32 = N3(l3)*N2(l2)
                    do l1 = max(0, 1 - first1), min(degree(1), nc(1) - first1)
                        idx1 = first1 + l1
                        idx = (idx3 - 1)*n12 + (idx2 - 1)*nc(1) + idx1
                        Tgc(idx) = b32*N1(l1)
                        dTgc(idx,1) = b32*dN1(l1)
                        dTgc(idx,2) = N3(l3)*dN2(l2)*N1(l1)
                        dTgc(idx,3) = dN3(l3)*N2(l2)*N1(l1)
                    end do
                end do
            end do
        else
            nloc = size(elem)
            allocate(dTgc(nloc, 3), Tgc(nloc), source=0.0_rk)
            call basis_bspline_der(map_parameter(Xt(1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
                knot1, nc(1), degree(1), 1, first1, N1, dN1)
            call basis_bspline_der(map_parameter(Xt(2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
                knot2, nc(2), degree(2), 1, first2, N2, dN2)
            call basis_bspline_der(map_parameter(Xt(3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
                knot3, nc(3), degree(3), 1, first3, N3, dN3)
            do a = 1, nloc
                if (elem(a) < 1 .or. elem(a) > nc(1)*nc(2)*nc(3)) cycle
                flat = elem(a) - 1
                i1 = mod(flat, nc(1)) + 1
                i2 = mod(flat/nc(1), nc(2)) + 1
                i3 = flat/n12 + 1
                l1 = i1 - first1
                l2 = i2 - first2
                l3 = i3 - first3
                if (l1 < 0 .or. l1 > degree(1) .or. l2 < 0 .or. l2 > degree(2) .or. &
                    l3 < 0 .or. l3 > degree(3)) cycle
                b32 = N3(l3)*N2(l2)
                Tgc(a) = b32*N1(l1)
                dTgc(a,1) = b32*dN1(l1)
                dTgc(a,2) = N3(l3)*dN2(l2)*N1(l1)
                dTgc(a,3) = dN3(l3)*N2(l2)*N1(l1)
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate rational volume basis values, gradients, and Hessians at many triples.
    !!
    !! For `ncp=product(nc)`, Hessian component `(a,b)` of basis `A` is stored
    !! at `d2Tgc(g,(a-1)*ncp+A,b)`. Quotient-rule normalization includes all
    !! first and second derivatives of the rational denominator.
    pure subroutine compute_d2Tgc_nurbs_3d_vector(Xt, knot1, knot2, knot3, degree, nc, ng, Wc, d2Tgc, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: ng(3)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:,:)
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3))
        real(rk) :: dN1(0:degree(1)), dN2(0:degree(2)), dN3(0:degree(3))
        real(rk) :: d2N1(0:degree(1)), d2N2(0:degree(2)), d2N3(0:degree(3))
        real(rk) :: denom, dden1, dden2, dden3, d2den11, d2den12, d2den13, d2den22, d2den23, d2den33
        real(rk) :: b, db1, db2, db3, d11, d12, d13, d22, d23, d33, wi, tval, dt1, dt2, dt3, wtol
        integer :: i, j1, j2, j3, idx1, idx2, idx3, idx, first1, first2, first3, n21, ncp, ng_
        integer :: j10, j11, j20, j21, j30, j31, weight_exponent
        logical :: wrap_parameters_(3)

        if (present(ng)) then
            ng_ = product(ng)
        else
            ng_ = size(Xt, 1)
        end if

        n21 = nc(1)*nc(2)
        ncp = n21*nc(3)
        allocate(d2Tgc(ng_, 3*ncp, 3), source=0.0_rk)
        allocate(dTgc(ng_, ncp, 3), Tgc(ng_, ncp), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        wtol = 0.0_rk
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N1, N2, N3, dN1, dN2, dN3, d2N1, d2N2, d2N3, denom, &
            dden1, dden2, dden3, d2den11, d2den12, d2den13, d2den22, d2den23, d2den33, &
            b, db1, db2, db3, d11, d12, d13, d22, d23, d33, wi, tval, dt1, dt2, dt3, &
            j1, j2, j3, idx1, idx2, idx3, idx, first1, first2, first3, &
            j10, j11, j20, j21, j30, j31, weight_exponent)
#endif
            call basis_bspline_der(map_parameter(Xt(i,1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
                knot1, nc(1), degree(1), 2, first1, N1, dN1, d2N1)
            call basis_bspline_der(map_parameter(Xt(i,2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
                knot2, nc(2), degree(2), 2, first2, N2, dN2, d2N2)
            call basis_bspline_der(map_parameter(Xt(i,3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
                knot3, nc(3), degree(3), 2, first3, N3, dN3, d2N3)
            j10 = max(0, 1 - first1)
            j11 = min(degree(1), nc(1) - first1)
            j20 = max(0, 1 - first2)
            j21 = min(degree(2), nc(2) - first2)
            j30 = max(0, 1 - first3)
            j31 = min(degree(3), nc(3) - first3)
            idx = ((first3+j30-1)*nc(2) + first2+j20-1)*nc(1) + first1+j10
            weight_exponent = exponent(Wc(idx))
            do j3 = j30, j31
                do j2 = j20, j21
                    do j1 = j10, j11
                        idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1+j1
                        weight_exponent = max(weight_exponent, exponent(Wc(idx)))
                    end do
                end do
            end do
            denom = 0.0_rk
            dden1 = 0.0_rk
            dden2 = 0.0_rk
            dden3 = 0.0_rk
            d2den11 = 0.0_rk
            d2den12 = 0.0_rk
            d2den13 = 0.0_rk
            d2den22 = 0.0_rk
            d2den23 = 0.0_rk
            d2den33 = 0.0_rk
            do j3 = j30, j31
                idx3 = first3 + j3
                do j2 = j20, j21
                    idx2 = first2 + j2
                    do j1 = j10, j11
                        idx1 = first1 + j1
                        idx = (idx3 - 1)*n21 + (idx2 - 1)*nc(1) + idx1
                        wi = scale(Wc(idx), -weight_exponent)
                        b = N3(j3)*N2(j2)*N1(j1)
                        db1 = N3(j3)*N2(j2)*dN1(j1)
                        db2 = N3(j3)*dN2(j2)*N1(j1)
                        db3 = dN3(j3)*N2(j2)*N1(j1)
                        d11 = N3(j3)*N2(j2)*d2N1(j1)
                        d12 = N3(j3)*dN2(j2)*dN1(j1)
                        d13 = dN3(j3)*N2(j2)*dN1(j1)
                        d22 = N3(j3)*d2N2(j2)*N1(j1)
                        d23 = dN3(j3)*dN2(j2)*N1(j1)
                        d33 = d2N3(j3)*N2(j2)*N1(j1)
                        denom = denom + b*wi
                        dden1 = dden1 + db1*wi
                        dden2 = dden2 + db2*wi
                        dden3 = dden3 + db3*wi
                        d2den11 = d2den11 + d11*wi
                        d2den12 = d2den12 + d12*wi
                        d2den13 = d2den13 + d13*wi
                        d2den22 = d2den22 + d22*wi
                        d2den23 = d2den23 + d23*wi
                        d2den33 = d2den33 + d33*wi
                    end do
                end do
            end do
            if (abs(denom) > wtol) then
                do j3 = j30, j31
                    idx3 = first3 + j3
                    do j2 = j20, j21
                        idx2 = first2 + j2
                        do j1 = j10, j11
                            idx1 = first1 + j1
                            idx = (idx3 - 1)*n21 + (idx2 - 1)*nc(1) + idx1
                            wi = scale(Wc(idx), -weight_exponent)
                            b = N3(j3)*N2(j2)*N1(j1)
                            db1 = N3(j3)*N2(j2)*dN1(j1)
                            db2 = N3(j3)*dN2(j2)*N1(j1)
                            db3 = dN3(j3)*N2(j2)*N1(j1)
                            d11 = N3(j3)*N2(j2)*d2N1(j1)
                            d12 = N3(j3)*dN2(j2)*dN1(j1)
                            d13 = dN3(j3)*N2(j2)*dN1(j1)
                            d22 = N3(j3)*d2N2(j2)*N1(j1)
                            d23 = dN3(j3)*dN2(j2)*N1(j1)
                            d33 = d2N3(j3)*N2(j2)*N1(j1)
                            tval = b*wi/denom
                            dt1 = (db1*wi - tval*dden1)/denom
                            dt2 = (db2*wi - tval*dden2)/denom
                            dt3 = (db3*wi - tval*dden3)/denom
                            Tgc(i,idx) = tval
                            dTgc(i,idx,1) = dt1
                            dTgc(i,idx,2) = dt2
                            dTgc(i,idx,3) = dt3
                            d2Tgc(i,idx,1) = (d11*wi - 2.0_rk*dt1*dden1 - tval*d2den11)/denom
                            d2Tgc(i,ncp+idx,1) = (d12*wi - dt1*dden2 - dt2*dden1 - tval*d2den12)/denom
                            d2Tgc(i,2*ncp+idx,1) = (d13*wi - dt1*dden3 - dt3*dden1 - tval*d2den13)/denom
                            d2Tgc(i,idx,2) = d2Tgc(i,ncp+idx,1)
                            d2Tgc(i,ncp+idx,2) = (d22*wi - 2.0_rk*dt2*dden2 - tval*d2den22)/denom
                            d2Tgc(i,2*ncp+idx,2) = (d23*wi - dt2*dden3 - dt3*dden2 - tval*d2den23)/denom
                            d2Tgc(i,idx,3) = d2Tgc(i,2*ncp+idx,1)
                            d2Tgc(i,ncp+idx,3) = d2Tgc(i,2*ncp+idx,2)
                            d2Tgc(i,2*ncp+idx,3) = (d33*wi - 2.0_rk*dt3*dden3 - tval*d2den33)/denom
                        end do
                    end do
                end do
            end if
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate one dense rational volume basis, gradient, and Hessian.
    !!
    !! `Tgc` has length `ncp`, `dTgc` shape `[ncp,3]`, and `d2Tgc` shape
    !! `[3*ncp,3]`. Packed Hessian indexing is
    !! `d2Tgc((a-1)*ncp+A,b)` for parameter directions `a,b`.
    pure subroutine compute_d2Tgc_nurbs_3d_scalar(Xt, knot1, knot2, knot3, degree, nc, Wc, d2Tgc, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        real(rk), intent(in), contiguous :: Wc(:)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3))
        real(rk) :: dN1(0:degree(1)), dN2(0:degree(2)), dN3(0:degree(3))
        real(rk) :: d2N1(0:degree(1)), d2N2(0:degree(2)), d2N3(0:degree(3))
        real(rk) :: denom, dden1, dden2, dden3, d2den11, d2den12, d2den13, d2den22, d2den23, d2den33
        real(rk) :: b, db1, db2, db3, d11, d12, d13, d22, d23, d33, wi, tval, dt1, dt2, dt3, wtol
        integer :: j1, j2, j3, idx1, idx2, idx3, idx, first1, first2, first3, n21, ncp
        integer :: j10, j11, j20, j21, j30, j31, weight_exponent
        logical :: wrap_parameters_(3)

        n21 = nc(1)*nc(2)
        ncp = n21*nc(3)
        allocate(d2Tgc(3*ncp, 3), dTgc(ncp, 3), Tgc(ncp), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        wtol = 0.0_rk
        call basis_bspline_der(map_parameter(Xt(1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
            knot1, nc(1), degree(1), 2, first1, N1, dN1, d2N1)
        call basis_bspline_der(map_parameter(Xt(2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
            knot2, nc(2), degree(2), 2, first2, N2, dN2, d2N2)
        call basis_bspline_der(map_parameter(Xt(3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
            knot3, nc(3), degree(3), 2, first3, N3, dN3, d2N3)
        j10 = max(0, 1 - first1)
        j11 = min(degree(1), nc(1) - first1)
        j20 = max(0, 1 - first2)
        j21 = min(degree(2), nc(2) - first2)
        j30 = max(0, 1 - first3)
        j31 = min(degree(3), nc(3) - first3)
        idx = ((first3+j30-1)*nc(2) + first2+j20-1)*nc(1) + first1+j10
        weight_exponent = exponent(Wc(idx))
        do j3 = j30, j31
            do j2 = j20, j21
                do j1 = j10, j11
                    idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1+j1
                    weight_exponent = max(weight_exponent, exponent(Wc(idx)))
                end do
            end do
        end do
        denom = 0.0_rk
        dden1 = 0.0_rk
        dden2 = 0.0_rk
        dden3 = 0.0_rk
        d2den11 = 0.0_rk
        d2den12 = 0.0_rk
        d2den13 = 0.0_rk
        d2den22 = 0.0_rk
        d2den23 = 0.0_rk
        d2den33 = 0.0_rk
        do j3 = j30, j31
            idx3 = first3 + j3
            do j2 = j20, j21
                idx2 = first2 + j2
                do j1 = j10, j11
                    idx1 = first1 + j1
                    idx = (idx3 - 1)*n21 + (idx2 - 1)*nc(1) + idx1
                    wi = scale(Wc(idx), -weight_exponent)
                    b = N3(j3)*N2(j2)*N1(j1)
                    db1 = N3(j3)*N2(j2)*dN1(j1)
                    db2 = N3(j3)*dN2(j2)*N1(j1)
                    db3 = dN3(j3)*N2(j2)*N1(j1)
                    d11 = N3(j3)*N2(j2)*d2N1(j1)
                    d12 = N3(j3)*dN2(j2)*dN1(j1)
                    d13 = dN3(j3)*N2(j2)*dN1(j1)
                    d22 = N3(j3)*d2N2(j2)*N1(j1)
                    d23 = dN3(j3)*dN2(j2)*N1(j1)
                    d33 = d2N3(j3)*N2(j2)*N1(j1)
                    denom = denom + b*wi
                    dden1 = dden1 + db1*wi
                    dden2 = dden2 + db2*wi
                    dden3 = dden3 + db3*wi
                    d2den11 = d2den11 + d11*wi
                    d2den12 = d2den12 + d12*wi
                    d2den13 = d2den13 + d13*wi
                    d2den22 = d2den22 + d22*wi
                    d2den23 = d2den23 + d23*wi
                    d2den33 = d2den33 + d33*wi
                end do
            end do
        end do
        if (abs(denom) > wtol) then
            do j3 = j30, j31
                idx3 = first3 + j3
                do j2 = j20, j21
                    idx2 = first2 + j2
                    do j1 = j10, j11
                        idx1 = first1 + j1
                        idx = (idx3 - 1)*n21 + (idx2 - 1)*nc(1) + idx1
                        wi = scale(Wc(idx), -weight_exponent)
                        b = N3(j3)*N2(j2)*N1(j1)
                        db1 = N3(j3)*N2(j2)*dN1(j1)
                        db2 = N3(j3)*dN2(j2)*N1(j1)
                        db3 = dN3(j3)*N2(j2)*N1(j1)
                        d11 = N3(j3)*N2(j2)*d2N1(j1)
                        d12 = N3(j3)*dN2(j2)*dN1(j1)
                        d13 = dN3(j3)*N2(j2)*dN1(j1)
                        d22 = N3(j3)*d2N2(j2)*N1(j1)
                        d23 = dN3(j3)*dN2(j2)*N1(j1)
                        d33 = d2N3(j3)*N2(j2)*N1(j1)
                        tval = b*wi/denom
                        dt1 = (db1*wi - tval*dden1)/denom
                        dt2 = (db2*wi - tval*dden2)/denom
                        dt3 = (db3*wi - tval*dden3)/denom
                        Tgc(idx) = tval
                        dTgc(idx,1) = dt1
                        dTgc(idx,2) = dt2
                        dTgc(idx,3) = dt3
                        d2Tgc(idx,1) = (d11*wi - 2.0_rk*dt1*dden1 - tval*d2den11)/denom
                        d2Tgc(ncp+idx,1) = (d12*wi - dt1*dden2 - dt2*dden1 - tval*d2den12)/denom
                        d2Tgc(2*ncp+idx,1) = (d13*wi - dt1*dden3 - dt3*dden1 - tval*d2den13)/denom
                        d2Tgc(idx,2) = d2Tgc(ncp+idx,1)
                        d2Tgc(ncp+idx,2) = (d22*wi - 2.0_rk*dt2*dden2 - tval*d2den22)/denom
                        d2Tgc(2*ncp+idx,2) = (d23*wi - dt2*dden3 - dt3*dden2 - tval*d2den23)/denom
                        d2Tgc(idx,3) = d2Tgc(2*ncp+idx,1)
                        d2Tgc(ncp+idx,3) = d2Tgc(2*ncp+idx,2)
                        d2Tgc(2*ncp+idx,3) = (d33*wi - 2.0_rk*dt3*dden3 - tval*d2den33)/denom
                    end do
                end do
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate polynomial volume basis values, gradients, and Hessians at many triples.
    !!
    !! The Hessian uses shape `[ng_,3*ncp,3]` and packed indexing
    !! `d2Tgc(g,(a-1)*ncp+A,b)`. Mixed derivatives are written to both symmetric
    !! positions; values outside local support remain zero.
    pure subroutine compute_d2Tgc_bspline_3d_vector(Xt, knot1, knot2, knot3, degree, nc, ng, d2Tgc, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: ng(3)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:,:)
        real(rk), allocatable, intent(out) :: dTgc(:,:,:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        integer :: i, j1, j2, j3, idx1, idx2, idx3, idx, first1, first2, first3, n21, ncp, ng_
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3))
        real(rk) :: dN1(0:degree(1)), dN2(0:degree(2)), dN3(0:degree(3))
        real(rk) :: d2N1(0:degree(1)), d2N2(0:degree(2)), d2N3(0:degree(3))
        real(rk) :: n32, dn21, dn31, dn32
        logical :: wrap_parameters_(3)

        if (present(ng)) then
            ng_ = product(ng)
        else
            ng_ = size(Xt, 1)
        end if

        n21 = nc(1)*nc(2)
        ncp = n21*nc(3)
        allocate(d2Tgc(ng_, 3*ncp, 3), source=0.0_rk)
        allocate(dTgc(ng_, ncp, 3), source=0.0_rk)
        allocate(Tgc(ng_, ncp), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N1, N2, N3, dN1, dN2, dN3, d2N1, d2N2, d2N3, &
            n32, dn21, dn31, dn32, j1, j2, j3, idx1, idx2, idx3, idx, first1, first2, first3)
#endif
            call basis_bspline_der(map_parameter(Xt(i,1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
                knot1, nc(1), degree(1), 2, first1, N1, dN1, d2N1)
            call basis_bspline_der(map_parameter(Xt(i,2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
                knot2, nc(2), degree(2), 2, first2, N2, dN2, d2N2)
            call basis_bspline_der(map_parameter(Xt(i,3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
                knot3, nc(3), degree(3), 2, first3, N3, dN3, d2N3)
            do j3 = max(0, 1 - first3), min(degree(3), nc(3) - first3)
                idx3 = first3 + j3
                do j2 = max(0, 1 - first2), min(degree(2), nc(2) - first2)
                    idx2 = first2 + j2
                    n32 = N3(j3)*N2(j2)
                    do j1 = max(0, 1 - first1), min(degree(1), nc(1) - first1)
                        idx1 = first1 + j1
                        idx = (idx3 - 1)*n21 + (idx2 - 1)*nc(1) + idx1
                        dn21 = N3(j3)*dN2(j2)*dN1(j1)
                        dn31 = dN3(j3)*N2(j2)*dN1(j1)
                        dn32 = dN3(j3)*dN2(j2)*N1(j1)
                        Tgc(i,idx) = n32*N1(j1)
                        dTgc(i,idx,1) = n32*dN1(j1)
                        dTgc(i,idx,2) = N3(j3)*dN2(j2)*N1(j1)
                        dTgc(i,idx,3) = dN3(j3)*N2(j2)*N1(j1)
                        d2Tgc(i,idx,1) = n32*d2N1(j1)
                        d2Tgc(i,ncp+idx,1) = dn21
                        d2Tgc(i,2*ncp+idx,1) = dn31
                        d2Tgc(i,idx,2) = dn21
                        d2Tgc(i,ncp+idx,2) = N3(j3)*d2N2(j2)*N1(j1)
                        d2Tgc(i,2*ncp+idx,2) = dn32
                        d2Tgc(i,idx,3) = dn31
                        d2Tgc(i,ncp+idx,3) = dn32
                        d2Tgc(i,2*ncp+idx,3) = d2N3(j3)*N2(j2)*N1(j1)
                    end do
                end do
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate one dense polynomial volume basis, gradient, and Hessian.
    !!
    !! `Tgc` has length `ncp`, `dTgc` shape `[ncp,3]`, and the complete packed
    !! Hessian has shape `[3*ncp,3]`. Derivatives above a directional degree are
    !! zero.
    pure subroutine compute_d2Tgc_bspline_3d_scalar(Xt, knot1, knot2, knot3, degree, nc, d2Tgc, dTgc, Tgc, wrap_parameters)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        real(rk), allocatable, intent(out) :: d2Tgc(:,:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable, intent(out) :: Tgc(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        integer :: j1, j2, j3, idx1, idx2, idx3, idx, first1, first2, first3, n21, ncp
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3))
        real(rk) :: dN1(0:degree(1)), dN2(0:degree(2)), dN3(0:degree(3))
        real(rk) :: d2N1(0:degree(1)), d2N2(0:degree(2)), d2N3(0:degree(3))
        real(rk) :: n32, dn21, dn31, dn32
        logical :: wrap_parameters_(3)

        n21 = nc(1)*nc(2)
        ncp = n21*nc(3)
        allocate(d2Tgc(3*ncp, 3), dTgc(ncp, 3), Tgc(ncp), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        call basis_bspline_der(map_parameter(Xt(1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
            knot1, nc(1), degree(1), 2, first1, N1, dN1, d2N1)
        call basis_bspline_der(map_parameter(Xt(2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
            knot2, nc(2), degree(2), 2, first2, N2, dN2, d2N2)
        call basis_bspline_der(map_parameter(Xt(3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
            knot3, nc(3), degree(3), 2, first3, N3, dN3, d2N3)
        do j3 = max(0, 1 - first3), min(degree(3), nc(3) - first3)
            idx3 = first3 + j3
            do j2 = max(0, 1 - first2), min(degree(2), nc(2) - first2)
                idx2 = first2 + j2
                n32 = N3(j3)*N2(j2)
                do j1 = max(0, 1 - first1), min(degree(1), nc(1) - first1)
                    idx1 = first1 + j1
                    idx = (idx3 - 1)*n21 + (idx2 - 1)*nc(1) + idx1
                    dn21 = N3(j3)*dN2(j2)*dN1(j1)
                    dn31 = dN3(j3)*N2(j2)*dN1(j1)
                    dn32 = dN3(j3)*dN2(j2)*N1(j1)
                    Tgc(idx) = n32*N1(j1)
                    dTgc(idx,1) = n32*dN1(j1)
                    dTgc(idx,2) = N3(j3)*dN2(j2)*N1(j1)
                    dTgc(idx,3) = dN3(j3)*N2(j2)*N1(j1)
                    d2Tgc(idx,1) = n32*d2N1(j1)
                    d2Tgc(ncp+idx,1) = dn21
                    d2Tgc(2*ncp+idx,1) = dn31
                    d2Tgc(idx,2) = dn21
                    d2Tgc(ncp+idx,2) = N3(j3)*d2N2(j2)*N1(j1)
                    d2Tgc(2*ncp+idx,2) = dn32
                    d2Tgc(idx,3) = dn31
                    d2Tgc(ncp+idx,3) = dn32
                    d2Tgc(2*ncp+idx,3) = d2N3(j3)*N2(j2)*N1(j1)
                end do
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the dense rational volume basis at many parameter triples.
    !!
    !! The allocated result has shape `[ng_,product(nc)]`, with direction 1
    !! varying fastest in the basis columns. Valid rows form a partition of
    !! unity after positive-weight normalization.
    pure function compute_Tgc_nurbs_3d_vector(Xt, knot1, knot2, knot3, degree, nc, ng, Wc, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: ng(3)
        real(rk), intent(in), contiguous :: Wc(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Tgc(:,:)
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3))
        real(rk) :: wi, wsum, n32, wtol
        integer :: i, j1, j2, j3, idx, first1, first2, first3, ng_
        integer :: j10, j11, j20, j21, j30, j31, weight_exponent
        logical :: wrap_parameters_(3)

        if (present(ng)) then
            ng_ = product(ng)
        else
            ng_ = size(Xt, 1)
        end if

        allocate(Tgc(ng_, nc(1)*nc(2)*nc(3)), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        wtol = 0.0_rk
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N1, N2, N3, wi, wsum, n32, j1, j2, j3, idx, &
            first1, first2, first3, j10, j11, j20, j21, j30, j31, weight_exponent)
#endif
            call basis_bspline_der(map_parameter(Xt(i,1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
                knot1, nc(1), degree(1), 0, first1, N1)
            call basis_bspline_der(map_parameter(Xt(i,2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
                knot2, nc(2), degree(2), 0, first2, N2)
            call basis_bspline_der(map_parameter(Xt(i,3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
                knot3, nc(3), degree(3), 0, first3, N3)
            j10 = max(0, 1 - first1)
            j11 = min(degree(1), nc(1) - first1)
            j20 = max(0, 1 - first2)
            j21 = min(degree(2), nc(2) - first2)
            j30 = max(0, 1 - first3)
            j31 = min(degree(3), nc(3) - first3)
            idx = ((first3+j30-1)*nc(2) + first2+j20-1)*nc(1) + first1+j10
            weight_exponent = exponent(Wc(idx))
            do j3 = j30, j31
                do j2 = j20, j21
                    do j1 = j10, j11
                        idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1+j1
                        weight_exponent = max(weight_exponent, exponent(Wc(idx)))
                    end do
                end do
            end do
            wsum = 0.0_rk
            do j3 = j30, j31
                do j2 = j20, j21
                    n32 = N3(j3)*N2(j2)
                    do j1 = j10, j11
                        idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1 + j1
                        wi = scale(Wc(idx), -weight_exponent)
                        wsum = wsum + n32*N1(j1)*wi
                    end do
                end do
            end do
            if (abs(wsum) > wtol) then
                do j3 = j30, j31
                    do j2 = j20, j21
                        n32 = N3(j3)*N2(j2)
                        do j1 = j10, j11
                            idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1 + j1
                            wi = scale(Wc(idx), -weight_exponent)
                            Tgc(i,idx) = n32*N1(j1)*wi/wsum
                        end do
                    end do
                end do
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the rational volume basis at one parameter triple.
    !!
    !! Dense output has length `product(nc)`. With `elem`, result order follows
    !! the supplied flattened global indices and `Wc` may be global or local;
    !! inactive and invalid listed controls remain zero.
    pure function compute_Tgc_nurbs_3d_scalar(Xt, knot1, knot2, knot3, degree, nc, Wc, elem, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        real(rk), intent(in), contiguous :: Wc(:)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Tgc(:)
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3))
        real(rk) :: wi, wsum, n32, wtol
        integer :: j1, j2, j3, idx, idx1, idx2, idx3, a1, a2, a3
        integer :: first1, first2, first3, n21, ncp, nloc, rem
        integer :: j10, j11, j20, j21, j30, j31, weight_exponent
        logical :: local_weights
        logical :: wrap_parameters_(3)

        n21 = nc(1)*nc(2)
        ncp = n21*nc(3)
        if (present(elem)) then
            nloc = size(elem)
            allocate(Tgc(nloc), source=0.0_rk)
        else
            allocate(Tgc(ncp), source=0.0_rk)
        end if
        wtol = 0.0_rk
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        call basis_bspline_der(map_parameter(Xt(1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
            knot1, nc(1), degree(1), 0, first1, N1)
        call basis_bspline_der(map_parameter(Xt(2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
            knot2, nc(2), degree(2), 0, first2, N2)
        call basis_bspline_der(map_parameter(Xt(3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
            knot3, nc(3), degree(3), 0, first3, N3)
        wsum = 0.0_rk
        if (present(elem)) then
            local_weights = size(Wc) == nloc .and. size(Wc) /= ncp
            if (local_weights) then
                weight_exponent = exponent(maxval(Wc))
                do j1 = 1, nloc
                    idx = elem(j1)
                    if (idx < 1 .or. idx > ncp) cycle
                    idx3 = (idx - 1)/n21 + 1
                    rem = idx - (idx3 - 1)*n21
                    idx2 = (rem - 1)/nc(1) + 1
                    idx1 = rem - (idx2 - 1)*nc(1)
                    a1 = idx1 - first1
                    a2 = idx2 - first2
                    a3 = idx3 - first3
                    if (a1 < 0 .or. a1 > degree(1) .or. a2 < 0 .or. a2 > degree(2) .or. &
                        a3 < 0 .or. a3 > degree(3)) cycle
                    wi = scale(Wc(j1), -weight_exponent)
                    wsum = wsum + N3(a3)*N2(a2)*N1(a1)*wi
                end do
            else
                j10 = max(0, 1 - first1)
                j11 = min(degree(1), nc(1) - first1)
                j20 = max(0, 1 - first2)
                j21 = min(degree(2), nc(2) - first2)
                j30 = max(0, 1 - first3)
                j31 = min(degree(3), nc(3) - first3)
                idx = ((first3+j30-1)*nc(2) + first2+j20-1)*nc(1) + first1+j10
                weight_exponent = exponent(Wc(idx))
                do j3 = j30, j31
                    do j2 = j20, j21
                        do j1 = j10, j11
                            idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1+j1
                            weight_exponent = max(weight_exponent, exponent(Wc(idx)))
                        end do
                    end do
                end do
                do j3 = j30, j31
                    do j2 = j20, j21
                        n32 = N3(j3)*N2(j2)
                        do j1 = j10, j11
                            idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1 + j1
                            wi = scale(Wc(idx), -weight_exponent)
                            wsum = wsum + n32*N1(j1)*wi
                        end do
                    end do
                end do
            end if
            if (abs(wsum) > wtol) then
                do j1 = 1, nloc
                    idx = elem(j1)
                    if (idx < 1 .or. idx > ncp) cycle
                    idx3 = (idx - 1)/n21 + 1
                    rem = idx - (idx3 - 1)*n21
                    idx2 = (rem - 1)/nc(1) + 1
                    idx1 = rem - (idx2 - 1)*nc(1)
                    a1 = idx1 - first1
                    a2 = idx2 - first2
                    a3 = idx3 - first3
                    if (a1 < 0 .or. a1 > degree(1) .or. a2 < 0 .or. a2 > degree(2) .or. &
                        a3 < 0 .or. a3 > degree(3)) cycle
                    if (local_weights) then
                        wi = scale(Wc(j1), -weight_exponent)
                    else
                        wi = scale(Wc(idx), -weight_exponent)
                    end if
                    Tgc(j1) = N3(a3)*N2(a2)*N1(a1)*wi/wsum
                end do
            end if
        else
            j10 = max(0, 1 - first1)
            j11 = min(degree(1), nc(1) - first1)
            j20 = max(0, 1 - first2)
            j21 = min(degree(2), nc(2) - first2)
            j30 = max(0, 1 - first3)
            j31 = min(degree(3), nc(3) - first3)
            idx = ((first3+j30-1)*nc(2) + first2+j20-1)*nc(1) + first1+j10
            weight_exponent = exponent(Wc(idx))
            do j3 = j30, j31
                do j2 = j20, j21
                    do j1 = j10, j11
                        idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1+j1
                        weight_exponent = max(weight_exponent, exponent(Wc(idx)))
                    end do
                end do
            end do
            do j3 = j30, j31
                do j2 = j20, j21
                    n32 = N3(j3)*N2(j2)
                    do j1 = j10, j11
                        idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1 + j1
                        wi = scale(Wc(idx), -weight_exponent)
                        wsum = wsum + n32*N1(j1)*wi
                    end do
                end do
            end do
            if (abs(wsum) > wtol) then
                do j3 = j30, j31
                    do j2 = j20, j21
                        n32 = N3(j3)*N2(j2)
                        do j1 = j10, j11
                            idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1 + j1
                            wi = scale(Wc(idx), -weight_exponent)
                            Tgc(idx) = n32*N1(j1)*wi/wsum
                        end do
                    end do
                end do
            end if
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the dense polynomial volume basis at many parameter triples.
    !!
    !! The result has shape `[ng_,product(nc)]`. Direction 1 varies fastest and
    !! every row contains only the active tensor-product support.
    pure function compute_Tgc_bspline_3d_vector(Xt, knot1, knot2, knot3, degree, nc, ng, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:,:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: ng(3)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Tgc(:,:)
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3)), n32
        integer :: i, j1, j2, j3, idx, first1, first2, first3, ng_
        logical :: wrap_parameters_(3)

        if (present(ng)) then
            ng_ = product(ng)
        else
            ng_ = size(Xt, 1)
        end if

        allocate(Tgc(ng_, nc(1)*nc(2)*nc(3)), source=0.0_rk)
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
#if defined(__NVCOMPILER) || defined(__GFORTRAN__)
        do i = 1, ng_
#else
        do concurrent (i = 1: ng_) local(N1, N2, N3, n32, j1, j2, j3, idx, first1, first2, first3)
#endif
            call basis_bspline_der(map_parameter(Xt(i,1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
                knot1, nc(1), degree(1), 0, first1, N1)
            call basis_bspline_der(map_parameter(Xt(i,2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
                knot2, nc(2), degree(2), 0, first2, N2)
            call basis_bspline_der(map_parameter(Xt(i,3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
                knot3, nc(3), degree(3), 0, first3, N3)
            do j3 = max(0, 1 - first3), min(degree(3), nc(3) - first3)
                do j2 = max(0, 1 - first2), min(degree(2), nc(2) - first2)
                    n32 = N3(j3)*N2(j2)
                    do j1 = max(0, 1 - first1), min(degree(1), nc(1) - first1)
                        idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1 + j1
                        Tgc(i,idx) = n32*N1(j1)
                    end do
                end do
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the polynomial volume basis at one parameter triple.
    !!
    !! Without `elem`, the result is dense in direction-1-fastest control order.
    !! With `elem`, it follows the supplied one-based control-index list and
    !! inactive or invalid entries remain zero.
    pure function compute_Tgc_bspline_3d_scalar(Xt, knot1, knot2, knot3, degree, nc, elem, wrap_parameters) result(Tgc)
        real(rk), intent(in), contiguous :: Xt(:)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
        integer, intent(in) :: degree(3)
        integer, intent(in) :: nc(3)
        integer, intent(in), optional :: elem(:)
        logical, intent(in), contiguous, optional :: wrap_parameters(:)
        real(rk), allocatable :: Tgc(:)
        real(rk) :: N1(0:degree(1)), N2(0:degree(2)), N3(0:degree(3)), n32
        integer :: j1, j2, j3, idx, idx1, idx2, idx3, a1, a2, a3
        integer :: first1, first2, first3, n21, ncp, nloc, rem
        logical :: wrap_parameters_(3)

        n21 = nc(1)*nc(2)
        ncp = n21*nc(3)
        if (present(elem)) then
            nloc = size(elem)
            allocate(Tgc(nloc), source=0.0_rk)
        else
            allocate(Tgc(ncp), source=0.0_rk)
        end if
        wrap_parameters_ = .false.
        if (present(wrap_parameters)) then
            wrap_parameters_(1:min(3,size(wrap_parameters))) = &
                wrap_parameters(1:min(3,size(wrap_parameters)))
        end if
        call basis_bspline_der(map_parameter(Xt(1), knot1, nc(1), degree(1), wrap_parameters_(1)), &
            knot1, nc(1), degree(1), 0, first1, N1)
        call basis_bspline_der(map_parameter(Xt(2), knot2, nc(2), degree(2), wrap_parameters_(2)), &
            knot2, nc(2), degree(2), 0, first2, N2)
        call basis_bspline_der(map_parameter(Xt(3), knot3, nc(3), degree(3), wrap_parameters_(3)), &
            knot3, nc(3), degree(3), 0, first3, N3)
        if (present(elem)) then
            do j1 = 1, nloc
                idx = elem(j1)
                if (idx < 1 .or. idx > ncp) cycle
                idx3 = (idx - 1)/n21 + 1
                rem = idx - (idx3 - 1)*n21
                idx2 = (rem - 1)/nc(1) + 1
                idx1 = rem - (idx2 - 1)*nc(1)
                a1 = idx1 - first1
                a2 = idx2 - first2
                a3 = idx3 - first3
                if (a1 < 0 .or. a1 > degree(1) .or. a2 < 0 .or. a2 > degree(2) .or. &
                    a3 < 0 .or. a3 > degree(3)) cycle
                Tgc(j1) = N3(a3)*N2(a2)*N1(a1)
            end do
        else
            do j3 = max(0, 1 - first3), min(degree(3), nc(3) - first3)
                do j2 = max(0, 1 - first2), min(degree(2), nc(2) - first2)
                    n32 = N3(j3)*N2(j2)
                    do j1 = max(0, 1 - first1), min(degree(1), nc(1) - first1)
                        idx = ((first3+j3-1)*nc(2) + first2+j2-1)*nc(1) + first1 + j1
                        Tgc(idx) = n32*N1(j1)
                    end do
                end do
            end do
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Fit polynomial tensor-product B-spline controls to parameterized volume data.
    !!
    !! `ndata=[n1,n2,n3]` defines the rows used. A recognized Cartesian grid is
    !! solved by three directional banded normal-equation factorizations. The
    !! grid test requires exactly `product(ndata)` parameter and data rows and
    !! compares repeated directional coordinates with tolerance
    !! `32*epsilon(rk)*max(1,maxval(abs(Xt)))`; otherwise the first
    !! `product(ndata)` rows use one banded tensor-product normal system.
    !!
    !! The fit uses the polynomial tensor basis even if explicit weights are
    !! stored; it neither reads nor clears those weights. Thus
    !! `this%is_rational()` must be false if later object evaluation is to match
    !! the fitted objective. `Xt` and `Xdata` must provide enough finite rows,
    !! `Xt` must have at least three columns, `Xdata` at least one, and each
    !! directional parameter set must give full-rank active-domain coverage.
    !! These shape and finiteness preconditions are not completely diagnosed.
    !! Only `Xc` is replaced; cached `Xg` is not recomputed.
    !!
    !! @warning Normal equations require full-rank, well-scaled parameter
    !! coverage and square the design matrix condition number. No QR/SVD
    !! fallback is used.
    !! @endwarning
    pure subroutine lsq_fit_bspline(this, Xt, Xdata, ndata)
        use forcad_utils, only: solve_spd_banded
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in), contiguous :: Xt(:,:)
            !! Parameter triples `[at least product(ndata),3]`.
        real(rk), intent(in), contiguous :: Xdata(:,:)
            !! Data coordinates `[at least product(ndata),ncoord]`.
        integer, intent(in) :: ndata(3)
            !! Directional data counts `[3]`, each at least the corresponding `nc`.
        real(rk), allocatable :: TtT(:,:), TtX(:,:), Xsol(:,:)
        real(rk), allocatable :: G1band(:,:), G2band(:,:), G3band(:,:), RHS(:,:)
        real(rk), allocatable :: rhs1(:,:), rhs2(:,:), rhs3(:,:), sol1(:,:), sol2(:,:), sol3(:,:)
        real(rk), allocatable :: Phi1(:,:), Phi2(:,:), Phi3(:,:)
        real(rk) :: N1(0:this%degree(1)), N2(0:this%degree(2)), N3(0:this%degree(3))
        real(rk) :: ba, bb, n32a, n32b
        integer, allocatable :: F1(:), F2(:), F3(:)
        real(rk) :: grid_tol
        integer :: i, n, ncp, dim_, row, col, i1, i2, i3, bandwidth
        integer :: a1, a2, a3, b1, b2, b3, k, ia, ib, first1, first2, first3
        logical :: tensor_grid

        if (.not. this%err%ok) return

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
            grid_tol = 32.0_rk*epsilon(1.0_rk)*max(1.0_rk, maxval(abs(Xt)))
            do i3 = 1, ndata(3)
                do i2 = 1, ndata(2)
                    do i1 = 1, ndata(1)
                        row = ((i3-1)*ndata(2) + i2-1)*ndata(1) + i1
                        if (abs(Xt(row,1) - Xt(i1,1)) > grid_tol .or. &
                            abs(Xt(row,2) - Xt((i2-1)*ndata(1)+1,2)) > grid_tol .or. &
                            abs(Xt(row,3) - Xt(((i3-1)*ndata(2))*ndata(1)+1,3)) > grid_tol) then
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

            if (allocated(this%Xc)) then
                if (size(this%Xc,1) /= ncp .or. size(this%Xc,2) /= dim_) deallocate(this%Xc)
            end if
            if (.not. allocated(this%Xc)) allocate(this%Xc(ncp, dim_))
            allocate(rhs1(this%nc(1), this%nc(2)*this%nc(3)))
            allocate(rhs2(this%nc(2), this%nc(1)*this%nc(3)))
            allocate(rhs3(this%nc(3), this%nc(1)*this%nc(2)))

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
                    call this%err%set(&
                        code       = 108,&
                        severity   = 1,&
                        category   = "forcad_nurbs_volume",&
                        message    = "Directional banded least-squares solve failed.",&
                        location   = "lsq_fit_bspline",&
                        suggestion = "Check parameter coverage in the first volume direction.")
                    return
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
                    call this%err%set(&
                        code       = 108,&
                        severity   = 1,&
                        category   = "forcad_nurbs_volume",&
                        message    = "Directional banded least-squares solve failed.",&
                        location   = "lsq_fit_bspline",&
                        suggestion = "Check parameter coverage in the second volume direction.")
                    return
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
                    call this%err%set(&
                        code       = 108,&
                        severity   = 1,&
                        category   = "forcad_nurbs_volume",&
                        message    = "Directional banded least-squares solve failed.",&
                        location   = "lsq_fit_bspline",&
                        suggestion = "Check parameter coverage in the third volume direction.")
                    return
                end if
                do i3 = 1, this%nc(3)
                    do i2 = 1, this%nc(2)
                        do i1 = 1, this%nc(1)
                            this%Xc(((i3-1)*this%nc(2)+i2-1)*this%nc(1)+i1,k) = sol3(i3,(i2-1)*this%nc(1)+i1)
                        end do
                    end do
                end do
            end do
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
            call this%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_nurbs_volume",&
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
        this%rational = max(maxval(this%Wc)-this%Wc(1),this%Wc(1)-minval(this%Wc)) > &
            32.0_rk*epsilon(1.0_rk)*maxval(this%Wc)

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
            this%rational = max(maxval(this%Wc)-this%Wc(1),this%Wc(1)-minval(this%Wc)) > &
                32.0_rk*epsilon(1.0_rk)*maxval(this%Wc)
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
