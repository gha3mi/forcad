!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Topology, continuity constraints, and assembly data for surface patches.
!!
!! Patches remain complete [[nurbs_surface]] objects. Connections identify
!! oriented boundary curves and require affine-equivalent tangential trace
!! spaces. [[nurbs_multipatch_surface:cmp_dof_map]] merges controls paired by
!! recorded \(C^0\) topology only where both normal endpoint bases are
!! interpolatory, while
!! [[nurbs_multipatch_surface:cmp_dof_constraint]] supplies sparse equations for
!! normal derivatives through any order
!! \(0\le n\le\min(p_{n,A},p_{n,B})\). A compact geometric connection uses the
!! separable map \(\Phi(\tau,\eta)=(T\tau,\phi(\eta))\). A general connection
!! supplies Bernstein fields for every normal derivative of both mapped
!! coordinates and, for rational patches, every derivative of a nonvanishing
!! projective factor. The resulting tangentially varying multivariate
!! \(G^n\) equations include all mixed chain-rule and projective-product terms.
!!
!! Global element order is patch insertion order followed by each patch's
!! direction-1-fastest active-span order. The container does not assemble PDE
!! operators; it provides the topology and algebraic data required to do so.
!!
!! Unmerged scalar control variables are concatenated by patch,
!! \(\mathbf q=[(\mathbf q^{(1)})^T,\ldots]^T\). An order-zero interface may be
!! represented by the compact map from
!! [[nurbs_multipatch_surface:cmp_dof_map]] only when its normal endpoints are
!! interpolatory. General unclamped traces and higher continuity are returned as
!! \(\mathbf A\mathbf q=\mathbf0\) in one-based CSR form and are not silently
!! folded into element connectivity. Raw rows act on polynomial coefficients
!! or homogeneous rational controls. Rational scalar fields require weighted
!! numerator controls \(w_iq_i\), together with a compatible denominator.
module forcad_multipatch_surface

    use forcad_kinds, only: rk
    use forcad_multipatch, only: multipatch_chain_rule_coefficients, multipatch_composition_coefficients, &
        multipatch_connection, multipatch_compatible_trace_space, multipatch_growth_capacity, &
        multipatch_projective_weight_scale, multipatch_valid_reparameterization
    use forcad_nurbs_surface, only: nurbs_surface
    use forcad_utils, only: active_knot_multiplicity, active_span_count, basis_bernstein, &
        basis_bspline_der_order_active, boundary_index, boundary_layer_index, disjoint_set_map, disjoint_set_union, &
        interpolatory_boundary_layer, knot_end, knot_start, show_pyvista_multipatch, structural_nonzero
    use fordebug, only: debug

    implicit none

    private
    public :: nurbs_multipatch_surface

    integer, parameter :: u_min_ = 1
    integer, parameter :: u_max_ = 2
    integer, parameter :: v_min_ = 3
    integer, parameter :: v_max_ = 4

    !===============================================================================
    !> Ordered collection of connected NURBS surface patches.
    !!
    !! Valid side labels are `u_min`, `u_max`, `v_min`, and `v_max`, with the
    !! documented aliases accepted by `connect`. `reverse` reverses the single
    !! tangential boundary coordinate. `continuity=-1` records topology without
    !! constraints. Nonnegative continuity contributes sparse equations, while
    !! interpolatory normal endpoints additionally permit direct control sharing.
    !!
    !! Diagnostic codes are: 201 invalid patch identifier; 202 invalid side;
    !! 203 invalid continuity metadata; 204 incompatible interface geometry or
    !! trace space; and 205 invalid patch state.
    type nurbs_multipatch_surface
        type(nurbs_surface), allocatable, private :: patches(:) !! Owned patch copies with spare append capacity.
        type(multipatch_connection), allocatable, private :: connections(:) !! Validated interface metadata.
        integer, private :: npatch = 0 !! Number of active patches.
        integer, private :: nconnection = 0 !! Number of active connections.

        type(debug) :: err !! Recoverable diagnostic state for the most recent failed operation.
    contains
        procedure :: add_patch => add_patch_surface !!> Append a validated patch copy and optionally return its identifier.
        procedure :: connect => connect_surface !!> Validate metadata and trace-space shape, then record an edge interface.
        procedure :: create => create_surface !!> Sample every patch using common directional parameters.
        procedure :: finalize => clear_surface !!> Release patches and connections and reset counts.
        procedure :: get_npatch => get_npatch_surface !!> Return the active patch count.
        procedure :: get_nconnection => get_nconnection_surface !!> Return the active connection count.
        procedure :: get_patch => get_patch_surface !!> Return a copy of one patch.
        procedure :: get_connection => get_connection_surface !!> Return a copy of one connection.
        procedure :: is_valid => is_valid_surface !!> Validate patch states, traces, orientations, and continuity.
        procedure :: is_rational => is_rational_surface !!> Report whether any patch has nonuniform rational weights.
        procedure :: cmp_dof_offsets => cmp_dof_offsets_surface !!> Return prefix offsets for unmerged patch controls.
        procedure :: cmp_dof_map => cmp_dof_map_surface !!> Compact directly shareable \(C^0\) controls.
        procedure :: cmp_dof_constraint => cmp_dof_constraint_surface !!> Assemble derivative-continuity equations in CSR form.
        procedure :: cmp_elem => cmp_elem_surface !!> Concatenate patch element connectivity with optional shared numbering.
        procedure :: get_elem => cmp_elem_surface !!> Alias for [[nurbs_multipatch_surface:cmp_elem]].
        procedure :: cmp_elem_patch => cmp_elem_patch_surface !!> Return the owning patch identifier for every global element.
        procedure :: cmp_elem_local => cmp_elem_local_surface !!> Return the patch-local identifier for every global element.
        procedure :: export_Xc => export_Xc_surface !!> Export each patch control net to a separate VTK file.
        procedure :: export_Xg => export_Xg_surface !!> Export each patch's cached geometry to a separate VTK file.
        procedure :: export_Xth_in_Xg => export_Xth_in_Xg_surface !!> Export embedded parameter grids for every patch.
        procedure :: rotate_Xc => rotate_Xc_surface !!> Rotate all patch control points.
        procedure :: rotate_Xg => rotate_Xg_surface !!> Rotate all cached patch geometry points.
        procedure :: translate_Xc => translate_Xc_surface !!> Translate all patch control points.
        procedure :: translate_Xg => translate_Xg_surface !!> Translate all cached patch geometry points.
        procedure :: show => show_surface !!> Display selected patch representations together with PyVista.
    end type
    !===============================================================================

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Append a copy of a valid surface patch.
    pure subroutine add_patch_surface(this, patch, id)
        class(nurbs_multipatch_surface), intent(inout) :: this
            !! Multipatch container.
        type(nurbs_surface), intent(in) :: patch
            !! Complete, internally consistent surface geometry with a clear diagnostic state.
        integer, intent(out), optional :: id
            !! Assigned one-based patch identifier; zero is returned on failure.
        type(nurbs_surface), allocatable :: tmp(:)
        integer :: new_capacity

        if (present(id)) id = 0
        if (.not. this%err%ok) return

        if (.not. patch%err%ok) then
            call this%err%set(&
                code       = 205,&
                severity   = 1,&
                category   = "forcad_multipatch_surface",&
                message    = "Patch has an active error state.",&
                location   = "add_patch",&
                suggestion = "Reset or fix the patch error state before adding it to a multipatch surface.")
            return
        end if
        if (.not. patch%is_valid()) then
            call this%err%set(&
                code       = 205,&
                severity   = 1,&
                category   = "forcad_multipatch_surface",&
                message    = "Patch geometry is incomplete or internally inconsistent.",&
                location   = "add_patch",&
                suggestion = "Initialize consistent knots, controls, degree, and optional weights before adding the patch.")
            return
        end if
        if (.not. allocated(this%patches)) then
            allocate(this%patches(16))
        else if (this%npatch == size(this%patches)) then
            new_capacity = multipatch_growth_capacity(size(this%patches))
            allocate(tmp(new_capacity))
            if (this%npatch > 0) tmp(1:this%npatch) = this%patches(1:this%npatch)
            call move_alloc(tmp, this%patches)
        end if
        this%npatch = this%npatch + 1
        this%patches(this%npatch) = patch
        if (present(id)) id = this%npatch
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Record two surface edges with parametric \(C^n\) or geometric metadata.
    !! `geometric` defaults to `.false.` for C^n continuity.
    !! For the compact G^n representation, `reparameterization` stores the
    !! scalar normal transition derivatives. The general representation uses
    !! `transition_jet(i,c,r)`, the Bernstein coefficient \(i\) of
    !! \(\partial_{\eta_A}^r\Phi_c(\tau,0)\). Component `c=1` is the mapped
    !! patch-B tangential coordinate and `c=2` its normal coordinate.
    !! `projective_jet(i,r+1)` similarly stores
    !! \(\partial_{\eta_A}^r\lambda(\tau,0)\) for homogeneous rational
    !! continuity. Tangential orientation is controlled independently by
    !! `reverse`. Do not pass both compact and full transition jets.
    !!
    !! The connected boundary curves must have equal control counts, equal
    !! tangential degree, and affine-equivalent normalized knot vectors. Their
    !! physical traces, derivative residuals and, for rational patches,
    !! projective weight scales are checked by
    !! [[nurbs_multipatch_surface:is_valid]]. Call that method before assembly or
    !! analysis; successful `connect` only establishes valid metadata and trace
    !! space shape.
    !!
    pure subroutine connect_surface(this, patch_a, side_a, patch_b, side_b, continuity, reverse, &
        geometric, reparameterization, transition_jet, projective_jet)
        class(nurbs_multipatch_surface), intent(inout) :: this
            !! Multipatch container.
        integer, intent(in) :: patch_a
            !! One-based first patch identifier.
        integer, intent(in) :: patch_b
            !! One-based second patch identifier.
        character(len=*), intent(in) :: side_a
            !! First side label (`u_min`, `u_max`, `v_min`, or `v_max`).
        character(len=*), intent(in) :: side_b
            !! Second side label.
        integer, intent(in), optional :: continuity
            !! Requested order, default zero; `-1` disables constraints.
        logical, intent(in), optional :: reverse
            !! Reverse the tangential coordinate of patch B.
        logical, intent(in), optional :: geometric
            !! Select geometric rather than parametric continuity.
        real(rk), intent(in), contiguous, optional :: reparameterization(:)
            !! Compact normalized normal transition derivatives for \(G^n\).
        real(rk), intent(in), contiguous, optional :: transition_jet(:,:,:)
            !! Full Bernstein transition coefficients `[q+1,2,n]`.
        real(rk), intent(in), contiguous, optional :: projective_jet(:,:)
            !! Bernstein projective-factor coefficients `[q+1,n+1]`.
        type(multipatch_connection) :: conn
        type(multipatch_connection), allocatable :: tmp(:)
        real(rk), allocatable :: transition_jet_(:,:,:,:), projective_jet_(:,:,:)
        integer, allocatable :: ia(:), ib(:)
        integer :: side_a_id, side_b_id, continuity_, new_capacity, dir_a, dir_b, normal_sign
        logical :: reverse_, geometric_

        if (.not. this%err%ok) return

        if (patch_a < 1 .or. patch_a > this%get_npatch()) then
            call this%err%set(&
                code       = 201,&
                severity   = 1,&
                category   = "forcad_multipatch_surface",&
                message    = "Invalid first patch id.",&
                location   = "connect",&
                suggestion = "Use a patch id from add_patch's optional id argument.")
            return
        end if
        if (patch_b < 1 .or. patch_b > this%get_npatch()) then
            call this%err%set(&
                code       = 201,&
                severity   = 1,&
                category   = "forcad_multipatch_surface",&
                message    = "Invalid second patch id.",&
                location   = "connect",&
                suggestion = "Use a patch id from add_patch's optional id argument.")
            return
        end if
        side_a_id = 0
        side_b_id = 0
        select case(trim(side_a))
        case("u_min", "umin", "u-", "U_MIN", "UMIN", "U-")
            side_a_id = u_min_
        case("u_max", "umax", "u+", "U_MAX", "UMAX", "U+")
            side_a_id = u_max_
        case("v_min", "vmin", "v-", "V_MIN", "VMIN", "V-")
            side_a_id = v_min_
        case("v_max", "vmax", "v+", "V_MAX", "VMAX", "V+")
            side_a_id = v_max_
        case default
            continue
        end select
        select case(trim(side_b))
        case("u_min", "umin", "u-", "U_MIN", "UMIN", "U-")
            side_b_id = u_min_
        case("u_max", "umax", "u+", "U_MAX", "UMAX", "U+")
            side_b_id = u_max_
        case("v_min", "vmin", "v-", "V_MIN", "VMIN", "V-")
            side_b_id = v_min_
        case("v_max", "vmax", "v+", "V_MAX", "VMAX", "V+")
            side_b_id = v_max_
        case default
            continue
        end select
        if (side_a_id == 0 .or. side_b_id == 0) then
            call this%err%set(&
                code       = 202,&
                severity   = 1,&
                category   = "forcad_multipatch_surface",&
                message    = "Invalid surface side label.",&
                location   = "connect",&
                suggestion = "Use u_min/u_max/v_min/v_max or u-/u+/v-/v+.")
            return
        end if
        continuity_ = 0
        if (present(continuity)) continuity_ = continuity
        geometric_ = .false.
        if (present(geometric)) geometric_ = geometric
        if (continuity_ < -1) then
            call this%err%set(&
                code       = 203,&
                severity   = 1,&
                category   = "forcad_multipatch_surface",&
                message    = "Invalid continuity value.",&
                location   = "connect",&
                suggestion = "Use -1 for discontinuous or Cn with n >= 0.")
            return
        end if
        dir_a = merge(1, 2, side_a_id <= u_max_)
        dir_b = merge(1, 2, side_b_id <= u_max_)
        if (continuity_ > min(this%patches(patch_a)%get_degree(dir_a), this%patches(patch_b)%get_degree(dir_b))) then
            call this%err%set(&
                code       = 203,&
                severity   = 1,&
                category   = "forcad_multipatch_surface",&
                message    = "Continuity is higher than the connected patch degree.",&
                location   = "connect",&
                suggestion = "Use a continuity not larger than the minimum degree normal to the connected sides.")
            return
        end if
        reverse_ = .false.
        if (present(reverse)) reverse_ = reverse
        if (present(transition_jet)) then
            allocate(transition_jet_(size(transition_jet,1),1,size(transition_jet,2),size(transition_jet,3)))
            transition_jet_(:,1,:,:) = transition_jet
        else
            allocate(transition_jet_(0,0,0,0))
        end if
        if (present(projective_jet)) then
            allocate(projective_jet_(size(projective_jet,1),1,size(projective_jet,2)))
            projective_jet_(:,1,:) = projective_jet
        else
            allocate(projective_jet_(0,0,0))
        end if
        call conn%set(&
            patch_a            = patch_a,&
            side_a             = side_a_id,&
            patch_b            = patch_b,&
            side_b             = side_b_id,&
            continuity         = continuity_,&
            reverse            = reverse_,&
            geometric          = geometric_,&
            reparameterization = reparameterization,&
            transition_jet     = transition_jet_,&
            projective_jet     = projective_jet_)
        normal_sign = merge(-1, 1, mod(side_a_id, 2) == mod(side_b_id, 2))
        if (.not. multipatch_valid_reparameterization(conn, normal_sign, 2)) then
            call this%err%set(&
                code       = 203,&
                severity   = 1,&
                category   = "forcad_multipatch_surface",&
                message    = "Invalid continuity reparameterization.",&
                location   = "connect",&
                suggestion = "Pass either n compact derivatives or a finite `[q+1,2,n]` full jet "//&
                    "whose first normal derivative has the side-consistent sign.")
            return
        end if
        if (continuity_ >= 0) then
            call boundary_index(this%patches(patch_a)%get_nc(), side_a_id, .false., ia)
            call boundary_index(this%patches(patch_b)%get_nc(), side_b_id, reverse_, ib)
            if (size(ia) /= size(ib)) then
                call this%err%set(&
                    code       = 204,&
                    severity   = 1,&
                    category   = "forcad_multipatch_surface",&
                    message    = "Connected surface sides have incompatible control point counts.",&
                    location   = "connect",&
                    suggestion = "Use matching boundary discretizations or refine one side before connecting patches.")
                return
            end if
            if (.not. connection_trace_space_surface(this, conn, sqrt(epsilon(1.0_rk)))) then
                call this%err%set(&
                    code       = 204,&
                    severity   = 1,&
                    category   = "forcad_multipatch_surface",&
                    message    = "Connected surface sides have incompatible tangential trace spaces.",&
                    location   = "connect",&
                    suggestion = "Use matching tangential degree and affine-equivalent knot vectors before connecting patches.")
                return
            end if
        end if
        if (.not. allocated(this%connections)) then
            allocate(this%connections(16))
        else if (this%nconnection == size(this%connections)) then
            new_capacity = multipatch_growth_capacity(size(this%connections))
            allocate(tmp(new_capacity))
            if (this%nconnection > 0) tmp(1:this%nconnection) = this%connections(1:this%nconnection)
            call move_alloc(tmp, this%connections)
        end if
        this%nconnection = this%nconnection + 1
        this%connections(this%nconnection) = conn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Release all owned surface patches and connections.
    !!
    !! Active counts are reset to zero and the multipatch diagnostic is
    !! preserved.
    pure elemental subroutine clear_surface(this)
        class(nurbs_multipatch_surface), intent(inout) :: this

        if (allocated(this%patches)) deallocate(this%patches)
        if (allocated(this%connections)) deallocate(this%connections)
        this%npatch = 0
        this%nconnection = 0
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the number of active surface patches.
    pure elemental integer function get_npatch_surface(this) result(n)
        class(nurbs_multipatch_surface), intent(in) :: this

        n = this%npatch
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return exact allocation bounds for one general geometric edge constraint.
    pure subroutine general_constraint_shape_surface(this, conn, nrow, nnz)
        class(nurbs_multipatch_surface), intent(in) :: this
        type(multipatch_connection), intent(in) :: conn
        integer, intent(out) :: nrow, nnz
        real(rk), allocatable :: knot(:), transition_jet(:,:,:,:), projective_jet(:,:,:)
        integer :: r, degree_row, transition_degree, projective_degree, nspan
        integer :: pa, pb, normal_a, normal_b, tangent_a, tangent_b
        integer :: nc_a(2), nc_b(2), degree_a(2), degree_b(2), support

        pa = conn%get_patch_a()
        pb = conn%get_patch_b()
        normal_a = (conn%get_side_a() + 1)/2
        normal_b = (conn%get_side_b() + 1)/2
        tangent_a = 3 - normal_a
        tangent_b = 3 - normal_b
        nc_a = this%patches(pa)%get_nc()
        nc_b = this%patches(pb)%get_nc()
        degree_a = this%patches(pa)%get_degree()
        degree_b = this%patches(pb)%get_degree()
        knot = this%patches(pa)%get_knot(tangent_a)
        nspan = active_span_count(knot, nc_a(tangent_a), degree_a(tangent_a))
        transition_degree = 0
        projective_degree = 0
        if (conn%has_transition_jet()) then
            transition_jet = conn%get_transition_jet()
            transition_degree = size(transition_jet,1) - 1
        end if
        if (conn%has_projective_jet()) then
            projective_jet = conn%get_projective_jet()
            projective_degree = size(projective_jet,1) - 1
        end if

        nrow = 0
        do r = 0, conn%get_continuity()
            degree_row = max(degree_a(tangent_a), &
                degree_b(tangent_b) + r*transition_degree + projective_degree)
            nrow = nrow + nspan*(degree_row + 1)
        end do
        support = (degree_a(tangent_a) + 1)*(degree_a(normal_a) + 1) + &
            (degree_b(tangent_b) + 1)*(degree_b(normal_b) + 1)
        nnz = nrow*support
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Append one full multivariate surface \(G^n\) constraint to CSR storage.
    !!
    !! On each nonzero tangential knot span, an order-\(r\) residual has degree
    !! at most
    !! \(d_r=\max(p_A,p_B+r q_\Phi+q_\lambda)\). Evaluation at `d_r+1`
    !! distinct interior points is therefore equivalent to coefficient
    !! equality for the polynomial residual on that span. Mixed patch-B
    !! derivatives are weighted by the multivariate Faa di Bruno coefficients;
    !! projective derivatives use the Leibniz rule.
    pure subroutine append_general_constraint_surface(this, conn, offsets, rowptr, col, val, row, pos)
        class(nurbs_multipatch_surface), intent(in) :: this
        type(multipatch_connection), intent(in) :: conn
        integer, intent(in), contiguous :: offsets(:)
        integer, intent(inout), contiguous :: rowptr(:), col(:)
        real(rk), intent(inout), contiguous :: val(:)
        integer, intent(inout) :: row, pos
        real(rk), allocatable :: knot_na(:), knot_nb(:), knot_ta(:), knot_tb(:)
        real(rk), allocatable :: basis_na(:,:), basis_nb(:,:), basis_ta(:), basis_tb(:,:)
        real(rk), allocatable :: bernstein_transition(:), bernstein_projective(:)
        real(rk), allocatable :: transition_jet(:,:,:,:), projective_jet(:,:,:)
        real(rk), allocatable :: transition_value(:,:), projective_value(:)
        real(rk), allocatable :: composition(:,:,:,:), reparameterization(:)
        real(rk) :: start_na, start_nb, start_ta, start_tb
        real(rk) :: range_na, range_nb, range_ta, range_tb
        real(rk) :: parameter_na, parameter_nb, parameter_ta, parameter_tb
        real(rk) :: span_left, span_right, tau, coefficient_b, mixed, choose
        real(rk) :: weight_scale
        integer, allocatable :: first_na(:), first_nb(:), first_tb(:)
        integer :: pa, pb, side_a, side_b, normal_a, normal_b, tangent_a, tangent_b
        integer :: nc_a(2), nc_b(2), degree_a(2), degree_b(2), index_a(2), index_b(2)
        integer :: n, r, j, alpha_t, alpha_n, local_t, local_n, span, point
        integer :: first_ta, degree_row, transition_degree, projective_degree
        logical :: compatible_weights

        pa = conn%get_patch_a()
        pb = conn%get_patch_b()
        side_a = conn%get_side_a()
        side_b = conn%get_side_b()
        normal_a = (side_a + 1)/2
        normal_b = (side_b + 1)/2
        tangent_a = 3 - normal_a
        tangent_b = 3 - normal_b
        nc_a = this%patches(pa)%get_nc()
        nc_b = this%patches(pb)%get_nc()
        degree_a = this%patches(pa)%get_degree()
        degree_b = this%patches(pb)%get_degree()
        knot_na = this%patches(pa)%get_knot(normal_a)
        knot_nb = this%patches(pb)%get_knot(normal_b)
        knot_ta = this%patches(pa)%get_knot(tangent_a)
        knot_tb = this%patches(pb)%get_knot(tangent_b)
        start_na = knot_start(knot_na, nc_a(normal_a), degree_a(normal_a))
        start_nb = knot_start(knot_nb, nc_b(normal_b), degree_b(normal_b))
        start_ta = knot_start(knot_ta, nc_a(tangent_a), degree_a(tangent_a))
        start_tb = knot_start(knot_tb, nc_b(tangent_b), degree_b(tangent_b))
        range_na = knot_end(knot_na, nc_a(normal_a), degree_a(normal_a)) - start_na
        range_nb = knot_end(knot_nb, nc_b(normal_b), degree_b(normal_b)) - start_nb
        range_ta = knot_end(knot_ta, nc_a(tangent_a), degree_a(tangent_a)) - start_ta
        range_tb = knot_end(knot_tb, nc_b(tangent_b), degree_b(tangent_b)) - start_tb
        parameter_na = merge(start_na, start_na + range_na, mod(side_a,2) == 1)
        parameter_nb = merge(start_nb, start_nb + range_nb, mod(side_b,2) == 1)
        n = conn%get_continuity()

        transition_degree = 0
        projective_degree = 0
        if (conn%has_transition_jet()) then
            transition_jet = conn%get_transition_jet()
            transition_degree = size(transition_jet,1) - 1
        else
            reparameterization = conn%get_reparameterization()
        end if
        if (conn%has_projective_jet()) then
            projective_jet = conn%get_projective_jet()
            projective_degree = size(projective_jet,1) - 1
            weight_scale = 1.0_rk
        else
            call connection_weight_scale_surface(&
                this,&
                conn,&
                sqrt(epsilon(1.0_rk)),&
                weight_scale,&
                compatible_weights)
            if (.not. compatible_weights) weight_scale = 1.0_rk
        end if

        allocate(&
            basis_na(0:degree_a(normal_a),0:n),&
            basis_nb(0:degree_b(normal_b),0:n),&
            basis_ta(0:degree_a(tangent_a)),&
            basis_tb(0:degree_b(tangent_b),0:n),&
            first_na(0:n),&
            first_nb(0:n),&
            first_tb(0:n),&
            transition_value(2,n),&
            projective_value(0:n))
        do j = 0, n
            call basis_bspline_der_order_active(&
                parameter_na,knot_na,nc_a(normal_a),degree_a(normal_a),j,first_na(j),basis_na(:,j))
            call basis_bspline_der_order_active(&
                parameter_nb,knot_nb,nc_b(normal_b),degree_b(normal_b),j,first_nb(j),basis_nb(:,j))
        end do

        do r = 0, n
            degree_row = max(degree_a(tangent_a), &
                degree_b(tangent_b) + r*transition_degree + projective_degree)
            do span = degree_a(tangent_a) + 1, nc_a(tangent_a)
                if (knot_ta(span+1) <= knot_ta(span)) cycle
                span_left = (knot_ta(span) - start_ta)/range_ta
                span_right = (knot_ta(span+1) - start_ta)/range_ta
                do point = 1, degree_row + 1
                    tau = span_left + (span_right - span_left)*&
                        (real(point, rk) - 0.5_rk)/real(degree_row + 1, rk)
                    parameter_ta = start_ta + tau*range_ta
                    parameter_tb = start_tb + merge(1.0_rk - tau, tau, conn%is_reversed())*range_tb
                    call basis_bspline_der_order_active(&
                        parameter_ta,knot_ta,nc_a(tangent_a),degree_a(tangent_a),0,first_ta,basis_ta)
                    do j = 0, n
                        call basis_bspline_der_order_active(&
                            parameter_tb,knot_tb,nc_b(tangent_b),degree_b(tangent_b),j,first_tb(j),basis_tb(:,j))
                    end do

                    transition_value = 0.0_rk
                    if (conn%has_transition_jet()) then
                        bernstein_transition = basis_bernstein(tau, size(transition_jet,1))
                        do j = 1, n
                            transition_value(1,j) = dot_product(&
                                bernstein_transition,transition_jet(:,1,1,j))
                            transition_value(2,j) = dot_product(&
                                bernstein_transition,transition_jet(:,1,2,j))
                        end do
                    else if (n > 0) then
                        transition_value(2,:) = reparameterization
                    end if
                    call multipatch_composition_coefficients(transition_value, composition)

                    projective_value = 0.0_rk
                    if (conn%has_projective_jet()) then
                        bernstein_projective = basis_bernstein(tau, size(projective_jet,1))
                        do j = 0, n
                            projective_value(j) = dot_product(&
                                bernstein_projective,projective_jet(:,1,j+1))
                        end do
                    else
                        projective_value(0) = weight_scale
                    end if

                    row = row + 1
                    do local_t = 0, degree_a(tangent_a)
                        index_a(tangent_a) = first_ta + local_t
                        do local_n = 0, degree_a(normal_a)
                            index_a(normal_a) = first_na(r) + local_n
                            coefficient_b = basis_ta(local_t)*range_na**r*basis_na(local_n,r)
                            if (.not. structural_nonzero(coefficient_b)) cycle
                            col(pos) = offsets(pa) + index_a(1) + (index_a(2)-1)*nc_a(1)
                            val(pos) = coefficient_b
                            pos = pos + 1
                        end do
                    end do
                    do local_t = 0, degree_b(tangent_b)
                        index_b(tangent_b) = first_tb(0) + local_t
                        do local_n = 0, degree_b(normal_b)
                            index_b(normal_b) = first_nb(0) + local_n
                            coefficient_b = 0.0_rk
                            choose = 1.0_rk
                            do j = 0, r
                                mixed = 0.0_rk
                                do alpha_n = 0, j
                                    do alpha_t = 0, j - alpha_n
                                        mixed = mixed + composition(j,alpha_t,alpha_n,0)*&
                                            range_tb**alpha_t*basis_tb(local_t,alpha_t)*&
                                            range_nb**alpha_n*basis_nb(local_n,alpha_n)
                                    end do
                                end do
                                coefficient_b = coefficient_b + &
                                    choose*projective_value(r-j)*mixed
                                if (j < r) choose = choose*real(r-j, rk)/real(j+1, rk)
                            end do
                            if (.not. structural_nonzero(coefficient_b)) cycle
                            col(pos) = offsets(pb) + index_b(1) + (index_b(2)-1)*nc_b(1)
                            val(pos) = -coefficient_b
                            pos = pos + 1
                        end do
                    end do
                    rowptr(row+1) = pos
                end do
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the number of recorded edge connections.
    pure elemental integer function get_nconnection_surface(this) result(n)
        class(nurbs_multipatch_surface), intent(in) :: this

        n = this%nconnection
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return a value copy of patch `i`.
    !!
    !! An out-of-range one-based identifier returns a default empty patch.
    pure function get_patch_surface(this, i) result(patch)
        class(nurbs_multipatch_surface), intent(in) :: this
        integer, intent(in) :: i
        type(nurbs_surface) :: patch

        if (i < 1 .or. i > this%npatch) return
        patch = this%patches(i)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return a value copy of connection `i`.
    !!
    !! An out-of-range one-based identifier returns a default connection.
    pure function get_connection_surface(this, i) result(conn)
        class(nurbs_multipatch_surface), intent(in) :: this
        integer, intent(in) :: i
        type(multipatch_connection) :: conn

        if (i < 1 .or. i > this%nconnection) return
        conn = this%connections(i)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return prefix offsets into patch-concatenated control-variable numbering.
    !! The result starts at zero and its final entry is the total unmerged
    !! control-point count across all patches.
    pure function cmp_dof_offsets_surface(this) result(offsets)
        class(nurbs_multipatch_surface), intent(in) :: this
        integer, allocatable :: offsets(:)
        integer :: i, nc(2)

        allocate(offsets(this%get_npatch()+1))
        offsets(1) = 0
        do i = 1, this%get_npatch()
            nc = this%patches(i)%get_nc()
            offsets(i+1) = offsets(i) + product(nc)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build compact global identifiers for directly shareable edge controls.
    !! Direct identification is performed only when both normal endpoint basis
    !! vectors are exactly one-hot. For a rational trace,
    !!
    !! \[
    !! R_j^{(P)}(\tau)=
    !! \frac{N_j(\tau)\widehat w_j^{(P)}}
    !! {\sum_k N_k(\tau)\widehat w_k^{(P)}} ,
    !! \]
    !!
    !! corresponding scalar controls are shareable only when the oriented trace
    !! weights satisfy
    !! \(\widehat w_j^{(A)}=s\widehat w_j^{(B)}\) for one constant nonzero
    !! \(s\). A tangentially varying projective factor changes the rational basis
    !! and therefore remains in
    !! [[nurbs_multipatch_surface:cmp_dof_constraint]]. General unclamped
    !! \(C^0\) traces likewise remain constraints. The map trusts recorded
    !! topology; call `is_valid` first.
    !!
    !! Disjoint-set implementation reference: R. E. Tarjan, *JACM* 22 (1975),
    !! 215--225.
    !! [doi:10.1145/321879.321884](https://doi.org/10.1145/321879.321884).
    pure function cmp_dof_map_surface(this) result(map)
        class(nurbs_multipatch_surface), intent(in) :: this
        integer, allocatable :: map(:)
        integer, allocatable :: parent(:), offsets(:), ia(:), ib(:)
        real(rk), allocatable :: knot_a(:), knot_b(:)
        real(rk) :: weight_scale
        integer :: i, j, total, pa, pb, dir_a, dir_b, dega, degb, layer_a, layer_b
        integer :: nca(2), ncb(2)
        logical :: compatible_weights

        offsets = this%cmp_dof_offsets()
        total = offsets(size(offsets))
        allocate(parent(total))
        do i = 1, total
            parent(i) = i
        end do

        do i = 1, this%get_nconnection()
            associate(conn => this%connections(i))
                if (conn%get_continuity() < 0) cycle
                pa = conn%get_patch_a()
                pb = conn%get_patch_b()
                nca = this%patches(pa)%get_nc()
                ncb = this%patches(pb)%get_nc()
                dir_a = (conn%get_side_a() + 1)/2
                dir_b = (conn%get_side_b() + 1)/2
                dega = this%patches(pa)%get_degree(dir_a)
                degb = this%patches(pb)%get_degree(dir_b)
                knot_a = this%patches(pa)%get_knot(dir_a)
                knot_b = this%patches(pb)%get_knot(dir_b)
                layer_a = interpolatory_boundary_layer(knot_a, nca(dir_a), dega, conn%get_side_a())
                layer_b = interpolatory_boundary_layer(knot_b, ncb(dir_b), degb, conn%get_side_b())
                if (layer_a < 0 .or. layer_b < 0) cycle
                call connection_weight_scale_surface(&
                    this,&
                    conn,&
                    128.0_rk*epsilon(1.0_rk),&
                    weight_scale,&
                    compatible_weights)
                if (.not. compatible_weights) cycle
                call boundary_layer_index(nca, conn%get_side_a(), layer_a, .false., ia)
                call boundary_layer_index(ncb, conn%get_side_b(), layer_b, conn%is_reversed(), ib)
                do j = 1, size(ia)
                    call disjoint_set_union(parent, offsets(pa)+ia(j), offsets(pb)+ib(j))
                end do
            end associate
        end do
        call disjoint_set_map(parent, map)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Check projective compatibility of rational edge traces.
    !!
    !! Normal endpoint basis values are applied before comparing the oriented
    !! tangential denominator coefficients. Compatible traces satisfy
    !! \(\widehat w_j^{(A)}=\lambda\widehat w_j^{(B)}\) within `tol`;
    !! `weight_scale` returns \(\lambda\). Polynomial traces use unit scale.
    pure subroutine connection_weight_scale_surface(this, conn, tol, weight_scale, ok)
        class(nurbs_multipatch_surface), intent(in) :: this
        type(multipatch_connection), intent(in) :: conn
        real(rk), intent(in) :: tol
        real(rk), intent(out) :: weight_scale
        logical, intent(out) :: ok
        real(rk), allocatable :: basis_a(:), basis_b(:), knot_a(:), knot_b(:), trace_a(:), trace_b(:)
        real(rk) :: xa, xb
        integer, allocatable :: ia(:), ib(:)
        integer :: j, layer, pa, pb, dir_a, dir_b, dega, degb, first_a, first_b
        integer :: normal_a, normal_b, active_a, active_b
        integer :: nca(2), ncb(2)
        logical :: rational_a, rational_b

        weight_scale = 1.0_rk
        ok = .true.
        pa = conn%get_patch_a()
        pb = conn%get_patch_b()
        rational_a = this%patches(pa)%is_rational()
        rational_b = this%patches(pb)%is_rational()
        if (.not. rational_a .and. .not. rational_b) return
        if (rational_a .neqv. rational_b) then
            ok = .false.
            return
        end if

        nca = this%patches(pa)%get_nc()
        ncb = this%patches(pb)%get_nc()
        dir_a = (conn%get_side_a() + 1)/2
        dir_b = (conn%get_side_b() + 1)/2
        dega = this%patches(pa)%get_degree(dir_a)
        degb = this%patches(pb)%get_degree(dir_b)
        knot_a = this%patches(pa)%get_knot(dir_a)
        knot_b = this%patches(pb)%get_knot(dir_b)
        xa = merge(knot_start(knot_a, nca(dir_a), dega), knot_end(knot_a, nca(dir_a), dega), &
            mod(conn%get_side_a(), 2) == 1)
        xb = merge(knot_start(knot_b, ncb(dir_b), degb), knot_end(knot_b, ncb(dir_b), degb), &
            mod(conn%get_side_b(), 2) == 1)
        allocate(basis_a(0:dega), basis_b(0:degb))
        call basis_bspline_der_order_active(xa, knot_a, nca(dir_a), dega, 0, first_a, basis_a)
        call basis_bspline_der_order_active(xb, knot_b, ncb(dir_b), degb, 0, first_b, basis_b)

        allocate(trace_a(nca(3-dir_a)), source=0.0_rk)
        do active_a = 0, dega
            normal_a = first_a + active_a
            if (normal_a < 1 .or. normal_a > nca(dir_a)) cycle
            if (.not. structural_nonzero(basis_a(active_a))) cycle
            layer = merge(normal_a - 1, nca(dir_a) - normal_a, mod(conn%get_side_a(), 2) == 1)
            call boundary_layer_index(nca, conn%get_side_a(), layer, .false., ia)
            do j = 1, size(ia)
                trace_a(j) = trace_a(j) + basis_a(active_a)*this%patches(pa)%get_Wc(ia(j))
            end do
        end do
        allocate(trace_b(ncb(3-dir_b)), source=0.0_rk)
        do active_b = 0, degb
            normal_b = first_b + active_b
            if (normal_b < 1 .or. normal_b > ncb(dir_b)) cycle
            if (.not. structural_nonzero(basis_b(active_b))) cycle
            layer = merge(normal_b - 1, ncb(dir_b) - normal_b, mod(conn%get_side_b(), 2) == 1)
            call boundary_layer_index(ncb, conn%get_side_b(), layer, conn%is_reversed(), ib)
            do j = 1, size(ib)
                trace_b(j) = trace_b(j) + basis_b(active_b)*this%patches(pb)%get_Wc(ib(j))
            end do
        end do
        call multipatch_projective_weight_scale(trace_a, trace_b, tol, weight_scale, ok)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Test whether connected edges have compatible tangential spline spaces.
    !!
    !! The mapped traces must have equal degree, equal control count, and
    !! affine-equivalent normalized knot vectors. Reversal of patch B is
    !! included in the comparison.
    pure logical function connection_trace_space_surface(this, conn, tol) result(ok)
        class(nurbs_multipatch_surface), intent(in) :: this
        type(multipatch_connection), intent(in) :: conn
        real(rk), intent(in) :: tol
        real(rk), allocatable :: knot_a(:), knot_b(:)
        integer :: pa, pb, tan_a, tan_b

        pa = conn%get_patch_a()
        pb = conn%get_patch_b()
        tan_a = merge(2, 1, conn%get_side_a() <= u_max_)
        tan_b = merge(2, 1, conn%get_side_b() <= u_max_)
        knot_a = this%patches(pa)%get_knot(tan_a)
        knot_b = this%patches(pb)%get_knot(tan_b)
        ok = multipatch_compatible_trace_space(knot_a, this%patches(pa)%get_nc(tan_a), &
            this%patches(pa)%get_degree(tan_a), knot_b, this%patches(pb)%get_nc(tan_b), &
            this%patches(pb)%get_degree(tan_b), conn%is_reversed(), tol)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Assemble sparse parametric and geometric constraints for connected edges.
    !!
    !! The one-based CSR arrays define \(Aq=0\) over unmerged,
    !! patch-concatenated control variables. One row is emitted for every
    !! tangential boundary control and derivative order. `geometric`, when
    !! present, selects only parametric (false) or geometric (true) connections;
    !! when absent, both kinds are included. Compact geometric rows use the
    !! scalar normal transition. General rows include tangentially varying
    !! multivariate transition and projective jets. Rational projective scaling
    !! is incorporated for compatible boundary weights.
    !! Apply the raw matrix to polynomial coefficients or homogeneous geometry
    !! columns. For a rational scalar field with coefficients \(q_i\), apply it
    !! to weighted numerator controls \(w_iq_i\), not directly to \(q_i\).
    !!
    !! For every mapped tangential control and derivative order \(r\), one row
    !! discretizes
    !!
    !! \[
    !! \partial_{\eta_A}^{\,r}f_A(\tau)
    !! -\sum_{k=0}^{r}B_{r,k}(\phi',\ldots)
    !!  \partial_{\eta_B}^{\,k}f_B(\tau)=0.
    !! \]
    !!
    !! Tangential orientation is applied before row assembly. Parameter-range
    !! powers convert stored-knot derivatives to normalized normal derivatives.
    !! Call `is_valid` first; this routine assumes compatible connections.
    !! Coefficients with absolute value at most `128*epsilon(rk)` are omitted
    !! from CSR storage, so the returned matrix is numerically sparsified rather
    !! than symbolically exact.
    !!
    pure subroutine cmp_dof_constraint_surface(this, rowptr, col, val, geometric)
        class(nurbs_multipatch_surface), intent(in) :: this
            !! Valid multipatch surface.
        integer, allocatable, intent(out) :: rowptr(:)
            !! CSR row pointers of size `nrow+1`, starting at one.
        integer, allocatable, intent(out) :: col(:)
            !! One-based unmerged global column identifiers.
        real(rk), allocatable, intent(out) :: val(:)
            !! Constraint coefficients corresponding to `col`.
        logical, intent(in), optional :: geometric
            !! Optional filter: false for parametric, true for geometric connections.
        real(rk), allocatable :: ca(:,:), cb(:,:), chain(:,:), knot_a(:), knot_b(:)
        real(rk) :: xa, xb, weight_scale, range_a, range_b, der_scale_a, coefficient_b, eps
        integer, allocatable :: first_a(:), first_b(:), offsets(:), ia(:), ib(:), layer_a(:,:), layer_b(:,:)
        integer :: i, j, k, r, layer, row, pos, nrow, nnz, pa, pb, side_a, side_b
        integer :: dir_a, dir_b, normal_a, normal_b, normal_sign, ntan
        integer :: active_a, active_b, general_nrow, general_nnz
        integer :: nca(2), ncb(2), dega, degb, max_layer_a, max_layer_b
        logical :: compatible_weights

        offsets = this%cmp_dof_offsets()
        nrow = 0
        nnz = 0
        do i = 1, this%get_nconnection()
            associate(conn => this%connections(i))
                if (conn%get_continuity() < 0) cycle
                if (present(geometric)) then
                    if (conn%is_geometric() .neqv. geometric) cycle
                end if
                if (conn%has_transition_jet() .or. conn%has_projective_jet()) then
                    call general_constraint_shape_surface(this, conn, general_nrow, general_nnz)
                    nrow = nrow + general_nrow
                    nnz = nnz + general_nnz
                    cycle
                end if
                pa = conn%get_patch_a()
                pb = conn%get_patch_b()
                call boundary_index(this%patches(pa)%get_nc(), conn%get_side_a(), .false., ia)
                ntan = size(ia)
                dir_a = merge(1, 2, conn%get_side_a() <= u_max_)
                dir_b = merge(1, 2, conn%get_side_b() <= u_max_)
                dega = this%patches(pa)%get_degree(dir_a)
                degb = this%patches(pb)%get_degree(dir_b)
                nrow = nrow + ntan*(conn%get_continuity() + 1)
                nnz = nnz + ntan*(conn%get_continuity() + 1)*(&
                    dega + degb + 2)
            end associate
        end do

        allocate(rowptr(nrow+1))
        if (nrow == 0) then
            rowptr(1) = 1
            allocate(col(0), val(0))
            return
        end if
        allocate(col(nnz), val(nnz))
        row = 0
        pos = 1
        rowptr(1) = 1
        eps = 128.0_rk*epsilon(1.0_rk)
        do i = 1, this%get_nconnection()
            associate(conn => this%connections(i))
                if (conn%get_continuity() < 0) cycle
                if (present(geometric)) then
                    if (conn%is_geometric() .neqv. geometric) cycle
                end if
                if (conn%has_transition_jet() .or. conn%has_projective_jet()) then
                    call append_general_constraint_surface(&
                        this,&
                        conn,&
                        offsets,&
                        rowptr,&
                        col,&
                        val,&
                        row,&
                        pos)
                    cycle
                end if
                pa = conn%get_patch_a()
                pb = conn%get_patch_b()
                side_a = conn%get_side_a()
                side_b = conn%get_side_b()
                nca = this%patches(pa)%get_nc()
                ncb = this%patches(pb)%get_nc()
                dir_a = merge(1, 2, side_a <= u_max_)
                dir_b = merge(1, 2, side_b <= u_max_)
                dega = this%patches(pa)%get_degree(dir_a)
                degb = this%patches(pb)%get_degree(dir_b)
                normal_sign = merge(-1, 1, mod(side_a, 2) == mod(side_b, 2))
                knot_a = this%patches(pa)%get_knot(dir_a)
                knot_b = this%patches(pb)%get_knot(dir_b)
                xa = merge(knot_start(knot_a, nca(dir_a), dega), &
                    knot_end(knot_a, nca(dir_a), dega), mod(side_a, 2) == 1)
                xb = merge(knot_start(knot_b, ncb(dir_b), degb), &
                    knot_end(knot_b, ncb(dir_b), degb), mod(side_b, 2) == 1)
                range_a = knot_end(knot_a, nca(dir_a), dega) - knot_start(knot_a, nca(dir_a), dega)
                range_b = knot_end(knot_b, ncb(dir_b), degb) - knot_start(knot_b, ncb(dir_b), degb)
                call connection_weight_scale_surface(this, conn, sqrt(epsilon(1.0_rk)), weight_scale, compatible_weights)
                if (.not. compatible_weights) weight_scale = 1.0_rk
                allocate(ca(0:dega,0:conn%get_continuity()), cb(0:degb,0:conn%get_continuity()))
                allocate(first_a(0:conn%get_continuity()), first_b(0:conn%get_continuity()))
                do r = 0, conn%get_continuity()
                    call basis_bspline_der_order_active(xa, knot_a, nca(dir_a), dega, r, first_a(r), ca(:,r))
                    call basis_bspline_der_order_active(xb, knot_b, ncb(dir_b), degb, r, first_b(r), cb(:,r))
                end do
                call multipatch_chain_rule_coefficients(conn, normal_sign, chain)
                call boundary_index(nca, side_a, .false., ia)
                ntan = size(ia)
                max_layer_a = min(nca(dir_a) - 1, dega)
                max_layer_b = min(ncb(dir_b) - 1, degb)
                allocate(layer_a(ntan, max_layer_a + 1), layer_b(ntan, max_layer_b + 1))
                do layer = 0, max_layer_a
                    call boundary_layer_index(nca, side_a, layer, .false., ia)
                    layer_a(:, layer + 1) = ia
                end do
                do layer = 0, max_layer_b
                    call boundary_layer_index(ncb, side_b, layer, conn%is_reversed(), ib)
                    layer_b(:, layer + 1) = ib
                end do
                do r = 0, conn%get_continuity()
                    der_scale_a = range_a**r
                    do j = 1, ntan
                        row = row + 1
                        do layer = 0, max_layer_a
                            normal_a = merge(1 + layer, nca(dir_a) - layer, mod(side_a, 2) == 1)
                            active_a = normal_a - first_a(r)
                            if (active_a >= 0 .and. active_a <= dega .and. abs(der_scale_a*ca(active_a,r)) > eps) then
                                col(pos) = offsets(pa) + layer_a(j, layer + 1)
                                val(pos) = der_scale_a*ca(active_a,r)
                                pos = pos + 1
                            end if
                        end do
                        do layer = 0, max_layer_b
                            normal_b = merge(1 + layer, ncb(dir_b) - layer, mod(side_b, 2) == 1)
                            coefficient_b = 0.0_rk
                            do k = 0, r
                                active_b = normal_b - first_b(k)
                                if (active_b >= 0 .and. active_b <= degb) then
                                    coefficient_b = coefficient_b + chain(r,k)*range_b**k*cb(active_b,k)
                                end if
                            end do
                            if (abs(weight_scale*coefficient_b) > eps) then
                                col(pos) = offsets(pb) + layer_b(j, layer + 1)
                                val(pos) = -weight_scale*coefficient_b
                                pos = pos + 1
                            end if
                        end do
                        rowptr(row+1) = pos
                    end do
                end do
                deallocate(ca, cb, chain, first_a, first_b, layer_a, layer_b)
            end associate
        end do
        col = col(:pos-1)
        val = val(:pos-1)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Concatenate active tensor-product elements from all surface patches.
    !! Rows follow patch insertion order and direction-1-fastest local element
    !! order. Different local basis sizes are zero-padded to the largest row.
    !! `shared=.true.` applies direct sharing from the compact map; the default
    !! retains unmerged patch-offset numbering. General unclamped traces and
    !! rational traces with nonconstant projective scale still require
    !! `cmp_dof_constraint`.
    pure function cmp_elem_surface(this, shared) result(elem)
        class(nurbs_multipatch_surface), intent(in) :: this
        logical, intent(in), optional :: shared
        integer, allocatable :: elem(:,:)
        integer, allocatable :: offsets(:), map(:), mul1(:), mul2(:), nspan_patch(:,:), row_start(:)
        real(rk), allocatable :: knot1(:), knot2(:)
        integer :: a, i, row, ne, maxnen, nc(2), degree(2), nspan1, nspan2
        integer :: e1, e2, i1, i2, pref1, pref2, node, np
        logical :: use_shared

        use_shared = .false.
        if (present(shared)) use_shared = shared
        np = this%get_npatch()
        offsets = this%cmp_dof_offsets()
        ne = 0
        maxnen = 0
        allocate(nspan_patch(2, np), row_start(np + 1))
        row_start(1) = 0
        do i = 1, np
            nc = this%patches(i)%get_nc()
            degree = this%patches(i)%get_degree()
            knot1 = this%patches(i)%get_knot(1)
            knot2 = this%patches(i)%get_knot(2)
            nspan_patch(1, i) = active_span_count(knot1, nc(1), degree(1))
            nspan_patch(2, i) = active_span_count(knot2, nc(2), degree(2))
            row_start(i + 1) = row_start(i) + nspan_patch(1, i)*nspan_patch(2, i)
            maxnen = max(maxnen, (degree(1) + 1)*(degree(2) + 1))
        end do
        ne = row_start(np + 1)
        allocate(elem(ne, maxnen), source=0)
        if (use_shared) map = this%cmp_dof_map()
        do i = 1, np
            nc = this%patches(i)%get_nc()
            degree = this%patches(i)%get_degree()
            knot1 = this%patches(i)%get_knot(1)
            knot2 = this%patches(i)%get_knot(2)
            nspan1 = nspan_patch(1, i)
            nspan2 = nspan_patch(2, i)
            if (nspan1 < 1 .or. nspan2 < 1) cycle
            mul1 = active_knot_multiplicity(knot1, nc(1), degree(1))
            mul2 = active_knot_multiplicity(knot2, nc(2), degree(2))
            if (size(mul1) < nspan1 + 1 .or. size(mul2) < nspan2 + 1) cycle
            row = row_start(i)
            pref2 = degree(2) + 1
            do e2 = 1, nspan2
                pref1 = degree(1) + 1
                do e1 = 1, nspan1
                    row = row + 1
                    do i2 = 0, degree(2)
                        do i1 = 0, degree(1)
                            a = i2*(degree(1) + 1) + i1 + 1
                            node = ((pref2 - degree(2) + i2) - 1)*nc(1) + pref1 - degree(1) + i1
                            if (use_shared) then
                                elem(row,a) = map(offsets(i) + node)
                            else
                                elem(row,a) = offsets(i) + node
                            end if
                        end do
                    end do
                    if (e1 < nspan1) pref1 = pref1 + mul1(e1+1)
                end do
                if (e2 < nspan2) pref2 = pref2 + mul2(e2+1)
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the owning patch identifier for every concatenated element.
    !!
    !! Result ordering is identical to [[nurbs_multipatch_surface:cmp_elem]].
    pure function cmp_elem_patch_surface(this) result(patch_id)
        class(nurbs_multipatch_surface), intent(in) :: this
        integer, allocatable :: patch_id(:)
        integer, allocatable :: nloc(:), row_start(:)
        integer :: i, e, ne, np

        np = this%get_npatch()
        allocate(nloc(np), row_start(np + 1))
        row_start(1) = 0
        do i = 1, np
            nloc(i) = active_span_count(this%patches(i)%get_knot(1), this%patches(i)%get_nc(1), &
                this%patches(i)%get_degree(1)) * active_span_count(this%patches(i)%get_knot(2), &
                this%patches(i)%get_nc(2), this%patches(i)%get_degree(2))
            row_start(i + 1) = row_start(i) + nloc(i)
        end do
        ne = row_start(np + 1)
        allocate(patch_id(ne))
        do i = 1, np
            do e = row_start(i) + 1, row_start(i + 1)
                patch_id(e) = i
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return each concatenated element's one-based patch-local identifier.
    !!
    !! Result ordering is identical to [[nurbs_multipatch_surface:cmp_elem]].
    pure function cmp_elem_local_surface(this) result(local_id)
        class(nurbs_multipatch_surface), intent(in) :: this
        integer, allocatable :: local_id(:)
        integer, allocatable :: nloc(:), row_start(:)
        integer :: i, e, ne, np

        np = this%get_npatch()
        allocate(nloc(np), row_start(np + 1))
        row_start(1) = 0
        do i = 1, np
            nloc(i) = active_span_count(this%patches(i)%get_knot(1), this%patches(i)%get_nc(1), &
                this%patches(i)%get_degree(1)) * active_span_count(this%patches(i)%get_knot(2), &
                this%patches(i)%get_nc(2), this%patches(i)%get_degree(2))
            row_start(i + 1) = row_start(i) + nloc(i)
        end do
        ne = row_start(np + 1)
        allocate(local_id(ne))
        do i = 1, np
            do e = 1, nloc(i)
                local_id(row_start(i) + e) = e
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Validate patch states and every oriented surface interface.
    !!
    !! A valid multipatch contains at least one complete, internally consistent
    !! surface patch. Checks also include identifiers, sides, continuity order, transition jets,
    !! tangential trace spaces, physical endpoint traces, homogeneous
    !! normal-derivative residuals in the compact or multivariate model,
    !! embedding dimensions, and rational projective compatibility. `tol` is one
    !! absolute threshold for assembled homogeneous residuals; the default is
    !! `sqrt(epsilon(rk))`. Each component of \(A\mathbf H\) is compared
    !! directly with `tol`. Weight proportionality
    !! uses `max(tol,128*epsilon(rk))*max(abs(wA),abs(lambda*wB))` per pair.
    !! `tol>=0` is an unchecked caller precondition. The test conservatively
    !! requires equal
    !! whole-patch `is_rational` classifications when no explicit projective
    !! jet is supplied. A full projective jet validates the complete selected
    !! homogeneous-lift relation instead.
    pure function is_valid_surface(this, tol) result(ok)
        class(nurbs_multipatch_surface), intent(in) :: this
        real(rk), intent(in), optional :: tol
            !! Optional nonnegative absolute coordinate/residual tolerance; default `sqrt(epsilon(rk))`.
        logical :: ok
        real(rk) :: eps
        integer :: i, j, d, dega, degb, dir_a, dir_b
        integer :: pa, pb, dim_a, dim_b, normal_sign
        integer :: row, a, c, total, ncomp
        integer, allocatable :: offsets(:), rowptr(:), col(:)
        real(rk), allocatable :: csrval(:), dof(:,:), Xc(:,:), Wc(:)
        real(rk) :: residual, weight_scale
        logical :: rational_a, rational_b, compatible_weights, use_homogeneous

        if (.not. this%err%ok) then
            ok = .false.
            return
        end if
        if (this%npatch < 1 .or. .not. allocated(this%patches)) then
            ok = .false.
            return
        end if
        if (this%npatch > size(this%patches) .or. this%nconnection < 0) then
            ok = .false.
            return
        end if
        if (this%nconnection > 0) then
            if (.not. allocated(this%connections)) then
                ok = .false.
                return
            end if
            if (this%nconnection > size(this%connections)) then
                ok = .false.
                return
            end if
        end if
        do i = 1, this%npatch
            if (.not. this%patches(i)%is_valid()) then
                ok = .false.
                return
            end if
        end do

        eps = sqrt(epsilon(1.0_rk))
        if (present(tol)) eps = tol
        ok = .true.
        do i = 1, this%get_nconnection()
            associate(conn => this%connections(i))
                pa = conn%get_patch_a()
                pb = conn%get_patch_b()
                if (pa < 1 .or. pa > this%get_npatch()) ok = .false.
                if (pb < 1 .or. pb > this%get_npatch()) ok = .false.
                if (.not. ok) return
                if (conn%get_side_a() < u_min_ .or. conn%get_side_a() > v_max_ .or. &
                    conn%get_side_b() < u_min_ .or. conn%get_side_b() > v_max_) then
                    ok = .false.; return
                end if
                if (conn%get_continuity() < -1) then
                    ok = .false.; return
                end if
                normal_sign = merge(-1, 1, mod(conn%get_side_a(), 2) == mod(conn%get_side_b(), 2))
                if (.not. multipatch_valid_reparameterization(conn, normal_sign, 2)) then
                    ok = .false.; return
                end if
                dir_a = merge(1, 2, conn%get_side_a() <= u_max_)
                dir_b = merge(1, 2, conn%get_side_b() <= u_max_)
                dega = this%patches(pa)%get_degree(dir_a)
                degb = this%patches(pb)%get_degree(dir_b)
                if (conn%get_continuity() > min(dega, degb)) then
                    ok = .false.; return
                end if
                if (conn%get_continuity() >= 0) then
                    if (.not. connection_trace_space_surface(this, conn, eps)) then
                        ok = .false.; return
                    end if
                    dim_a = size(this%patches(pa)%get_Xc(1))
                    dim_b = size(this%patches(pb)%get_Xc(1))
                    if (dim_a /= dim_b) then
                        ok = .false.; return
                    end if
                    if (.not. conn%has_projective_jet()) then
                        rational_a = this%patches(pa)%is_rational()
                        rational_b = this%patches(pb)%is_rational()
                        if (rational_a .neqv. rational_b) then
                            ok = .false.; return
                        end if
                        call connection_weight_scale_surface(this, conn, eps, weight_scale, compatible_weights)
                        if (.not. compatible_weights) then
                            ok = .false.; return
                        end if
                    end if
                end if
            end associate
        end do

        offsets = this%cmp_dof_offsets()
        total = offsets(size(offsets))
        if (total < 1) then
            ok = .false.
            return
        end if
        ncomp = 0
        use_homogeneous = this%is_rational()
        do i = 1, this%get_nconnection()
            if (this%connections(i)%has_projective_jet()) use_homogeneous = .true.
        end do
        do i = 1, this%get_npatch()
            ncomp = max(ncomp, size(this%patches(i)%get_Xc(1)))
        end do
        if (use_homogeneous) ncomp = ncomp + 1
        allocate(dof(total, ncomp), source=0.0_rk)
        do i = 1, this%get_npatch()
            Xc = this%patches(i)%get_Xc()
            if (use_homogeneous .and. this%patches(i)%is_rational()) Wc = this%patches(i)%get_Wc()
            do a = 1, size(Xc, 1)
                if (use_homogeneous .and. this%patches(i)%is_rational()) then
                    do d = 1, size(Xc, 2)
                        dof(offsets(i)+a, d) = Xc(a, d)*Wc(a)
                    end do
                    dof(offsets(i)+a, ncomp) = Wc(a)
                else
                    do d = 1, size(Xc, 2)
                        dof(offsets(i)+a, d) = Xc(a, d)
                    end do
                    if (use_homogeneous) dof(offsets(i)+a, ncomp) = 1.0_rk
                end if
            end do
        end do
        call this%cmp_dof_constraint(rowptr, col, csrval)
        do row = 1, size(rowptr)-1
            do c = 1, ncomp
                residual = 0.0_rk
                do j = rowptr(row), rowptr(row+1)-1
                    residual = residual + csrval(j)*dof(col(j), c)
                end do
                if (abs(residual) > eps) then
                    ok = .false.
                    return
                end if
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Sample every surface patch with common directional sampling arguments.
    !! Arguments are forwarded unchanged to [[nurbs_surface:create]]. Any patch
    !! failure is converted to multipatch diagnostic code 205.
    pure subroutine create_surface(this, res1, res2, Xt1, Xt2, Xt)
        class(nurbs_multipatch_surface), intent(inout) :: this
        integer, intent(in), optional :: res1, res2
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:)
        real(rk), intent(in), contiguous, optional :: Xt(:,:)
        integer :: i

        if (.not. this%err%ok) return

        if (this%get_npatch() < 1) then
            call this%err%set(&
                code       = 205,&
                severity   = 1,&
                category   = "forcad_multipatch_surface",&
                message    = "No patches are available.",&
                location   = "create",&
                suggestion = "Add at least one patch before calling create().")
            return
        end if
        do i = 1, this%get_npatch()
            call this%patches(i)%create(res1=res1, res2=res2, Xt1=Xt1, Xt2=Xt2, Xt=Xt)
            if (.not. this%patches(i)%err%ok) then
                call this%err%set(&
                    code       = 205,&
                    severity   = 1,&
                    category   = "forcad_multipatch_surface",&
                    message    = "Patch create failed.",&
                    location   = "create",&
                    suggestion = "Inspect the patch error state and fix the patch inputs.")
                return
            end if
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report whether at least one active patch uses nonuniform rational weights.
    pure elemental logical function is_rational_surface(this) result(r)
        class(nurbs_multipatch_surface), intent(in) :: this
        integer :: i

        r = .false.
        do i = 1, this%get_npatch()
            if (this%patches(i)%is_rational()) then
                r = .true.
                return
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export every patch control net to a separate legacy VTK file.
    !!
    !! Patch `i` is written as `<prefix>_Xc_<i>.vtk`.
    impure subroutine export_Xc_surface(this, prefix, encoding)
        class(nurbs_multipatch_surface), intent(inout) :: this
        character(len=*), intent(in) :: prefix
        character(len=*), intent(in), optional :: encoding
        character(len=512) :: vtkfile
        integer :: i

        if (.not. this%err%ok) return

        do i = 1, this%get_npatch()
            write(vtkfile, "(a,a,i0,a)") trim(prefix), "_Xc_", i, ".vtk"
            call this%patches(i)%export_Xc(vtkfile, encoding=encoding)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export every patch's cached physical samples to a separate VTK file.
    !!
    !! Patch `i` is written as `<prefix>_Xg_<i>.vtk`.
    impure subroutine export_Xg_surface(this, prefix, encoding)
        class(nurbs_multipatch_surface), intent(inout) :: this
        character(len=*), intent(in) :: prefix
        character(len=*), intent(in), optional :: encoding
        character(len=512) :: vtkfile
        integer :: i

        if (.not. this%err%ok) return

        do i = 1, this%get_npatch()
            write(vtkfile, "(a,a,i0,a)") trim(prefix), "_Xg_", i, ".vtk"
            call this%patches(i)%export_Xg(vtkfile, encoding=encoding)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export parameter lines mapped into every surface patch.
    impure subroutine export_Xth_in_Xg_surface(this, prefix, res, encoding)
        class(nurbs_multipatch_surface), intent(inout) :: this
        character(len=*), intent(in) :: prefix
        integer, intent(in), optional :: res
        character(len=*), intent(in), optional :: encoding
        character(len=512) :: vtkfile
        integer :: i

        if (.not. this%err%ok) return

        do i = 1, this%get_npatch()
            write(vtkfile, "(a,a,i0,a)") trim(prefix), "_Xth_", i, ".vtk"
            call this%patches(i)%export_Xth_in_Xg(vtkfile, res=res, encoding=encoding)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Apply one x-y-z Euler rotation in degrees to every control net.
    !!
    !! Cached physical samples are not changed.
    pure subroutine rotate_Xc_surface(this, alpha, beta, theta)
        class(nurbs_multipatch_surface), intent(inout) :: this
        real(rk), intent(in) :: alpha, beta, theta
        integer :: i

        if (.not. this%err%ok) return

        do i = 1, this%get_npatch()
            call this%patches(i)%rotate_Xc(alpha, beta, theta)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Apply one x-y-z Euler rotation in degrees to every sampled patch.
    !!
    !! Control points are not changed.
    pure subroutine rotate_Xg_surface(this, alpha, beta, theta)
        class(nurbs_multipatch_surface), intent(inout) :: this
        real(rk), intent(in) :: alpha, beta, theta
        integer :: i

        if (.not. this%err%ok) return

        do i = 1, this%get_npatch()
            call this%patches(i)%rotate_Xg(alpha, beta, theta)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Translate every patch control net by physical vector `vec`.
    !!
    !! Cached physical samples are not changed.
    pure subroutine translate_Xc_surface(this, vec)
        class(nurbs_multipatch_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: vec(:)
        integer :: i

        if (.not. this%err%ok) return

        do i = 1, this%get_npatch()
            call this%patches(i)%translate_Xc(vec)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Translate every patch's cached samples by physical vector `vec`.
    !!
    !! Control points are not changed.
    pure subroutine translate_Xg_surface(this, vec)
        class(nurbs_multipatch_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: vec(:)
        integer :: i

        if (.not. this%err%ok) return

        do i = 1, this%get_npatch()
            call this%patches(i)%translate_Xg(vec)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Display previously exported multipatch VTK files using PyVista.
    !! File arguments may be exact paths or glob patterns and are never created
    !! or overwritten by this procedure. Common scalar `Xg` point fields and
    !! colormaps are selected in the viewer.
    impure subroutine show_surface(this, vtkfile_Xc, vtkfile_Xg, vtkfile_Xth_in_Xg)
        class(nurbs_multipatch_surface), intent(inout) :: this
        character(len=*), intent(in) :: vtkfile_Xc
            !! Existing control-geometry VTK path or glob pattern.
        character(len=*), intent(in) :: vtkfile_Xg
            !! Existing sampled-geometry VTK path or glob pattern.
        character(len=*), intent(in), optional :: vtkfile_Xth_in_Xg
            !! Optional existing parameter-grid VTK path or glob pattern.

        if (.not. this%err%ok) return
#ifndef NOSHOW_PYVISTA
        call show_pyvista_multipatch(&
            vtkfile_Xc        = vtkfile_Xc,&
            vtkfile_Xg        = vtkfile_Xg,&
            vtkfile_Xth_in_Xg = vtkfile_Xth_in_Xg,&
            rank_name         = "multipatch surface")
#endif
    end subroutine
    !===============================================================================
end module forcad_multipatch_surface
