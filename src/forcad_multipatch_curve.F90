!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Topology, continuity constraints, and assembly data for curve patches.
!!
!! Patches remain complete [[nurbs_curve]] objects. The container records
!! oriented endpoint connections; [[nurbs_multipatch_curve:is_valid]] verifies
!! their physical and derivative traces. Two complementary global algebraic
!! objects are exposed:
!!
!! * [[nurbs_multipatch_curve:cmp_dof_map]] merges corresponding \(C^0\)
!!   interface controls only where both endpoint bases are interpolatory;
!! * [[nurbs_multipatch_curve:cmp_dof_constraint]] returns sparse equations for
!!   derivative continuity through any requested order
!!   \(0\le n\le\min(p_A,p_B)\).
!!
!! Element rows are concatenated in patch insertion order. Within a patch they
!! retain the single-patch active-span ordering. This module stores topology and
!! algebraic constraints; it does not assemble a PDE matrix or choose boundary
!! conditions for an analysis.
!!
!! If patch \(k\) owns \(n_k\) control variables, unmerged numbering is
!!
!! \[
!! \mathbf q=
!! [(\mathbf q^{(1)})^T,\ldots,(\mathbf q^{(n_p)})^T]^T,\qquad
!! o_k=1+\sum_{\ell<k}n_\ell .
!! \]
!!
!! The degree-of-freedom map identifies paired order-zero controls only when
!! each endpoint basis contains one unit coefficient. A general unclamped trace
!! is a linear combination of controls and remains in the sparse homogeneous
!! system \(\mathbf A\mathbf q=\mathbf0\), which can be eliminated, constrained
!! with multipliers, or supplied to another analysis package. Call `is_valid`
!! before using either object. For polynomial fields, \(\mathbf q\) contains
!! ordinary B-spline coefficients. For rational geometry, apply the matrix to
!! homogeneous controls \((w\mathbf P,w)\). A rational scalar field with
!! coefficients \(q_i\) uses weighted numerator controls \(w_iq_i\), not the
!! unweighted coefficient vector, while the denominator must satisfy the same
!! projective interface compatibility.
module forcad_multipatch_curve

    use forcad_kinds, only: rk
    use forcad_multipatch, only: multipatch_chain_rule_coefficients, multipatch_connection, &
        multipatch_growth_capacity, multipatch_projective_weight_scale, multipatch_valid_reparameterization
    use forcad_nurbs_curve, only: nurbs_curve
    use forcad_utils, only: active_knot_multiplicity, active_span_count, basis_bspline_der_order_active, &
        boundary_layer_index, disjoint_set_map, disjoint_set_union, interpolatory_boundary_layer, &
        knot_end, knot_start, show_pyvista_multipatch
    use fordebug, only: debug

    implicit none

    private
    public :: nurbs_multipatch_curve

    integer, parameter :: left_  = 1
    integer, parameter :: right_ = 2

    !===============================================================================
    !> Ordered collection of connected NURBS curve patches.
    !!
    !! Valid side labels accepted by `connect` are `left`, `right`, `min`, and
    !! `max` (case variants shown in the implementation are accepted).
    !! `continuity=-1` retains an interface without continuity equations.
    !! Nonnegative orders appear in the sparse constraint matrix; direct
    !! endpoint sharing is additionally available for interpolatory endpoints.
    !!
    !! Diagnostic codes are: 201 invalid patch identifier; 202 invalid side;
    !! 203 invalid continuity metadata; 204 incompatible interface geometry or
    !! trace space; and 205 invalid patch state.
    type nurbs_multipatch_curve
        type(nurbs_curve), allocatable, private :: patches(:) !! Owned patch copies with spare append capacity.
        type(multipatch_connection), allocatable, private :: connections(:) !! Validated interface metadata.
        integer, private :: npatch = 0 !! Number of active patches.
        integer, private :: nconnection = 0 !! Number of active connections.

        type(debug) :: err !! Recoverable diagnostic state for the most recent failed operation.
    contains
        procedure :: add_patch => add_patch_curve !!> Append a validated patch copy and optionally return its identifier.
        procedure :: connect => connect_curve !!> Validate connection metadata and record an oriented endpoint interface.
        procedure :: create => create_curve !!> Sample every patch using common or patch-specific parameters.
        procedure :: finalize => clear_curve !!> Release patches and connections and reset counts.
        procedure :: get_npatch => get_npatch_curve !!> Return the active patch count.
        procedure :: get_nconnection => get_nconnection_curve !!> Return the active connection count.
        procedure :: get_patch => get_patch_curve !!> Return a copy of one patch.
        procedure :: get_connection => get_connection_curve !!> Return a copy of one connection.
        procedure :: is_valid => is_valid_curve !!> Validate patch states, interfaces, orientations, and requested continuity.
        procedure :: is_rational => is_rational_curve !!> Report whether any patch has nonuniform rational weights.
        procedure :: cmp_dof_offsets => cmp_dof_offsets_curve !!> Return prefix offsets for unmerged patch control variables.
        procedure :: cmp_dof_map => cmp_dof_map_curve !!> Compact directly shareable \(C^0\) controls.
        procedure :: cmp_dof_constraint => cmp_dof_constraint_curve !!> Assemble derivative-continuity equations in CSR form.
        procedure :: cmp_elem => cmp_elem_curve !!> Concatenate patch element connectivity with optional shared numbering.
        procedure :: get_elem => cmp_elem_curve !!> Alias for [[nurbs_multipatch_curve:cmp_elem]].
        procedure :: cmp_elem_patch => cmp_elem_patch_curve !!> Return the owning patch identifier for every global element.
        procedure :: cmp_elem_local => cmp_elem_local_curve !!> Return the patch-local identifier for every global element.
        procedure :: export_Xc => export_Xc_curve !!> Export each patch control polygon to a separate VTK file.
        procedure :: export_Xg => export_Xg_curve !!> Export each patch's cached geometry to a separate VTK file.
        procedure :: export_Xth_in_Xg => export_Xth_in_Xg_curve !!> Export embedded parameter lines for every patch.
        procedure :: rotate_Xc => rotate_Xc_curve !!> Rotate all patch control points.
        procedure :: rotate_Xg => rotate_Xg_curve !!> Rotate all cached patch geometry points.
        procedure :: translate_Xc => translate_Xc_curve !!> Translate all patch control points.
        procedure :: translate_Xg => translate_Xg_curve !!> Translate all cached patch geometry points.
        procedure :: show => show_curve !!> Display selected patch representations together with PyVista.
    end type
    !===============================================================================

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Append a copy of a valid curve patch.
    pure subroutine add_patch_curve(this, patch, id)
        class(nurbs_multipatch_curve), intent(inout) :: this
            !! Multipatch container.
        type(nurbs_curve), intent(in) :: patch
            !! Complete, internally consistent curve geometry with a clear diagnostic state.
        integer, intent(out), optional :: id
            !! Assigned one-based patch identifier; zero is returned on failure.
        type(nurbs_curve), allocatable :: tmp(:)
        integer :: new_capacity

        if (present(id)) id = 0
        if (.not. this%err%ok) return

        if (.not. patch%err%ok) then
            call this%err%set(&
                code       = 205,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
                message    = "Patch has an active error state.",&
                location   = "add_patch",&
                suggestion = "Reset or fix the patch error state before adding it to a multipatch curve.")
            return
        end if
        if (.not. patch%is_valid()) then
            call this%err%set(&
                code       = 205,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
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
    !> Record two curve endpoints with parametric \(C^n\) or geometric \(G^n\) metadata.
    !! `continuity` defaults to zero and `geometric` defaults to `.false.`.
    !! A G^n connection requires `reparameterization(k)` to
    !! contain the kth derivative of patch B's normalized parameter with
    !! respect to patch A's normalized local parameter at the interface. The
    !! scalar transition is the complete local reparameterization model for a
    !! curve endpoint. `reverse` has no mathematical effect because an endpoint
    !! has no tangential coordinate; it is stored only for API symmetry.
    !!
    !! This routine validates identifiers, side labels, continuity order, and
    !! transition-jet shape before recording the connection. Endpoint
    !! coincidence, derivative residuals, coordinate dimensions, and rational
    !! projective compatibility are checked later by
    !! [[nurbs_multipatch_curve:is_valid]].
    !!
    !! @warning A recorded connection is not proof of a valid geometric join;
    !! call [[nurbs_multipatch_curve:is_valid]] before assembly or analysis.
    !! @endwarning
    pure subroutine connect_curve(this, patch_a, side_a, patch_b, side_b, continuity, reverse, &
        geometric, reparameterization)
        class(nurbs_multipatch_curve), intent(inout) :: this
            !! Multipatch container.
        integer, intent(in) :: patch_a
            !! One-based first patch identifier.
        integer, intent(in) :: patch_b
            !! One-based second patch identifier.
        character(len=*), intent(in) :: side_a
            !! First endpoint label: `left`/`min` or `right`/`max`.
        character(len=*), intent(in) :: side_b
            !! Second endpoint label.
        integer, intent(in), optional :: continuity
            !! Requested order, default zero; `-1` disables constraints.
        logical, intent(in), optional :: reverse
            !! Rank-consistent metadata flag; it has no effect for a single endpoint.
        logical, intent(in), optional :: geometric
            !! Select geometric rather than parametric continuity.
        real(rk), intent(in), contiguous, optional :: reparameterization(:)
            !! Normalized transition derivatives for \(G^n\).
        type(multipatch_connection) :: conn
        type(multipatch_connection), allocatable :: tmp(:)
        integer :: side_a_id, side_b_id, continuity_, new_capacity, normal_sign
        logical :: reverse_, geometric_

        if (.not. this%err%ok) return

        if (patch_a < 1 .or. patch_a > this%get_npatch()) then
            call this%err%set(&
                code       = 201,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
                message    = "Invalid first patch id.",&
                location   = "connect",&
                suggestion = "Use a patch id from add_patch's optional id argument.")
            return
        end if
        if (patch_b < 1 .or. patch_b > this%get_npatch()) then
            call this%err%set(&
                code       = 201,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
                message    = "Invalid second patch id.",&
                location   = "connect",&
                suggestion = "Use a patch id from add_patch's optional id argument.")
            return
        end if
        side_a_id = 0
        side_b_id = 0
        select case(trim(side_a))
        case("left", "LEFT", "Left", "min", "MIN", "Min")
            side_a_id = left_
        case("right", "RIGHT", "Right", "max", "MAX", "Max")
            side_a_id = right_
        case default
            continue
        end select
        select case(trim(side_b))
        case("left", "LEFT", "Left", "min", "MIN", "Min")
            side_b_id = left_
        case("right", "RIGHT", "Right", "max", "MAX", "Max")
            side_b_id = right_
        case default
            continue
        end select
        if (side_a_id == 0 .or. side_b_id == 0) then
            call this%err%set(&
                code       = 202,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
                message    = "Invalid curve side label.",&
                location   = "connect",&
                suggestion = "Use left/right or min/max.")
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
                category   = "forcad_multipatch_curve",&
                message    = "Invalid continuity value.",&
                location   = "connect",&
                suggestion = "Use -1 for discontinuous or Cn with n >= 0.")
            return
        end if
        if (continuity_ > min(this%patches(patch_a)%get_degree(), this%patches(patch_b)%get_degree())) then
            call this%err%set(&
                code       = 203,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
                message    = "Continuity is higher than the connected patch degree.",&
                location   = "connect",&
                suggestion = "Use a continuity not larger than the minimum degree of the two curve patches.")
            return
        end if
        reverse_ = .false.
        if (present(reverse)) reverse_ = reverse
        call conn%set(&
            patch_a           = patch_a,&
            side_a            = side_a_id,&
            patch_b           = patch_b,&
            side_b            = side_b_id,&
            continuity        = continuity_,&
            reverse           = reverse_,&
            geometric          = geometric_,&
            reparameterization = reparameterization)
        normal_sign = merge(-1, 1, mod(side_a_id, 2) == mod(side_b_id, 2))
        if (.not. multipatch_valid_reparameterization(conn, normal_sign)) then
            call this%err%set(&
                code       = 203,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
                message    = "Invalid continuity reparameterization.",&
                location   = "connect",&
                suggestion = "For G^n, pass n finite normalized derivatives with a "//&
                    "first-derivative sign consistent with the sides.")
            return
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
    !> Release all owned curve patches and connections.
    !!
    !! Active counts are reset to zero and the multipatch diagnostic is
    !! preserved.
    pure elemental subroutine clear_curve(this)
        class(nurbs_multipatch_curve), intent(inout) :: this

        if (allocated(this%patches)) deallocate(this%patches)
        if (allocated(this%connections)) deallocate(this%connections)
        this%npatch = 0
        this%nconnection = 0
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the number of active curve patches.
    pure elemental integer function get_npatch_curve(this) result(n)
        class(nurbs_multipatch_curve), intent(in) :: this

        n = this%npatch
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the number of recorded endpoint connections.
    pure elemental integer function get_nconnection_curve(this) result(n)
        class(nurbs_multipatch_curve), intent(in) :: this

        n = this%nconnection
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return a value copy of patch `i`.
    !!
    !! An out-of-range one-based identifier returns a default empty patch.
    pure function get_patch_curve(this, i) result(patch)
        class(nurbs_multipatch_curve), intent(in) :: this
        integer, intent(in) :: i
        type(nurbs_curve) :: patch

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
    pure function get_connection_curve(this, i) result(conn)
        class(nurbs_multipatch_curve), intent(in) :: this
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
    !! The result has size `get_npatch()+1`, starts at zero, and the final entry
    !! is the total number of unmerged patch control points.
    pure function cmp_dof_offsets_curve(this) result(offsets)
        class(nurbs_multipatch_curve), intent(in) :: this
        integer, allocatable :: offsets(:)
        integer :: i

        allocate(offsets(this%get_npatch()+1))
        offsets(1) = 0
        do i = 1, this%get_npatch()
            offsets(i+1) = offsets(i) + this%patches(i)%get_nc()
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build compact global identifiers for directly shareable endpoints.
    !!
    !! The input position is the unmerged index defined by
    !! [[nurbs_multipatch_curve:cmp_dof_offsets]]. A connection with continuity
    !! at least zero merges one control from each side only when both endpoint
    !! basis vectors are exactly one-hot. General unclamped \(C^0\) traces are
    !! linear combinations and remain unmerged; obtain their equations from
    !! `cmp_dof_constraint`. This method trusts recorded topology and does not
    !! call `is_valid`.
    !!
    !! Disjoint-set implementation reference: R. E. Tarjan, "Efficiency of a
    !! Good But Not Linear Set Union Algorithm," *JACM* 22 (1975), 215--225.
    !! [doi:10.1145/321879.321884](https://doi.org/10.1145/321879.321884).
    pure function cmp_dof_map_curve(this) result(map)
        class(nurbs_multipatch_curve), intent(in) :: this
        integer, allocatable :: map(:)
        integer, allocatable :: parent(:), offsets(:), ia(:), ib(:)
        real(rk), allocatable :: knot_a(:), knot_b(:)
        integer :: i, j, total, pa, pb, nca, ncb, dega, degb, layer_a, layer_b

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
                dega = this%patches(pa)%get_degree()
                degb = this%patches(pb)%get_degree()
                knot_a = this%patches(pa)%get_knot()
                knot_b = this%patches(pb)%get_knot()
                layer_a = interpolatory_boundary_layer(knot_a, nca, dega, conn%get_side_a())
                layer_b = interpolatory_boundary_layer(knot_b, ncb, degb, conn%get_side_b())
                if (layer_a < 0 .or. layer_b < 0) cycle
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
    !> Check projective compatibility of rational endpoint traces.
    !!
    !! The endpoint denominator is evaluated from all nonzero normal basis
    !! values. Compatible traces satisfy
    !! \(\widehat w_A=\lambda\widehat w_B\) within `tol`;
    !! `weight_scale` returns \(\lambda\). Polynomial traces use unit scale.
    pure subroutine connection_weight_scale_curve(this, conn, tol, weight_scale, ok)
        class(nurbs_multipatch_curve), intent(in) :: this
        type(multipatch_connection), intent(in) :: conn
        real(rk), intent(in) :: tol
        real(rk), intent(out) :: weight_scale
        logical, intent(out) :: ok
        real(rk), allocatable :: basis_a(:), basis_b(:), knot_a(:), knot_b(:)
        real(rk) :: trace_a(1), trace_b(1), xa, xb
        integer :: i, first_a, first_b, pa, pb, nca, ncb, dega, degb
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
        dega = this%patches(pa)%get_degree()
        degb = this%patches(pb)%get_degree()
        knot_a = this%patches(pa)%get_knot()
        knot_b = this%patches(pb)%get_knot()
        xa = merge(knot_start(knot_a, nca, dega), knot_end(knot_a, nca, dega), &
            mod(conn%get_side_a(), 2) == 1)
        xb = merge(knot_start(knot_b, ncb, degb), knot_end(knot_b, ncb, degb), &
            mod(conn%get_side_b(), 2) == 1)
        allocate(basis_a(0:dega), basis_b(0:degb))
        call basis_bspline_der_order_active(xa, knot_a, nca, dega, 0, first_a, basis_a)
        call basis_bspline_der_order_active(xb, knot_b, ncb, degb, 0, first_b, basis_b)

        trace_a(1) = 0.0_rk
        do i = 0, dega
            if (first_a+i >= 1 .and. first_a+i <= nca) then
                trace_a(1) = trace_a(1) + basis_a(i)*this%patches(pa)%get_Wc(first_a+i)
            end if
        end do
        trace_b(1) = 0.0_rk
        do i = 0, degb
            if (first_b+i >= 1 .and. first_b+i <= ncb) then
                trace_b(1) = trace_b(1) + basis_b(i)*this%patches(pb)%get_Wc(first_b+i)
            end if
        end do
        call multipatch_projective_weight_scale(trace_a, trace_b, tol, weight_scale, ok)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Assemble sparse \(C^n\) and \(G^n\) interface constraints.
    !! When `geometric` is present, only connections of that kind are
    !! included. \(G^n\) rows use the complete arbitrary-order Faa di Bruno chain
    !! rule for the stored normalized normal reparameterization jet.
    !!
    !! The returned arrays encode a one-based compressed sparse row matrix
    !! \(A\): nonzeros of row `i` occupy
    !! `rowptr(i):rowptr(i+1)-1`, and continuity is imposed as \(Aq=0\).
    !! Columns use unmerged patch-concatenated numbering. An absent filter
    !! includes both \(C^n\) and \(G^n\) connections. Apply the raw matrix to
    !! polynomial B-spline coefficients or to each homogeneous geometry column.
    !! For a rational scalar field with control coefficients \(q_i\), apply it
    !! to \(w_iq_i\) (equivalently, use \(A\operatorname{diag}(w)\mathbf q\));
    !! applying \(A\) directly to the unweighted coefficients is generally not
    !! sufficient.
    !!
    !! In normalized boundary coordinates the row for derivative order \(r\)
    !! represents
    !!
    !! \[
    !! \partial_{\eta_A}^{\,r}f_A
    !! -\sum_{k=0}^{r}B_{r,k}(\phi',\ldots)
    !!  \partial_{\eta_B}^{\,k}f_B=0.
    !! \]
    !!
    !! Parameter-range powers convert derivatives of the stored knot
    !! coordinates to normalized derivatives. Compatible rational traces also
    !! include their single projective boundary-weight scale. Call `is_valid`
    !! first; this assembly routine assumes compatible recorded connections.
    !! Coefficients with absolute value at most `128*epsilon(rk)` are omitted
    !! from the returned CSR storage. The result is therefore a numerically
    !! sparsified floating-point matrix, not a symbolic exact matrix.
    !!
    pure subroutine cmp_dof_constraint_curve(this, rowptr, col, val, geometric)
        class(nurbs_multipatch_curve), intent(in) :: this
            !! Valid multipatch curve.
        integer, allocatable, intent(out) :: rowptr(:)
            !! CSR row pointers of size `nrow+1`, starting at one.
        integer, allocatable, intent(out) :: col(:)
            !! One-based unmerged global column identifiers.
        real(rk), allocatable, intent(out) :: val(:)
            !! Constraint coefficients corresponding to `col`.
        logical, intent(in), optional :: geometric
            !! Optional connection-kind filter: false for \(C^n\), true for \(G^n\).
        real(rk), allocatable :: ca(:,:), cb(:,:), chain(:,:), knot_a(:), knot_b(:)
        real(rk) :: xa, xb, weight_scale, range_a, range_b, der_scale_a, coefficient_b, eps
        integer, allocatable :: first_a(:), first_b(:), offsets(:)
        integer :: i, k, r, row, pos, nrow, nnz, pa, pb, side_a, side_b, normal_sign
        integer :: nca, ncb, dega, degb, layer, idxa, idxb, active_b
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
                pa = conn%get_patch_a()
                pb = conn%get_patch_b()
                dega = this%patches(pa)%get_degree()
                degb = this%patches(pb)%get_degree()
                nrow = nrow + conn%get_continuity() + 1
                nnz = nnz + (conn%get_continuity() + 1)*(&
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
        do i = 1, this%get_nconnection()
            associate(conn => this%connections(i))
                if (conn%get_continuity() < 0) cycle
                if (present(geometric)) then
                    if (conn%is_geometric() .neqv. geometric) cycle
                end if
                pa = conn%get_patch_a()
                pb = conn%get_patch_b()
                side_a = conn%get_side_a()
                side_b = conn%get_side_b()
                normal_sign = merge(-1, 1, mod(side_a, 2) == mod(side_b, 2))
                nca = this%patches(pa)%get_nc()
                ncb = this%patches(pb)%get_nc()
                dega = this%patches(pa)%get_degree()
                degb = this%patches(pb)%get_degree()
                knot_a = this%patches(pa)%get_knot()
                knot_b = this%patches(pb)%get_knot()
                xa = merge(knot_start(knot_a, nca, dega), knot_end(knot_a, nca, dega), side_a == left_)
                xb = merge(knot_start(knot_b, ncb, degb), knot_end(knot_b, ncb, degb), side_b == left_)
                range_a = knot_end(knot_a, nca, dega) - knot_start(knot_a, nca, dega)
                range_b = knot_end(knot_b, ncb, degb) - knot_start(knot_b, ncb, degb)
                call connection_weight_scale_curve(this, conn, sqrt(epsilon(1.0_rk)), weight_scale, compatible_weights)
                if (.not. compatible_weights) weight_scale = 1.0_rk
                allocate(ca(0:dega,0:conn%get_continuity()), cb(0:degb,0:conn%get_continuity()))
                allocate(first_a(0:conn%get_continuity()), first_b(0:conn%get_continuity()))
                do r = 0, conn%get_continuity()
                    call basis_bspline_der_order_active(xa, knot_a, nca, dega, r, first_a(r), ca(:,r))
                    call basis_bspline_der_order_active(xb, knot_b, ncb, degb, r, first_b(r), cb(:,r))
                end do
                call multipatch_chain_rule_coefficients(conn, normal_sign, chain)
                eps = 128.0_rk*epsilon(1.0_rk)
                do r = 0, conn%get_continuity()
                    der_scale_a = range_a**r
                    row = row + 1
                    do layer = 0, dega
                        idxa = first_a(r) + layer
                        if (idxa >= 1 .and. idxa <= nca .and. abs(der_scale_a*ca(layer,r)) > eps) then
                            col(pos) = offsets(pa) + idxa
                            val(pos) = der_scale_a*ca(layer,r)
                            pos = pos + 1
                        end if
                    end do
                    do layer = 0, degb
                        idxb = first_b(0) + layer
                        coefficient_b = 0.0_rk
                        do k = 0, r
                            active_b = idxb - first_b(k)
                            if (active_b >= 0 .and. active_b <= degb) then
                                coefficient_b = coefficient_b + chain(r,k)*range_b**k*cb(active_b,k)
                            end if
                        end do
                        if (idxb >= 1 .and. idxb <= ncb .and. abs(weight_scale*coefficient_b) > eps) then
                            col(pos) = offsets(pb) + idxb
                            val(pos) = -weight_scale*coefficient_b
                            pos = pos + 1
                        end if
                    end do
                    rowptr(row+1) = pos
                end do
                deallocate(ca, cb, chain, first_a, first_b)
            end associate
        end do
        col = col(:pos-1)
        val = val(:pos-1)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Concatenate active-span element connectivity from all patches.
    !!
    !! Rows follow patch insertion order and then each patch's span order.
    !! Columns are padded with zero when patch degrees differ. By default node
    !! identifiers use unmerged patch offsets; with `shared=.true.`, the result
    !! uses identifiers from [[nurbs_multipatch_curve:cmp_dof_map]]. Unclamped
    !! traces that are not directly shareable still require the returned sparse
    !! constraint matrix.
    !!
    pure function cmp_elem_curve(this, shared) result(elem)
        class(nurbs_multipatch_curve), intent(in) :: this
            !! Multipatch curve.
        logical, intent(in), optional :: shared
            !! Compact directly shareable \(C^0\) controls; default false.
        integer, allocatable :: elem(:,:)
        integer, allocatable :: offsets(:), map(:), mul(:), nspan_patch(:), row_start(:)
        real(rk), allocatable :: knot(:)
        integer :: a, i, e, row, ne, maxnen, nen, nc, degree, nspan, pref, node, np
        logical :: use_shared

        use_shared = .false.
        if (present(shared)) use_shared = shared
        np = this%get_npatch()
        offsets = this%cmp_dof_offsets()
        ne = 0
        maxnen = 0
        allocate(nspan_patch(np), row_start(np + 1))
        row_start(1) = 0
        do i = 1, np
            nc = this%patches(i)%get_nc()
            degree = this%patches(i)%get_degree()
            knot = this%patches(i)%get_knot()
            nspan_patch(i) = active_span_count(knot, nc, degree)
            row_start(i + 1) = row_start(i) + nspan_patch(i)
            maxnen = max(maxnen, degree + 1)
        end do
        ne = row_start(np + 1)
        allocate(elem(ne, maxnen), source=0)
        if (use_shared) map = this%cmp_dof_map()
        do i = 1, np
            nc = this%patches(i)%get_nc()
            degree = this%patches(i)%get_degree()
            knot = this%patches(i)%get_knot()
            nspan = nspan_patch(i)
            if (nspan < 1) cycle
            mul = active_knot_multiplicity(knot, nc, degree)
            if (size(mul) < nspan + 1) cycle
            nen = degree + 1
            pref = degree + 1
            row = row_start(i)
            do e = 1, nspan
                row = row + 1
                do a = 1, nen
                    node = pref - degree + a - 1
                    if (use_shared) then
                        elem(row,a) = map(offsets(i) + node)
                    else
                        elem(row,a) = offsets(i) + node
                    end if
                end do
                if (e < nspan) pref = pref + mul(e+1)
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the owning patch identifier for every concatenated element.
    !!
    !! Result ordering is identical to [[nurbs_multipatch_curve:cmp_elem]].
    pure function cmp_elem_patch_curve(this) result(patch_id)
        class(nurbs_multipatch_curve), intent(in) :: this
        integer, allocatable :: patch_id(:)
        integer, allocatable :: nloc(:), row_start(:)
        integer :: i, e, ne, np

        np = this%get_npatch()
        allocate(nloc(np), row_start(np + 1))
        row_start(1) = 0
        do i = 1, np
            nloc(i) = active_span_count(this%patches(i)%get_knot(), this%patches(i)%get_nc(), &
                this%patches(i)%get_degree())
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
    !! Result ordering is identical to [[nurbs_multipatch_curve:cmp_elem]].
    pure function cmp_elem_local_curve(this) result(local_id)
        class(nurbs_multipatch_curve), intent(in) :: this
        integer, allocatable :: local_id(:)
        integer, allocatable :: nloc(:), row_start(:)
        integer :: i, e, ne, np

        np = this%get_npatch()
        allocate(nloc(np), row_start(np + 1))
        row_start(1) = 0
        do i = 1, np
            nloc(i) = active_span_count(this%patches(i)%get_knot(), this%patches(i)%get_nc(), &
                this%patches(i)%get_degree())
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
    !> Validate all patches and connection invariants without changing state.
    !!
    !! A valid multipatch contains at least one complete, internally consistent
    !! curve patch. Validation also covers patch/side ranges, continuity order,
    !! endpoint-trace coincidence, transition jets, homogeneous derivative
    !! residuals, and rational projective weight compatibility. `tol` is used as
    !! one absolute threshold for assembled homogeneous residuals; the default
    !! is `sqrt(epsilon(rk))`. Each component of \(A\mathbf H\) is compared
    !! directly with `tol`. Weight proportionality uses
    !! `max(tol,128*epsilon(rk))*max(abs(wA),abs(lambda*wB))`.
    !! `tol>=0` is a caller precondition and is not separately diagnosed. The
    !! test conservatively requires both
    !! patches to have the same whole-patch `is_rational` classification and, for
    !! rational patches, enforces the constant-projective-scale homogeneous model
    !! documented by [[forcad_multipatch]]. It can therefore reject a Euclidean
    !! join admitted by a more general projective gauge.
    pure function is_valid_curve(this, tol) result(ok)
        class(nurbs_multipatch_curve), intent(in) :: this
        real(rk), intent(in), optional :: tol
            !! Optional nonnegative absolute coordinate/residual tolerance; default `sqrt(epsilon(rk))`.
        logical :: ok
        real(rk) :: eps
        integer :: i, j, d, pa, pb, dim_a, dim_b, row, a, c, total, ncomp, normal_sign
        integer, allocatable :: offsets(:), rowptr(:), col(:)
        real(rk), allocatable :: csrval(:), dof(:,:), Xc(:,:), Wc(:)
        real(rk) :: residual, weight_scale
        logical :: rational_a, rational_b, compatible_weights

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
                if (conn%get_side_a() < left_ .or. conn%get_side_a() > right_ .or. &
                    conn%get_side_b() < left_ .or. conn%get_side_b() > right_) then
                    ok = .false.; return
                end if
                if (conn%get_continuity() < -1) then
                    ok = .false.; return
                end if
                normal_sign = merge(-1, 1, mod(conn%get_side_a(), 2) == mod(conn%get_side_b(), 2))
                if (.not. multipatch_valid_reparameterization(conn, normal_sign)) then
                    ok = .false.; return
                end if
                if (conn%get_continuity() > min(&
                    this%patches(pa)%get_degree(), &
                    this%patches(pb)%get_degree())) then
                    ok = .false.; return
                end if
                if (conn%get_continuity() >= 0) then
                    dim_a = size(this%patches(pa)%get_Xc(1))
                    dim_b = size(this%patches(pb)%get_Xc(1))
                    if (dim_a /= dim_b) then
                        ok = .false.; return
                    end if
                    rational_a = this%patches(pa)%is_rational()
                    rational_b = this%patches(pb)%is_rational()
                    if (rational_a .neqv. rational_b) then
                        ok = .false.; return
                    end if
                    call connection_weight_scale_curve(this, conn, eps, weight_scale, compatible_weights)
                    if (.not. compatible_weights) then
                        ok = .false.; return
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
        do i = 1, this%get_npatch()
            ncomp = max(ncomp, size(this%patches(i)%get_Xc(1)))
        end do
        if (this%is_rational()) ncomp = ncomp + 1
        allocate(dof(total, ncomp), source=0.0_rk)
        do i = 1, this%get_npatch()
            Xc = this%patches(i)%get_Xc()
            if (this%is_rational() .and. this%patches(i)%is_rational()) Wc = this%patches(i)%get_Wc()
            do a = 1, size(Xc, 1)
                if (this%is_rational() .and. this%patches(i)%is_rational()) then
                    do d = 1, size(Xc, 2)
                        dof(offsets(i)+a, d) = Xc(a, d)*Wc(a)
                    end do
                    dof(offsets(i)+a, ncomp) = Wc(a)
                else
                    do d = 1, size(Xc, 2)
                        dof(offsets(i)+a, d) = Xc(a, d)
                    end do
                    if (this%is_rational()) dof(offsets(i)+a, ncomp) = 1.0_rk
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
    !> Sample every patch and cache its single-patch geometry data.
    !! `res` is passed to every patch. When `Xt` is present, the same explicit
    !! parameter vector is used for every patch and `res` is ignored by the
    !! corresponding single-patch overload.
    pure subroutine create_curve(this, res, Xt)
        class(nurbs_multipatch_curve), intent(inout) :: this
        integer, intent(in), optional :: res
        real(rk), intent(in), contiguous, optional :: Xt(:)
        integer :: i

        if (.not. this%err%ok) return

        if (this%get_npatch() < 1) then
            call this%err%set(&
                code       = 205,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
                message    = "No patches are available.",&
                location   = "create",&
                suggestion = "Add at least one patch before calling create().")
            return
        end if
        do i = 1, this%get_npatch()
            call this%patches(i)%create(res=res, Xt=Xt)
            if (.not. this%patches(i)%err%ok) then
                call this%err%set(&
                    code       = 205,&
                    severity   = 1,&
                    category   = "forcad_multipatch_curve",&
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
    pure elemental logical function is_rational_curve(this) result(r)
        class(nurbs_multipatch_curve), intent(in) :: this
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
    !> Export every patch control polygon to a separate legacy VTK file.
    !!
    !! Patch `i` is written as `<prefix>_Xc_<i>.vtk`.
    impure subroutine export_Xc_curve(this, prefix, encoding)
        class(nurbs_multipatch_curve), intent(inout) :: this
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
    impure subroutine export_Xg_curve(this, prefix, encoding)
        class(nurbs_multipatch_curve), intent(inout) :: this
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
    !> Export parameter lines mapped into every curve patch.
    impure subroutine export_Xth_in_Xg_curve(this, prefix, res, encoding)
        class(nurbs_multipatch_curve), intent(inout) :: this
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
    !> Apply one x-y-z Euler rotation in degrees to every control polygon.
    !!
    !! Cached physical samples are not changed.
    pure subroutine rotate_Xc_curve(this, alpha, beta, theta)
        class(nurbs_multipatch_curve), intent(inout) :: this
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
    pure subroutine rotate_Xg_curve(this, alpha, beta, theta)
        class(nurbs_multipatch_curve), intent(inout) :: this
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
    !> Translate every patch control polygon by physical vector `vec`.
    !!
    !! Cached physical samples are not changed.
    pure subroutine translate_Xc_curve(this, vec)
        class(nurbs_multipatch_curve), intent(inout) :: this
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
    pure subroutine translate_Xg_curve(this, vec)
        class(nurbs_multipatch_curve), intent(inout) :: this
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
    impure subroutine show_curve(this, vtkfile_Xc, vtkfile_Xg, vtkfile_Xth_in_Xg)
        class(nurbs_multipatch_curve), intent(inout) :: this
        character(len=*), intent(in) :: vtkfile_Xc
            !! Existing control-geometry VTK path or glob pattern.
        character(len=*), intent(in) :: vtkfile_Xg
            !! Existing sampled-geometry VTK path or glob pattern.
        character(len=*), intent(in), optional :: vtkfile_Xth_in_Xg
            !! Optional existing parameter-line VTK path or glob pattern.

        if (.not. this%err%ok) return
#ifndef NOSHOW_PYVISTA
        call show_pyvista_multipatch(&
            vtkfile_Xc        = vtkfile_Xc,&
            vtkfile_Xg        = vtkfile_Xg,&
            vtkfile_Xth_in_Xg = vtkfile_Xth_in_Xg,&
            rank_name         = "multipatch curve")
#endif
    end subroutine
    !===============================================================================
end module forcad_multipatch_curve
