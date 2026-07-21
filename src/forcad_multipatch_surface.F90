!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Topology, continuity constraints, and assembly data for surface patches.
!!
!! Patches remain complete [[nurbs_surface]] objects. Connections identify
!! oriented boundary curves and require affine-equivalent tangential trace
!! spaces. [[nurbs_multipatch_surface:cmp_dof_map]] merges controls paired by
!! recorded \(C^0\) topology (geometric coincidence is established separately
!! by `is_valid`), while
!! [[nurbs_multipatch_surface:cmp_dof_constraint]] supplies sparse equations for
!! normal derivatives through any order
!! \(0\le n\le\min(p_{n,A},p_{n,B})\). Geometric rows
!! use the restricted separable map
!! \(\Phi(\tau,\eta)=(T\tau,\phi(\eta))\); they do not represent a fully general
!! tangentially coupled surface \(G^n\) reparameterization.
!!
!! Global element order is patch insertion order followed by each patch's
!! direction-1-fastest active-span order. The container does not assemble PDE
!! operators; it provides the topology and algebraic data required to do so.
!!
!! Unmerged scalar control variables are concatenated by patch,
!! \(\mathbf q=[(\mathbf q^{(1)})^T,\ldots]^T\). Order-zero interfaces may be
!! represented by the compact map from
!! [[nurbs_multipatch_surface:cmp_dof_map]]. Higher continuity is returned as
!! \(\mathbf A\mathbf q=\mathbf0\) in one-based CSR form and is not silently
!! folded into element connectivity. Raw rows act on polynomial coefficients or
!! homogeneous rational controls. Rational scalar fields require weighted
!! numerator controls \(w_iq_i\), together with a compatible denominator.
module forcad_multipatch_surface

    use forcad_kinds, only: rk
    use forcad_multipatch, only: multipatch_chain_rule_coefficients, multipatch_connection, &
        multipatch_compatible_trace_space, multipatch_growth_capacity, &
        multipatch_valid_reparameterization
    use forcad_nurbs_surface, only: nurbs_surface
    use forcad_utils, only: active_knot_multiplicity, active_span_count, basis_bspline_der_order_active, &
        boundary_index, boundary_layer_index, &
        disjoint_set_map, disjoint_set_union, knot_end, knot_start, show_pyvista_multipatch
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
    !! constraints; order zero merges matching boundary controls.
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
        procedure :: cmp_dof_map => cmp_dof_map_surface !!> Map unmerged controls to compact shared \(C^0\) identifiers.
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
            !! Surface patch whose diagnostic state must be clear.
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
    !> Record two surface edges with parametric \(C^n\) or restricted geometric metadata.
    !! `geometric` defaults to `.false.` for C^n continuity.
    !! For G^n, `reparameterization` stores the arbitrary-order transition jet
    !! between normalized coordinates normal to the connected sides. Tangential
    !! orientation is controlled independently by `reverse`. The geometric map
    !! is restricted to \(\Phi(\tau,\eta)=(T\tau,\phi(\eta))\); no transition
    !! derivative may vary along the edge or couple normal motion into the
    !! tangential coordinate.
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
        geometric, reparameterization)
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
            !! Select the restricted normal-only geometric model rather than parametric continuity.
        real(rk), intent(in), contiguous, optional :: reparameterization(:)
            !! Normalized normal transition derivatives for \(G^n\).
        type(multipatch_connection) :: conn
        type(multipatch_connection), allocatable :: tmp(:)
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
        call conn%set(&
            patch_a            = patch_a,&
            side_a             = side_a_id,&
            patch_b            = patch_b,&
            side_b             = side_b_id,&
            continuity         = continuity_,&
            reverse            = reverse_,&
            geometric          = geometric_,&
            reparameterization = reparameterization)
        normal_sign = merge(-1, 1, mod(side_a_id, 2) == mod(side_b_id, 2))
        if (.not. multipatch_valid_reparameterization(conn, normal_sign)) then
            call this%err%set(&
                code       = 203,&
                severity   = 1,&
                category   = "forcad_multipatch_surface",&
                message    = "Invalid continuity reparameterization.",&
                location   = "connect",&
                suggestion = "For G^n, pass n finite normalized derivatives with a "//&
                    "first-derivative sign consistent with the sides.")
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
    !> Build compact global identifiers for paired \(C^0\) boundary controls.
    !! Higher normal-derivative conditions remain explicit in
    !! [[nurbs_multipatch_surface:cmp_dof_constraint]]. The map trusts recorded
    !! topology and does not prove that paired controls are geometrically
    !! coincident; call `is_valid` first.
    !!
    !! Disjoint-set implementation reference: R. E. Tarjan, *JACM* 22 (1975),
    !! 215--225.
    !! [doi:10.1145/321879.321884](https://doi.org/10.1145/321879.321884).
    pure function cmp_dof_map_surface(this) result(map)
        class(nurbs_multipatch_surface), intent(in) :: this
        integer, allocatable :: map(:)
        integer, allocatable :: parent(:), offsets(:), ia(:), ib(:)
        integer :: i, j, total, nca(2), ncb(2)

        offsets = this%cmp_dof_offsets()
        total = offsets(size(offsets))
        allocate(parent(total))
        do i = 1, total
            parent(i) = i
        end do

        do i = 1, this%get_nconnection()
            associate(conn => this%connections(i))
                if (conn%get_continuity() < 0) cycle
                nca = this%patches(conn%get_patch_a())%get_nc()
                ncb = this%patches(conn%get_patch_b())%get_nc()
                call boundary_index(nca, conn%get_side_a(), .false., ia)
                call boundary_index(ncb, conn%get_side_b(), conn%is_reversed(), ib)
                do j = 1, size(ia)
                    call disjoint_set_union(parent, offsets(conn%get_patch_a())+ia(j), offsets(conn%get_patch_b())+ib(j))
                end do
            end associate
        end do
        call disjoint_set_map(parent, map)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Check projective compatibility of rational edge weights.
    !!
    !! Compatible mapped traces satisfy \(w_j^{(A)}=\lambda w_j^{(B)}\)
    !! within `tol`; `weight_scale` returns \(\lambda\).
    !! Polynomial traces use unit scale.
    pure subroutine connection_weight_scale_surface(this, conn, tol, weight_scale, ok)
        class(nurbs_multipatch_surface), intent(in) :: this
        type(multipatch_connection), intent(in) :: conn
        real(rk), intent(in) :: tol
        real(rk), intent(out) :: weight_scale
        logical, intent(out) :: ok
        real(rk) :: ref, wa, wb, wtol
        integer, allocatable :: ia(:), ib(:)
        integer :: j, pa, pb, nca(2), ncb(2)
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
        call boundary_index(nca, conn%get_side_a(), .false., ia)
        call boundary_index(ncb, conn%get_side_b(), conn%is_reversed(), ib)
        if (size(ia) /= size(ib) .or. size(ia) < 1) then
            ok = .false.
            return
        end if

        wb = this%patches(pb)%get_Wc(ib(1))
        if (abs(wb) <= tiny(1.0_rk)) then
            ok = .false.
            return
        end if

        weight_scale = this%patches(pa)%get_Wc(ia(1))/wb
        do j = 1, size(ia)
            wa = this%patches(pa)%get_Wc(ia(j))
            wb = this%patches(pb)%get_Wc(ib(j))
            ref = max(abs(wa), abs(weight_scale*wb))
            wtol = max(tol*ref, 128.0_rk*epsilon(1.0_rk)*ref)
            if (abs(wa - weight_scale*wb) > wtol) then
                weight_scale = 1.0_rk
                ok = .false.
                return
            end if
        end do
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
    !> Assemble sparse parametric and restricted geometric constraints for connected edges.
    !!
    !! The one-based CSR arrays define \(Aq=0\) over unmerged,
    !! patch-concatenated control variables. One row is emitted for every
    !! tangential boundary control and derivative order. `geometric`, when
    !! present, selects only parametric (false) or geometric (true) connections;
    !! when absent, both kinds are included. Geometric rows use only the scalar
    !! normal transition jet documented by [[forcad_multipatch]]. Rational
    !! projective scaling is incorporated for compatible boundary weights.
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
            !! Optional filter: false for parametric, true for restricted geometric connections.
        real(rk), allocatable :: ca(:,:), cb(:,:), chain(:,:), knot_a(:), knot_b(:)
        real(rk) :: xa, xb, weight_scale, range_a, range_b, der_scale_a, coefficient_b, eps
        integer, allocatable :: first_a(:), first_b(:), offsets(:), ia(:), ib(:), layer_a(:,:), layer_b(:,:)
        integer :: i, j, k, r, layer, row, pos, nrow, nnz, pa, pb, side_a, side_b
        integer :: dir_a, dir_b, normal_a, normal_b, normal_sign, ntan
        integer :: active_a, active_b
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
    !! `shared=.true.` applies the compact \(C^0\) map; the default retains
    !! unmerged patch-offset numbering.
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
    !! Checks include identifiers, sides, continuity order, transition jets,
    !! tangential trace spaces, paired boundary-control geometry, homogeneous
    !! normal-derivative residuals in the supported separable model, embedding
    !! dimensions, and rational projective weight compatibility. `tol` is one
    !! absolute threshold for coordinates and assembled homogeneous residuals;
    !! the default is `sqrt(epsilon(rk))`. Paired coordinates and each component
    !! of \(A\mathbf H\) are compared directly with `tol`. Weight proportionality
    !! uses `max(tol,128*epsilon(rk))*max(abs(wA),abs(lambda*wB))` per pair.
    !! `tol>=0` is an unchecked caller precondition. The test conservatively
    !! requires equal
    !! whole-patch `is_rational` classifications and the constant-projective-
    !! scale homogeneous model documented by [[forcad_multipatch]]. It is
    !! sufficient but not exhaustive for Euclidean rational continuity.
    pure function is_valid_surface(this, tol) result(ok)
        class(nurbs_multipatch_surface), intent(in) :: this
        real(rk), intent(in), optional :: tol
            !! Optional nonnegative absolute coordinate/residual tolerance; default `sqrt(epsilon(rk))`.
        logical :: ok
        real(rk) :: eps
        integer :: i, j, d, nca(2), ncb(2), dega, degb, dir_a, dir_b
        integer :: pa, pb, dim_a, dim_b, normal_sign
        integer :: row, a, c, total, ncomp
        integer, allocatable :: ia(:), ib(:), offsets(:), rowptr(:), col(:)
        real(rk), allocatable :: csrval(:), dof(:,:), Xc(:,:), Wc(:)
        real(rk) :: residual, weight_scale
        logical :: rational_a, rational_b, compatible_weights

        if (.not. this%err%ok) then
            ok = .false.
            return
        end if

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
                if (.not. multipatch_valid_reparameterization(conn, normal_sign)) then
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
                    nca = this%patches(pa)%get_nc()
                    ncb = this%patches(pb)%get_nc()
                    call boundary_index(nca, conn%get_side_a(), .false., ia)
                    call boundary_index(ncb, conn%get_side_b(), conn%is_reversed(), ib)
                    if (size(ia) /= size(ib)) then
                        ok = .false.; return
                    end if
                    dim_a = size(this%patches(pa)%get_Xc(ia(1)))
                    dim_b = size(this%patches(pb)%get_Xc(ib(1)))
                    if (dim_a /= dim_b) then
                        ok = .false.; return
                    end if
                    rational_a = this%patches(pa)%is_rational()
                    rational_b = this%patches(pb)%is_rational()
                    if (rational_a .neqv. rational_b) then
                        ok = .false.; return
                    end if
                    call connection_weight_scale_surface(this, conn, eps, weight_scale, compatible_weights)
                    if (.not. compatible_weights) then
                        ok = .false.; return
                    end if
                    do j = 1, size(ia)
                        do d = 1, dim_a
                            if (abs(this%patches(pa)%get_Xc(ia(j), d) - this%patches(pb)%get_Xc(ib(j), d)) > eps) then
                                ok = .false.; return
                            end if
                        end do
                    end do
                end if
            end associate
        end do

        offsets = this%cmp_dof_offsets()
        total = offsets(size(offsets))
        if (total < 1) return
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
