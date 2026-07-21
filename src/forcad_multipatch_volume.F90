!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Topology, continuity constraints, and assembly data for volume patches.
!!
!! Patches remain complete [[nurbs_volume]] mappings. Connections identify
!! oriented boundary surfaces and validate both tangential trace spaces after
!! the requested reversal, axis exchange, and flips.
!! [[nurbs_multipatch_volume:cmp_dof_map]] merges controls paired by recorded
!! \(C^0\) face topology only where both normal endpoint bases are
!! interpolatory; [[nurbs_multipatch_volume:cmp_dof_constraint]] returns sparse
!! normal-derivative equations for any order
!! \(0\le n\le\min(p_{n,A},p_{n,B})\). A compact geometric connection uses the
!! separable map
!! \(\Phi(\boldsymbol\tau,\eta)=(T\boldsymbol\tau,\phi(\eta))\). A general
!! connection supplies tensor-product Bernstein fields for every normal
!! derivative of all three mapped coordinates and, for rational patches, every
!! derivative of a nonvanishing projective factor. The assembled \(G^n\)
!! equations therefore include tangentially varying normal jets, mixed
!! derivatives, and the complete projective product rule.
!!
!! Global elements are ordered by patch insertion and then by each volume's
!! direction-1-fastest active-span ordering. The type supplies topology and
!! assembly metadata but does not assemble or solve a physical model.
!!
!! Unmerged scalar control variables are concatenated by patch,
!! \(\mathbf q=[(\mathbf q^{(1)})^T,\ldots]^T\). An order-zero interface may be
!! represented by the compact map from
!! [[nurbs_multipatch_volume:cmp_dof_map]] only when its normal endpoints are
!! interpolatory. General unclamped traces and higher continuity are returned as
!! \(\mathbf A\mathbf q=\mathbf0\) in one-based CSR form. Face reversal,
!! tangential-axis exchange, and flips are applied before equations are formed.
!! Raw rows act on polynomial coefficients or homogeneous rational controls.
!! Rational scalar fields require weighted numerator controls \(w_iq_i\),
!! together with a compatible denominator.
module forcad_multipatch_volume

    use forcad_kinds, only: rk
    use forcad_multipatch, only: multipatch_chain_rule_coefficients, multipatch_composition_coefficients, &
        multipatch_connection, multipatch_compatible_trace_space, multipatch_growth_capacity, &
        multipatch_projective_weight_scale, multipatch_valid_reparameterization
    use forcad_nurbs_volume, only: nurbs_volume
    use forcad_utils, only: active_knot_multiplicity, active_span_count, basis_bernstein, &
        basis_bspline_der_order_active, boundary_index, boundary_layer_index, disjoint_set_map, disjoint_set_union, &
        interpolatory_boundary_layer, knot_end, knot_start, show_pyvista_multipatch, structural_nonzero
    use fordebug, only: debug

    implicit none

    private
    public :: nurbs_multipatch_volume

    integer, parameter :: u_min_ = 1
    integer, parameter :: u_max_ = 2
    integer, parameter :: v_min_ = 3
    integer, parameter :: v_max_ = 4
    integer, parameter :: w_min_ = 5
    integer, parameter :: w_max_ = 6
    logical, parameter :: no_flip_(2) = [.false., .false.]

    !===============================================================================
    !> Ordered collection of connected trivariate NURBS patches.
    !!
    !! Valid side labels are `u_min/u_max`, `v_min/v_max`, and `w_min/w_max`,
    !! including the short aliases accepted by `connect`. The patch-B face map
    !! can reverse orientation, swap its two tangential axes, and flip either
    !! mapped axis. `continuity=-1` stores topology without constraints.
    !!
    !! Diagnostic codes are: 201 invalid patch identifier; 202 invalid side;
    !! 203 invalid continuity metadata; 204 incompatible interface geometry or
    !! trace space; and 205 invalid patch state.
    type nurbs_multipatch_volume
        type(nurbs_volume), allocatable, private :: patches(:) !! Owned patch copies with spare append capacity.
        type(multipatch_connection), allocatable, private :: connections(:) !! Validated interface metadata.
        integer, private :: npatch = 0 !! Number of active patches.
        integer, private :: nconnection = 0 !! Number of active connections.

        type(debug) :: err !! Recoverable diagnostic state for the most recent failed operation.
    contains
        procedure :: add_patch => add_patch_volume !!> Append a validated patch copy and optionally return its identifier.
        procedure :: connect => connect_volume !!> Validate metadata and trace-space shape, then record a face interface.
        procedure :: create => create_volume !!> Sample every patch using common directional parameters.
        procedure :: finalize => clear_volume !!> Release patches and connections and reset counts.
        procedure :: get_npatch => get_npatch_volume !!> Return the active patch count.
        procedure :: get_nconnection => get_nconnection_volume !!> Return the active connection count.
        procedure :: get_patch => get_patch_volume !!> Return a copy of one patch.
        procedure :: get_connection => get_connection_volume !!> Return a copy of one connection.
        procedure :: is_valid => is_valid_volume !!> Validate patch states, face maps, traces, and continuity.
        procedure :: is_rational => is_rational_volume !!> Report whether any patch has nonuniform rational weights.
        procedure :: cmp_dof_offsets => cmp_dof_offsets_volume !!> Return prefix offsets for unmerged patch controls.
        procedure :: cmp_dof_map => cmp_dof_map_volume !!> Compact directly shareable \(C^0\) controls.
        procedure :: cmp_dof_constraint => cmp_dof_constraint_volume !!> Assemble derivative-continuity equations in CSR form.
        procedure :: cmp_elem => cmp_elem_volume !!> Concatenate patch element connectivity with optional shared numbering.
        procedure :: get_elem => cmp_elem_volume !!> Alias for [[nurbs_multipatch_volume:cmp_elem]].
        procedure :: cmp_elem_patch => cmp_elem_patch_volume !!> Return the owning patch identifier for every global element.
        procedure :: cmp_elem_local => cmp_elem_local_volume !!> Return the patch-local identifier for every global element.
        procedure :: export_Xc => export_Xc_volume !!> Export each patch control lattice to a separate VTK file.
        procedure :: export_Xg => export_Xg_volume !!> Export each patch's cached geometry to a separate VTK file.
        procedure :: export_Xth_in_Xg => export_Xth_in_Xg_volume !!> Export embedded parameter grids for every patch.
        procedure :: rotate_Xc => rotate_Xc_volume !!> Rotate all patch control points.
        procedure :: rotate_Xg => rotate_Xg_volume !!> Rotate all cached patch geometry points.
        procedure :: translate_Xc => translate_Xc_volume !!> Translate all patch control points.
        procedure :: translate_Xg => translate_Xg_volume !!> Translate all cached patch geometry points.
        procedure :: show => show_volume !!> Display selected patch representations together with PyVista.
    end type
    !===============================================================================

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Append a copy of a valid volume patch.
    pure subroutine add_patch_volume(this, patch, id)
        class(nurbs_multipatch_volume), intent(inout) :: this
            !! Multipatch container.
        type(nurbs_volume), intent(in) :: patch
            !! Complete, internally consistent volume geometry with a clear diagnostic state.
        integer, intent(out), optional :: id
            !! Assigned one-based patch identifier; zero is returned on failure.
        type(nurbs_volume), allocatable :: tmp(:)
        integer :: new_capacity

        if (present(id)) id = 0
        if (.not. this%err%ok) return

        if (.not. patch%err%ok) then
            call this%err%set(&
                code       = 205,&
                severity   = 1,&
                category   = "forcad_multipatch_volume",&
                message    = "Patch has an active error state.",&
                location   = "add_patch",&
                suggestion = "Reset or fix the patch error state before adding it to a multipatch volume.")
            return
        end if
        if (.not. patch%is_valid()) then
            call this%err%set(&
                code       = 205,&
                severity   = 1,&
                category   = "forcad_multipatch_volume",&
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
    !> Record two volume faces with parametric \(C^n\) or geometric metadata.
    !! `geometric` defaults to `.false.` for C^n continuity.
    !! For the compact G^n representation, `reparameterization` stores the
    !! scalar normal transition derivatives. The general representation uses
    !! `transition_jet(i,j,c,r)`, the tensor-product Bernstein coefficient
    !! `(i,j)` of \(\partial_{\eta_A}^r\Phi_c(\tau_1,\tau_2,0)\).
    !! Components `c=1:2` follow the two mapped patch-B tangential directions
    !! after orientation and `c=3` is the patch-B normal coordinate.
    !! `projective_jet(i,j,r+1)` stores the corresponding homogeneous
    !! projective-factor derivatives. Do not pass both compact and full
    !! transition jets.
    !!
    !! Both mapped tangential trace spaces must have equal degrees and
    !! affine-equivalent normalized knots. Boundary control counts must match.
    !! Geometry, derivative residuals, and projective rational weights are
    !! checked by `is_valid`. Call it before assembly or analysis; successful
    !! `connect` establishes only valid metadata and trace-space shape.
    !!
    pure subroutine connect_volume(this, patch_a, side_a, patch_b, side_b, continuity, reverse, swap, flip, &
        geometric, reparameterization, transition_jet, projective_jet)
        class(nurbs_multipatch_volume), intent(inout) :: this
            !! Multipatch container.
        integer, intent(in) :: patch_a
            !! One-based first patch identifier.
        integer, intent(in) :: patch_b
            !! One-based second patch identifier.
        character(len=*), intent(in) :: side_a
            !! First face label (`u_min` through `w_max`).
        character(len=*), intent(in) :: side_b
            !! Second face label.
        integer, intent(in), optional :: continuity
            !! Requested order, default zero; `-1` disables constraints.
        logical, intent(in), optional :: reverse
            !! Reverse both mapped tangential orientations.
        logical, intent(in), optional :: swap
            !! Exchange patch-B tangential axes before matching.
        logical, intent(in), optional :: geometric
            !! Select geometric rather than parametric continuity.
        logical, intent(in), optional, contiguous :: flip(:)
            !! Reverse either mapped patch-B tangential axis.
        real(rk), intent(in), contiguous, optional :: reparameterization(:)
            !! Compact normalized normal transition derivatives for \(G^n\).
        real(rk), intent(in), contiguous, optional :: transition_jet(:,:,:,:)
            !! Full Bernstein transition coefficients `[q1+1,q2+1,3,n]`.
        real(rk), intent(in), contiguous, optional :: projective_jet(:,:,:)
            !! Bernstein projective-factor coefficients `[q1+1,q2+1,n+1]`.
        type(multipatch_connection) :: conn
        type(multipatch_connection), allocatable :: tmp(:)
        integer, allocatable :: ia(:), ib(:)
        logical :: flip_(2), reverse_, swap_, geometric_
        integer :: n, side_a_id, side_b_id, continuity_, new_capacity, dir_a, dir_b, normal_sign

        if (.not. this%err%ok) return

        if (patch_a < 1 .or. patch_a > this%get_npatch()) then
            call this%err%set(&
                code       = 201,&
                severity   = 1,&
                category   = "forcad_multipatch_volume",&
                message    = "Invalid first patch id.",&
                location   = "connect",&
                suggestion = "Use a patch id from add_patch's optional id argument.")
            return
        end if
        if (patch_b < 1 .or. patch_b > this%get_npatch()) then
            call this%err%set(&
                code       = 201,&
                severity   = 1,&
                category   = "forcad_multipatch_volume",&
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
        case("w_min", "wmin", "w-", "W_MIN", "WMIN", "W-")
            side_a_id = w_min_
        case("w_max", "wmax", "w+", "W_MAX", "WMAX", "W+")
            side_a_id = w_max_
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
        case("w_min", "wmin", "w-", "W_MIN", "WMIN", "W-")
            side_b_id = w_min_
        case("w_max", "wmax", "w+", "W_MAX", "WMAX", "W+")
            side_b_id = w_max_
        case default
            continue
        end select
        if (side_a_id == 0 .or. side_b_id == 0) then
            call this%err%set(&
                code       = 202,&
                severity   = 1,&
                category   = "forcad_multipatch_volume",&
                message    = "Invalid volume side label.",&
                location   = "connect",&
                suggestion = "Use u/v/w min/max side labels, such as u_min or w+.")
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
                category   = "forcad_multipatch_volume",&
                message    = "Invalid continuity value.",&
                location   = "connect",&
                suggestion = "Use -1 for discontinuous or Cn with n >= 0.")
            return
        end if
        dir_a = (side_a_id + 1)/2
        dir_b = (side_b_id + 1)/2
        if (continuity_ > min(this%patches(patch_a)%get_degree(dir_a), this%patches(patch_b)%get_degree(dir_b))) then
            call this%err%set(&
                code       = 203,&
                severity   = 1,&
                category   = "forcad_multipatch_volume",&
                message    = "Continuity is higher than the connected patch degree.",&
                location   = "connect",&
                suggestion = "Use a continuity not larger than the minimum degree normal to the connected faces.")
            return
        end if
        reverse_ = .false.
        swap_ = .false.
        flip_ = .false.
        if (present(reverse)) reverse_ = reverse
        if (present(swap)) swap_ = swap
        if (present(flip)) then
            n = min(2, size(flip))
            flip_(1:n) = flip(1:n)
        end if
        call conn%set(&
            patch_a            = patch_a,&
            side_a             = side_a_id,&
            patch_b            = patch_b,&
            side_b             = side_b_id,&
            continuity         = continuity_,&
            reverse            = reverse_,&
            swap               = swap_,&
            flip               = flip_,&
            geometric          = geometric_,&
            reparameterization = reparameterization,&
            transition_jet     = transition_jet,&
            projective_jet     = projective_jet)
        normal_sign = merge(-1, 1, mod(side_a_id, 2) == mod(side_b_id, 2))
        if (.not. multipatch_valid_reparameterization(conn, normal_sign, 3)) then
            call this%err%set(&
                code       = 203,&
                severity   = 1,&
                category   = "forcad_multipatch_volume",&
                message    = "Invalid continuity reparameterization.",&
                location   = "connect",&
                suggestion = "Pass either n compact derivatives or a finite `[q1+1,q2+1,3,n]` "//&
                    "full jet whose first normal derivative has the side-consistent sign.")
            return
        end if
        if (continuity_ >= 0) then
            call boundary_index(this%patches(patch_a)%get_nc(), side_a_id, .false., .false., no_flip_, ia)
            call boundary_index(this%patches(patch_b)%get_nc(), side_b_id, reverse_, swap_, flip_, ib)
            if (size(ia) /= size(ib)) then
                call this%err%set(&
                    code       = 204,&
                    severity   = 1,&
                    category   = "forcad_multipatch_volume",&
                    message    = "Connected volume faces have incompatible control point counts.",&
                    location   = "connect",&
                    suggestion = "Use matching boundary discretizations or refine one face before connecting patches.")
                return
            end if
            if (.not. connection_trace_space_volume(this, conn, sqrt(epsilon(1.0_rk)))) then
                call this%err%set(&
                    code       = 204,&
                    severity   = 1,&
                    category   = "forcad_multipatch_volume",&
                    message    = "Connected volume faces have incompatible tangential trace spaces.",&
                    location   = "connect",&
                    suggestion = "Use matching tangential degrees and affine-equivalent knot vectors before connecting patches.")
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
    !> Release all owned volume patches and connections.
    !!
    !! Active counts are reset to zero and the multipatch diagnostic is
    !! preserved.
    pure elemental subroutine clear_volume(this)
        class(nurbs_multipatch_volume), intent(inout) :: this

        if (allocated(this%patches)) deallocate(this%patches)
        if (allocated(this%connections)) deallocate(this%connections)
        this%npatch = 0
        this%nconnection = 0
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the number of active volume patches.
    pure elemental integer function get_npatch_volume(this) result(n)
        class(nurbs_multipatch_volume), intent(in) :: this

        n = this%npatch
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return exact allocation bounds for one general geometric face constraint.
    pure subroutine general_constraint_shape_volume(this, conn, nrow, nnz)
        class(nurbs_multipatch_volume), intent(in) :: this
        type(multipatch_connection), intent(in) :: conn
        integer, intent(out) :: nrow, nnz
        real(rk), allocatable :: knot1(:), knot2(:), transition_jet(:,:,:,:), projective_jet(:,:,:)
        integer :: r, nspan1, nspan2, support, pa, pb, normal_a, normal_b
        integer :: tangent_a(2), tangent_b(2), mapped_b(2), degree_row(2)
        integer :: transition_degree(2), projective_degree(2)
        integer :: nc_a(3), nc_b(3), degree_a(3), degree_b(3)

        pa = conn%get_patch_a()
        pb = conn%get_patch_b()
        normal_a = (conn%get_side_a() + 1)/2
        normal_b = (conn%get_side_b() + 1)/2
        call volume_tangent_dirs(conn%get_side_a(), tangent_a)
        call volume_tangent_dirs(conn%get_side_b(), tangent_b)
        if (conn%is_swapped()) then
            mapped_b(1) = tangent_b(2)
            mapped_b(2) = tangent_b(1)
        else
            mapped_b = tangent_b
        end if
        nc_a = this%patches(pa)%get_nc()
        nc_b = this%patches(pb)%get_nc()
        degree_a = this%patches(pa)%get_degree()
        degree_b = this%patches(pb)%get_degree()
        knot1 = this%patches(pa)%get_knot(tangent_a(1))
        knot2 = this%patches(pa)%get_knot(tangent_a(2))
        nspan1 = active_span_count(knot1, nc_a(tangent_a(1)), degree_a(tangent_a(1)))
        nspan2 = active_span_count(knot2, nc_a(tangent_a(2)), degree_a(tangent_a(2)))
        transition_degree = 0
        projective_degree = 0
        if (conn%has_transition_jet()) then
            transition_jet = conn%get_transition_jet()
            transition_degree(1) = size(transition_jet,1) - 1
            transition_degree(2) = size(transition_jet,2) - 1
        end if
        if (conn%has_projective_jet()) then
            projective_jet = conn%get_projective_jet()
            projective_degree(1) = size(projective_jet,1) - 1
            projective_degree(2) = size(projective_jet,2) - 1
        end if

        nrow = 0
        do r = 0, conn%get_continuity()
            degree_row(1) = max(degree_a(tangent_a(1)), &
                degree_b(mapped_b(1)) + r*transition_degree(1) + projective_degree(1))
            degree_row(2) = max(degree_a(tangent_a(2)), &
                degree_b(mapped_b(2)) + r*transition_degree(2) + projective_degree(2))
            nrow = nrow + nspan1*nspan2*(degree_row(1) + 1)*(degree_row(2) + 1)
        end do
        support = (degree_a(1) + 1)*(degree_a(2) + 1)*(degree_a(3) + 1) + &
            (degree_b(1) + 1)*(degree_b(2) + 1)*(degree_b(3) + 1)
        nnz = nrow*support
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Append one full multivariate volume \(G^n\) constraint to CSR storage.
    !!
    !! Tensor-product, degree-complete interior collocation makes every emitted
    !! row set equivalent to the piecewise-polynomial homogeneous jet identity
    !! on its face knot cell. The composition includes all mixed derivatives in
    !! the two mapped tangential coordinates and the mapped normal coordinate.
    pure subroutine append_general_constraint_volume(this, conn, offsets, rowptr, col, val, row, pos)
        class(nurbs_multipatch_volume), intent(in) :: this
        type(multipatch_connection), intent(in) :: conn
        integer, intent(in), contiguous :: offsets(:)
        integer, intent(inout), contiguous :: rowptr(:), col(:)
        real(rk), intent(inout), contiguous :: val(:)
        integer, intent(inout) :: row, pos
        real(rk), allocatable :: knot_na(:), knot_nb(:), knot_ta1(:), knot_ta2(:), knot_tb1(:), knot_tb2(:)
        real(rk), allocatable :: basis_na(:,:), basis_nb(:,:), basis_ta1(:), basis_ta2(:)
        real(rk), allocatable :: basis_tb1(:,:), basis_tb2(:,:)
        real(rk), allocatable :: bernstein_t1(:), bernstein_t2(:)
        real(rk), allocatable :: transition_jet(:,:,:,:), projective_jet(:,:,:)
        real(rk), allocatable :: transition_value(:,:), projective_value(:)
        real(rk), allocatable :: composition(:,:,:,:), reparameterization(:)
        real(rk) :: start_na, start_nb, start_ta1, start_ta2, start_tb1, start_tb2
        real(rk) :: range_na, range_nb, range_ta1, range_ta2, range_tb1, range_tb2
        real(rk) :: coordinate_na, coordinate_nb, coordinate_ta1, coordinate_ta2
        real(rk) :: coordinate_tb1, coordinate_tb2, left1, right1, left2, right2
        real(rk) :: tau1, tau2, coefficient_b, mixed, choose, weight_scale
        integer, allocatable :: first_na(:), first_nb(:), first_tb1(:), first_tb2(:)
        integer :: pa, pb, side_a, side_b, normal_a, normal_b, tangent_a(2), tangent_b(2), mapped_b(2)
        integer :: nc_a(3), nc_b(3), degree_a(3), degree_b(3), index_a(3), index_b(3)
        integer :: transition_degree(2), projective_degree(2), degree_row(2)
        integer :: n, r, j, alpha1, alpha2, alpha_n, local1, local2, local_n
        integer :: span1, span2, point1, point2, coefficient1, coefficient2, component
        integer :: first_ta1, first_ta2
        logical :: reverse_dir(2), compatible_weights

        pa = conn%get_patch_a()
        pb = conn%get_patch_b()
        side_a = conn%get_side_a()
        side_b = conn%get_side_b()
        normal_a = (side_a + 1)/2
        normal_b = (side_b + 1)/2
        call volume_tangent_dirs(side_a, tangent_a)
        call volume_tangent_dirs(side_b, tangent_b)
        if (conn%is_swapped()) then
            mapped_b(1) = tangent_b(2)
            mapped_b(2) = tangent_b(1)
            reverse_dir(1) = conn%is_flipped(2)
            reverse_dir(2) = conn%is_flipped(1)
        else
            mapped_b = tangent_b
            reverse_dir(1) = conn%is_flipped(1)
            reverse_dir(2) = conn%is_flipped(2)
        end if
        if (conn%is_reversed()) reverse_dir = .not. reverse_dir

        nc_a = this%patches(pa)%get_nc()
        nc_b = this%patches(pb)%get_nc()
        degree_a = this%patches(pa)%get_degree()
        degree_b = this%patches(pb)%get_degree()
        knot_na = this%patches(pa)%get_knot(normal_a)
        knot_nb = this%patches(pb)%get_knot(normal_b)
        knot_ta1 = this%patches(pa)%get_knot(tangent_a(1))
        knot_ta2 = this%patches(pa)%get_knot(tangent_a(2))
        knot_tb1 = this%patches(pb)%get_knot(mapped_b(1))
        knot_tb2 = this%patches(pb)%get_knot(mapped_b(2))
        start_na = knot_start(knot_na, nc_a(normal_a), degree_a(normal_a))
        start_nb = knot_start(knot_nb, nc_b(normal_b), degree_b(normal_b))
        start_ta1 = knot_start(knot_ta1, nc_a(tangent_a(1)), degree_a(tangent_a(1)))
        start_ta2 = knot_start(knot_ta2, nc_a(tangent_a(2)), degree_a(tangent_a(2)))
        start_tb1 = knot_start(knot_tb1, nc_b(mapped_b(1)), degree_b(mapped_b(1)))
        start_tb2 = knot_start(knot_tb2, nc_b(mapped_b(2)), degree_b(mapped_b(2)))
        range_na = knot_end(knot_na, nc_a(normal_a), degree_a(normal_a)) - start_na
        range_nb = knot_end(knot_nb, nc_b(normal_b), degree_b(normal_b)) - start_nb
        range_ta1 = knot_end(knot_ta1, nc_a(tangent_a(1)), degree_a(tangent_a(1))) - start_ta1
        range_ta2 = knot_end(knot_ta2, nc_a(tangent_a(2)), degree_a(tangent_a(2))) - start_ta2
        range_tb1 = knot_end(knot_tb1, nc_b(mapped_b(1)), degree_b(mapped_b(1))) - start_tb1
        range_tb2 = knot_end(knot_tb2, nc_b(mapped_b(2)), degree_b(mapped_b(2))) - start_tb2
        coordinate_na = merge(start_na, start_na + range_na, mod(side_a,2) == 1)
        coordinate_nb = merge(start_nb, start_nb + range_nb, mod(side_b,2) == 1)
        n = conn%get_continuity()

        transition_degree = 0
        projective_degree = 0
        if (conn%has_transition_jet()) then
            transition_jet = conn%get_transition_jet()
            transition_degree(1) = size(transition_jet,1) - 1
            transition_degree(2) = size(transition_jet,2) - 1
        else
            reparameterization = conn%get_reparameterization()
        end if
        if (conn%has_projective_jet()) then
            projective_jet = conn%get_projective_jet()
            projective_degree(1) = size(projective_jet,1) - 1
            projective_degree(2) = size(projective_jet,2) - 1
            weight_scale = 1.0_rk
        else
            call connection_weight_scale_volume(&
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
            basis_ta1(0:degree_a(tangent_a(1))),&
            basis_ta2(0:degree_a(tangent_a(2))),&
            basis_tb1(0:degree_b(mapped_b(1)),0:n),&
            basis_tb2(0:degree_b(mapped_b(2)),0:n),&
            first_na(0:n),&
            first_nb(0:n),&
            first_tb1(0:n),&
            first_tb2(0:n),&
            transition_value(3,n),&
            projective_value(0:n))
        do j = 0, n
            call basis_bspline_der_order_active(&
                coordinate_na,knot_na,nc_a(normal_a),degree_a(normal_a),j,first_na(j),basis_na(:,j))
            call basis_bspline_der_order_active(&
                coordinate_nb,knot_nb,nc_b(normal_b),degree_b(normal_b),j,first_nb(j),basis_nb(:,j))
        end do

        do r = 0, n
            degree_row(1) = max(degree_a(tangent_a(1)), &
                degree_b(mapped_b(1)) + r*transition_degree(1) + projective_degree(1))
            degree_row(2) = max(degree_a(tangent_a(2)), &
                degree_b(mapped_b(2)) + r*transition_degree(2) + projective_degree(2))
            do span2 = degree_a(tangent_a(2)) + 1, nc_a(tangent_a(2))
                if (knot_ta2(span2+1) <= knot_ta2(span2)) cycle
                left2 = (knot_ta2(span2) - start_ta2)/range_ta2
                right2 = (knot_ta2(span2+1) - start_ta2)/range_ta2
                do span1 = degree_a(tangent_a(1)) + 1, nc_a(tangent_a(1))
                    if (knot_ta1(span1+1) <= knot_ta1(span1)) cycle
                    left1 = (knot_ta1(span1) - start_ta1)/range_ta1
                    right1 = (knot_ta1(span1+1) - start_ta1)/range_ta1
                    do point2 = 1, degree_row(2) + 1
                        tau2 = left2 + (right2 - left2)*&
                            (real(point2, rk) - 0.5_rk)/real(degree_row(2) + 1, rk)
                        coordinate_ta2 = start_ta2 + tau2*range_ta2
                        coordinate_tb2 = start_tb2 + merge(1.0_rk - tau2, tau2, reverse_dir(2))*range_tb2
                        do point1 = 1, degree_row(1) + 1
                            tau1 = left1 + (right1 - left1)*&
                                (real(point1, rk) - 0.5_rk)/real(degree_row(1) + 1, rk)
                            coordinate_ta1 = start_ta1 + tau1*range_ta1
                            coordinate_tb1 = start_tb1 + merge(1.0_rk - tau1, tau1, reverse_dir(1))*range_tb1
                            call basis_bspline_der_order_active(&
                                coordinate_ta1,knot_ta1,nc_a(tangent_a(1)),degree_a(tangent_a(1)),&
                                0,first_ta1,basis_ta1)
                            call basis_bspline_der_order_active(&
                                coordinate_ta2,knot_ta2,nc_a(tangent_a(2)),degree_a(tangent_a(2)),&
                                0,first_ta2,basis_ta2)
                            do j = 0, n
                                call basis_bspline_der_order_active(&
                                    coordinate_tb1,knot_tb1,nc_b(mapped_b(1)),degree_b(mapped_b(1)),&
                                    j,first_tb1(j),basis_tb1(:,j))
                                call basis_bspline_der_order_active(&
                                    coordinate_tb2,knot_tb2,nc_b(mapped_b(2)),degree_b(mapped_b(2)),&
                                    j,first_tb2(j),basis_tb2(:,j))
                            end do

                            transition_value = 0.0_rk
                            if (conn%has_transition_jet()) then
                                bernstein_t1 = basis_bernstein(tau1, size(transition_jet,1))
                                bernstein_t2 = basis_bernstein(tau2, size(transition_jet,2))
                                do j = 1, n
                                    do component = 1, 3
                                        do coefficient2 = 1, size(transition_jet,2)
                                            do coefficient1 = 1, size(transition_jet,1)
                                                transition_value(component,j) = transition_value(component,j) + &
                                                    bernstein_t1(coefficient1)*bernstein_t2(coefficient2)*&
                                                    transition_jet(coefficient1,coefficient2,component,j)
                                            end do
                                        end do
                                    end do
                                end do
                            else if (n > 0) then
                                transition_value(3,:) = reparameterization
                            end if
                            call multipatch_composition_coefficients(transition_value, composition)

                            projective_value = 0.0_rk
                            if (conn%has_projective_jet()) then
                                bernstein_t1 = basis_bernstein(tau1, size(projective_jet,1))
                                bernstein_t2 = basis_bernstein(tau2, size(projective_jet,2))
                                do j = 0, n
                                    do coefficient2 = 1, size(projective_jet,2)
                                        do coefficient1 = 1, size(projective_jet,1)
                                            projective_value(j) = projective_value(j) + &
                                                bernstein_t1(coefficient1)*bernstein_t2(coefficient2)*&
                                                projective_jet(coefficient1,coefficient2,j+1)
                                        end do
                                    end do
                                end do
                            else
                                projective_value(0) = weight_scale
                            end if

                            row = row + 1
                            do local2 = 0, degree_a(tangent_a(2))
                                index_a(tangent_a(2)) = first_ta2 + local2
                                do local1 = 0, degree_a(tangent_a(1))
                                    index_a(tangent_a(1)) = first_ta1 + local1
                                    do local_n = 0, degree_a(normal_a)
                                        index_a(normal_a) = first_na(r) + local_n
                                        coefficient_b = basis_ta1(local1)*basis_ta2(local2)*&
                                            range_na**r*basis_na(local_n,r)
                                        if (.not. structural_nonzero(coefficient_b)) cycle
                                        col(pos) = offsets(pa) + index_a(1) + (index_a(2)-1)*nc_a(1) + &
                                            (index_a(3)-1)*nc_a(1)*nc_a(2)
                                        val(pos) = coefficient_b
                                        pos = pos + 1
                                    end do
                                end do
                            end do
                            do local2 = 0, degree_b(mapped_b(2))
                                index_b(mapped_b(2)) = first_tb2(0) + local2
                                do local1 = 0, degree_b(mapped_b(1))
                                    index_b(mapped_b(1)) = first_tb1(0) + local1
                                    do local_n = 0, degree_b(normal_b)
                                        index_b(normal_b) = first_nb(0) + local_n
                                        coefficient_b = 0.0_rk
                                        choose = 1.0_rk
                                        do j = 0, r
                                            mixed = 0.0_rk
                                            do alpha_n = 0, j
                                                do alpha2 = 0, j - alpha_n
                                                    do alpha1 = 0, j - alpha_n - alpha2
                                                        mixed = mixed + composition(j,alpha1,alpha2,alpha_n)*&
                                                            range_tb1**alpha1*basis_tb1(local1,alpha1)*&
                                                            range_tb2**alpha2*basis_tb2(local2,alpha2)*&
                                                            range_nb**alpha_n*basis_nb(local_n,alpha_n)
                                                    end do
                                                end do
                                            end do
                                            coefficient_b = coefficient_b + &
                                                choose*projective_value(r-j)*mixed
                                            if (j < r) choose = choose*real(r-j, rk)/real(j+1, rk)
                                        end do
                                        if (.not. structural_nonzero(coefficient_b)) cycle
                                        col(pos) = offsets(pb) + index_b(1) + (index_b(2)-1)*nc_b(1) + &
                                            (index_b(3)-1)*nc_b(1)*nc_b(2)
                                        val(pos) = -coefficient_b
                                        pos = pos + 1
                                    end do
                                end do
                            end do
                            rowptr(row+1) = pos
                        end do
                    end do
                end do
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the number of recorded face connections.
    pure elemental integer function get_nconnection_volume(this) result(n)
        class(nurbs_multipatch_volume), intent(in) :: this

        n = this%nconnection
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return a value copy of patch `i`.
    !!
    !! An out-of-range one-based identifier returns a default empty patch.
    pure function get_patch_volume(this, i) result(patch)
        class(nurbs_multipatch_volume), intent(in) :: this
        integer, intent(in) :: i
        type(nurbs_volume) :: patch

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
    pure function get_connection_volume(this, i) result(conn)
        class(nurbs_multipatch_volume), intent(in) :: this
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
    !! control-point count across all volume patches.
    pure function cmp_dof_offsets_volume(this) result(offsets)
        class(nurbs_multipatch_volume), intent(in) :: this
        integer, allocatable :: offsets(:)
        integer :: i, nc(3)

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
    !> Build compact global identifiers for directly shareable face controls.
    !! The complete face orientation map is applied only when both normal
    !! endpoint basis vectors are exactly one-hot. For a rational tensor-product
    !! trace,
    !!
    !! \[
    !! R_{jk}^{(P)}(\tau_1,\tau_2)=
    !! \frac{N_j(\tau_1)M_k(\tau_2)\widehat w_{jk}^{(P)}}
    !! {\sum_{a,b}N_a(\tau_1)M_b(\tau_2)\widehat w_{ab}^{(P)}} ,
    !! \]
    !!
    !! corresponding scalar controls are shareable only when all oriented trace
    !! weights differ by one constant nonzero scale. A tangentially varying
    !! projective factor changes these rational basis functions and remains in
    !! `cmp_dof_constraint`. General unclamped \(C^0\) traces likewise remain
    !! constraints. The map trusts recorded topology; call `is_valid` first.
    !!
    !! Disjoint-set implementation reference: R. E. Tarjan, *JACM* 22 (1975),
    !! 215--225.
    !! [doi:10.1145/321879.321884](https://doi.org/10.1145/321879.321884).
    pure function cmp_dof_map_volume(this) result(map)
        class(nurbs_multipatch_volume), intent(in) :: this
        integer, allocatable :: map(:)
        integer, allocatable :: parent(:), offsets(:), ia(:), ib(:)
        real(rk), allocatable :: knot_a(:), knot_b(:)
        real(rk) :: weight_scale
        integer :: i, j, total, pa, pb, dir_a, dir_b, dega, degb, layer_a, layer_b
        integer :: nca(3), ncb(3)
        logical :: compatible_weights, flip_b(2)

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
                call connection_weight_scale_volume(&
                    this,&
                    conn,&
                    128.0_rk*epsilon(1.0_rk),&
                    weight_scale,&
                    compatible_weights)
                if (.not. compatible_weights) cycle
                flip_b(1) = conn%is_flipped(1)
                flip_b(2) = conn%is_flipped(2)
                call boundary_layer_index(nca, conn%get_side_a(), layer_a, .false., .false., no_flip_, ia)
                call boundary_layer_index(ncb, conn%get_side_b(), layer_b, &
                    conn%is_reversed(), conn%is_swapped(), flip_b, ib)
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
    !> Check projective compatibility of rational face traces.
    !!
    !! Normal endpoint basis values are applied before comparing the fully
    !! oriented tangential denominator coefficients. Compatible traces satisfy
    !! \(\widehat w_j^{(A)}=\lambda\widehat w_j^{(B)}\) within `tol`;
    !! `weight_scale` returns \(\lambda\). Polynomial traces use unit scale.
    pure subroutine connection_weight_scale_volume(this, conn, tol, weight_scale, ok)
        class(nurbs_multipatch_volume), intent(in) :: this
        type(multipatch_connection), intent(in) :: conn
        real(rk), intent(in) :: tol
        real(rk), intent(out) :: weight_scale
        logical, intent(out) :: ok
        real(rk), allocatable :: basis_a(:), basis_b(:), knot_a(:), knot_b(:), trace_a(:), trace_b(:)
        real(rk) :: xa, xb
        integer, allocatable :: ia(:), ib(:)
        integer :: j, layer, pa, pb, dir_a, dir_b, dega, degb, first_a, first_b
        integer :: normal_a, normal_b, active_a, active_b
        integer :: nca(3), ncb(3)
        logical :: rational_a, rational_b
        logical :: flip_b(2)

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

        flip_b(1) = conn%is_flipped(1)
        flip_b(2) = conn%is_flipped(2)
        allocate(trace_a(product(nca)/nca(dir_a)), source=0.0_rk)
        do active_a = 0, dega
            normal_a = first_a + active_a
            if (normal_a < 1 .or. normal_a > nca(dir_a)) cycle
            if (.not. structural_nonzero(basis_a(active_a))) cycle
            layer = merge(normal_a - 1, nca(dir_a) - normal_a, mod(conn%get_side_a(), 2) == 1)
            call boundary_layer_index(nca, conn%get_side_a(), layer, .false., .false., no_flip_, ia)
            do j = 1, size(ia)
                trace_a(j) = trace_a(j) + basis_a(active_a)*this%patches(pa)%get_Wc(ia(j))
            end do
        end do
        allocate(trace_b(product(ncb)/ncb(dir_b)), source=0.0_rk)
        do active_b = 0, degb
            normal_b = first_b + active_b
            if (normal_b < 1 .or. normal_b > ncb(dir_b)) cycle
            if (.not. structural_nonzero(basis_b(active_b))) cycle
            layer = merge(normal_b - 1, ncb(dir_b) - normal_b, mod(conn%get_side_b(), 2) == 1)
            call boundary_layer_index(ncb, conn%get_side_b(), layer, &
                conn%is_reversed(), conn%is_swapped(), flip_b, ib)
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
    !> Return the two tangential parameter directions for a face identifier.
    !!
    !! Faces 1--2 are normal to direction 1, faces 3--4 to direction 2, and
    !! faces 5--6 to direction 3. Invalid identifiers return `[0,0]`.
    pure subroutine volume_tangent_dirs(side, dirs)
        integer, intent(in) :: side
        integer, intent(out) :: dirs(2)

        select case ((side + 1)/2)
        case (1)
            dirs(1) = 2
            dirs(2) = 3
        case (2)
            dirs(1) = 1
            dirs(2) = 3
        case (3)
            dirs(1) = 1
            dirs(2) = 2
        case default
            dirs(1) = 0
            dirs(2) = 0
        end select
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Test whether connected faces have compatible tangential spline spaces.
    !!
    !! After applying swap and flip metadata, both tangential traces must have
    !! equal degrees, equal control counts, and affine-equivalent normalized
    !! knot vectors.
    pure logical function connection_trace_space_volume(this, conn, tol) result(ok)
        class(nurbs_multipatch_volume), intent(in) :: this
        type(multipatch_connection), intent(in) :: conn
        real(rk), intent(in) :: tol
        real(rk), allocatable :: knot_a(:), knot_b(:)
        integer :: pa, pb, tan_a(2), tan_b(2), mapped_b(2), a
        logical :: reverse_dir(2)

        pa = conn%get_patch_a()
        pb = conn%get_patch_b()
        call volume_tangent_dirs(conn%get_side_a(), tan_a)
        call volume_tangent_dirs(conn%get_side_b(), tan_b)
        if (any(tan_a == 0) .or. any(tan_b == 0)) then
            ok = .false.
            return
        end if

        if (conn%is_swapped()) then
            mapped_b(1) = tan_b(2)
            mapped_b(2) = tan_b(1)
            reverse_dir(1) = conn%is_flipped(2)
            reverse_dir(2) = conn%is_flipped(1)
        else
            mapped_b(1) = tan_b(1)
            mapped_b(2) = tan_b(2)
            reverse_dir(1) = conn%is_flipped(1)
            reverse_dir(2) = conn%is_flipped(2)
        end if
        if (conn%is_reversed()) then
            reverse_dir(1) = .not. reverse_dir(1)
            reverse_dir(2) = .not. reverse_dir(2)
        end if

        ok = .true.
        do a = 1, 2
            knot_a = this%patches(pa)%get_knot(tan_a(a))
            knot_b = this%patches(pb)%get_knot(mapped_b(a))
            if (.not. multipatch_compatible_trace_space(knot_a, this%patches(pa)%get_nc(tan_a(a)), &
                this%patches(pa)%get_degree(tan_a(a)), knot_b, this%patches(pb)%get_nc(mapped_b(a)), &
                this%patches(pb)%get_degree(mapped_b(a)), reverse_dir(a), tol)) then
                ok = .false.
                return
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Assemble sparse parametric and geometric constraints for connected faces.
    !!
    !! The one-based CSR arrays encode \(Aq=0\) over unmerged,
    !! patch-concatenated control variables. For each mapped tangential face
    !! control, rows enforce normalized normal derivatives from order zero
    !! through the requested order. `geometric` optionally filters parametric
    !! (false) or geometric (true) connections. Compact geometric rows use the
    !! scalar normal transition; general rows include all tangentially varying
    !! multivariate transition and projective jet terms. Rational projective
    !! scaling is included for compatible boundary weights.
    !! Apply the raw matrix to polynomial coefficients or homogeneous geometry
    !! columns. For a rational scalar field with coefficients \(q_i\), apply it
    !! to weighted numerator controls \(w_iq_i\), not directly to \(q_i\).
    !!
    !! For every mapped pair of tangential control indices and derivative order
    !! \(r\), one row discretizes
    !!
    !! \[
    !! \partial_{\eta_A}^{\,r}f_A(\tau_1,\tau_2)
    !! -\sum_{k=0}^{r}B_{r,k}(\phi',\ldots)
    !!  \partial_{\eta_B}^{\,k}f_B(\tau_1,\tau_2)=0.
    !! \]
    !!
    !! Parameter-range powers convert stored-knot derivatives to normalized
    !! normal derivatives. Call `is_valid` first; this routine assumes
    !! compatible connections. Coefficients with absolute value at most
    !! `128*epsilon(rk)` are omitted from CSR storage, so the returned matrix is
    !! numerically sparsified rather than symbolically exact.
    !!
    pure subroutine cmp_dof_constraint_volume(this, rowptr, col, val, geometric)
        class(nurbs_multipatch_volume), intent(in) :: this
            !! Valid multipatch volume.
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
        integer :: nca(3), ncb(3), dega, degb, max_layer_a, max_layer_b
        logical :: compatible_weights, flip_b(2)

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
                    call general_constraint_shape_volume(this, conn, general_nrow, general_nnz)
                    nrow = nrow + general_nrow
                    nnz = nnz + general_nnz
                    cycle
                end if
                pa = conn%get_patch_a()
                pb = conn%get_patch_b()
                call boundary_index(this%patches(pa)%get_nc(), conn%get_side_a(), .false., .false., no_flip_, ia)
                ntan = size(ia)
                dir_a = (conn%get_side_a() + 1)/2
                dir_b = (conn%get_side_b() + 1)/2
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
                    call append_general_constraint_volume(&
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
                flip_b(1) = conn%is_flipped(1)
                flip_b(2) = conn%is_flipped(2)
                nca = this%patches(pa)%get_nc()
                ncb = this%patches(pb)%get_nc()
                dir_a = (side_a + 1)/2
                dir_b = (side_b + 1)/2
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
                call connection_weight_scale_volume(this, conn, sqrt(epsilon(1.0_rk)), weight_scale, compatible_weights)
                if (.not. compatible_weights) weight_scale = 1.0_rk
                allocate(ca(0:dega,0:conn%get_continuity()), cb(0:degb,0:conn%get_continuity()))
                allocate(first_a(0:conn%get_continuity()), first_b(0:conn%get_continuity()))
                do r = 0, conn%get_continuity()
                    call basis_bspline_der_order_active(xa, knot_a, nca(dir_a), dega, r, first_a(r), ca(:,r))
                    call basis_bspline_der_order_active(xb, knot_b, ncb(dir_b), degb, r, first_b(r), cb(:,r))
                end do
                call multipatch_chain_rule_coefficients(conn, normal_sign, chain)
                call boundary_index(nca, side_a, .false., .false., no_flip_, ia)
                ntan = size(ia)
                max_layer_a = min(nca(dir_a) - 1, dega)
                max_layer_b = min(ncb(dir_b) - 1, degb)
                allocate(layer_a(ntan, max_layer_a + 1), layer_b(ntan, max_layer_b + 1))
                do layer = 0, max_layer_a
                    call boundary_layer_index(nca, side_a, layer, .false., .false., no_flip_, ia)
                    layer_a(:, layer + 1) = ia
                end do
                do layer = 0, max_layer_b
                    call boundary_layer_index(ncb, side_b, layer, conn%is_reversed(), conn%is_swapped(), flip_b, ib)
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
    !> Concatenate active tensor-product elements from all volume patches.
    !! Rows follow patch insertion order and direction-1-fastest local element
    !! order. Different local basis sizes are zero-padded to the largest row.
    !! `shared=.true.` applies direct sharing from the compact map; the default
    !! retains unmerged patch-offset numbering. General unclamped traces and
    !! rational traces with nonconstant projective scale still require
    !! `cmp_dof_constraint`.
    pure function cmp_elem_volume(this, shared) result(elem)
        class(nurbs_multipatch_volume), intent(in) :: this
        logical, intent(in), optional :: shared
        integer, allocatable :: elem(:,:)
        integer, allocatable :: offsets(:), map(:), mul1(:), mul2(:), mul3(:), nspan_patch(:,:), row_start(:)
        real(rk), allocatable :: knot1(:), knot2(:), knot3(:)
        integer :: a, i, row, ne, maxnen, nc(3), degree(3), nspan1, nspan2, nspan3
        integer :: e1, e2, e3, i1, i2, i3, pref1, pref2, pref3, node, nc12, np
        logical :: use_shared

        use_shared = .false.
        if (present(shared)) use_shared = shared
        np = this%get_npatch()
        offsets = this%cmp_dof_offsets()
        ne = 0
        maxnen = 0
        allocate(nspan_patch(3, np), row_start(np + 1))
        row_start(1) = 0
        do i = 1, np
            nc = this%patches(i)%get_nc()
            degree = this%patches(i)%get_degree()
            knot1 = this%patches(i)%get_knot(1)
            knot2 = this%patches(i)%get_knot(2)
            knot3 = this%patches(i)%get_knot(3)
            nspan_patch(1, i) = active_span_count(knot1, nc(1), degree(1))
            nspan_patch(2, i) = active_span_count(knot2, nc(2), degree(2))
            nspan_patch(3, i) = active_span_count(knot3, nc(3), degree(3))
            row_start(i + 1) = row_start(i) + nspan_patch(1, i)*nspan_patch(2, i)*nspan_patch(3, i)
            maxnen = max(maxnen, (degree(1) + 1)*(degree(2) + 1)*(degree(3) + 1))
        end do
        ne = row_start(np + 1)
        allocate(elem(ne, maxnen), source=0)
        if (use_shared) map = this%cmp_dof_map()
        do i = 1, np
            nc = this%patches(i)%get_nc()
            degree = this%patches(i)%get_degree()
            knot1 = this%patches(i)%get_knot(1)
            knot2 = this%patches(i)%get_knot(2)
            knot3 = this%patches(i)%get_knot(3)
            nspan1 = nspan_patch(1, i)
            nspan2 = nspan_patch(2, i)
            nspan3 = nspan_patch(3, i)
            if (nspan1 < 1 .or. nspan2 < 1 .or. nspan3 < 1) cycle
            mul1 = active_knot_multiplicity(knot1, nc(1), degree(1))
            mul2 = active_knot_multiplicity(knot2, nc(2), degree(2))
            mul3 = active_knot_multiplicity(knot3, nc(3), degree(3))
            if (size(mul1) < nspan1 + 1 .or. size(mul2) < nspan2 + 1 .or. size(mul3) < nspan3 + 1) cycle
            nc12 = nc(1)*nc(2)
            row = row_start(i)
            pref3 = degree(3) + 1
            do e3 = 1, nspan3
                pref2 = degree(2) + 1
                do e2 = 1, nspan2
                    pref1 = degree(1) + 1
                    do e1 = 1, nspan1
                        row = row + 1
                        do i3 = 0, degree(3)
                            do i2 = 0, degree(2)
                                do i1 = 0, degree(1)
                                    a = i3*(degree(2) + 1)*(degree(1) + 1) + i2*(degree(1) + 1) + i1 + 1
                                    node = ((pref3 - degree(3) + i3) - 1)*nc12 + &
                                        ((pref2 - degree(2) + i2) - 1)*nc(1) + pref1 - degree(1) + i1
                                    if (use_shared) then
                                        elem(row,a) = map(offsets(i) + node)
                                    else
                                        elem(row,a) = offsets(i) + node
                                    end if
                                end do
                            end do
                        end do
                        if (e1 < nspan1) pref1 = pref1 + mul1(e1+1)
                    end do
                    if (e2 < nspan2) pref2 = pref2 + mul2(e2+1)
                end do
                if (e3 < nspan3) pref3 = pref3 + mul3(e3+1)
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the owning patch identifier for every concatenated element.
    !!
    !! Result ordering is identical to [[nurbs_multipatch_volume:cmp_elem]].
    pure function cmp_elem_patch_volume(this) result(patch_id)
        class(nurbs_multipatch_volume), intent(in) :: this
        integer, allocatable :: patch_id(:)
        integer, allocatable :: nloc(:), row_start(:)
        integer :: i, e, ne, np

        np = this%get_npatch()
        allocate(nloc(np), row_start(np + 1))
        row_start(1) = 0
        do i = 1, np
            nloc(i) = active_span_count(this%patches(i)%get_knot(1), this%patches(i)%get_nc(1), &
                this%patches(i)%get_degree(1)) * active_span_count(this%patches(i)%get_knot(2), &
                this%patches(i)%get_nc(2), this%patches(i)%get_degree(2)) * &
                active_span_count(this%patches(i)%get_knot(3), this%patches(i)%get_nc(3), &
                this%patches(i)%get_degree(3))
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
    !! Result ordering is identical to [[nurbs_multipatch_volume:cmp_elem]].
    pure function cmp_elem_local_volume(this) result(local_id)
        class(nurbs_multipatch_volume), intent(in) :: this
        integer, allocatable :: local_id(:)
        integer, allocatable :: nloc(:), row_start(:)
        integer :: i, e, ne, np

        np = this%get_npatch()
        allocate(nloc(np), row_start(np + 1))
        row_start(1) = 0
        do i = 1, np
            nloc(i) = active_span_count(this%patches(i)%get_knot(1), this%patches(i)%get_nc(1), &
                this%patches(i)%get_degree(1)) * active_span_count(this%patches(i)%get_knot(2), &
                this%patches(i)%get_nc(2), this%patches(i)%get_degree(2)) * &
                active_span_count(this%patches(i)%get_knot(3), this%patches(i)%get_nc(3), &
                this%patches(i)%get_degree(3))
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
    !> Validate patch states and every oriented volume interface.
    !!
    !! A valid multipatch contains at least one complete, internally consistent
    !! volume patch. Checks also include identifiers, face maps, continuity order, transition jets,
    !! both tangential trace spaces, physical endpoint traces,
    !! homogeneous normal-derivative residuals in the compact or multivariate model,
    !! physical dimensions, and rational projective weight compatibility. `tol`
    !! is one absolute threshold for assembled homogeneous residuals; the
    !! default is `sqrt(epsilon(rk))`. Each component of \(A\mathbf H\) is
    !! compared directly with `tol`.
    !! Weight proportionality uses
    !! `max(tol,128*epsilon(rk))*max(abs(wA),abs(lambda*wB))` per pair.
    !! `tol>=0` is an unchecked caller precondition. The test conservatively
    !! requires equal whole-patch `is_rational` classifications when no
    !! explicit projective jet is supplied. A full projective jet validates the
    !! complete selected homogeneous-lift relation instead.
    pure function is_valid_volume(this, tol) result(ok)
        class(nurbs_multipatch_volume), intent(in) :: this
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
                if (conn%get_side_a() < u_min_ .or. conn%get_side_a() > w_max_ .or. &
                    conn%get_side_b() < u_min_ .or. conn%get_side_b() > w_max_) then
                    ok = .false.; return
                end if
                if (conn%get_continuity() < -1) then
                    ok = .false.; return
                end if
                normal_sign = merge(-1, 1, mod(conn%get_side_a(), 2) == mod(conn%get_side_b(), 2))
                if (.not. multipatch_valid_reparameterization(conn, normal_sign, 3)) then
                    ok = .false.; return
                end if
                dir_a = (conn%get_side_a() + 1)/2
                dir_b = (conn%get_side_b() + 1)/2
                dega = this%patches(pa)%get_degree(dir_a)
                degb = this%patches(pb)%get_degree(dir_b)
                if (conn%get_continuity() > min(dega, degb)) then
                    ok = .false.; return
                end if
                if (conn%get_continuity() >= 0) then
                    if (.not. connection_trace_space_volume(this, conn, eps)) then
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
                        call connection_weight_scale_volume(this, conn, eps, weight_scale, compatible_weights)
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
    !> Sample every volume patch with common directional sampling arguments.
    !! Arguments are forwarded unchanged to [[nurbs_volume:create]]. Any patch
    !! failure is converted to multipatch diagnostic code 205.
    pure subroutine create_volume(this, res1, res2, res3, Xt1, Xt2, Xt3, Xt)
        class(nurbs_multipatch_volume), intent(inout) :: this
        integer, intent(in), optional :: res1, res2, res3
        real(rk), intent(in), contiguous, optional :: Xt1(:), Xt2(:), Xt3(:)
        real(rk), intent(in), contiguous, optional :: Xt(:,:)
        integer :: i

        if (.not. this%err%ok) return

        if (this%get_npatch() < 1) then
            call this%err%set(&
                code       = 205,&
                severity   = 1,&
                category   = "forcad_multipatch_volume",&
                message    = "No patches are available.",&
                location   = "create",&
                suggestion = "Add at least one patch before calling create().")
            return
        end if
        do i = 1, this%get_npatch()
            call this%patches(i)%create(res1=res1, res2=res2, res3=res3, Xt1=Xt1, Xt2=Xt2, Xt3=Xt3, Xt=Xt)
            if (.not. this%patches(i)%err%ok) then
                call this%err%set(&
                    code       = 205,&
                    severity   = 1,&
                    category   = "forcad_multipatch_volume",&
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
    pure elemental logical function is_rational_volume(this) result(r)
        class(nurbs_multipatch_volume), intent(in) :: this
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
    !> Export every patch control lattice to a separate legacy VTK file.
    !!
    !! Patch `i` is written as `<prefix>_Xc_<i>.vtk`.
    impure subroutine export_Xc_volume(this, prefix, encoding)
        class(nurbs_multipatch_volume), intent(inout) :: this
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
    impure subroutine export_Xg_volume(this, prefix, encoding)
        class(nurbs_multipatch_volume), intent(inout) :: this
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
    !> Export parameter lines mapped into every volume patch.
    impure subroutine export_Xth_in_Xg_volume(this, prefix, res, encoding)
        class(nurbs_multipatch_volume), intent(inout) :: this
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
    !> Apply one x-y-z Euler rotation in degrees to every control lattice.
    !!
    !! Cached physical samples are not changed.
    pure subroutine rotate_Xc_volume(this, alpha, beta, theta)
        class(nurbs_multipatch_volume), intent(inout) :: this
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
    pure subroutine rotate_Xg_volume(this, alpha, beta, theta)
        class(nurbs_multipatch_volume), intent(inout) :: this
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
    !> Translate every patch control lattice by physical vector `vec`.
    !!
    !! Cached physical samples are not changed.
    pure subroutine translate_Xc_volume(this, vec)
        class(nurbs_multipatch_volume), intent(inout) :: this
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
    pure subroutine translate_Xg_volume(this, vec)
        class(nurbs_multipatch_volume), intent(inout) :: this
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
    impure subroutine show_volume(this, vtkfile_Xc, vtkfile_Xg, vtkfile_Xth_in_Xg)
        class(nurbs_multipatch_volume), intent(inout) :: this
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
            rank_name         = "multipatch volume")
#endif
    end subroutine
    !===============================================================================
end module forcad_multipatch_volume
