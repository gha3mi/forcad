! author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Topology, continuity constraints, and assembly data for curve patches.
!!
!! Patches remain complete [[nurbs_curve]] objects. The container records
!! oriented endpoint connections; [[nurbs_multipatch_curve:is_valid]] verifies
!! their physical and derivative traces. Two complementary global algebraic
!! objects are exposed:
!!
!! * [[nurbs_multipatch_curve:cmp_dof_map]] composes intrinsic periodic equivalences and corresponding \(C^0\)
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

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    use forcad_kinds, only: rk
    use forcad_multipatch, only: multipatch_connection, multipatch_growth_capacity, check_new_patch, check_patch_transform
    use forcad_multipatch, only: constraint_rows, valid_multipatch, sample_count
    use forcad_multipatch, only: resolve_patch_interface, record_connection, require_patches
    use forcad_multipatch, only: multipatch_dof_map, multipatch_dof_offsets, assemble_design_gradient
    use forcad_multipatch, only: multipatch_elements, element_layout, multipatch_element_ids, requested
    use forcad_nurbs_curve, only: nurbs_curve
    use forcad_utils, only: active_span_count, export_vtk_legacy
    use fordebug, only: debug

    implicit none

    private
    public :: nurbs_multipatch_curve


    !===============================================================================
    !> Ordered collection of connected NURBS curve patches.
    !!
    !! Valid side labels accepted by `connect` are `left`, `right`, `min`, and
    !! `max` (case variants shown in the implementation are accepted).
    !! `continuity=-1` retains an interface without continuity equations.
    !! Nonnegative orders appear in the sparse constraint matrix; direct
    !! endpoint sharing is additionally available for interpolatory endpoints.
    !!
    !! Sampling is independent of the IGA knot discretization. In
    !! `create(res=n)`, `n` is broadcast to every patch; `create(res=r)` with
    !! rank-one integer `r` uses `r(k)` for patch `k`; and
    !! `create(target_size=h)` estimates every patch count from its physical
    !! control-polygon length.
    !!
    !! Diagnostic codes are: 201 invalid patch identifier; 202 invalid side;
    !! 203 invalid continuity metadata; 204 incompatible interface geometry or
    !! trace space; 205 invalid patch state; 206 invalid sampling request; and
    !! 207 invalid transformation data and 208 invalid or unrepresentable sensitivity data.
    type nurbs_multipatch_curve
        type(nurbs_curve), allocatable, private :: patches(:) !! Owned patch copies with spare append capacity.
        type(multipatch_connection), allocatable, private :: connections(:) !! Validated interface metadata.
        integer, private :: npatch = 0 !! Number of active patches.
        integer, private :: nconnection = 0 !! Number of active connections.

        type(debug) :: err !! Recoverable diagnostic state for the most recent failed operation.
    contains
        procedure :: add_patch => add_patch_curve !!> Append a validated patch copy and optionally return its identifier.
        procedure :: connect => connect_curve !!> Validate connection metadata and record an oriented endpoint interface.
        procedure, private :: create_common => create_curve
        procedure, private :: create_patches => create_curve_patches
        procedure, private :: create_automatic => create_curve_automatic
        generic :: create => create_common, create_patches, create_automatic
            !!> Sample patches with common, patch-specific, or physically scaled resolutions.
        procedure :: finalize => clear_curve !!> Release patches and connections and reset counts.
        procedure :: get_npatch => get_npatch_curve !!> Return the active patch count.
        procedure :: get_nconnection => get_nconnection_curve !!> Return the active connection count.
        procedure :: get_patch => get_patch_curve !!> Return a copy of one patch.
        procedure :: get_connection => get_connection_curve !!> Return a copy of one connection.
        procedure :: is_valid => is_valid_curve !!> Validate patch states, interfaces, orientations, and requested continuity.
        procedure :: is_rational => is_rational_curve !!> Report whether any patch has nonuniform rational weights.
        procedure :: cmp_dof_offsets => cmp_dof_offsets_curve !!> Return prefix offsets for unmerged patch control variables.
        procedure :: assemble_sensitivities => assemble_sensitivities_curve !!> Assemble independent control and log-weight gradients.
        procedure :: cmp_dof_map => cmp_dof_map_curve !!> Compact directly shareable \(C^0\) controls.
        procedure :: cmp_dof_constraint => cmp_dof_constraint_curve !!> Assemble derivative-continuity equations in CSR form.
        procedure :: cmp_elem => cmp_elem_curve !!> Concatenate patch element connectivity with optional shared numbering.
        procedure :: get_elem => cmp_elem_curve !!> Alias for [[nurbs_multipatch_curve:cmp_elem]].
        procedure :: cmp_elem_patch => cmp_elem_patch_curve !!> Return the owning patch identifier for every global element.
        procedure :: cmp_elem_local => cmp_elem_local_curve !!> Return the patch-local identifier for every global element.
        procedure :: export_Xc => export_Xc_curve !!> Export separate or combined patch control polygons.
        procedure :: export_Xg => export_Xg_curve !!> Export separate or combined cached patch geometry.
        procedure :: export_Xth_in_Xg => export_Xth_in_Xg_curve !!> Export separate or combined embedded parameter lines.
        procedure :: transform => transform_curve !!> Apply one affine map to every patch and cache.
        procedure :: rotate_Xc => rotate_Xc_curve !!> Rotate all patch control points.
        procedure :: rotate_Xg => rotate_Xg_curve !!> Rotate all cached patch geometry points.
        procedure :: translate_Xc => translate_Xc_curve !!> Translate all patch control points.
        procedure :: translate_Xg => translate_Xg_curve !!> Translate all cached patch geometry points.
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

        call check_new_patch(1,patch,patch%err%ok,this%err,id)
        if (.not.this%err%ok) return
        if (.not. allocated(this%patches)) then
            allocate(this%patches(multipatch_growth_capacity(0)))
        else if (this%npatch == size(this%patches)) then
            allocate(tmp(multipatch_growth_capacity(this%npatch)))
            tmp(:this%npatch) = this%patches
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
    !> Apply one common affine transformation to all curve patches.
    !!
    !! Every control polygon and existing sampled-geometry cache receives
    !! \(\widetilde{\mathbf x}=\mathbf c+\mathbf A(\mathbf x-\mathbf c)+
    !! \mathbf b\). A common affine map preserves matching interface jets and
    !! homogeneous projective relations, so recorded \(C^n\) and \(G^n\)
    !! connections remain valid whenever the original multipatch was valid.
    pure subroutine transform_curve(this, linear, translation, center)
        class(nurbs_multipatch_curve), intent(inout) :: this
        real(rk), intent(in), contiguous :: linear(:,:)
            !! Square linear map shared by every patch.
        real(rk), intent(in), contiguous, optional :: translation(:)
            !! Shared translation vector.
        real(rk), intent(in), contiguous, optional :: center(:)
            !! Shared transformation center.
        integer :: i

        call check_patch_transform(this%patches,this%npatch,1,linear,translation,center,this%err)
        if (.not.this%err%ok) return
        do i=1,this%npatch
            call this%patches(i)%transform(linear,translation,center)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export one combined visualization dataset without merging patch points.
    !!
    !! Keeping patch-local points distinct preserves discontinuous point data.
    !! `patch_id` records the one-based owner of every output cell.
    impure subroutine export_combined_curve(this, filename, source, res, encoding)
        class(nurbs_multipatch_curve), intent(inout) :: this
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: source
        integer, intent(in), optional :: res
        character(len=*), intent(in), optional :: encoding

        character(len=8), parameter :: cell_field_names(1) = ["patch_id"]
        real(rk), allocatable :: patch_points(:,:), points(:,:), cell_data(:,:)
        real(rk), allocatable :: knot(:)
        integer, allocatable :: patch_elem(:,:), elemConn(:,:)
        integer :: cell, cell_offset, i, line, node, point_offset
        integer :: ncell, ncoord, nnode, npoint, patch_npoint, res_min, span_count, vtk_cell_type

        npoint = 0
        ncell = 0
        nnode = 0
        do i = 1, this%get_npatch()
            select case (source)
            case ("Xc")
                patch_npoint = this%patches(i)%get_nc()
                ncoord = size(this%patches(i)%get_Xc(1))
                patch_elem = this%patches(i)%get_elem_Xc_vis()
                if (size(patch_elem) == 0) then
                    patch_elem = this%patches(i)%cmp_elem_Xc_vis()
                end if
            case ("Xg")
                patch_npoint = this%patches(i)%get_ng()
                ncoord = size(this%patches(i)%get_Xg(1))
                patch_elem = this%patches(i)%get_elem_Xg_vis()
                if (size(patch_elem) == 0) then
                    patch_elem = this%patches(i)%cmp_elem_Xg_vis()
                end if
            case ("Xth")
                knot = this%patches(i)%get_knot()
                span_count = active_span_count(&
                    knot,&
                    this%patches(i)%get_nc(),&
                    this%patches(i)%get_degree())
                ncoord = size(this%patches(i)%get_Xc(1))
                if (.not. this%patches(i)%err%ok .or. span_count < 1 .or. ncoord < 2 .or. ncoord > 3) then
                    call this%err%set(&
                        code       = 205,&
                        severity   = 1,&
                        category   = "forcad_multipatch_curve",&
                        message    = "A patch visualization dataset is unavailable.",&
                        location   = "export_combined",&
                        suggestion = "Initialize and sample every patch before combined export.")
                    return
                end if
                res_min = 10
                if (present(res)) res_min = max(2,res)
                patch_npoint = span_count*(res_min-1) + 1
                npoint = npoint+patch_npoint
                ncell = ncell+patch_npoint-1
                nnode = 2
                cycle
            case default
                call this%err%set(&
                    code       = 205,&
                    severity   = 1,&
                    category   = "forcad_multipatch_curve",&
                    message    = "Unknown visualization representation.",&
                    location   = "export_combined",&
                    suggestion = "Use a supported multipatch export procedure.")
                return
            end select
            if (.not. allocated(patch_elem)) allocate(patch_elem(0,0))
            if (.not. this%patches(i)%err%ok .or. patch_npoint == 0 .or. &
                size(patch_elem,1) == 0) then
                call this%err%set(&
                    code       = 205,&
                    severity   = 1,&
                    category   = "forcad_multipatch_curve",&
                    message    = "A patch visualization dataset is unavailable.",&
                    location   = "export_combined",&
                    suggestion = "Initialize and sample every patch before combined export.")
                return
            end if
            if (ncoord < 1 .or. ncoord > 3) then
                call this%err%set(&
                    code       = 205,&
                    severity   = 1,&
                    category   = "forcad_multipatch_curve",&
                    message    = "A patch has an unsupported visualization dimension.",&
                    location   = "export_combined",&
                    suggestion = "Use one-, two-, or three-coordinate patch geometry.")
                return
            end if
            npoint = npoint+patch_npoint
            if (nnode == 0) nnode = size(patch_elem,2)
            if (size(patch_elem,2) /= nnode) then
                call this%err%set(&
                    code       = 205,&
                    severity   = 1,&
                    category   = "forcad_multipatch_curve",&
                    message    = "Patch visualization cell orders differ.",&
                    location   = "export_combined",&
                    suggestion = "Use one visualization cell order for all patches.")
                return
            end if
            ncell = ncell+size(patch_elem,1)
        end do

        vtk_cell_type = 3
        allocate(points(npoint,3), elemConn(ncell,nnode), cell_data(ncell,1))
        points = 0.0_rk
        point_offset = 0
        cell_offset = 0
        do i = 1, this%get_npatch()
            select case (source)
            case ("Xc")
                patch_points = this%patches(i)%get_Xc()
                patch_elem = this%patches(i)%get_elem_Xc_vis()
                if (size(patch_elem) == 0) patch_elem = this%patches(i)%cmp_elem_Xc_vis()
            case ("Xg")
                patch_points = this%patches(i)%get_Xg()
                patch_elem = this%patches(i)%get_elem_Xg_vis()
                if (size(patch_elem) == 0) patch_elem = this%patches(i)%cmp_elem_Xg_vis()
            case ("Xth")
                call this%patches(i)%cmp_Xth_in_Xg(&
                    Xg       = patch_points,&
                    elemConn = patch_elem,&
                    res      = res)
            case default
                return
            end select
            points(point_offset+1:point_offset+size(patch_points,1),1:size(patch_points,2)) = patch_points
            if (source == "Xth") then
                do line = 1, size(patch_elem,1)
                    do node = 1, size(patch_elem,2) - 1
                        cell_offset = cell_offset+1
                        elemConn(cell_offset,:) = [&
                            patch_elem(line,node)+point_offset,&
                            patch_elem(line,node+1)+point_offset]
                        cell_data(cell_offset,1) = real(i,rk)
                    end do
                end do
            else
                do cell = 1, size(patch_elem,1)
                    cell_offset = cell_offset+1
                    elemConn(cell_offset,:) = patch_elem(cell,:)+point_offset
                    cell_data(cell_offset,1) = real(i,rk)
                end do
            end if
            point_offset = point_offset+size(patch_points,1)
        end do
        if (point_offset /= npoint .or. cell_offset /= ncell) then
            call this%err%set(&
                code       = 205,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
                message    = "Combined visualization sizes are inconsistent.",&
                location   = "export_combined",&
                suggestion = "Report the patch data and requested representation.")
            return
        end if

        call export_vtk_legacy(&
            filename         = filename,&
            points           = points,&
            elemConn         = elemConn,&
            vtkCellType      = vtk_cell_type,&
            encoding         = encoding,&
            cell_data        = cell_data,&
            cell_field_names = cell_field_names)
    end subroutine export_combined_curve
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
        integer :: side_a_id, side_b_id, continuity_

        call resolve_patch_interface(this%patches,this%npatch,1,patch_a,side_a,patch_b,side_b,continuity, &
            side_a_id,side_b_id,continuity_,this%err)
        if (.not.this%err%ok) return
        call conn%set(patch_a, side_a_id, patch_b, side_b_id, continuity_, reverse, geometric=geometric, &
            reparameterization=reparameterization)
        call record_connection(this%patches,this%npatch,1,this%connections,this%nconnection,conn,this%err)
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
        offsets = multipatch_dof_offsets(this%patches,this%npatch)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Number scalar controls across periodic seams and directly shareable endpoints.
    !! The input index follows patch storage order. A directly shared trace requires
    !! an interpolatory endpoint on each side and proportional oriented trace weights.
    !! Call `is_valid` to check the geometry. General traces use
    !! [[nurbs_multipatch_curve:cmp_dof_constraint]].
    pure function cmp_dof_map_curve(this, active_only) result(map)
        class(nurbs_multipatch_curve), intent(in) :: this
        logical, optional, intent(in) :: active_only
            !! Number participating controls in first-active storage order and return zero elsewhere.
        integer, allocatable :: map(:)

        call multipatch_dof_map(this%patches,this%connections,this%nconnection,this%cmp_dof_offsets(),active_only,map)
    end function

    !> Assemble patch-local partials in independent design variables. See `assemble_design_gradient`.
    pure subroutine assemble_sensitivities_curve(this,local_dXc,local_dlogWc,dXc,dlogWc,err,map)
        class(nurbs_multipatch_curve), intent(in) :: this
        real(rk), optional, intent(in) :: local_dXc(:,:),local_dlogWc(:)
        real(rk), optional, allocatable, intent(out) :: dXc(:,:),dlogWc(:)
        type(debug), optional, intent(out) :: err
        integer, optional, allocatable, intent(out) :: map(:)
        call assemble_design_gradient(this%patches,this%connections,this%nconnection,this%cmp_dof_offsets(), &
            this%is_valid(),local_dXc,local_dlogWc,dXc,dlogWc,err,map)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Assemble one-based sparse continuity equations for connected endpoints.
    !! Call `is_valid` to verify the geometry. See `constraint_rows`
    !! for the field, solver-numbering and connection-filter options.
    pure subroutine cmp_dof_constraint_curve(this, rowptr, col, val, geometric, field, err, active_only)
        class(nurbs_multipatch_curve), intent(in) :: this
        integer, allocatable, intent(out) :: rowptr(:),col(:)
        real(rk), allocatable, intent(out) :: val(:)
        logical, optional, intent(in) :: geometric,field,active_only
        type(debug), optional, intent(out) :: err
        call constraint_rows(this%patches,this%connections,this%nconnection,this%cmp_dof_offsets(), &
            rowptr,col,val,geometric,field,err,active_only=active_only)
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
    pure function cmp_elem_curve(this, shared, active_only) result(elem)
        class(nurbs_multipatch_curve), intent(in) :: this
            !! Multipatch curve.
        logical, intent(in), optional :: shared
            !! Compact directly shareable \(C^0\) controls; default false.
        logical, intent(in), optional :: active_only
            !! Use the shared active solver IDs, with zero padding for shorter element rows.
        integer, allocatable :: elem(:,:), map(:)

        if (requested(shared) .or. requested(active_only)) map=this%cmp_dof_map(active_only)
        call multipatch_elements(this%patches,this%cmp_dof_offsets(),elem,map)
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
        integer, allocatable :: rows(:)

        call element_layout(this%patches, this%get_npatch(), rows)
        call multipatch_element_ids(rows, local=.false., ids=patch_id)
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
        integer, allocatable :: rows(:)

        call element_layout(this%patches, this%get_npatch(), rows)
        call multipatch_element_ids(rows, local=.true., ids=local_id)
    end function

    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Validate the stored curve patches, their oriented interfaces and homogeneous continuity residuals.
    !! See `valid_multipatch` for the residual test used by geometry sensitivity assembly.
    pure function is_valid_curve(this, rtol, atol) result(ok)
        class(nurbs_multipatch_curve), intent(in) :: this
        real(rk), intent(in), optional :: rtol
            !! Relative trace and homogeneous-residual tolerance; default `sqrt(epsilon(rk))`.
        real(rk), intent(in), optional :: atol
            !! Absolute homogeneous-residual tolerance; default zero.
        logical :: ok

        ok=this%err%ok
        if (ok) ok=valid_multipatch(this%patches,this%npatch,this%connections,this%nconnection,1,rtol,atol)
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

        call require_patches(1,this%npatch,this%err)
        if (.not.this%err%ok) return
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
    !> Sample curve patches with independent resolutions.
    !!
    !! `res(i)` is the number of cached parameter samples for patch `i`.
    !! Curve interfaces are points and therefore impose no tangential sampling
    !! compatibility between neighboring patches.
    pure subroutine create_curve_patches(this, res)
        class(nurbs_multipatch_curve), intent(inout) :: this
        integer, intent(in), contiguous :: res(:)
            !! Positive sample count for every active patch.
        integer :: i

        call require_patches(1,this%npatch,this%err)
        if (.not.this%err%ok) return
        if (size(res) /= this%get_npatch()) then
            call this%err%set(&
                code       = 206,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
                message    = "Patch-resolution count does not match the patch count.",&
                location   = "create",&
                suggestion = "Pass res with one positive entry per active patch.")
            return
        end if
        if (any(res < 1)) then
            call this%err%set(&
                code       = 206,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
                message    = "Patch resolutions must be positive.",&
                location   = "create",&
                suggestion = "Set every entry of res to at least one.")
            return
        end if

        do i = 1, this%get_npatch()
            call this%patches(i)%create(res=res(i))
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
    !> Choose patch resolutions from a physical target segment size.
    !!
    !! For patch \(k\), let \(L_k\) be the length of its control polygon. The
    !! generated sample count is
    !!
    !! \[
    !! n_k=\min\!\left(n_{\max},
    !!     \max\!\left(n_{\min},1+\left\lceil L_k/h\right\rceil\right)\right),
    !! \]
    !!
    !! where \(h=\mathtt{target\_size}\). Positive rational weights preserve
    !! the convex-hull and variation-diminishing construction underlying this
    !! inexpensive physical-space density estimate. It controls visualization
    !! density; it is not a certified Hausdorff-error bound.
    pure subroutine create_curve_automatic(this, target_size, min_res, max_res)
        class(nurbs_multipatch_curve), intent(inout) :: this
        real(rk), intent(in) :: target_size
            !! Positive finite target physical segment size.
        integer, intent(in), optional :: min_res
            !! Minimum samples per patch; default is two.
        integer, intent(in), optional :: max_res
            !! Maximum samples per patch; default is the integer range limit.
        integer, allocatable :: res(:)
        real(rk) :: delta, length, segment
        integer :: c, i, j, maximum, minimum, nc, ncoord

        call require_patches(1,this%npatch,this%err)
        if (.not.this%err%ok) return
        if (.not. ieee_is_finite(target_size) .or. target_size <= 0.0_rk) then
            call this%err%set(&
                code       = 206,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
                message    = "Automatic target size must be positive and finite.",&
                location   = "create",&
                suggestion = "Pass target_size in the physical coordinate units.")
            return
        end if

        minimum = 2
        maximum = huge(0)
        if (present(min_res)) minimum = min_res
        if (present(max_res)) maximum = max_res
        if (minimum < 1 .or. maximum < minimum) then
            call this%err%set(&
                code       = 206,&
                severity   = 1,&
                category   = "forcad_multipatch_curve",&
                message    = "Automatic resolution bounds are inconsistent.",&
                location   = "create",&
                suggestion = "Require 1 <= min_res <= max_res.")
            return
        end if

        allocate(res(this%get_npatch()))
        do i = 1, this%get_npatch()
            nc = this%patches(i)%get_nc()
            ncoord = size(this%patches(i)%get_Xc(1))
            length = 0.0_rk
            do j = 2, nc
                segment = 0.0_rk
                do c = 1, ncoord
                    delta = this%patches(i)%get_Xc(j,c) - this%patches(i)%get_Xc(j-1,c)
                    segment = hypot(segment, delta)
                end do
                length = length + segment
            end do
            res(i) = sample_count(length,target_size,minimum,maximum)
        end do

        call create_curve_patches(this, res)
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
    !> Export patch control polygons in legacy VTK format.
    !!
    !! By default patch `i` is written as `<prefix>_Xc_<i>.vtk`. With
    !! `combined=.true.`, all cells are written to `<prefix>_Xc.vtk` and retain
    !! their one-based owner in the cell field `patch_id`.
    impure subroutine export_Xc_curve(this, prefix, encoding, combined)
        class(nurbs_multipatch_curve), intent(inout) :: this
        character(len=*), intent(in) :: prefix
        character(len=*), intent(in), optional :: encoding
        logical, intent(in), optional :: combined
        character(len=512) :: vtkfile
        integer :: i

        if (.not. this%err%ok) return

        if (present(combined)) then
            if (combined) then
                vtkfile = trim(prefix)//"_Xc.vtk"
                call export_combined_curve(&
                    this     = this,&
                    filename = vtkfile,&
                    source   = "Xc",&
                    encoding = encoding)
                return
            end if
        end if
        do i = 1, this%get_npatch()
            write(vtkfile, "(a,a,i0,a)") trim(prefix), "_Xc_", i, ".vtk"
            call this%patches(i)%export_Xc(vtkfile, encoding=encoding)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Export cached patch geometry in legacy VTK format.
    !!
    !! By default patch `i` is written as `<prefix>_Xg_<i>.vtk`. With
    !! `combined=.true.`, all cells are written to `<prefix>_Xg.vtk` and retain
    !! their one-based owner in the cell field `patch_id`.
    impure subroutine export_Xg_curve(this, prefix, encoding, combined)
        class(nurbs_multipatch_curve), intent(inout) :: this
        character(len=*), intent(in) :: prefix
        character(len=*), intent(in), optional :: encoding
        logical, intent(in), optional :: combined
        character(len=512) :: vtkfile
        integer :: i

        if (.not. this%err%ok) return

        if (present(combined)) then
            if (combined) then
                vtkfile = trim(prefix)//"_Xg.vtk"
                call export_combined_curve(&
                    this     = this,&
                    filename = vtkfile,&
                    source   = "Xg",&
                    encoding = encoding)
                return
            end if
        end if
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
    !!
    !! By default patch `i` is written as `<prefix>_Xth_<i>.vtk`. With
    !! `combined=.true.`, all line segments are written to
    !! `<prefix>_Xth.vtk` with cell field `patch_id`.
    impure subroutine export_Xth_in_Xg_curve(this, prefix, res, encoding, combined)
        class(nurbs_multipatch_curve), intent(inout) :: this
        character(len=*), intent(in) :: prefix
        integer, intent(in), optional :: res
        character(len=*), intent(in), optional :: encoding
        logical, intent(in), optional :: combined
        character(len=512) :: vtkfile
        integer :: i

        if (.not. this%err%ok) return

        if (present(combined)) then
            if (combined) then
                vtkfile = trim(prefix)//"_Xth.vtk"
                call export_combined_curve(&
                    this     = this,&
                    filename = vtkfile,&
                    source   = "Xth",&
                    res      = res,&
                    encoding = encoding)
                return
            end if
        end if
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
end module forcad_multipatch_curve
