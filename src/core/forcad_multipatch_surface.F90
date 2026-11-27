!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Topology, continuity constraints, and assembly data for surface patches.
!!
!! Patches remain complete [[nurbs_surface]] objects. Connections identify
!! oriented boundary curves. [[nurbs_multipatch_surface:cmp_dof_map]] merges
!! directly matching controls. [[nurbs_multipatch_surface:cmp_dof_constraint]]
!! couples nested clamped edges by knot insertion at order zero and supplies
!! sparse equations on matching trace spaces for normal derivatives through order
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
    use forcad_multipatch, only: multipatch_connection, multipatch_growth_capacity, check_new_patch, check_patch_transform
    use forcad_multipatch, only: constraint_rows, valid_multipatch, interface_sampling, sampling_bounds, sample_count
    use forcad_multipatch, only: resolve_patch_interface, record_connection, require_patches
    use forcad_multipatch, only: multipatch_dof_map, multipatch_dof_offsets, assemble_design_gradient, interface_quadrature
    use forcad_multipatch, only: multipatch_elements, element_layout, multipatch_element_ids, requested
    use forcad_nurbs_surface, only: nurbs_surface
    use forcad_utils, only: active_span_count, export_vtk_legacy
    use fordebug, only: debug

    implicit none

    private
    public :: nurbs_multipatch_surface

    !===============================================================================
    !> Ordered collection of connected NURBS surface patches.
    !!
    !! Valid side labels are `u_min`, `u_max`, `v_min`, and `v_max`, with the
    !! documented aliases accepted by `connect`. `reverse` reverses the single
    !! tangential boundary coordinate. `continuity=-1` records topology without
    !! constraints. Nonnegative continuity contributes sparse equations, while
    !! interpolatory normal endpoints additionally permit direct control sharing.
    !!
    !! Sampling is independent of the IGA knot discretization. The legacy
    !! `res1/res2` form broadcasts two counts. A rank-one `res` contains common
    !! directional counts, while rank-two `res(d,k)` specifies direction `d` of
    !! patch `k`. `create(target_size=h)` estimates counts in physical space.
    !! New resolution overloads synchronize connected tangential directions by
    !! default; `synchronize=.false.` retains exact independent requests.
    !!
    !! Diagnostic codes are: 201 invalid patch identifier; 202 invalid side;
    !! 203 invalid continuity metadata; 204 incompatible interface geometry or
    !! trace space; 205 invalid patch state; 206 invalid sampling request; and
    !! 207 invalid transformation data and 208 invalid or unrepresentable sensitivity data.
    type nurbs_multipatch_surface
        type(nurbs_surface), allocatable, private :: patches(:) !! Owned patch copies with spare append capacity.
        type(multipatch_connection), allocatable, private :: connections(:) !! Validated interface metadata.
        integer, private :: npatch = 0 !! Number of active patches.
        integer, private :: nconnection = 0 !! Number of active connections.

        type(debug) :: err !! Recoverable diagnostic state for the most recent failed operation.
    contains
        procedure :: add_patch => add_patch_surface !!> Append a validated patch copy and optionally return its identifier.
        procedure :: connect => connect_surface !!> Validate metadata and trace-space shape, then record an edge interface.
        procedure, private :: create_common => create_surface
        procedure, private :: create_uniform => create_surface_uniform
        procedure, private :: create_patches => create_surface_patches
        procedure, private :: create_automatic => create_surface_automatic
        generic :: create => create_common, create_uniform, create_patches, create_automatic
            !!> Sample patches with common, patch-specific, or physically scaled resolutions.
        procedure :: finalize => clear_surface !!> Release patches and connections and reset counts.
        procedure :: get_npatch => get_npatch_surface !!> Return the active patch count.
        procedure :: get_nconnection => get_nconnection_surface !!> Return the active connection count.
        procedure :: get_patch => get_patch_surface !!> Return a copy of one patch.
        procedure :: get_connection => get_connection_surface !!> Return a copy of one connection.
        procedure :: is_valid => is_valid_surface !!> Validate patch states, traces, orientations, and continuity.
        procedure :: is_rational => is_rational_surface !!> Report whether any patch has nonuniform rational weights.
        procedure :: cmp_dof_offsets => cmp_dof_offsets_surface !!> Return prefix offsets for unmerged patch controls.
        procedure :: assemble_sensitivities => assemble_sensitivities_surface !!> Assemble independent control and log-weight gradients.
        procedure :: cmp_dof_map => cmp_dof_map_surface !!> Compact directly shareable \(C^0\) controls.
        procedure :: interface_quadrature => interface_quadrature_surface !!> Paired integration on a connected boundary.
        procedure :: cmp_dof_constraint => cmp_dof_constraint_surface !!> Assemble derivative-continuity equations in CSR form.
        procedure :: cmp_elem => cmp_elem_surface !!> Concatenate patch element connectivity with optional shared numbering.
        procedure :: get_elem => cmp_elem_surface !!> Alias for [[nurbs_multipatch_surface:cmp_elem]].
        procedure :: cmp_elem_patch => cmp_elem_patch_surface !!> Return the owning patch identifier for every global element.
        procedure :: cmp_elem_local => cmp_elem_local_surface !!> Return the patch-local identifier for every global element.
        procedure :: export_Xc => export_Xc_surface !!> Export separate or combined patch control nets.
        procedure :: export_Xg => export_Xg_surface !!> Export separate or combined cached patch geometry.
        procedure :: export_Xth_in_Xg => export_Xth_in_Xg_surface !!> Export separate or combined embedded parameter grids.
        procedure :: transform => transform_surface !!> Apply one affine map to every patch and cache.
        procedure :: rotate_Xc => rotate_Xc_surface !!> Rotate all patch control points.
        procedure :: rotate_Xg => rotate_Xg_surface !!> Rotate all cached patch geometry points.
        procedure :: translate_Xc => translate_Xc_surface !!> Translate all patch control points.
        procedure :: translate_Xg => translate_Xg_surface !!> Translate all cached patch geometry points.
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

        call check_new_patch(2,patch,patch%err%ok,this%err,id)
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
    !> Apply one common affine transformation to all surface patches.
    !!
    !! Every control net and existing sampled-geometry cache receives
    !! \(\widetilde{\mathbf x}=\mathbf c+\mathbf A(\mathbf x-\mathbf c)+
    !! \mathbf b\). A common affine map preserves matching interface jets and
    !! homogeneous projective relations, so recorded \(C^n\) and \(G^n\)
    !! connections remain valid whenever the original multipatch was valid.
    pure subroutine transform_surface(this, linear, translation, center)
        class(nurbs_multipatch_surface), intent(inout) :: this
        real(rk), intent(in), contiguous :: linear(:,:)
            !! Square linear map shared by every patch.
        real(rk), intent(in), contiguous, optional :: translation(:)
            !! Shared translation vector.
        real(rk), intent(in), contiguous, optional :: center(:)
            !! Shared transformation center.
        integer :: i

        call check_patch_transform(this%patches,this%npatch,2,linear,translation,center,this%err)
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
    impure subroutine export_combined_surface(this, filename, source, res, encoding)
        class(nurbs_multipatch_surface), intent(inout) :: this
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: source
        integer, intent(in), optional :: res
        character(len=*), intent(in), optional :: encoding

        character(len=8), parameter :: cell_field_names(1) = ["patch_id"]
        real(rk), allocatable :: patch_points(:,:), points(:,:), cell_data(:,:)
        real(rk), allocatable :: knot1(:), knot2(:)
        integer, allocatable :: patch_elem(:,:), elemConn(:,:)
        integer :: a, b, cell, cell_offset, i, line, line_count, node, nodes_per_line, point_offset, t
        integer :: ncell, ncoord, nnode, npoint, patch_npoint, res_min, span_count1, span_count2, vtk_cell_type

        npoint = 0
        ncell = 0
        nnode = 0
        do i = 1, this%get_npatch()
            select case (source)
            case ("Xc")
                patch_npoint = product(this%patches(i)%get_nc())
                ncoord = size(this%patches(i)%get_Xc(1))
                patch_elem = this%patches(i)%get_elem_Xc_vis()
                if (size(patch_elem) == 0) then
                    patch_elem = this%patches(i)%cmp_elem_Xc_vis()
                end if
            case ("Xg")
                patch_npoint = product(this%patches(i)%get_ng())
                ncoord = size(this%patches(i)%get_Xg(1))
                patch_elem = this%patches(i)%get_elem_Xg_vis()
                if (size(patch_elem) == 0) then
                    patch_elem = this%patches(i)%cmp_elem_Xg_vis()
                end if
            case ("Xth")
                knot1 = this%patches(i)%get_knot(1)
                knot2 = this%patches(i)%get_knot(2)
                span_count1 = active_span_count(&
                    knot1,&
                    this%patches(i)%get_nc(1),&
                    this%patches(i)%get_degree(1))
                span_count2 = active_span_count(&
                    knot2,&
                    this%patches(i)%get_nc(2),&
                    this%patches(i)%get_degree(2))
                ncoord = size(this%patches(i)%get_Xc(1))
                if (.not. this%patches(i)%err%ok .or. span_count1 < 1 .or. span_count2 < 1 .or. &
                    ncoord < 2 .or. ncoord > 3) then
                    call this%err%set(&
                        code       = 205,&
                        severity   = 1,&
                        category   = "forcad_multipatch_surface",&
                        message    = "A patch visualization dataset is unavailable.",&
                        location   = "export_combined",&
                        suggestion = "Initialize and sample every patch before combined export.")
                    return
                end if
                a = span_count1
                b = span_count2
                do
                    t = mod(a,b)
                    if (t == 0) exit
                    a = b
                    b = t
                end do
                res_min = 10
                if (present(res)) res_min = max(2,res)
                nodes_per_line = (span_count1/b)*span_count2*max(1,res_min-1) + 1
                line_count = span_count1+span_count2+2
                npoint = npoint+line_count*nodes_per_line
                ncell = ncell+line_count*(nodes_per_line-1)
                nnode = 2
                cycle
            case default
                call this%err%set(&
                    code       = 205,&
                    severity   = 1,&
                    category   = "forcad_multipatch_surface",&
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
                    category   = "forcad_multipatch_surface",&
                    message    = "A patch visualization dataset is unavailable.",&
                    location   = "export_combined",&
                    suggestion = "Initialize and sample every patch before combined export.")
                return
            end if
            if (ncoord < 1 .or. ncoord > 3) then
                call this%err%set(&
                    code       = 205,&
                    severity   = 1,&
                    category   = "forcad_multipatch_surface",&
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
                    category   = "forcad_multipatch_surface",&
                    message    = "Patch visualization cell orders differ.",&
                    location   = "export_combined",&
                    suggestion = "Use one visualization cell order for all patches.")
                return
            end if
            ncell = ncell+size(patch_elem,1)
        end do

        if (source == "Xth") then
            nnode = 2
            vtk_cell_type = 3
        else
            vtk_cell_type = 9
        end if
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
                category   = "forcad_multipatch_surface",&
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
    end subroutine export_combined_surface
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
    !! Matching traces have equal tangential degree and affine-equivalent knot
    !! vectors. At order zero with an implicit weight scale, interpolating sides
    !! also accept clamped edge knots obtained by insertion on the coarser side.
    !! Call [[nurbs_multipatch_surface:is_valid]] to verify geometry and weights.
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
        real(rk), allocatable :: transition_jet_(:,:,:,:), projective_jet_(:,:,:)
        integer :: side_a_id, side_b_id, continuity_

        call resolve_patch_interface(this%patches,this%npatch,2,patch_a,side_a,patch_b,side_b,continuity, &
            side_a_id,side_b_id,continuity_,this%err)
        if (.not.this%err%ok) return
        if (present(transition_jet)) transition_jet_=spread(transition_jet,2,1)
        if (present(projective_jet)) projective_jet_=spread(projective_jet,2,1)
        call conn%set(patch_a, side_a_id, patch_b, side_b_id, continuity_, reverse, geometric=geometric, &
            reparameterization=reparameterization, transition_jet=transition_jet_, projective_jet=projective_jet_)
        call record_connection(this%patches,this%npatch,2,this%connections,this%nconnection,conn,this%err)
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
        offsets = multipatch_dof_offsets(this%patches,this%npatch)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Number scalar controls across periodic seams and directly shareable edges.
    !! The input index follows patch storage order. A directly shared trace requires
    !! an interpolatory endpoint on each side and proportional oriented trace weights.
    !! Call `is_valid` to check the geometry. General traces use
    !! [[nurbs_multipatch_surface:cmp_dof_constraint]].
    pure function cmp_dof_map_surface(this, active_only) result(map)
        class(nurbs_multipatch_surface), intent(in) :: this
        logical, optional, intent(in) :: active_only
            !! Number participating controls in first-active storage order and return zero elsewhere.
        integer, allocatable :: map(:)

        call multipatch_dof_map(this%patches,this%connections,this%nconnection,this%cmp_dof_offsets(),active_only,map)
    end function

    !> Quadrature on both trace meshes. See [[forcad_multipatch:interface_quadrature]].
    pure subroutine interface_quadrature_surface(this,connection,Xt_a,Xt_b,weight,ngauss,physical,err)
        class(nurbs_multipatch_surface), intent(in) :: this
        integer, intent(in) :: connection
        real(rk), allocatable, intent(out) :: Xt_a(:,:),Xt_b(:,:),weight(:)
        integer, optional, intent(in) :: ngauss(:)
        logical, optional, intent(in) :: physical
        type(debug), optional, intent(out) :: err
        call interface_quadrature(this%patches,this%npatch,this%connections,this%nconnection,2,this%err, &
            connection,Xt_a,Xt_b,weight,ngauss,physical,err)
    end subroutine

    !> Assemble patch-local partials in independent design variables. See [[forcad_multipatch:assemble_design_gradient]].
    pure subroutine assemble_sensitivities_surface(this,local_dXc,local_dlogWc,dXc,dlogWc,err,map)
        class(nurbs_multipatch_surface), intent(in) :: this
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
    !> Assemble one-based sparse continuity equations for connected edges.
    !! Call `is_valid` to verify the geometry. See [[forcad_multipatch:constraint_rows]]
    !! for the field, solver-numbering and connection-filter options.
    pure subroutine cmp_dof_constraint_surface(this, rowptr, col, val, geometric, field, err, active_only)
        class(nurbs_multipatch_surface), intent(in) :: this
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
    !> Concatenate active tensor-product elements from all surface patches.
    !! Rows follow patch insertion order and direction-1-fastest local element
    !! order. Different local basis sizes are zero-padded to the largest row.
    !! `shared=.true.` applies direct sharing from the compact map; the default
    !! retains unmerged patch-offset numbering. General unclamped traces and
    !! rational traces with nonconstant projective scale still require
    !! `cmp_dof_constraint`.
    pure function cmp_elem_surface(this, shared, active_only) result(elem)
        class(nurbs_multipatch_surface), intent(in) :: this
        logical, intent(in), optional :: shared
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
    !! Result ordering is identical to [[nurbs_multipatch_surface:cmp_elem]].
    pure function cmp_elem_patch_surface(this) result(patch_id)
        class(nurbs_multipatch_surface), intent(in) :: this
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
    !! Result ordering is identical to [[nurbs_multipatch_surface:cmp_elem]].
    pure function cmp_elem_local_surface(this) result(local_id)
        class(nurbs_multipatch_surface), intent(in) :: this
        integer, allocatable :: local_id(:)
        integer, allocatable :: rows(:)

        call element_layout(this%patches, this%get_npatch(), rows)
        call multipatch_element_ids(rows, local=.true., ids=local_id)
    end function

    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Validate the stored surface patches, their oriented interfaces and homogeneous continuity residuals.
    !! See [[valid_multipatch]] for the residual test used by geometry sensitivity assembly.
    pure function is_valid_surface(this, rtol, atol) result(ok)
        class(nurbs_multipatch_surface), intent(in) :: this
        real(rk), intent(in), optional :: rtol
            !! Relative trace and homogeneous-residual tolerance; default `sqrt(epsilon(rk))`.
        real(rk), intent(in), optional :: atol
            !! Absolute homogeneous-residual tolerance; default zero.
        logical :: ok

        ok=this%err%ok
        if (ok) ok=valid_multipatch(this%patches,this%npatch,this%connections,this%nconnection,2,rtol,atol)
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

        call require_patches(2,this%npatch,this%err)
        if (.not.this%err%ok) return
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
    !> Sample every surface patch with one common directional resolution.
    !!
    !! `res=[n_u,n_v]` is broadcast to all patches before connected tangential
    !! directions are synchronized. This rank-one overload complements the
    !! legacy `res1/res2` keywords without changing their behavior.
    pure subroutine create_surface_uniform(this, res, synchronize)
        class(nurbs_multipatch_surface), intent(inout) :: this
        integer, intent(in), contiguous :: res(:)
            !! Two positive directional sample counts.
        logical, intent(in), optional :: synchronize
            !! Synchronize connected tangential counts; default is true.
        if (.not.this%err%ok) return
        if (size(res)/=2) then
            call this%err%set(code=206,severity=1,category="forcad_multipatch_surface",location="create", &
                message="Invalid common surface resolution.",suggestion="Pass one count per parameter direction.")
            return
        end if
        call create_surface_patches(this,spread(res,2,this%npatch),synchronize)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Sample surface patches with independent directional resolutions.
    !!
    !! `res(d,k)` requests the sample count in parameter direction `d` of patch
    !! `k`. With synchronization enabled, requested counts are lower bounds:
    !! every connected edge receives the largest count requested in its
    !! transitive tangential-direction component. Reversal changes parameter
    !! order but not the count.
    pure subroutine create_surface_patches(this, res, synchronize)
        class(nurbs_multipatch_surface), intent(inout) :: this
        integer, intent(in), contiguous :: res(:,:)
            !! Positive resolution matrix with shape `[2,npatch]`.
        logical, intent(in), optional :: synchronize
            !! Synchronize connected tangential counts; default is true.
        integer, allocatable :: resolution(:,:)
        integer :: i

        call interface_sampling(2,this%npatch,this%connections,this%nconnection,res,synchronize,resolution,this%err)
        if (.not.this%err%ok) return
        do i = 1, this%get_npatch()
            call this%patches(i)%create(&
                res1 = resolution(1,i),&
                res2 = resolution(2,i))
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
    !> Choose surface-patch resolutions from a physical target edge size.
    !!
    !! In each parameter direction, the estimator takes the longest
    !! corresponding control-polygon line and applies
    !! \(n=\min(n_{\max},\max(n_{\min},1+\lceil L/h\rceil))\).
    !! Connected tangential directions are subsequently synchronized unless
    !! disabled. The control-net estimate is inexpensive and scale aware, but
    !! is not a certified surface-tessellation error bound.
    pure subroutine create_surface_automatic(this, target_size, min_res, max_res, synchronize)
        class(nurbs_multipatch_surface), intent(inout) :: this
        real(rk), intent(in) :: target_size
            !! Positive finite target physical edge size.
        integer, intent(in), contiguous, optional :: min_res(:)
            !! Two minimum directional counts; defaults to `[2,2]`.
        integer, intent(in), contiguous, optional :: max_res(:)
            !! Two maximum directional counts; defaults to the integer range limit.
        logical, intent(in), optional :: synchronize
            !! Synchronize connected tangential counts; default is true.
        integer, allocatable :: res(:,:)
        real(rk) :: delta, length(2), line_length, segment
        integer :: c, i, j, k, maximum(2), minimum(2), nc(2), ncoord

        call sampling_bounds(this%npatch,target_size,min_res,max_res,minimum,maximum,this%err)
        if (.not.this%err%ok) return
        allocate(res(2,this%get_npatch()))
        do k = 1, this%get_npatch()
            nc = this%patches(k)%get_nc()
            ncoord = size(this%patches(k)%get_Xc(1))
            length = 0.0_rk
            do j = 1, nc(2)
                line_length = 0.0_rk
                do i = 2, nc(1)
                    segment = 0.0_rk
                    do c = 1, ncoord
                        delta = this%patches(k)%get_Xc(i + nc(1)*(j-1),c) - &
                            this%patches(k)%get_Xc(i-1 + nc(1)*(j-1),c)
                        segment = hypot(segment, delta)
                    end do
                    line_length = line_length + segment
                end do
                length(1) = max(length(1),line_length)
            end do
            do i = 1, nc(1)
                line_length = 0.0_rk
                do j = 2, nc(2)
                    segment = 0.0_rk
                    do c = 1, ncoord
                        delta = this%patches(k)%get_Xc(i + nc(1)*(j-1),c) - &
                            this%patches(k)%get_Xc(i + nc(1)*(j-2),c)
                        segment = hypot(segment, delta)
                    end do
                    line_length = line_length + segment
                end do
                length(2) = max(length(2),line_length)
            end do
            res(:,k) = sample_count(length,target_size,minimum,maximum)
        end do

        call create_surface_patches(this, res, synchronize)
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
    !> Export patch control nets in legacy VTK format.
    !!
    !! By default patch `i` is written as `<prefix>_Xc_<i>.vtk`. With
    !! `combined=.true.`, all cells are written to `<prefix>_Xc.vtk` and retain
    !! their one-based owner in the cell field `patch_id`.
    impure subroutine export_Xc_surface(this, prefix, encoding, combined)
        class(nurbs_multipatch_surface), intent(inout) :: this
        character(len=*), intent(in) :: prefix
        character(len=*), intent(in), optional :: encoding
        logical, intent(in), optional :: combined
        character(len=512) :: vtkfile
        integer :: i

        if (.not. this%err%ok) return

        if (present(combined)) then
            if (combined) then
                vtkfile = trim(prefix)//"_Xc.vtk"
                call export_combined_surface(&
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
    impure subroutine export_Xg_surface(this, prefix, encoding, combined)
        class(nurbs_multipatch_surface), intent(inout) :: this
        character(len=*), intent(in) :: prefix
        character(len=*), intent(in), optional :: encoding
        logical, intent(in), optional :: combined
        character(len=512) :: vtkfile
        integer :: i

        if (.not. this%err%ok) return

        if (present(combined)) then
            if (combined) then
                vtkfile = trim(prefix)//"_Xg.vtk"
                call export_combined_surface(&
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
    !> Export parameter lines mapped into every surface patch.
    !!
    !! By default patch `i` is written as `<prefix>_Xth_<i>.vtk`. With
    !! `combined=.true.`, all line segments are written to
    !! `<prefix>_Xth.vtk` with cell field `patch_id`.
    impure subroutine export_Xth_in_Xg_surface(this, prefix, res, encoding, combined)
        class(nurbs_multipatch_surface), intent(inout) :: this
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
                call export_combined_surface(&
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
end module forcad_multipatch_surface
