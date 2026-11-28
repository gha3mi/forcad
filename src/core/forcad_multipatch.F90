! author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Shared topology and continuity metadata for NURBS multipatches.
!!
!! A [[multipatch_connection]] identifies two oriented patch boundaries and the
!! continuity requested across their interface. For a patch parameter \(\xi_P\)
!! with active interval \([a_P,b_P]\), define the normalized local increment
!! from its selected boundary by
!!
!! \[
!! \eta_P=\frac{\xi_P-\xi_{\Gamma,P}}{b_P-a_P},\qquad
!! \xi_{\Gamma,P}\in\{a_P,b_P\}.
!! \]
!!
!! Thus \(\eta_P=0\) at the interface and its positive direction follows the
!! stored parameter axis; at a maximum side it is not an inward coordinate.
!! ForCAD's parametric \(C^n\) convention equates derivatives with respect to
!! these normalized coordinates through order \(n\), after applying interface
!! orientation. Geometric curve connections compose the patch-B parameter with
!! a scalar transition map \(\phi\). The compact transition jet is
!!
!! \[
!! \boldsymbol\phi_\Gamma=
!! [\phi'(0),\phi''(0),\ldots,\phi^{(n)}(0)].
!! \]
!!
!! For the compact separable surface and volume representation, corresponding
!! mapped tangential coordinates \(\boldsymbol\tau\) remain fixed and the
!! order-\(r\) interface equation is
!!
!! \[
!! \left.\frac{\partial^r f_A}{\partial\eta_A^r}\right|_{\Gamma}
!! =
!! \left.\frac{\mathrm d^r}{\mathrm d\eta_A^r}
!! f_B\!\left(\boldsymbol\tau,\phi(\eta_A)\right)\right|_{\eta_A=0}
!! =
!! \sum_{k=0}^{r}
!! B_{r,k}\!\left(\phi'(0),\ldots,\phi^{(r-k+1)}(0)\right)
!! \left.\frac{\partial^k f_B}{\partial\eta_B^k}\right|_{\Gamma},
!! \]
!!
!! where \(B_{r,k}\) is a partial exponential Bell polynomial. Parametric
!! \(C^n\) uses the signed affine transition
!! \(\phi(\eta_A)=\sigma\eta_A\), so only the coefficient
!! \(B_{r,r}=\sigma^r\) remains. Here \(\sigma=+1\) when one connected side is
!! a minimum and the other a maximum, and \(\sigma=-1\) when both are the same
!! side type. Geometric \(G^n\) uses the complete supplied transition jet. Order
!! zero always enforces trace equality.
!!
!! For curves, this univariate transition is the complete local
!! parameter-transition model. Surface and volume connections additionally
!! accept a full interface jet of a multivariate transition
!!
!! \[
!! \Phi(\boldsymbol\tau,\eta_A)
!! =T\boldsymbol\tau+
!! \sum_{r=1}^{n}\frac{\boldsymbol a_r(\boldsymbol\tau)}{r!}\eta_A^r
!! +O(\eta_A^{n+1}),
!! \]
!!
!! where all mapped tangential coordinates and the mapped normal coordinate
!! have independently specified normal derivatives. Each component of
!! \(\boldsymbol a_r\) is stored as a univariate or tensor-product Bernstein
!! polynomial over the normalized patch-A interface. This permits tangentially
!! varying normal jets and normal-dependent tangential coordinates. The
!! resulting mixed derivatives are assembled with the multivariate
!! Faa di Bruno formula.
!!
!! Rational connections act on homogeneous lifts. They may supply the
!! Bernstein jet of a nonvanishing projective factor
!! \(\lambda(\boldsymbol\tau,\eta_A)\) in
!!
!! \[
!! \boldsymbol H_A
!! =\lambda\,(\boldsymbol H_B\circ\Phi).
!! \]
!!
!! Omitting the full transition and projective jets selects the compact,
!! faster separable model with one constant projective scale. The general
!! surface and volume paths enforce the resulting piecewise-polynomial
!! identities by degree-complete collocation on every interface knot span.
!!
!! Tangential orientation is represented by reversal, swap, and per-axis
!! flips. Rank-specific `connect` methods validate metadata and compatible trace
!! spaces; rank-specific `is_valid` methods additionally check physical traces,
!! rational projective weights, and assembled derivative residuals. A
!! continuity value of `-1` records topology without merging degrees of freedom;
!! values `0:n` request sparse continuity equations through that derivative
!! order. Direct integer-map sharing is possible only for one-hot endpoint
!! bases. Rational traces additionally require corresponding denominator
!! coefficients to differ by one constant nonzero scale, which makes their
!! rational trace basis functions identical. General unclamped traces and
!! tangentially varying projective factors remain sparse linear constraints.
!! Every floating-point coefficient that is structurally nonzero is retained
!! in those constraints, irrespective of magnitude. In particular, a small
!! projective factor is part of the homogeneous identity and is not a numerical
!! zero. Solver equilibration must treat row and unknown scaling separately:
!! for \(\mathbf q=D_q\widehat{\mathbf q}\), a nonsingular diagonal row scaling
!! \(D_r\) gives
!! \(D_r A D_q\widehat{\mathbf q}=\mathbf0\), followed by recovery
!! \(\mathbf q=D_q\widehat{\mathbf q}\).
!!
!! **Geometric-continuity reference:** P. Kiciak, "Conditions for
!! geometric continuity between polynomial and rational surface patches,"
!! *Computer Aided Geometric Design* 13 (1996), 709--741.
!! [doi:10.1016/0167-8396(96)00006-4](https://doi.org/10.1016/0167-8396(96)00006-4).
!!
!! **Multivariate-chain-rule reference:** G. M. Constantine and
!! T. H. Savits, "A multivariate Faa di Bruno formula with applications,"
!! *Transactions of the American Mathematical Society* 348 (1996),
!! 503--520.
!! [doi:10.1090/S0002-9947-96-01501-2](https://doi.org/10.1090/S0002-9947-96-01501-2).
module forcad_multipatch

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fordebug, only: debug
    use forcad_utils, only: bernstein_span, equilibrate_field_row, eval_bernstein_endpoint, &
        prepare_bernstein_pair, eval_bernstein_pair
    use forcad_kinds, only: rk
    use forcad_utils, only: nested_edge_transfer, axis_transfer, tensor_transfer
    use forcad_nurbs_curve, only: nurbs_curve
    use forcad_nurbs_surface, only: nurbs_surface
    use forcad_nurbs_volume, only: nurbs_volume
    use forcad_utils, only: boundary_axes, disjoint_set_union, disjoint_set_map, interpolatory_boundary_layer
    use forcad_utils, only: basis_bspline_der_all_active, boundary_layer_index, structural_nonzero
    use forcad_utils, only: scaled_difference, affine_knot_residual, interval_ratio, two_sum, add_scaled_term

    implicit none

    private
    public multipatch_connection
#ifndef FORCAD_PUBLIC_DOCS
    public check_constraint_residuals, store_constraint_points, finish_constraint_block, finish_constraint_matrix
    public valid_multipatch, record_connection, resolve_patch_interface, interface_sampling
    public check_new_patch, check_patch_transform, require_patches, sampling_bounds, sample_count
    public requested, constraint_rows, storage_weights, interface_quadrature
    public multipatch_compatible_trace_space, multipatch_projective_weight_scale, &
        multipatch_growth_capacity, multipatch_chain_rule_coefficients, multipatch_composition_coefficients, &
        multipatch_valid_reparameterization, multipatch_element_ids, multipatch_element_block, &
        append_separable_constraint, include_constraint

    public connection_weight_scale, trace_weight_stencil, weighted_trace_sum, scale_projective_jet, continuity_tolerances
    public add_separable_shape, seed_patch_dofs, interpolated_trace_ids, merge_trace_dofs
    public multipatch_dof_map, multipatch_dof_offsets, volume_interface_axes, assemble_design_gradient
    public append_bernstein_constraint, add_bernstein_shape
    public element_layout, multipatch_elements
    public resolve_connection_sides, append_connection, valid_connection_metadata
#endif

    !> Bernstein coefficients of the supplied transition and projective jets on one interface cell.
    type interface_jet
        private
        real(rk), allocatable :: transition(:,:,:), projective(:,:,:)
    end type

    !> Overlay breaks, scaled native A cell widths, and paired one-dimensional Gauss nodes.
    type paired_axis
        integer :: order,ncell
        real(rk), allocatable :: break(:,:),offset(:,:),width(:),node(:,:),w(:)
        integer, allocatable :: power(:)
    end type

    !> One fine boundary control per row; coarse columns use patch-local storage IDs.
    type trace_prolongation
        logical :: coarse_a
        integer, allocatable :: ptr(:),coarse(:),fine(:)
        real(rk), allocatable :: val(:)
    end type

    !===============================================================================
    !> Oriented interface and continuity metadata for two patches.
    !!
    !! Patch identifiers are one-based. Side identifiers use the rank-specific
    !! ordering: curve `left/right = 1/2`; surface
    !! `u_min/u_max/v_min/v_max = 1:4`; and volume adds
    !! `w_min/w_max = 5/6`.
    !!
    !! `reverse` reverses the single tangential direction. Curve endpoints have
    !! no tangential direction, so `reverse` has no effect for curve patches and
    !! is retained only for rank-consistent metadata. For volume faces,
    !! `swap` exchanges the two tangential directions and `flip(1:2)` reverses
    !! either direction after that mapping. For geometric continuity,
    !! `reparameterization(k)` is \(\phi^{(k)}(0)\), the \(k\)-th derivative of
    !! the normalized patch-B normal coordinate with respect to patch A's
    !! normalized local increment at the interface. Surface and volume
    !! connections may instead store all mapped-coordinate normal derivatives
    !! as Bernstein fields in `transition_jet`; `projective_jet` stores the
    !! Bernstein fields for the homogeneous scale and its normal derivatives.
    !!
    !! @warning [[multipatch_connection:set]] is a low-level metadata setter. It
    !! does not have access to patches and cannot validate topology or spline
    !! compatibility. Prefer the rank-specific `connect` binding.
    !! @endwarning
    type multipatch_connection
        type(trace_prolongation), allocatable, private :: prolongation !! Exact transfer from coarse to fine trace coefficients.
        integer, private :: patch_a = 0 !! One-based first patch identifier.
        integer, private :: side_a = 0 !! First patch boundary identifier.
        integer, private :: patch_b = 0 !! One-based second patch identifier.
        integer, private :: side_b = 0 !! Second patch boundary identifier.
        integer, private :: continuity = 0 !! Requested order, or `-1` for no continuity constraint.
        logical, private :: geometric = .false. !! False for normalized \(C^n\), true for geometric continuity.
        real(rk), allocatable, private :: reparameterization(:) !! Normal transition derivatives 1 through \(n\).
        real(rk), allocatable, private :: transition_jet(:,:,:,:) !! Bernstein coefficients `[q1+1,q2+1,ncoord,n]`.
        real(rk), allocatable, private :: projective_jet(:,:,:) !! Bernstein coefficients `[q1+1,q2+1,n+1]`.
        logical, private :: reverse = .false. !! Reverse the boundary orientation.
        logical, private :: swap = .false. !! Swap two tangential directions on a volume face.
        logical, private :: flip(2) = [.false., .false.] !! Flip mapped tangential directions.
    contains
        procedure :: set => set_connection !!> Replace all connection metadata without patch-level validation.
        procedure :: get_patch_a => get_connection_patch_a !!> Return the first patch identifier.
        procedure :: get_side_a => get_connection_side_a !!> Return the first boundary identifier.
        procedure :: get_patch_b => get_connection_patch_b !!> Return the second patch identifier.
        procedure :: get_side_b => get_connection_side_b !!> Return the second boundary identifier.
        procedure :: get_continuity => get_connection_continuity !!> Return the requested continuity order.
        procedure :: get_reparameterization => get_connection_reparameterization !!> Return a copy of the transition jet.
        procedure :: get_transition_jet => get_connection_transition_jet !!> Return the multivariate transition jet.
        procedure :: get_projective_jet => get_connection_projective_jet !!> Return the projective-factor jet.
        procedure :: has_transition_jet => connection_has_transition_jet !!> Report a full multivariate transition jet.
        procedure :: has_projective_jet => connection_has_projective_jet !!> Report a variable projective-factor jet.
        procedure :: is_geometric => get_connection_geometric !!> Report whether geometric continuity is selected.
        procedure :: is_reversed => get_connection_reverse !!> Report tangential reversal.
        procedure :: is_swapped => get_connection_swap !!> Report tangential-axis exchange.
        procedure :: is_flipped => get_connection_flip !!> Report a selected tangential-axis flip.
    end type
    !===============================================================================

contains

    !> Polynomial degree sufficient for the order-r residual along one tangential axis.
    pure elemental integer function jet_row_degree(conn,p_a,p_b,order,axis) result(p)
        type(multipatch_connection), intent(in) :: conn
        integer, intent(in) :: p_a,p_b,order,axis
        p = p_b
        if (allocated(conn%transition_jet)) p = p + order*(size(conn%transition_jet,axis)-1)
        if (allocated(conn%projective_jet)) p = p + size(conn%projective_jet,axis)-1
        p = max(p_a,p)
    end function

    !> Sum the tensor collocation rows and their maximum active support over all derivative orders.
    pure subroutine add_bernstein_shape(conn,degree_a,degree_b,ncells,nrow,nnz,value_rows)
        type(multipatch_connection), intent(in) :: conn
        integer, intent(in) :: degree_a(3),degree_b(3),ncells
        integer, intent(inout) :: nrow,nnz
        integer, intent(out), optional :: value_rows
        integer :: r,rows
        rows = 0
        do r=0,conn%continuity
            rows = rows + ncells*product(1+jet_row_degree(conn,degree_a(:2),degree_b(:2),r,[1,2]))
            if (r==0 .and. present(value_rows)) value_rows=rows
        end do
        nrow = nrow + rows
        nnz = nnz + rows*(product(1+degree_a) + product(1+degree_b))
    end subroutine

    !> Corresponding parameter axes on both faces, in the order of side A.
    pure subroutine volume_interface_axes(conn, tangent_a, mapped_b, reverse_dir)
        type(multipatch_connection), intent(in) :: conn
        integer, intent(out) :: tangent_a(3), mapped_b(3)
        logical, optional, intent(out) :: reverse_dir(2)
        tangent_a=boundary_axes(conn%side_a)
        mapped_b=boundary_axes(conn%side_b,conn%swap,conn%flip,conn%reverse)
        if (present(reverse_dir)) reverse_dir=mapped_b(:2)<0
        mapped_b=abs(mapped_b)
    end subroutine

    !> Pair Gauss nodes of the two whole traces using the stored affine orientation.
    !! Rows of Xt_a and Xt_b are native patch parameters. Weights integrate side A's
    !! tangential parameters, or physical arc length/area on A when physical is true.
    pure subroutine interface_quadrature(patches,npatch,connections,nconnection,ndim,object_err,connection, &
            Xt_a,Xt_b,weight,ngauss,physical,err)
        use forcad_utils, only: active_knots, scaled_measure
        class(*), optional, intent(in) :: patches(:)
        type(multipatch_connection), optional, intent(in) :: connections(:)
        integer, intent(in) :: npatch,nconnection,ndim,connection
        type(debug), intent(in) :: object_err
        real(rk), allocatable, intent(out) :: Xt_a(:,:),Xt_b(:,:),weight(:)
        integer, optional, intent(in) :: ngauss(:)
        logical, optional, intent(in) :: physical
        type(debug), optional, intent(out) :: err
        type(debug) :: status
        integer, allocatable :: weight_power(:)
        type(paired_axis) :: rule(2)
        real(rk), allocatable :: knot(:),ka(:),kb(:)
        real(rk) :: endpoint(2),value
        integer :: axes(3,2),nc(3,2),p(3,2),side(2),patch(2),d,j,s,q,nq,i,power
        logical :: reflect(2),ok

        status=object_err
        if (.not.status%ok) goto 900
        if (connection<1.or.connection>nconnection.or..not.present(connections).or..not.present(patches)) goto 800
        if (ndim/=2.and.ndim/=3) goto 800
        if (present(ngauss)) then
            if (size(ngauss)/=ndim-1.or.any(ngauss<1)) goto 800
        end if
        associate(conn=>connections(connection))
            if (.not.valid_connection_metadata(conn,npatch,ndim)) goto 800
            call volume_interface_axes(conn,axes(:,1),axes(:,2),reflect)
            side=[conn%side_a,conn%side_b]
            patch=[conn%patch_a,conn%patch_b]
        end associate
        do s=1,2
            if (.not.patch_is_valid(patches(patch(s)))) goto 800
            call patch_basis(patches(patch(s)),nc(:,s),p(:,s))
            j=axes(3,s)
            knot=patch_axis_knots(patches(patch(s)),j)
            endpoint(s)=knot(merge(p(j,s)+1,nc(j,s)+1,mod(side(s),2)==1))
        end do
        nq=1
        do d=1,ndim-1
            j=axes(d,1)
            ka=active_knots(patch_axis_knots(patches(patch(1)),j),nc(j,1),p(j,1))
            j=axes(d,2)
            kb=active_knots(patch_axis_knots(patches(patch(2)),j),nc(j,2),p(j,2))
            if (reflect(d)) kb=-kb(size(kb):1:-1)
            rule(d)%order=max(p(axes(d,1),1),p(j,2))+1
            if (present(ngauss)) rule(d)%order=ngauss(d)
            call merge_trace_breaks(ka,kb,rule(d))
            if (rule(d)%ncell>huge(nq)/rule(d)%order/nq/ndim) goto 810
            nq=nq*rule(d)%order*rule(d)%ncell
        end do
        do d=1,ndim-1
            call sample_trace_axis(rule(d),ok)
            if (.not.ok) goto 810
            if (reflect(d)) rule(d)%node(:,2)=-rule(d)%node(:,2)
        end do
        allocate(Xt_a(nq,ndim),Xt_b(nq,ndim),weight(nq),weight_power(nq))
        Xt_a(:,axes(3,1))=endpoint(1)
        Xt_b(:,axes(3,2))=endpoint(2)
        do q=1,nq
            value=1.0_rk
            power=0
            i=q-1
            do d=1,ndim-1
                j=1+mod(i,size(rule(d)%w))
                i=i/size(rule(d)%w)
                Xt_a(q,axes(d,1))=rule(d)%node(j,1)
                Xt_b(q,axes(d,2))=rule(d)%node(j,2)
                value=value*rule(d)%w(j)
                power=power+rule(d)%power(1+(j-1)/rule(d)%order)
            end do
            weight(q)=fraction(value)
            weight_power(q)=power+exponent(value)
        end do
        if (requested(physical)) then
            call weight_physical_trace(patches(patch(1)),side(1),Xt_a,weight,weight_power)
        else
            weight=scaled_measure(weight,weight_power)
        end if
        if (any(weight<=0.0_rk)) then
            call status%set(code=208,severity=1,category="forcad_multipatch",location="interface_quadrature", &
                message="A quadrature weight is singular or outside the working real range.", &
                suggestion="Check the tangents on side A and the requested integration measure.")
            goto 900
        end if
        if (present(err)) err=status
        return
800     call status%set(code=205,severity=1,category="forcad_multipatch",location="interface_quadrature", &
            message="Invalid interface, patch, or quadrature order.",suggestion="Use a stored interface and positive tangent orders.")
        goto 900
810     call status%set(code=208,severity=1,category="forcad_multipatch",location="interface_quadrature", &
            message="The paired Gauss rule exceeds the available index or parameter precision.", &
            suggestion="Use representable knot intervals and Gauss orders.")
900     if (allocated(weight)) deallocate(weight,Xt_a,Xt_b)
        allocate(weight(0),Xt_a(0,ndim),Xt_b(0,ndim))
        if (present(err)) err=status
    end subroutine

    !> Bound derivative workspace while evaluating the actual rational geometry on side A.
    pure subroutine weight_physical_trace(patch,side,points,weight,power)
        use forcad_utils, only: weighted_tangent_measure
        use forcad_nurbs_surface, only: surface_geometry_derivative
        use forcad_nurbs_volume, only: volume_geometry_derivative
        class(*), intent(in) :: patch
        integer, intent(in) :: side,power(:)
        real(rk), intent(in) :: points(:,:)
        real(rk), intent(inout) :: weight(:)
        real(rk), allocatable :: tangent(:,:,:),derivative(:,:)
        integer :: axes(3),order(size(points,2)),trace(size(points,2)),nt,ncoord,d,i,first,last,n
        integer, parameter :: batch=256
        integer :: measure_power(batch)
        integer, allocatable :: tangent_power(:)

        axes=boundary_axes(side)
        nt=size(points,2)-1
        ncoord=patch_ncoord(patch)
        trace=0
        trace(axes(3))=merge(1,-1,mod(side,2)==1)
        allocate(tangent(ncoord,nt,min(batch,size(weight))))
        do first=1,size(weight),batch
            last=min(first+batch-1,size(weight))
            n=last-first+1
            measure_power(:n)=power(first:last)
            do d=1,nt
                order=0
                order(axes(d))=1
                select type (patch)
                type is (nurbs_surface)
                    call surface_geometry_derivative(patch,points(first:last,:),order,derivative,trace,tangent_power)
                type is (nurbs_volume)
                    call volume_geometry_derivative(patch,points(first:last,:),order,derivative,trace,tangent_power)
                end select
                if (any(shape(derivative)/=[n,ncoord])) then
                    weight(first)=0.0_rk
                    return
                end if
                tangent(:,d,:n)=transpose(derivative)
                measure_power(:n)=measure_power(:n)+tangent_power
            end do
            do i=first,last
                weight(i)=weighted_tangent_measure(tangent(:,:,i-first+1),weight(i),measure_power(i-first+1))
            end do
        end do
    end subroutine

    !> Merge increasing break lists in their native coordinates. The determinant compares
    !! (a-a0)/(a1-a0) with (b-b0)/(b1-b0) before either normalized coordinate is rounded.
    pure subroutine merge_trace_breaks(a,b,rule)
        real(rk), intent(in) :: a(:),b(:)
        type(paired_axis), intent(inout) :: rule
        real(rk) :: start(2),finish(2),domain(2),point(2),det,delta
        integer :: ia,ib,last,source,n,previous,shift(2),power

        start=[a(1),b(1)]
        finish=[a(size(a)),b(size(b))]
        call scaled_difference(start,finish,domain,shift)
        shift=shift+exponent(domain)
        domain=fraction(domain)
        n=size(a)+size(b)-2
        allocate(rule%break(2,n),rule%width(n-1),rule%power(n-1))
        allocate(rule%offset(2,n),source=0.0_rk)
        rule%break(:,1)=start
        last=1
        previous=1
        ia=2
        ib=2
        do while (ia<=size(a))
            point=[a(ia),b(ib)]
            call affine_knot_residual(start,finish,point,det,power,signed_result=.true.)
            source=merge(1,2,det>=0.0_rk)
            if (det>=0.0_rk) ia=ia+1
            if (det<=0.0_rk) ib=ib+1
            if (det/=0.0_rk) then
                block
                    integer :: other
                    real(rk) :: x
                    other=3-source
                    x=point(source)
                    if (start(source)==start(other).and.finish(source)==finish(other)) then
                        point(other)=x
                    else
                        point(other)=affine_value_at(start(source),finish(source),x,start(other),finish(other))
                        ! Preserve the fractional ULP at a translated native origin.
                        call affine_knot_residual(start,finish,point,delta,power,signed_result=.true.)
                        delta=scale(merge(-delta,delta,source==1)/domain(source),power-shift(source))
                        x=point(other)
                        call two_sum(x,delta,point(other),rule%offset(other,last+1))
                    end if
                end block
            end if
            ! A same-list width is a native difference. A mixed width equals
            ! the determinant of its two defining breaks divided by B's domain width.
            if (source==previous) then
                call scaled_difference(rule%break(source,last),point(source),delta,power)
                if (source==2) then
                    power=power+exponent(delta)
                    delta=fraction(delta)*domain(1)/domain(2)
                    power=power+shift(1)-shift(2)
                end if
            else
                block
                    real(rk) :: ends(2)
                    ends=point
                    ends(previous)=rule%break(previous,last)
                    call affine_knot_residual(start,finish,ends,delta,power)
                end block
                delta=delta/domain(2)
                power=power-shift(2)
            end if
            rule%power(last)=power+exponent(delta)
            rule%width(last)=fraction(delta)
            last=last+1
            rule%break(:,last)=point
            previous=source
        end do
        rule%ncell=last-1
    end subroutine

    !> Sample each open common interval. A rounded endpoint or repeated Gauss node
    !! cannot be evaluated unambiguously through the native patch parameter API.
    pure subroutine sample_trace_axis(rule,ok)
        use forcad_utils, only: gauss_leg
        type(paired_axis), intent(inout) :: rule
        logical, intent(out) :: ok
        real(rk), allocatable :: z(:),w(:)
        real(rk) :: lo,hi,left,right
        integer :: i,s,first,last
        ok=.false.
        call gauss_leg([0.0_rk,1.0_rk],rule%order-1,z,w)
        if (.not.all(ieee_is_finite(z)).or..not.all(w>0.0_rk)) return
        allocate(rule%node(rule%order*rule%ncell,2),rule%w(rule%order*rule%ncell))
        do i=1,rule%ncell
            first=1+(i-1)*rule%order
            last=i*rule%order
            rule%w(first:last)=w*rule%width(i)
            do s=1,2
                lo=rule%break(s,i)
                hi=rule%break(s,i+1)
                left=rule%offset(s,i)
                right=rule%offset(s,i+1)
                associate(x=>rule%node(first:last,s))
                    x=(1.0_rk-z)*left+z*right
                    if (lo<0.0_rk.and.hi>0.0_rk) then
                        x=blossom_at(0.0_rk,1.0_rk,z,1.0_rk-z,z,lo,hi)+x
                    else
                        where (z<=0.5_rk)
                            x=lo+(z*(hi-lo)+x)
                        else where
                            x=hi+(x-(1.0_rk-z)*(hi-lo))
                        end where
                    end if
                    if (x(1)<lo.or.(x(1)==lo.and.left>=0.0_rk)) return
                    if (x(size(x))>hi.or.(x(size(x))==hi.and.right<=0.0_rk)) return
                    if (any(x(2:)<=x(:size(x)-1))) return
                end associate
            end do
        end do
        ok=.true.
    end subroutine

    pure subroutine require_patches(ndim,npatch,err)
        integer, intent(in) :: ndim,npatch
        type(debug), intent(inout) :: err
        if (err%ok .and. npatch<1) call err%set(code=205,severity=1,category=patch_category(ndim),location="create", &
            message="No patches are available.",suggestion="Add at least one patch before calling create().")
    end subroutine

    !> Directional limits for automatic sampling of connected edges and faces.
    pure subroutine sampling_bounds(npatch,h,min_res,max_res,minimum,maximum,err)
        integer, intent(in) :: npatch
        real(rk), intent(in) :: h
        integer, optional, intent(in) :: min_res(:),max_res(:)
        integer, intent(out) :: minimum(:),maximum(:)
        type(debug), intent(inout) :: err
        integer :: ndim
        logical :: valid
        ndim=size(minimum)
        call require_patches(ndim,npatch,err)
        if (.not.err%ok) return
        valid=ieee_is_finite(h).and.h>0
        minimum=2
        maximum=huge(0)
        if (present(min_res)) then
            valid=valid.and.size(min_res)==ndim
            if (valid) minimum=min_res
        end if
        if (present(max_res)) then
            valid=valid.and.size(max_res)==ndim
            if (valid) maximum=max_res
        end if
        if (.not.valid .or. any(minimum<1) .or. any(maximum<minimum)) then
            call err%set(code=206,severity=1,category=patch_category(ndim),location="create", &
                message="Invalid target size or directional bounds.", &
                suggestion="Use a positive finite target size and one pair of ordered positive bounds per direction.")
        end if
    end subroutine

    !> Clip the sample count before converting the length-to-size ratio to an integer.
    pure elemental integer function sample_count(length,h,minimum,maximum) result(n)
        real(rk), intent(in) :: length,h
        integer, intent(in) :: minimum,maximum
        real(rk) :: limit
        limit=real(maximum-1,rk)
        n=maximum
        if (h<1.0_rk) then
            if (length>=h*limit) return
        else
            if (length/h>=limit) return
        end if
        n=max(minimum,ceiling(length/h)+1)
    end function

    !> Connected tangential axes share the largest requested sampling count.
    pure subroutine interface_sampling(ndim,npatch,connections,nconnection,res,synchronize,resolution,err)
        integer, intent(in) :: ndim,npatch,nconnection,res(:,:)
        type(multipatch_connection), optional, intent(in) :: connections(:)
        logical, optional, intent(in) :: synchronize
        integer, allocatable, intent(out) :: resolution(:,:)
        type(debug), intent(inout) :: err
        integer, allocatable :: parent(:),component(:),counts(:)
        integer :: i,d,a(3),b(3),group

        call require_patches(ndim,npatch,err)
        if (.not.err%ok) return
        if (any(shape(res)/=[ndim,npatch]) .or. any(res<1)) then
            call err%set(code=206,severity=1,category=patch_category(ndim),location="create", &
                message="Invalid patch resolutions.",suggestion="Use positive counts with shape [ndim,npatch].")
        end if
        if (.not.err%ok) return
        resolution=res
        if (nconnection==0) return
        if (present(synchronize)) then
            if (.not.synchronize) return
        end if
        parent=[(i,i=1,size(res))]
        do i=1,nconnection
            associate(conn=>connections(i))
                call volume_interface_axes(conn,a,b)
                do d=1,ndim-1
                    call disjoint_set_union(parent,ndim*(conn%patch_a-1)+a(d),ndim*(conn%patch_b-1)+b(d))
                end do
            end associate
        end do
        call disjoint_set_map(parent,component)
        allocate(counts(maxval(component)),source=0)
        do i=1,npatch
            do d=1,ndim
                group=component(ndim*(i-1)+d)
                counts(group)=max(counts(group),res(d,i))
            end do
        end do
        do i=1,npatch
            resolution(:,i)=counts(component(ndim*(i-1)+1:ndim*i))
        end do
    end subroutine

    !> Tensor metadata in storage axes. Missing axes have degree zero and one constant basis.
    pure subroutine patch_basis(patch,nc,degree,stride)
        class(*), intent(in) :: patch
        integer, intent(out) :: nc(3),degree(3)
        integer, optional, intent(out) :: stride(3)
        nc=1
        degree=0
        select type (patch)
        type is (nurbs_curve)
            nc(1)=patch%get_nc()
            degree(1)=patch%get_degree()
            if (present(stride)) stride=[1,0,0]
        type is (nurbs_surface)
            nc(:2)=patch%get_nc()
            degree(:2)=patch%get_degree()
            if (present(stride)) stride=[1,nc(1),0]
        type is (nurbs_volume)
            nc=patch%get_nc()
            degree=patch%get_degree()
            if (present(stride)) stride=[1,nc(1),product(nc(:2))]
        end select
    end subroutine

    pure subroutine constraint_shape(patches,conn,nrow,nnz,value_rows)
        use forcad_utils, only: active_span_count
        class(*), intent(in) :: patches(:)
        type(multipatch_connection), intent(in) :: conn
        integer, intent(inout) :: nrow,nnz
        integer, intent(out) :: value_rows
        integer :: axes_a(3),axes_b(3),nc_a(3),nc_b(3),p_a(3),p_b(3),axis,d,ncells

        call patch_basis(patches(conn%patch_a),nc_a,p_a)
        call patch_basis(patches(conn%patch_b),nc_b,p_b)
        if (.not.(conn%has_transition_jet().or.conn%has_projective_jet())) then
            call add_separable_shape(conn,nc_a,p_a,p_b,nrow,nnz,value_rows)
            return
        end if
        call volume_interface_axes(conn,axes_a,axes_b)
        ncells=1
        do d=1,2
            axis=axes_a(d)
            ncells=ncells*active_span_count(patch_axis_knots(patches(conn%patch_a),axis),nc_a(axis),p_a(axis))
        end do
        call add_bernstein_shape(conn,p_a(axes_a),p_b(axes_b),ncells,nrow,nnz,value_rows)
    end subroutine

    !> Pair the two oriented tangents and the normal in the common Bernstein jet assembler.
    pure subroutine append_explicit_constraint(patches,conn,offsets,rowptr,split,col,val,row,pos,power,resolved)
        class(*), intent(in) :: patches(:)
        type(multipatch_connection), intent(in) :: conn
        integer, intent(in), contiguous :: offsets(:)
        integer, intent(inout), contiguous :: rowptr(:), split(:), col(:)
        real(rk), intent(inout), contiguous :: val(:)
        integer, intent(inout) :: row,pos,power
        logical, intent(out) :: resolved
        integer :: a(3),b(3),nc(3),p_a(3),p_b(3),stride_a(3),stride_b(3)
        logical :: reverse_dir(2)

        call volume_interface_axes(conn,a,b,reverse_dir)
        associate(pa=>patches(conn%patch_a), pb=>patches(conn%patch_b))
            call patch_basis(pa,nc,p_a,stride_a)
            call patch_basis(pb,nc,p_b,stride_b)
            if (stride_a(3)==0) reverse_dir(2)=.false.
            call append_bernstein_constraint(conn,offsets,p_a(a),p_b(b),stride_a(a),stride_b(b),reverse_dir, &
                patch_axis_knots(pa,a(3)),patch_axis_knots(pb,b(3)),patch_axis_knots(pa,a(1)),patch_axis_knots(pb,b(1)), &
                patch_axis_knots(pa,a(2)),patch_axis_knots(pb,b(2)),rowptr,split,col,val,row,pos,power,resolved)
        end associate
    end subroutine

    !> Restrict each coefficient field by de Casteljau subdivision along the patch-A tangents.
    pure subroutine restrict_interface_jet(transition, projective, axis1, jet, axis2)
        real(rk), allocatable, intent(in) :: transition(:,:,:,:), projective(:,:,:)
        type(bernstein_span), intent(in) :: axis1
        type(bernstein_span), optional, intent(in) :: axis2
        type(interface_jet), intent(out) :: jet

        if (allocated(transition)) then
            jet%transition = reshape(transition,[size(transition,1),size(transition,2),size(transition,3)*size(transition,4)])
            call restrict_field(jet%transition)
        end if
        if (allocated(projective)) then
            jet%projective = projective
            call restrict_field(jet%projective)
        end if
    contains
        pure subroutine restrict_field(field)
            real(rk), intent(inout) :: field(:,:,:)
            integer :: c, j
            do c = 1,size(field,3)
                do j = 1,size(field,2)
                    call restrict_line(field(:,j,c),axis1)
                end do
                if (.not.present(axis2)) cycle
                do j = 1,size(field,1)
                    call restrict_line(field(j,:,c),axis2)
                end do
            end do
        end subroutine
    end subroutine

    !> Keep the right subdivision at the span start, then the left subdivision at its relative width.
    pure subroutine restrict_line(line,axis)
        real(rk), intent(inout) :: line(:)
        type(bernstein_span), intent(in) :: axis
        real(rk) :: lo,hi,x,w0,w1
        integer :: j,n
        n=size(line)
        if (n==1) return
        lo=axis%domain(1)
        hi=axis%domain(2)
        x=axis%bounds_knot(1)
        if (x/=lo) then
            w0=interval_ratio(x,hi,lo,hi)
            w1=interval_ratio(lo,x,lo,hi)
            do j=n-1,1,-1
                line(:j)=blossom_at(lo,hi,x,w0,w1,line(:j),line(2:j+1))
            end do
        end if
        lo=x
        x=axis%bounds_knot(2)
        if (x/=hi) then
            w0=interval_ratio(x,hi,lo,hi)
            w1=interval_ratio(lo,x,lo,hi)
            do j=2,n
                line(j:)=blossom_at(lo,hi,x,w0,w1,line(j-1:n-1),line(j:))
            end do
        end if
    end subroutine

    !> Interpolate a jet coefficient. A determinant preserves cancellation when a global
    !! normalized breakpoint rounds onto a zero of the coefficient field.
    pure elemental real(rk) function blossom_at(lo,hi,x,w0,w1,a,b) result(value)
        real(rk), intent(in) :: lo,hi,x,w0,w1,a,b
        if (a==b) then
            value=a
        else if ((a<0.0_rk .eqv. b<0.0_rk) .and. min(w0,w1)>=tiny(x)) then
            value=w0*a+w1*b
        else
            value=affine_value_at(lo,hi,x,a,b)
        end if
    end function

    !> Evaluate the affine interpolant from (lo,a) to (hi,b) with a scaled determinant.
    pure elemental real(rk) function affine_value_at(lo,hi,x,a,b) result(value)
        real(rk), intent(in) :: lo,hi,x,a,b
        real(rk) :: residual,width
        integer :: power,shift
        call affine_knot_residual([lo,a],[hi,b],[x,0.0_rk],residual,power,signed_result=.true.)
        call scaled_difference(lo,hi,width,shift)
        value=scale(-residual/fraction(width),power-shift-exponent(width))
    end function

    pure subroutine evaluate_interface_jet(jet,xi,transition,projective,reparameterization)
        type(interface_jet), intent(in) :: jet
        real(rk), intent(in) :: xi(:)
        real(rk), intent(out) :: transition(:,:), projective(0:)
        real(rk), allocatable, intent(in) :: reparameterization(:)

        transition = 0.0_rk
        if (allocated(jet%transition)) then
            transition = reshape(field_value(jet%transition,xi),shape(transition))
        else if (allocated(reparameterization)) then
            transition(size(transition,1),:) = reparameterization
        end if
        projective = 0.0_rk
        projective(0) = 1.0_rk
        if (allocated(jet%projective)) projective = field_value(jet%projective,xi)
    end subroutine

    !> Contract one or two Bernstein axes; the final index enumerates the jet components.
    pure function field_value(field,xi) result(value)
        real(rk), intent(in) :: field(:,:,:), xi(:)
        real(rk) :: value(size(field,3)), column(size(field,2))
        integer :: c,j
        do c=1,size(value)
            do j=1,size(column)
                column(j) = bernstein_value(field(:,j,c),xi(1))
            end do
            value(c) = column(1)
            if (size(xi)==2) value(c) = bernstein_value(column,xi(2))
        end do
    end function

    pure real(rk) function bernstein_value(coefficients,x) result(value)
        real(rk), intent(in) :: coefficients(:),x
        real(rk) :: work(size(coefficients))
        integer :: j
        work=coefficients
        do j=size(work)-1,1,-1
            work(:j)=(1.0_rk-x)*work(:j)+x*work(2:j+1)
        end do
        value=work(1)
    end function

    !> Degree-complete Bernstein collocation on an oriented edge or face.
    !! Degrees and strides follow the two mapped tangents, then the normal.
    !! An edge uses a constant second tangent with zero storage stride.
    pure subroutine append_bernstein_constraint(conn,offsets,degree_a,degree_b,stride_a,stride_b,reverse_dir, &
        knot_na,knot_nb,knot_ta1,knot_tb1,knot_ta2,knot_tb2,rowptr,split,col,val,row,pos,power,resolved)
        type(multipatch_connection), intent(in) :: conn
        integer, intent(in), contiguous :: offsets(:)
        integer, intent(inout), contiguous :: rowptr(:), split(:), col(:)
        real(rk), intent(inout), contiguous :: val(:)
        integer, intent(inout) :: row, pos, power
        logical, intent(out) :: resolved
        integer, intent(in) :: degree_a(3),degree_b(3),stride_a(3),stride_b(3)
        logical, intent(in) :: reverse_dir(2)
        real(rk), contiguous, intent(in) :: knot_na(:),knot_nb(:),knot_ta1(:),knot_tb1(:),knot_ta2(:),knot_tb2(:)
        real(rk), allocatable :: basis_na(:,:), basis_nb(:,:), basis_ta1(:,:), basis_ta2(:,:), basis_tb1(:,:), basis_tb2(:,:)
        type(bernstein_span) :: axis_a(2), axis_b(2)
        type(interface_jet) :: jet
        real(rk), allocatable :: projective_jet(:,:,:)
        real(rk), allocatable :: transition_value(:,:), projective_value(:)
        real(rk), allocatable :: composition(:,:,:,:), weights(:,:,:)
        real(rk) :: xi(2), entry, choose
        integer :: first_na, first_nb, tangent_order(2)
        integer :: base_a,base_b,ndim,d
        integer :: degree_row(2)
        integer :: n, r, j, alpha1, alpha2, alpha_n, local1, local2, local_n
        integer :: span1, span2, point1, point2
        ndim = merge(2,3,stride_a(2)==0)
        n = conn%get_continuity()

        tangent_order = 0
        if (allocated(conn%transition_jet)) then
            do d = 1, ndim-1
                do j = 1,n
                    if (all(conn%transition_jet(:,:,d,j)==0.0_rk)) cycle
                    ! A displacement starting at order j contributes at most floor(n/j) tangent derivatives.
                    tangent_order(d) = min(degree_b(d),n/j)
                    exit
                end do
            end do
        end if
        if (allocated(conn%projective_jet)) then
            projective_jet = conn%projective_jet
            call scale_projective_jet(projective_jet,power)
        end if

        allocate(basis_na(0:degree_a(3),0:n), basis_nb(0:degree_b(3),0:n))
        allocate(basis_ta1(0:degree_a(1),0:0), basis_tb1(0:degree_b(1),0:tangent_order(1)))
        allocate(basis_ta2(0:degree_a(2),0:0), basis_tb2(0:degree_b(2),0:tangent_order(2)))
        allocate(transition_value(3,n), projective_value(0:n), weights(0:n,0:n,0:n))
        call eval_bernstein_endpoint(knot_na,degree_a(3),conn%side_a,first_na,basis_na,resolved)
        if (.not. resolved) return
        call eval_bernstein_endpoint(knot_nb,degree_b(3),conn%side_b,first_nb,basis_nb,resolved)
        if (.not. resolved) return

        do r = 0, n
            degree_row = jet_row_degree(conn,degree_a(:2),degree_b(:2),r,[1,2])
            do span2 = degree_a(2) + 1, size(knot_ta2)-degree_a(2)-1
                if (knot_ta2(span2+1) <= knot_ta2(span2)) cycle
                call prepare_bernstein_pair(knot_ta2,knot_tb2,degree_a(2),span2, &
                    reverse_dir(2),axis_a(2),axis_b(2),resolved,order=min(r,tangent_order(2)))
                if (.not. resolved) return
                do span1 = degree_a(1) + 1, size(knot_ta1)-degree_a(1)-1
                    if (knot_ta1(span1+1) <= knot_ta1(span1)) cycle
                    call prepare_bernstein_pair(knot_ta1,knot_tb1,degree_a(1),span1, &
                        reverse_dir(1),axis_a(1),axis_b(1),resolved,order=min(r,tangent_order(1)))
                    if (.not. resolved) return
                    base_a = offsets(conn%patch_a)+1+dot_product(stride_a,[axis_a(1)%first,axis_a(2)%first,first_na]-1)
                    base_b = offsets(conn%patch_b)+1+dot_product(stride_b,[axis_b(1)%first,axis_b(2)%first,first_nb]-1)
                    call restrict_interface_jet(conn%transition_jet,projective_jet,axis_a(1),jet,axis_a(2))
                    do point2 = 1, degree_row(2)+1
                        xi(2) = (real(point2,rk)-0.5_rk)/real(degree_row(2)+1,rk)
                        call eval_bernstein_pair(axis_a(2),axis_b(2),reverse_dir(2),xi(2), &
                                basis_ta2,basis_tb2(:,:min(r,tangent_order(2))),resolved)
                        if (.not. resolved) return
                        do point1 = 1, degree_row(1)+1
                            xi(1) = (real(point1,rk)-0.5_rk)/real(degree_row(1)+1,rk)
                            call eval_bernstein_pair(axis_a(1),axis_b(1),reverse_dir(1),xi(1), &
                                basis_ta1,basis_tb1(:,:min(r,tangent_order(1))),resolved)
                            if (.not. resolved) return
                            call evaluate_interface_jet(jet,xi,transition_value(:ndim,:),projective_value,conn%reparameterization)
                            if (ndim==2) then
                                transition_value(3,:) = transition_value(2,:)
                                transition_value(2,:) = 0.0_rk
                            end if
                            call multipatch_composition_coefficients(transition_value(:,:r), composition)

                            weights(:r,:r,:r) = 0.0_rk
                            choose = 1.0_rk
                            do j = 0, r
                                weights(:r,:r,:r) = weights(:r,:r,:r) + (choose*projective_value(r-j))*composition(j,:,:,:)
                                if (j<r) choose = choose*real(r-j,rk)/real(j+1,rk)
                            end do
                            row = row + 1
                            do local2 = 0, degree_a(2)
                                do local1 = 0, degree_a(1)
                                    do local_n = 0, degree_a(3)
                                        entry = basis_ta1(local1,0)*basis_ta2(local2,0)*&
                                            basis_na(local_n,r)
                                        if (.not. structural_nonzero(entry)) cycle
                                        col(pos) = base_a + dot_product(stride_a,[local1,local2,local_n])
                                        val(pos) = entry
                                        pos = pos + 1
                                    end do
                                end do
                            end do
                            split(row) = pos
                            do local2 = 0, degree_b(2)
                                do local1 = 0, degree_b(1)
                                    do local_n = 0, degree_b(3)
                                        entry = 0.0_rk
                                        do alpha_n = 0,r
                                            do alpha2 = 0,min(tangent_order(2),r-alpha_n)
                                                do alpha1 = 0,min(tangent_order(1),r-alpha_n-alpha2)
                                                    entry = entry + weights(alpha1,alpha2,alpha_n)*basis_tb1(local1,alpha1)* &
                                                        basis_tb2(local2,alpha2)*basis_nb(local_n,alpha_n)
                                                end do
                                            end do
                                        end do
                                        if (.not. structural_nonzero(entry)) cycle
                                        col(pos) = base_b + dot_product(stride_b,[local1,local2,local_n])
                                        val(pos) = -entry
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

    !> Check incident patches and decode both boundary labels before indexing the patch table.
    pure subroutine resolve_connection_sides(npatch, ndim, pa, label_a, pb, label_b, sa, sb, err)
        integer, intent(in) :: npatch, ndim, pa, pb
        character(len=*), intent(in) :: label_a, label_b
        integer, intent(out) :: sa, sb
        type(debug), intent(inout) :: err
        character(len=7), parameter :: names(3) = [character(len=7) :: "curve", "surface", "volume"]
        character(len=*), parameter :: category = "forcad_multipatch_"

        sa = connection_side(label_a, ndim)
        sb = connection_side(label_b, ndim)
        if (min(pa,pb) < 1 .or. max(pa,pb) > npatch) then
            call err%set(code=201, severity=1, category=category//trim(names(ndim)), location="connect", &
                message=merge("Invalid first patch id. ", "Invalid second patch id.", pa<1 .or. pa>npatch), &
                suggestion="Use a patch id from add_patch's optional id argument.")
        else if (min(sa,sb) == 0) then
            call err%set(code=202, severity=1, category=category//trim(names(ndim)), location="connect", &
                message="Invalid "//trim(names(ndim))//" side label.", &
                suggestion="Use left/right for curves or a patch axis with min/max, such as u_min, for surfaces and volumes.")
        end if
    end subroutine

    !> Resolve the sides and check the normal orders before constructing the interface jet.
    pure subroutine resolve_patch_interface(patches,npatch,ndim,pa,label_a,pb,label_b,continuity,sa,sb,order,err)
        class(*), optional, intent(in) :: patches(:)
        integer, intent(in) :: npatch,ndim,pa,pb
        character(*), intent(in) :: label_a,label_b
        integer, optional, intent(in) :: continuity
        integer, intent(out) :: sa,sb,order
        type(debug), intent(inout) :: err
        if (.not.err%ok) return
        call resolve_connection_sides(npatch,ndim,pa,label_a,pb,label_b,sa,sb,err)
        if (.not.err%ok) return
        order=0
        if (present(continuity)) order=continuity
        if (order < -1 .or. order > normal_degree(patches,pa,sa,pb,sb)) then
            call err%set(code=203,severity=1,category=patch_category(ndim),location="connect", &
                message="Invalid interface continuity order.", &
                suggestion="Use -1 for independent patches, or an order up to the lower of the two normal degrees.")
        end if
    end subroutine

    pure integer function normal_degree(patches,pa,sa,pb,sb) result(degree)
        class(*), intent(in) :: patches(:)
        integer, intent(in) :: pa,sa,pb,sb
        integer :: nc(3),p(3),q(3)
        call patch_basis(patches(pa),nc,p)
        call patch_basis(patches(pb),nc,q)
        degree=min(p((sa+1)/2),q((sb+1)/2))
    end function

    !> Commit the connection once its jets and oriented tangential spaces are compatible.
    pure subroutine record_connection(patches,npatch,ndim,connections,nconnection,conn,err)
        class(*), intent(in) :: patches(:)
        integer, intent(in) :: npatch,ndim
        type(multipatch_connection), allocatable, intent(inout) :: connections(:)
        integer, intent(inout) :: nconnection
        type(multipatch_connection), intent(inout) :: conn
        type(debug), intent(inout) :: err
        if (allocated(conn%prolongation)) deallocate(conn%prolongation)
        if (.not.valid_connection_metadata(conn,npatch,ndim)) then
            call err%set(code=203,severity=1,category=patch_category(ndim),location="connect", &
                message="Invalid continuity reparameterization.", &
                suggestion="Check the jet dimensions, finite derivatives and the side-consistent normal sign.")
            return
        end if
        if (conn%continuity>=0) then
            if (.not.matching_trace_space(patches,conn,ndim,sqrt(epsilon(1.0_rk)))) then
                if (conn%continuity==0 .and. .not.conn%has_projective_jet()) call set_nested_trace(patches,conn,ndim)
                if (.not.allocated(conn%prolongation)) then
                    call err%set(code=204,severity=1,category=patch_category(ndim),location="connect", &
                        message="The connected sides have incompatible tangential trace spaces.", &
                        suggestion="Use matching spaces, or clamped nested edges or faces for order zero.")
                    return
                end if
            end if
        end if
        call append_connection(connections,nconnection,conn)
    end subroutine

    !> Form the knot-insertion factors in side A's tangential axis order. Each factor
    !! reflects its own knot coordinate; the tensor product accounts for the face swap.
    pure subroutine set_nested_trace(patches,conn,ndim)
        class(*), intent(in) :: patches(:)
        type(multipatch_connection), intent(inout) :: conn
        integer, intent(in) :: ndim
        type boundary_indices
            integer, allocatable :: ids(:)
        end type
        type(boundary_indices) :: side(2)
        type(axis_transfer) :: factor(ndim-1)
        type(trace_prolongation), allocatable :: transfer
        integer, allocatable :: col(:)
        integer :: patch(2),axes(3,2),nc(3),degree(3),extent(ndim-1,2),p(ndim-1,2),i,d,c,f
        logical :: reflect(2),ok

        patch=[conn%patch_a,conn%patch_b]
        call volume_interface_axes(conn,axes(:,1),axes(:,2),reflect)
        do i=1,2
            call patch_basis(patches(patch(i)),nc,degree)
            extent(:,i)=nc(axes(:ndim-1,i))
            p(:,i)=degree(axes(:ndim-1,i))
            side(i)%ids=interpolated_trace_ids(nc(:ndim),patch_axis_knots(patches(patch(i)),axes(3,i)), &
                conn,i==2,native=.true.)
            if (size(side(i)%ids)==0) return
        end do
        if (any(p(:,1)/=p(:,2)) .or. all(extent(:,1)==extent(:,2))) return
        c=merge(1,2,all(extent(:,1)<=extent(:,2)))
        f=3-c
        if (any(extent(:,c)>extent(:,f))) return
        do d=1,ndim-1
            associate(t=>factor(d))
                call nested_edge_transfer(patch_axis_knots(patches(patch(c)),axes(d,c)), &
                    patch_axis_knots(patches(patch(f)),axes(d,f)),p(d,c),t%ptr,t%col,t%val,ok,reflect(d))
            end associate
            if (.not.ok) return
        end do
        allocate(transfer)
        transfer%coarse_a=c==1
        call tensor_transfer(factor,extent(:,c),conn%swap.and.c==2,conn%swap.and.f==2,transfer%ptr,col,transfer%val,ok)
        if (.not.ok) return
        transfer%coarse=side(c)%ids(col)
        call move_alloc(side(f)%ids,transfer%fine)
        call move_alloc(transfer,conn%prolongation)
    end subroutine

    pure function patch_category(ndim) result(category)
        integer, intent(in) :: ndim
        character(:), allocatable :: category
        character(7), parameter :: names(3)=[character(7)::"curve","surface","volume"]
        category="forcad_multipatch_"//trim(names(ndim))
    end function

    !> Resolve a boundary label using the side numbering of the patch dimension.
    pure integer function connection_side(label, ndim) result(side)
        character(len=*), intent(in) :: label
        integer, intent(in) :: ndim

        if (ndim == 1) then
            select case (label)
            case ("left", "LEFT", "Left", "min", "MIN", "Min")
                side = 1
            case ("right", "RIGHT", "Right", "max", "MAX", "Max")
                side = 2
            case default
                side = 0
            end select
        else
            select case (label)
            case ("u_min", "umin", "u-", "U_MIN", "UMIN", "U-")
                side = 1
            case ("u_max", "umax", "u+", "U_MAX", "UMAX", "U+")
                side = 2
            case ("v_min", "vmin", "v-", "V_MIN", "VMIN", "V-")
                side = 3
            case ("v_max", "vmax", "v+", "V_MAX", "VMAX", "V+")
                side = 4
            case ("w_min", "wmin", "w-", "W_MIN", "WMIN", "W-")
                side = merge(5, 0, ndim == 3)
            case ("w_max", "wmax", "w+", "W_MAX", "WMAX", "W+")
                side = merge(6, 0, ndim == 3)
            case default
                side = 0
            end select
        end if
    end function

    !> Store a validated connection with amortized growth of the interface table.
    pure subroutine append_connection(connections, count, conn)
        type(multipatch_connection), allocatable, intent(inout) :: connections(:)
        integer, intent(inout) :: count
        type(multipatch_connection), intent(in) :: conn
        type(multipatch_connection), allocatable :: larger(:)

        if (.not. allocated(connections)) then
            allocate(connections(multipatch_growth_capacity(0)))
        else if (count == size(connections)) then
            allocate(larger(multipatch_growth_capacity(count)))
            larger(:count) = connections
            call move_alloc(larger, connections)
        end if
        count = count + 1
        connections(count) = conn
    end subroutine

    !> Validate stored incidence, side numbering and the supplied change of parameters.
    pure logical function valid_connection_metadata(conn, npatch, ndim) result(valid)
        type(multipatch_connection), intent(in) :: conn
        integer, intent(in) :: npatch, ndim
        integer :: normal_sign

        valid = .false.
        if (min(conn%patch_a,conn%patch_b)<1 .or. max(conn%patch_a,conn%patch_b)>npatch) return
        if (min(conn%side_a,conn%side_b)<1 .or. max(conn%side_a,conn%side_b)>2*ndim) return
        if (conn%continuity < -1) return
        normal_sign = merge(-1, 1, mod(conn%side_a,2)==mod(conn%side_b,2))
        valid = multipatch_valid_reparameterization(conn, normal_sign, ndim)
    end function

    !> Dimension of the physical control coordinates on a stored patch.
    pure integer function patch_ncoord(patch) result(ncoord)
        class(*), intent(in) :: patch
        select type (patch)
        type is (nurbs_curve)
            ncoord=size(patch%get_Xc(1))
        type is (nurbs_surface)
            ncoord=size(patch%get_Xc(1))
        type is (nurbs_volume)
            ncoord=size(patch%get_Xc(1))
        class default
            ncoord=0
        end select
    end function

    !> Validate a common physical coordinate map before changing any connected patch.
    pure subroutine check_patch_transform(patches,npatch,ndim,linear,translation,center,err)
        use forcad_utils, only: check_affine_map
        class(*), optional, intent(in) :: patches(:)
        integer, intent(in) :: npatch,ndim
        real(rk), intent(in) :: linear(:,:)
        real(rk), optional, intent(in) :: translation(:),center(:)
        type(debug), intent(inout) :: err
        integer :: i,ncoord

        if (.not.err%ok.or.npatch<1) return
        ncoord=patch_ncoord(patches(1))
        call check_affine_map(ncoord,linear,translation,center,err,patch_category(ndim),207,207)
        if (.not.err%ok) return
        do i=2,npatch
            if (patch_ncoord(patches(i))==ncoord) cycle
            call err%set(code=207,severity=1,category=patch_category(ndim),location="transform", &
                message="All patches must use one coordinate dimension.",suggestion="Embed every patch in the same physical space.")
            return
        end do
    end subroutine

    !> Compare tensor-product trace bases in the oriented tangential directions.
    pure logical function matching_trace_space(patches,conn,ndim,tol) result(ok)
        class(*), intent(in) :: patches(:)
        type(multipatch_connection), intent(in) :: conn
        integer, intent(in) :: ndim
        real(rk), intent(in) :: tol
        integer :: a(3),b(3),na(3),nb(3),pa(3),pb(3),d
        logical :: reverse(2)
        ok=.true.
        if (ndim==1 .or. allocated(conn%prolongation)) return
        call volume_interface_axes(conn,a,b,reverse)
        associate(left=>patches(conn%patch_a),right=>patches(conn%patch_b))
            call patch_basis(left,na,pa)
            call patch_basis(right,nb,pb)
            do d=1,ndim-1
                ok=multipatch_compatible_trace_space(patch_axis_knots(left,a(d)),na(a(d)),pa(a(d)), &
                    patch_axis_knots(right,b(d)),nb(b(d)),pb(b(d)),reverse(d),tol)
                if (.not.ok) return
            end do
        end associate
    end function

    !> Validate geometry and the homogeneous interface equations shared by analysis and design.
    !! In coordinates centered at a common physical origin, each row and component must satisfy
    !! \( |(AH)_{rc}| \le a_{\rm tol} + r_{\rm tol}\sum_j |A_{rj}H_{jc}| \),
    !! including summation roundoff. The denominator is checked with the same jet.
    !! Relative tolerance also controls the comparison of oriented knots and projective scales.
    pure logical function valid_multipatch(patches,npatch,connections,nconnection,ndim,rtol,atol) result(ok)
        class(*), optional, intent(in) :: patches(:)
        type(multipatch_connection), optional, intent(in) :: connections(:)
        integer, intent(in) :: npatch,nconnection,ndim
        real(rk), optional, intent(in) :: rtol,atol
        real(rk) :: relative,absolute
        real(rk), allocatable :: values(:),points(:,:),weights(:)
        integer, allocatable :: offsets(:),rowptr(:),col(:),powers(:)
        integer :: i,ncoord
        logical :: homogeneous
        type(debug) :: status

        ok=.false.
        if (npatch<1 .or. .not.present(patches) .or. nconnection<0) return
        if (npatch>size(patches)) return
        if (nconnection>0) then
            if (.not.present(connections)) return
            if (nconnection>size(connections)) return
        end if
        do i=1,npatch
            if (.not.patch_is_valid(patches(i))) return
        end do
        call continuity_tolerances(rtol,atol,relative,absolute,ok)
        if (.not.ok) return
        if (nconnection==0) return
        homogeneous=.false.
        do i=1,nconnection
            associate(conn=>connections(i))
                ok=valid_connection_metadata(conn,npatch,ndim)
                if (.not.ok) return
                if (conn%continuity<0) cycle
                associate(a=>patches(conn%patch_a),b=>patches(conn%patch_b))
                    ok=conn%continuity<=normal_degree(patches,conn%patch_a,conn%side_a,conn%patch_b,conn%side_b) .and. &
                        patch_ncoord(a)==patch_ncoord(b)
                end associate
                if (.not.ok) return
                ok=matching_trace_space(patches,conn,ndim,relative)
                if (.not.ok) return
                homogeneous=homogeneous.or.conn%has_projective_jet()
            end associate
        end do
        offsets=multipatch_dof_offsets(patches,npatch)
        call constraint_rows(patches,connections,nconnection,offsets,rowptr,col,values,err=status,powers=powers,rtol=relative)
        ok=status%ok
        if (.not.ok .or. size(col)==0) return
        ncoord=0
        do i=1,npatch
            homogeneous=homogeneous.or.patch_is_rational(patches(i))
            ncoord=max(ncoord,patch_ncoord(patches(i)))
        end do
        allocate(points(offsets(npatch+1),ncoord),weights(0))
        points(1,:)=0.0_rk
        if (homogeneous) weights=storage_weights(patches,offsets)
        do i=1,npatch
            select type (patch=>patches(i))
            type is (nurbs_curve)
                call store_constraint_points(points,offsets(i),patch%get_Xc())
            type is (nurbs_surface)
                call store_constraint_points(points,offsets(i),patch%get_Xc())
            type is (nurbs_volume)
                call store_constraint_points(points,offsets(i),patch%get_Xc())
            end select
        end do
        call check_constraint_residuals(rowptr,col,values,powers,points,weights,relative,absolute,ok)
    end function

    !> Select connections with continuity equations and the requested optional kind.
    pure logical function include_constraint(conn, geometric) result(selected)
        type(multipatch_connection), intent(in) :: conn
        logical, optional, intent(in) :: geometric

        selected = conn%continuity >= 0
        if (present(geometric)) selected = selected .and. (conn%geometric .eqv. geometric)
    end function

    !> The rational factor and all its derivatives share a homogeneous scale.
    !! Keep it outside Bernstein evaluation so a small factor survives collocation.
    pure subroutine scale_projective_jet(jet, power)
        real(rk), intent(inout) :: jet(:,:,:)
        integer, intent(out) :: power
        power = exponent(maxval(abs(jet)))
        if (power>0) power = min(power,max(0,minval(exponent(jet),mask=jet/=0.0_rk)-minexponent(1.0_rk)))
        jet = scale(jet,-power)
    end subroutine

    pure subroutine continuity_tolerances(rtol, atol, relative, absolute, valid)
        real(rk), optional, intent(in) :: rtol, atol
        real(rk), intent(out) :: relative, absolute
        logical, intent(out) :: valid
        relative = sqrt(epsilon(1.0_rk))
        absolute = 0.0_rk
        if (present(rtol)) relative = rtol
        if (present(atol)) absolute = atol
        valid = ieee_is_finite(relative) .and. ieee_is_finite(absolute) .and. relative>=0.0_rk .and. absolute>=0.0_rk
    end subroutine

    !> Count one normal jet per oriented tangential control for compact connections.
    pure subroutine add_separable_shape(conn, nc_a, degree_a, degree_b, nrow, nnz, value_rows)
        type(multipatch_connection), intent(in) :: conn
        integer, intent(in) :: nc_a(:), degree_a(:), degree_b(:)
        integer, intent(inout) :: nrow, nnz
        integer, intent(out), optional :: value_rows
        integer :: da, db, rows

        da = (conn%side_a+1)/2
        db = (conn%side_b+1)/2
        rows = product(nc_a(:da-1))*product(nc_a(da+1:))
        if (present(value_rows)) value_rows=rows
        rows = rows*(conn%continuity+1)
        nrow = nrow+rows
        nnz = nnz+rows*(degree_a(da)+degree_b(db)+2)
    end subroutine

    !> Append one connection's trace equations in patch-concatenated control numbering.
    pure subroutine append_separable_constraint(conn, nc_a, nc_b, knot_a, knot_b, offsets, rowptr, split, col, val, row, pos)
        type(multipatch_connection), intent(in) :: conn
        integer, contiguous, intent(in) :: nc_a(:), nc_b(:)
        real(rk), contiguous, intent(in) :: knot_a(:), knot_b(:)
        integer, intent(in) :: offsets(:)
        integer, intent(inout) :: rowptr(:), split(:), col(:), row, pos
        real(rk), intent(inout) :: val(:)
        real(rk), allocatable :: a(:,:), b(:,:), chain(:,:)
        real(rk) :: length_a, length_b
        integer, allocatable :: plane_a(:), plane_b(:)
        integer :: stride_a, stride_b, normal_sign, r, j

        call endpoint_stencil(nc_a, knot_a, conn%side_a, conn%continuity, .false., .false., [.false.,.false.], &
            plane_a, stride_a, length_a, a)
        call endpoint_stencil(nc_b, knot_b, conn%side_b, conn%continuity, conn%reverse, conn%swap, conn%flip, &
            plane_b, stride_b, length_b, b)
        plane_a = plane_a + offsets(conn%patch_a)
        plane_b = plane_b + offsets(conn%patch_b)
        normal_sign = merge(-1, 1, mod(conn%side_a,2) == mod(conn%side_b,2))
        call multipatch_chain_rule_coefficients(conn, normal_sign, chain)
        ! Rows encode D^r H_A - scale * D^r(H_B composed with the transition).
        do r = 1, conn%continuity
            a(r,:) = length_a**r * a(r,:)
            ! Combine parameter scales before applying a potentially large derivative.
            chain(r:,r) = chain(r:,r) * length_b**r
        end do
        do r = conn%continuity, 0, -1
            b(r,:) = -matmul(chain(r,0:r), b(0:r,:))
        end do
        do r = 0, conn%continuity
            do j = 1, size(plane_a)
                call append_row_side(a, r, plane_a(j), stride_a, col, val, pos)
                split(row+1) = pos
                call append_row_side(b, r, plane_b(j), stride_b, col, val, pos)
                row = row + 1
                rowptr(row+1) = pos
            end do
        end do
    end subroutine

    !> Apply the projective factor to side B of each row after assembling the basis jets.
    pure subroutine finish_constraint_block(first, last, rowptr, split, col, val, factor, power, valid, weights, powers)
        integer, intent(in) :: first, last, rowptr(:), split(:), col(:), power
        real(rk), intent(inout) :: val(:)
        real(rk), intent(in) :: factor
        logical, intent(out) :: valid
        real(rk), allocatable, intent(in) :: weights(:)
        integer, allocatable, optional, intent(inout) :: powers(:)
        integer :: r, lo, hi, k

        valid = all(ieee_is_finite(val(rowptr(first):rowptr(last+1)-1)))
        if (.not. valid) return
        do r = first,last
            lo = rowptr(r)
            hi = rowptr(r+1)-1
            block
                integer :: shift(lo:hi)
                shift = 0
                if (factor /= 1.0_rk .or. power /= 0) then
                    do k=split(r),hi
                        shift(k) = exponent(val(k)) + exponent(factor) + power
                        val(k) = fraction(val(k))*fraction(factor)
                    end do
                end if
                if (present(powers)) then
                    powers(lo:hi) = shift
                else if (allocated(weights)) then
                    call equilibrate_field_row(val(lo:hi),weights(col(lo:hi)),valid,shift)
                    if (.not. valid) return
                else
                    do k=split(r),hi
                        if (shift(k)==0) cycle
                        valid = exponent(val(k))+shift(k) <= maxexponent(factor)
                        if (.not. valid) return
                        val(k) = scale(val(k),shift(k))
                        valid = val(k) /= 0.0_rk .and. ieee_is_finite(val(k))
                        if (.not. valid) return
                    end do
                end if
            end block
        end do
    end subroutine

    !> An omitted logical option has its default value of false.
    pure logical function requested(option)
        logical, optional, intent(in) :: option
        requested = .false.
        if (present(option)) requested = option
    end function

    !> Express both homogeneous traces in the fine basis. Each row contains one
    !! fine control and the coarse controls whose knot-insertion coefficients reach it.
    !! Side A is positive and side B negative, also when B is coarse.
    pure subroutine append_nested_constraint(conn,offsets,rowptr,split,col,val,row,pos)
        type(multipatch_connection), intent(in) :: conn
        integer, intent(in) :: offsets(:)
        integer, intent(inout) :: rowptr(:),split(:),col(:),row,pos
        real(rk), intent(inout) :: val(:)
        integer :: i,lo,hi,n,c,f,coarse,fine
        real(rk) :: s

        associate(transfer=>conn%prolongation)
            coarse=merge(conn%patch_a,conn%patch_b,transfer%coarse_a)
            fine=merge(conn%patch_b,conn%patch_a,transfer%coarse_a)
            s=merge(1.0_rk,-1.0_rk,transfer%coarse_a)
            do i=1,size(transfer%fine)
                lo=transfer%ptr(i)
                hi=transfer%ptr(i+1)-1
                n=hi-lo+1
                c=merge(pos,pos+1,transfer%coarse_a)
                f=merge(pos+n,pos,transfer%coarse_a)
                col(c:c+n-1)=offsets(coarse)+transfer%coarse(lo:hi)
                val(c:c+n-1)=s*transfer%val(lo:hi)
                col(f)=offsets(fine)+transfer%fine(i)
                val(f)=-s
                row=row+1
                split(row)=max(c,f)
                pos=pos+n+1
                rowptr(row+1)=pos
            end do
        end associate
    end subroutine

    !> Assemble one-based CSR equations \(Aq=0\) with row pointers `rowptr`, columns `col` and coefficients `val`.
    !! Raw columns index homogeneous controls in patch storage order. A nested edge or face
    !! contributes one row per fine control, with both traces expressed in the fine basis.
    !! Set `field=.true.` for equilibrated ordinary-field rows. With `active_only=.true.`,
    !! they use the same IDs as the elements and omit rows already enforced by sharing.
    !! The optional `geometric` filter selects geometric or parametric connections.
    !! Inspect `err` for incompatible denominators or unrepresentable coefficients.
    pure subroutine constraint_rows(patches,connections,nconnection,offsets,rowptr,col,val, &
            geometric,field,err,powers,rtol,active_only)
        class(*), optional, intent(in) :: patches(:)
        type(multipatch_connection), optional, intent(in) :: connections(:)
        integer, intent(in) :: nconnection
        integer, contiguous, intent(in) :: offsets(:)
        integer, allocatable, intent(out) :: rowptr(:),col(:)
        real(rk), allocatable, intent(out) :: val(:)
        logical, optional, intent(in) :: geometric,field,active_only
        type(debug), optional, intent(out) :: err
        integer, optional, allocatable, intent(out) :: powers(:)
        real(rk), optional, intent(in) :: rtol
        integer, allocatable :: split(:),map(:)
        logical, allocatable :: shareable(:)
        real(rk), allocatable :: weights(:)
        real(rk) :: factor,tol
        integer :: i,pa,pb,nrow,nnz,row,pos,first,power,value_rows(nconnection)
        logical :: solver,ordinary,valid,compatible,resolved

        if (present(err)) call err%reset()
        solver = requested(active_only)
        ordinary = solver
        if (present(field)) ordinary = field
        if (solver .and. .not.ordinary) then
            rowptr=[1]
            allocate(col(0),val(0))
            if (present(err)) call err%set(code=204,severity=1,category="forcad_multipatch",location="cmp_dof_constraint", &
                message="Active solver IDs identify ordinary scalar coefficients.", &
                suggestion="Omit field=.false. for solver rows, or omit active_only for raw homogeneous rows.")
            return
        end if
        if (ordinary) then
            do i=1,size(offsets)-1
                if (patch_is_valid(patches(i))) cycle
                rowptr=[1]
                allocate(col(0),val(0))
                if (present(err)) call err%set(code=205,severity=1,category="forcad_multipatch",location="cmp_dof_constraint", &
                    message="Field constraints require valid stored patches.",suggestion="Recover or replace invalid patches.")
                return
            end do
        end if
        tol=sqrt(epsilon(1.0_rk))
        if (present(rtol)) tol=rtol
        nrow=0
        nnz=0
        do i=1,nconnection
            if (.not.include_constraint(connections(i),geometric)) cycle
            if (allocated(connections(i)%prolongation)) then
                value_rows(i)=size(connections(i)%prolongation%fine)
                nrow=nrow+value_rows(i)
                nnz=nnz+value_rows(i)+size(connections(i)%prolongation%val)
            else
                call constraint_shape(patches,connections(i),nrow,nnz,value_rows(i))
            end if
        end do
        allocate(rowptr(nrow+1),split(nrow),col(nnz),val(nnz))
        if (present(powers)) allocate(powers(nnz))
        rowptr(1)=1
        if (nrow==0) return
        if (ordinary) weights=storage_weights(patches,offsets)
        if (solver) call multipatch_dof_map(patches,connections,nconnection,offsets,.true.,map,shareable)
        compatible=.true.
        valid=.true.
        resolved=.true.
        row=0
        pos=1
        do i=1,nconnection
            associate(conn=>connections(i))
                if (.not.include_constraint(conn,geometric)) cycle
                first=row+1
                pa=conn%patch_a
                pb=conn%patch_b
                factor=1.0_rk
                power=0
                if (.not.conn%has_projective_jet()) then
                    call connection_weight_scale(patches,conn,tol,factor,compatible,power)
                    if (.not.compatible) exit
                end if
                if (allocated(conn%prolongation)) then
                    call append_nested_constraint(conn,offsets,rowptr,split,col,val,row,pos)
                else if (conn%has_transition_jet() .or. conn%has_projective_jet()) then
                    call append_explicit_constraint(patches,conn,offsets,rowptr,split,col,val,row,pos,power,resolved)
                    if (.not.resolved) exit
                else
                    call append_separable_constraint(conn,patch_control_counts(patches(pa)),patch_control_counts(patches(pb)), &
                        patch_axis_knots(patches(pa),(conn%side_a+1)/2),patch_axis_knots(patches(pb),(conn%side_b+1)/2), &
                        offsets,rowptr,split,col,val,row,pos)
                end if
                call finish_constraint_block(first,row,rowptr,split,col,val,factor,power,valid,weights,powers)
                if (.not.valid) exit
                ! A directly shared trace already enforces its value equations exactly.
                if (solver) then
                    if (shareable(i)) val(rowptr(first):rowptr(first+value_rows(i))-1)=0.0_rk
                end if
            end associate
        end do
        call finish_constraint_matrix(rowptr,col,val,pos,powers,compatible,valid,ordinary,err,resolved)
        if (solver .and. size(col)>0) call compact_solver_rows(map,rowptr,col,val,err)
    end subroutine

    !> Validate the patch before assigning a storage ID or recording any of its faces.
    pure subroutine check_new_patch(ndim,patch,ready,err,id)
        integer, intent(in) :: ndim
        class(*), intent(in) :: patch
        logical, intent(in) :: ready
        type(debug), intent(inout) :: err
        integer, optional, intent(out) :: id

        if (present(id)) id=0
        if (.not.err%ok) return
        if (.not.ready) then
            call err%set(code=205,severity=1,category=patch_category(ndim),location="add_patch", &
                message="Patch has an active error state.",suggestion="Resolve the patch error before adding it.")
        else if (.not.patch_is_valid(patch)) then
            call err%set(code=205,severity=1,category=patch_category(ndim),location="add_patch", &
                message="Patch geometry is incomplete or internally inconsistent.", &
                suggestion="Initialize consistent knots, controls, degrees and optional weights.")
        end if
    end subroutine

    pure logical function patch_is_valid(patch) result(valid)
        class(*), intent(in) :: patch
        select type (patch)
        type is (nurbs_curve)
            valid=patch%is_valid()
        type is (nurbs_surface)
            valid=patch%is_valid()
        type is (nurbs_volume)
            valid=patch%is_valid()
        class default
            valid=.false.
        end select
    end function

    !> Effective control weights in patch-concatenated storage order.
    pure function storage_weights(patches,offsets) result(weights)
        class(*), intent(in) :: patches(:)
        integer, intent(in) :: offsets(:)
        real(rk), allocatable :: weights(:)
        integer :: i
        allocate(weights(offsets(size(offsets))),source=1.0_rk)
        do i=1,size(offsets)-1
            if (.not.patch_is_rational(patches(i))) cycle
            associate(w=>weights(offsets(i)+1:offsets(i+1)))
                select type (patch=>patches(i))
                type is (nurbs_curve)
                    w=patch%get_Wc()
                type is (nurbs_surface)
                    w=patch%get_Wc()
                type is (nurbs_volume)
                    w=patch%get_Wc()
                end select
            end associate
        end do
    end function

    !> Apply the Boolean prolongation q_i = q_solver(map(i)) to ordinary-field CSR rows.
    !! Reuse the storage map as a sparse accumulator index after translating the columns.
    !! Only the correction for each touched coefficient needs extra row-local workspace.
    pure subroutine compact_solver_rows(map,rowptr,col,val,err)
        integer, intent(inout) :: map(:)
        integer, allocatable, intent(inout) :: rowptr(:),col(:)
        real(rk), allocatable, intent(inout) :: val(:)
        type(debug), optional, intent(inout) :: err
        real(rk), allocatable :: correction(:)
        real(rk) :: total,remainder,peak
        integer :: r,i,j,id,slot,first,last,dest,nrow,lo,hi

        allocate(correction(maxval(rowptr(2:)-rowptr(:size(rowptr)-1))))
        do i=1,size(col)
            col(i)=map(col(i))
        end do
        map=0
        nrow=0
        first=1
        lo=1
        do r=1,size(rowptr)-1
            hi=rowptr(r+1)-1
            last=first-1
            do i=lo,hi
                id=col(i)
                if (id==0) cycle
                slot=map(id)
                if (slot==0) then
                    last=last+1
                    col(last)=id
                    val(last)=val(i)
                    map(id)=last
                    correction(last-first+1)=0.0_rk
                else
                    ! Keep the rounding error from each addition until all occurrences of this ID are folded.
                    call two_sum(val(slot),val(i),total,remainder)
                    j=slot-first+1
                    correction(j)=correction(j)+remainder
                    val(slot)=total
                end if
            end do
            dest=first
            do i=first,last
                map(col(i))=0
                total=val(i)+correction(i-first+1)
                if (total==0.0_rk) cycle
                col(dest)=col(i)
                val(dest)=total
                dest=dest+1
            end do
            if (dest>first) then
                peak=maxval(abs(val(first:dest-1)))
                val(first:dest-1)=val(first:dest-1)/peak
                if (any(val(first:dest-1)==0.0_rk)) then
                    call finish_constraint_matrix(rowptr,col,val,1,compatible=.true.,valid=.false.,field=.true.,err=err)
                    return
                end if
                nrow=nrow+1
                rowptr(nrow+1)=dest
            end if
            first=dest
            lo=hi+1
        end do
        rowptr=rowptr(:nrow+1)
        col=col(:first-1)
        val=val(:first-1)
    end subroutine

    !> Publish the populated CSR or an empty matrix with the reason assembly failed.
    pure subroutine finish_constraint_matrix(rowptr,col,val,pos,powers,compatible,valid,field,err,resolved)
        integer, allocatable, intent(inout) :: rowptr(:), col(:)
        real(rk), allocatable, intent(inout) :: val(:)
        integer, intent(in) :: pos
        integer, optional, allocatable, intent(inout) :: powers(:)
        logical, intent(in) :: compatible, valid, field
        type(debug), optional, intent(inout) :: err
        logical, optional, intent(in) :: resolved
        character(:), allocatable :: message
        logical :: samples_ok

        samples_ok = .true.
        if (present(resolved)) samples_ok = resolved
        if (compatible .and. valid .and. samples_ok) then
            col = col(:pos-1)
            val = val(:pos-1)
            if (present(powers)) powers = powers(:pos-1)
            return
        end if
        rowptr = [1]
        col = [integer ::]
        val = [real(rk) ::]
        if (present(powers)) powers = [integer ::]
        if (.not. present(err)) return
        if (.not. samples_ok) then
            call err%set(code=204,severity=1,category="forcad_multipatch",location="cmp_dof_constraint", &
                message="The element-local Bernstein evaluation of the interface jet is unrepresentable.", &
                suggestion="Check the knot scaling, interface jets and requested derivative order.")
            return
        else if (.not. compatible) then
            message = "Boundary denominator coefficients are not proportional."
        else if (field) then
            message = "The ordinary-field row is unrepresentable or its denominator fails continuity."
        else
            message = "A homogeneous constraint coefficient is outside the working real range."
        end if
        call err%set(code=204,severity=1,category="forcad_multipatch",message=message,location="cmp_dof_constraint", &
            suggestion="Check interface weights and jets. Use field=.true. for ordinary scalar coefficients.")
    end subroutine

    !> Normal basis jets and oriented control IDs of the first participating plane.
    pure subroutine endpoint_stencil(nc, knot, side, order, reverse, swap, flip, plane, stride, length, ders)
        integer, contiguous, intent(in) :: nc(:)
        real(rk), contiguous, intent(in) :: knot(:)
        integer, intent(in) :: side, order
        logical, intent(in) :: reverse, swap, flip(2)
        integer, allocatable, intent(out) :: plane(:)
        integer, intent(out) :: stride
        real(rk), intent(out) :: length
        real(rk), allocatable, intent(out) :: ders(:,:)
        integer :: dir, degree, n, first, layer
        real(rk) :: lo, hi
        logical :: lower

        dir = (side+1)/2
        n = nc(dir)
        degree = size(knot)-n-1
        lo = knot(degree+1)
        hi = knot(n+1)
        length = 1.0_rk
        if (order > 0) length = hi-lo
        lower = mod(side,2) == 1
        allocate(ders(0:order,0:degree))
        call basis_bspline_der_all_active(merge(lo,hi,lower), knot, n, degree, order, first, ders)
        stride = product(nc(:dir-1))
        if (size(nc) == 1) then
            plane = [first]
            return
        end if
        ! On upper edges/faces the public CSR enumerates normal layers inward.
        layer = first-1
        if (.not. lower) then
            layer = n-(first+degree)
            ders = ders(:,degree:0:-1)
            stride = -stride
        end if
        if (size(nc) == 2) then
            call boundary_layer_index(nc, side, layer, reverse, plane)
        else
            call boundary_layer_index(nc, side, layer, reverse, swap, flip, plane)
        end if
    end subroutine

    pure subroutine append_row_side(ders, r, first, stride, col, val, pos)
        real(rk), intent(in) :: ders(0:,0:)
        integer, intent(in) :: r, first, stride
        integer, intent(inout) :: col(:), pos
        real(rk), intent(inout) :: val(:)
        integer :: a

        do a = 0, ubound(ders,2)
            if (.not. structural_nonzero(ders(r,a))) cycle
            col(pos) = first + a*stride
            val(pos) = ders(r,a)
            pos = pos + 1
        end do
    end subroutine

    !> Count positive tensor spans and the maximum number of controls per element.
    pure subroutine element_layout(patches, npatch, rows, width)
        use forcad_utils, only: active_span_count
        class(*), optional, intent(in) :: patches(:)
        integer, intent(in) :: npatch
        integer, allocatable, intent(out) :: rows(:)
        integer, optional, intent(out) :: width
        integer :: i, d, nc(3), degree(3), stride(3), cells

        allocate(rows(npatch+1))
        rows(1)=0
        if (present(width)) width=0
        do i=1,npatch
            call patch_basis(patches(i),nc,degree,stride)
            cells=1
            do d=1,3
                if (stride(d)==0) cycle
                cells=cells*active_span_count(patch_axis_knots(patches(i),d),nc(d),degree(d))
            end do
            rows(i+1)=rows(i)+cells
            if (present(width)) width=max(width,product(degree+1))
        end do
    end subroutine

    !> Assemble all patch elements using the supplied storage-to-solver map, if present.
    pure subroutine multipatch_elements(patches, offsets, elem, map)
        class(*), optional, intent(in) :: patches(:)
        integer, intent(in) :: offsets(:)
        integer, allocatable, intent(out) :: elem(:,:)
        integer, optional, intent(in) :: map(:)
        integer, allocatable :: rows(:)
        integer :: i, width

        call element_layout(patches,size(offsets)-1,rows,width)
        allocate(elem(rows(size(rows)),width),source=0)
        do i=1,size(offsets)-1
            if (rows(i)==rows(i+1)) cycle
            select type (patch=>patches(i))
            type is (nurbs_curve)
                call multipatch_element_block(elem,rows(i)+1,patch%cmp_elem(),offsets(i),map)
            type is (nurbs_surface)
                call multipatch_element_block(elem,rows(i)+1,patch%cmp_elem(),offsets(i),map)
            type is (nurbs_volume)
                call multipatch_element_block(elem,rows(i)+1,patch%cmp_elem(),offsets(i),map)
            end select
        end do
    end subroutine

    !> Expand cumulative element counts into patch owners or one-based patch-local rows.
    pure subroutine multipatch_element_ids(offsets, local, ids)
        integer, intent(in) :: offsets(:)
        logical, intent(in) :: local
        integer, allocatable, intent(out) :: ids(:)
        integer :: patch, first, last, e

        allocate(ids(offsets(size(offsets))))
        do patch=1,size(offsets)-1
            first=offsets(patch)+1
            last=offsets(patch+1)
            if (local) then
                do e=first,last
                    ids(e)=e-offsets(patch)
                end do
            else
                ids(first:last)=patch
            end if
        end do
    end subroutine multipatch_element_ids

    !> Place one patch's exact active supports in the global table; padding stays zero.
    pure subroutine multipatch_element_block(elem,first,local,offset,map)
        integer, intent(inout) :: elem(:,:)
        integer, intent(in) :: first,local(:,:),offset
        integer, optional, intent(in) :: map(:)
        integer :: a,last

        last=first+size(local,1)-1
        do a=1,size(local,2)
            if (present(map)) then
                elem(first:last,a)=map(offset+local(:,a))
            else
                elem(first:last,a)=offset+local(:,a)
            end if
        end do
    end subroutine multipatch_element_block

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the next geometric-growth capacity for patch metadata arrays.
    !! The result is at least 16 and doubles established capacities, giving
    !! amortized constant-time append without changing mathematical data.
    pure elemental integer function multipatch_growth_capacity(current) result(next)
        integer, intent(in) :: current
            !! Current allocation capacity; nonpositive values are allowed.

        next = max(16, 2*max(1, current))
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Store low-level connection metadata and reset omitted options to defaults.
    !!
    !! @warning This routine performs no patch-level validation. Rank-specific
    !! `connect` methods should be used for normal construction.
    !! @endwarning
    pure subroutine set_connection(this, patch_a, side_a, patch_b, side_b, continuity, reverse, swap, flip, &
        geometric, reparameterization, transition_jet, projective_jet)
        class(multipatch_connection), intent(inout) :: this
            !! Connection to replace.
        integer, intent(in) :: patch_a
            !! One-based first patch identifier.
        integer, intent(in) :: side_a
            !! First boundary identifier using rank-specific numbering.
        integer, intent(in) :: patch_b
            !! One-based second patch identifier.
        integer, intent(in) :: side_b
            !! Second boundary identifier using rank-specific numbering.
        integer, intent(in), optional :: continuity
            !! Requested order; defaults to zero.
        logical, intent(in), optional :: reverse
            !! Reverse the mapped boundary orientation; defaults to false.
        logical, intent(in), optional :: swap
            !! Exchange volume-face tangential directions; defaults to false.
        logical, intent(in), optional :: geometric
            !! Select the scalar normal-transition model instead of normalized \(C^n\); defaults to false.
        logical, intent(in), contiguous, optional :: flip(:)
            !! Up to two volume-face tangential flips; omitted entries are false.
        real(rk), intent(in), contiguous, optional :: reparameterization(:)
            !! Normal transition derivatives 1 through \(n\).
        real(rk), intent(in), contiguous, optional :: transition_jet(:,:,:,:)
            !! Full transition-jet Bernstein coefficients.
        real(rk), intent(in), contiguous, optional :: projective_jet(:,:,:)
            !! Projective-factor Bernstein coefficients for orders zero through \(n\).
        integer :: n

        this%patch_a = patch_a
        if (allocated(this%prolongation)) deallocate(this%prolongation)
        this%side_a = side_a
        this%patch_b = patch_b
        this%side_b = side_b
        this%continuity = 0
        this%geometric = .false.
        this%reverse = .false.
        this%swap = .false.
        this%flip = .false.
        if (allocated(this%reparameterization)) deallocate(this%reparameterization)
        if (allocated(this%transition_jet)) deallocate(this%transition_jet)
        if (allocated(this%projective_jet)) deallocate(this%projective_jet)
        if (present(continuity)) this%continuity = continuity
        if (present(geometric)) this%geometric = geometric
        if (present(reparameterization)) this%reparameterization = reparameterization
        if (present(transition_jet)) then
            if (size(transition_jet)>0) this%transition_jet = transition_jet
        end if
        if (present(projective_jet)) then
            if (size(projective_jet)>0) this%projective_jet = projective_jet
        end if
        if (present(reverse)) this%reverse = reverse
        if (present(swap)) this%swap = swap
        if (present(flip)) then
            n = min(2, size(flip))
            this%flip(1:n) = flip(1:n)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute multivariate Faa di Bruno coefficients for one transition jet.
    !!
    !! `transition(c,j)` is the \(j\)-th ordinary derivative of mapped
    !! coordinate \(c\) with respect to the patch-A normal coordinate.
    !! `coefficient(r,a1,a2,a3)` multiplies
    !! \(\partial_1^{a_1}\partial_2^{a_2}\partial_3^{a_3}f\) in the \(r\)-th
    !! derivative of \(f\circ\Phi\). For a surface, pass two coordinate rows;
    !! entries with `a3>0` remain zero. The recurrence is
    !!
    !! \[
    !! C_{r,\boldsymbol\alpha}=
    !! \sum_{c:\alpha_c>0}\sum_{j=1}^{r}
    !! {r-1\choose j-1}\Phi_c^{(j)}
    !! C_{r-j,\boldsymbol\alpha-\boldsymbol e_c},
    !! \qquad C_{0,\boldsymbol0}=1.
    !! \]
    !!
    !! This is the coefficient form of the multivariate Faa di Bruno formula
    !! used by the general surface and volume interface assemblers.
    pure subroutine multipatch_composition_coefficients(transition, coefficient)
        real(rk), intent(in), contiguous :: transition(:,:)
            !! Transition derivatives `[mapped coordinate,normal order]`.
        real(rk), allocatable, intent(out) :: coefficient(:,:,:,:)
            !! Composition coefficients indexed from zero.
        real(rk) :: choose, previous
        integer :: a1,a2,a3,j,r,c,next(3),continuity

        continuity=size(transition,2)
        allocate(coefficient(0:continuity,0:continuity,0:continuity,0:continuity), source=0.0_rk)
        coefficient(0,0,0,0)=1.0_rk
        do r=1,continuity
            choose=1.0_rk
            do j=1,r
                do a3=0,r-j
                    do a2=0,r-j-a3
                        do a1=0,r-j-a3-a2
                            previous=coefficient(r-j,a1,a2,a3)
                            if (previous==0.0_rk) cycle
                            do c=1,min(3,size(transition,1))
                                next=[a1,a2,a3]
                                next(c)=next(c)+1
                                associate(value => coefficient(r,next(1),next(2),next(3)))
                                    value=value+choose*transition(c,j)*previous
                                end associate
                            end do
                        end do
                    end do
                end do
                choose=choose*real(r-j,rk)/real(j,rk)
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the one-based identifier of the first connected patch.
    pure elemental integer function get_connection_patch_a(this) result(patch_id)
        class(multipatch_connection), intent(in) :: this

        patch_id = this%patch_a
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the rank-specific side identifier on the first patch.
    pure elemental integer function get_connection_side_a(this) result(side)
        class(multipatch_connection), intent(in) :: this

        side = this%side_a
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the one-based identifier of the second connected patch.
    pure elemental integer function get_connection_patch_b(this) result(patch_id)
        class(multipatch_connection), intent(in) :: this

        patch_id = this%patch_b
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the rank-specific side identifier on the second patch.
    pure elemental integer function get_connection_side_b(this) result(side)
        class(multipatch_connection), intent(in) :: this

        side = this%side_b
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the requested continuity order.
    !! The value `-1` records topology without continuity equations.
    pure elemental integer function get_connection_continuity(this) result(continuity)
        class(multipatch_connection), intent(in) :: this

        continuity = this%continuity
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return a copy of the normalized normal-coordinate transition derivatives.
    !! An unallocated transition jet is represented by a zero-length result.
    pure function get_connection_reparameterization(this) result(reparameterization)
        class(multipatch_connection), intent(in) :: this
        real(rk), allocatable :: reparameterization(:)

        if (allocated(this%reparameterization)) then
            allocate(reparameterization(size(this%reparameterization)), source=this%reparameterization)
        else
            allocate(reparameterization(0))
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the full multivariate transition-jet Bernstein coefficients.
    !! A surface jet has shape `[q+1,1,2,n]`; a volume jet has shape
    !! `[q1+1,q2+1,3,n]`. An absent jet returns a zero-sized array.
    pure function get_connection_transition_jet(this) result(transition_jet)
        class(multipatch_connection), intent(in) :: this
        real(rk), allocatable :: transition_jet(:,:,:,:)

        if (allocated(this%transition_jet)) then
            transition_jet = this%transition_jet
        else
            allocate(transition_jet(0,0,0,0))
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the projective-factor jet Bernstein coefficients.
    !! The final index stores derivative orders zero through \(n\). An absent
    !! jet returns a zero-sized array.
    pure function get_connection_projective_jet(this) result(projective_jet)
        class(multipatch_connection), intent(in) :: this
        real(rk), allocatable :: projective_jet(:,:,:)

        if (allocated(this%projective_jet)) then
            projective_jet = this%projective_jet
        else
            allocate(projective_jet(0,0,0))
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report whether full multivariate transition data are stored.
    pure elemental logical function connection_has_transition_jet(this) result(has_jet)
        class(multipatch_connection), intent(in) :: this

        has_jet = allocated(this%transition_jet)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report whether a variable projective-factor jet is stored.
    pure elemental logical function connection_has_projective_jet(this) result(has_jet)
        class(multipatch_connection), intent(in) :: this

        has_jet = allocated(this%projective_jet)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report whether the connection requests geometric continuity.
    pure elemental logical function get_connection_geometric(this) result(geometric)
        class(multipatch_connection), intent(in) :: this

        geometric = this%geometric
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report reversal of the mapped tangential ordering.
    pure elemental logical function get_connection_reverse(this) result(reverse)
        class(multipatch_connection), intent(in) :: this

        reverse = this%reverse
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report exchange of the two mapped tangential face directions.
    pure elemental logical function get_connection_swap(this) result(swap)
        class(multipatch_connection), intent(in) :: this

        swap = this%swap
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Report reversal of mapped tangential direction `dir`.
    !! Invalid directions return `.false.`.
    pure elemental logical function get_connection_flip(this, dir) result(flip)
        class(multipatch_connection), intent(in) :: this
        integer, intent(in) :: dir

        flip = .false.
        if (dir < 1 .or. dir > size(this%flip)) return
        flip = this%flip(dir)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Validate parametric, compact geometric, or multivariate geometric metadata.
    !! `parameter_dimension` is two for a surface and three for a volume. When
    !! it is absent, only the scalar curve/compact transition is accepted. A
    !! full transition stores Bernstein coefficients of every mapped coordinate
    !! derivative. Strictly signed first-normal-derivative coefficients provide
    !! a conservative certificate that the interface Jacobian is nonsingular.
    !! The same Bernstein convex-hull certificate is used for a supplied
    !! nonvanishing projective factor.
    pure logical function multipatch_valid_reparameterization(connection, normal_sign, parameter_dimension) result(ok)
        type(multipatch_connection), intent(in) :: connection
            !! Connection metadata to validate.
        integer, intent(in) :: normal_sign
            !! Expected sign of the first normalized derivative.
        integer, intent(in), optional :: parameter_dimension
            !! Number of mapped coordinates for a full transition jet.
        real(rk) :: threshold
        integer :: continuity, dimension
        logical :: full_transition, projective

        continuity = connection%continuity
        full_transition = allocated(connection%transition_jet)
        projective = allocated(connection%projective_jet)
        if (.not. connection%geometric) then
            ok = .not. allocated(connection%reparameterization) .and. .not. full_transition .and. .not. projective
            return
        end if
        if (continuity < 0) then
            ok = .false.
            return
        end if
        if (allocated(connection%reparameterization) .and. full_transition) then
            ok = .false.
            return
        end if

        dimension = 0
        if (present(parameter_dimension)) dimension = parameter_dimension
        if (full_transition) then
            if (dimension < 2 .or. dimension > 3) then
                ok = .false.
                return
            end if
            if (size(connection%transition_jet,1) < 1 .or. size(connection%transition_jet,2) < 1 .or. &
                size(connection%transition_jet,3) /= dimension .or. &
                size(connection%transition_jet,4) /= continuity) then
                ok = .false.
                return
            end if
            if (dimension == 2 .and. size(connection%transition_jet,2) /= 1) then
                ok = .false.
                return
            end if
            if (.not. all(ieee_is_finite(connection%transition_jet))) then
                ok = .false.
                return
            end if
        end if
        if (projective) then
            if (dimension < 2 .or. dimension > 3) then
                ok = .false.
                return
            end if
            if (size(connection%projective_jet,1) < 1 .or. size(connection%projective_jet,2) < 1 .or. &
                size(connection%projective_jet,3) /= continuity + 1) then
                ok = .false.
                return
            end if
            if (dimension == 2 .and. size(connection%projective_jet,2) /= 1) then
                ok = .false.
                return
            end if
            if (.not. all(ieee_is_finite(connection%projective_jet))) then
                ok = .false.
                return
            end if
            if (.not. all(sign(1.0_rk,connection%projective_jet(1,1,1))*&
                connection%projective_jet(:,:,1) > 0.0_rk)) then
                ok = .false.
                return
            end if
        end if

        if (continuity == 0) then
            ok = .not. full_transition
            if (allocated(connection%reparameterization)) ok = ok .and. size(connection%reparameterization) == 0
            return
        end if
        if (.not. allocated(connection%reparameterization) .and. .not. full_transition) then
            ok = .false.
            return
        end if
        threshold = sqrt(tiny(1.0_rk))
        if (full_transition) then
            ok = all(real(normal_sign, rk)*connection%transition_jet(:,:,dimension,1) > threshold)
        else
            ok = size(connection%reparameterization) == continuity
            if (.not. ok) return
            ok = all(ieee_is_finite(connection%reparameterization))
            if (.not. ok) return
            ok = real(normal_sign, rk)*connection%reparameterization(1) > threshold
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute arbitrary-order univariate Faa di Bruno coefficients.
    !! `coefficient(r,k)` multiplies the kth normalized normal derivative on
    !! patch B when forming the rth derivative of its composition with the
    !! connection reparameterization. For C^n, this reduces exactly to the
    !! signed affine transition implied by the connected sides.
    !!
    !! The returned lower-triangular matrix satisfies
    !!
    !! \[
    !! C_{0,0}=1,\qquad
    !! C_{r,k}=\sum_{j=1}^{r-k+1}
    !! {r-1\choose j-1}\phi^{(j)}(0)C_{r-j,k-1}.
    !! \]
    !!
    !! Thus \(C_{r,k}=B_{r,k}\), the partial exponential Bell polynomial used
    !! by the univariate Faa di Bruno formula. This compact routine is used
    !! when a surface or volume connection omits the full multivariate
    !! transition jet; `multipatch_composition_coefficients` handles the
    !! general mapped-coordinate jet.
    pure subroutine multipatch_chain_rule_coefficients(connection, normal_sign, coefficient)
        type(multipatch_connection), intent(in) :: connection
            !! Validated connection whose transition is applied.
        integer, intent(in) :: normal_sign
            !! Orientation sign for a parametric transition.
        real(rk), allocatable, intent(out) :: coefficient(:,:)
            !! Lower-triangular coefficient matrix indexed from zero.
        real(rk) :: choose
        integer :: j, k, r, continuity

        continuity = max(0, connection%continuity)
        allocate(coefficient(0:continuity,0:continuity), source=0.0_rk)
        coefficient(0,0) = 1.0_rk
        if (continuity == 0) return

        if (.not. connection%geometric) then
            do r = 1, continuity
                coefficient(r,r) = real(normal_sign, rk)**r
            end do
            return
        end if

        do r = 1, continuity
            do k = 1, r
                choose = 1.0_rk
                do j = 1, r - k + 1
                    coefficient(r,k) = coefficient(r,k) + &
                        choose*connection%reparameterization(j)*coefficient(r-j,k-1)
                    choose = choose*real(r-j, rk)/real(j, rk)
                end do
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Test whether two normalized univariate trace spaces coincide.
    !!
    !! Degrees and control counts must match. Knot arrays are assumed to be valid
    !! for the supplied degrees and counts. Values are compared after an
    !! affine map from each active interval to \([0,1]\); `reverse` compares the
    !! second trace in opposite orientation. Multiplicities agree exactly. At an
    !! interior break the displacement is measured against each adjacent positive
    !! span. A break of multiplicity p+1 must have the same mapped location.
    !!
    pure logical function multipatch_compatible_trace_space(knot_a, nc_a, degree_a, &
        knot_b, nc_b, degree_b, reverse, tol) result(ok)
        real(rk), intent(in), contiguous :: knot_a(:)
            !! First trace knot vector.
        real(rk), intent(in), contiguous :: knot_b(:)
            !! Second trace knot vector.
        integer, intent(in) :: nc_a
            !! First trace control-point count.
        integer, intent(in) :: nc_b
            !! Second trace control-point count.
        integer, intent(in) :: degree_a
            !! First trace polynomial degree.
        integer, intent(in) :: degree_b
            !! Second trace polynomial degree.
        logical, intent(in) :: reverse
            !! Reverse the normalized second trace before comparison.
        real(rk), intent(in) :: tol
            !! Tolerance relative to the active interval and to incident interior spans.
        real(rk) :: start(2), finish(2), width(2), point(2), residual, eps
        integer :: n, p, nc, first, last, j, step, width_power(2), residual_power

        ok = .false.
        n = size(knot_a)
        p = degree_a
        nc = nc_a
        if (p/=degree_b .or. nc/=nc_b .or. p<0 .or. nc<=p) return
        if (n/=nc+p+1 .or. size(knot_b)/=n) return
        start = [knot_a(p+1), knot_b(merge(nc+1,p+1,reverse))]
        finish = [knot_a(nc+1), knot_b(merge(p+1,nc+1,reverse))]
        if (start(1)>=finish(1) .or. knot_b(p+1)>=knot_b(nc+1)) return
        if (.not.all(ieee_is_finite(knot_a)) .or. .not.all(ieee_is_finite(knot_b)) .or. .not.ieee_is_finite(tol)) return
        if (.not.reverse) then
            ok = all(knot_a==knot_b)
            if (ok) return
        end if
        call scaled_difference(start,finish,width,width_power)
        eps = max(tol,128.0_rk*epsilon(tol))
        step = merge(-1,1,reverse)
        first = 1
        do while (first<=n)
            j = merge(n-first+1,first,reverse)
            last = first
            do while (last<n)
                if ((knot_a(last)==knot_a(last+1)) .neqv. (knot_b(j)==knot_b(j+step))) return
                if (knot_a(last)/=knot_a(last+1)) exit
                last = last+1
                j = j+step
            end do
            point = [knot_a(last),knot_b(j)]
            if (.not.(all(point==start) .or. all(point==finish))) then
                call affine_knot_residual(start,finish,point,residual,residual_power)
                if (residual/=0.0_rk) then
                    if (point(1)>start(1) .and. point(1)<finish(1)) then
                        if (last-first>=p) return
                        if (.not.within_span(point,[knot_a(first-1),knot_b(j-(last-first+1)*step)])) return
                        if (.not.within_span(point,[knot_a(last+1),knot_b(j+step)])) return
                    else if (.not.within_span(start,finish)) then
                        return
                    end if
                end if
            end if
            first = last+1
        end do
        ok = .true.
    contains
        ! For each side require |det| <= tol * (local span) * (other domain width).
        pure logical function within_span(left,right) result(matches)
            real(rk), intent(in) :: left(2),right(2)
            real(rk) :: span(2),limit(2)
            integer :: power(2)
            call scaled_difference(left,right,span,power)
            limit = abs(fraction(span)*fraction(width(2:1:-1)))*fraction(eps)
            power = power+exponent(span)+width_power(2:1:-1)+exponent(width(2:1:-1))+exponent(eps)+exponent(limit)
            matches = all(residual_power<power .or. (residual_power==power .and. residual<=fraction(limit)))
        end function
    end function
    !===============================================================================


    !> Prefix counts in patch-storage numbering, including exterior controls.
    pure function multipatch_dof_offsets(patches,npatch) result(offsets)
        class(*), optional, intent(in) :: patches(:)
        integer, intent(in) :: npatch
        integer, allocatable :: offsets(:)
        integer :: i,nc(3),degree(3)
        allocate(offsets(npatch+1))
        offsets(1)=0
        do i=1,npatch
            call patch_basis(patches(i),nc,degree)
            offsets(i+1)=offsets(i)+product(nc)
        end do
    end function

    pure function patch_control_counts(patch) result(nc)
        class(*), intent(in) :: patch
        integer, allocatable :: nc(:)
        select type (patch)
        type is (nurbs_curve)
            nc = [patch%get_nc()]
        type is (nurbs_surface)
            nc = patch%get_nc()
        type is (nurbs_volume)
            nc = patch%get_nc()
        end select
    end function

    pure function patch_axis_knots(patch, dir) result(knot)
        class(*), intent(in) :: patch
        integer, intent(in) :: dir
        real(rk), allocatable :: knot(:)
        select type (patch)
        type is (nurbs_curve)
            knot = patch%get_knot()
        type is (nurbs_surface)
            if (dir==3) then
                knot = [0.0_rk,1.0_rk] ! Constant second tangent in the edge Bernstein assembler.
            else
                knot = patch%get_knot(dir)
            end if
        type is (nurbs_volume)
            knot = patch%get_knot(dir)
        end select
    end function

    pure logical function patch_is_rational(patch) result(rational)
        class(*), intent(in) :: patch
        rational = .false.
        select type (patch)
        type is (nurbs_curve)
            rational = patch%is_rational()
        type is (nurbs_surface)
            rational = patch%is_rational()
        type is (nurbs_volume)
            rational = patch%is_rational()
        end select
    end function

    pure function patch_weights(patch, ids) result(weights)
        class(*), intent(in) :: patch
        integer, intent(in) :: ids(:)
        real(rk) :: weights(size(ids))
        weights = 1.0_rk
        select type (patch)
        type is (nurbs_curve)
            if (patch%is_rational()) weights = patch%get_Wc(ids)
        type is (nurbs_surface)
            if (patch%is_rational()) weights = patch%get_Wc(ids)
        type is (nurbs_volume)
            if (patch%is_rational()) weights = patch%get_Wc(ids)
        end select
    end function

    !> Number periodic and directly shareable controls in patch-storage order.
    !! The disjoint sets use path compression.
    !! [Tarjan (1975)](https://doi.org/10.1145/321879.321884).
    pure subroutine multipatch_dof_map(patches, connections, nconnection, offsets, active_only, map, shareable, include_aliases)
        class(*), optional, intent(in) :: patches(:)
        type(multipatch_connection), optional, intent(in) :: connections(:)
        integer, intent(in) :: nconnection, offsets(:)
        logical, optional, intent(in) :: active_only, include_aliases
        integer, allocatable, intent(out) :: map(:)
        logical, allocatable, optional, intent(out) :: shareable(:)
        integer, allocatable :: parent(:), a(:), b(:), nc(:)
        logical :: merged
        logical, allocatable :: active(:)
        integer :: i, d, pa, pb

        if (present(shareable)) allocate(shareable(nconnection),source=.false.)
        allocate(parent(offsets(size(offsets))))
        do i=1,size(parent)
            parent(i)=i
        end do
        if (requested(active_only)) allocate(active(size(parent)),source=.false.)
        do i=1,size(offsets)-1
            if (.not.patch_is_valid(patches(i))) cycle
            select type (patch => patches(i))
            type is (nurbs_curve)
                call seed_patch_dofs(parent,offsets(i),patch%cmp_dof_map())
            type is (nurbs_surface)
                call seed_patch_dofs(parent,offsets(i),patch%cmp_dof_map())
            type is (nurbs_volume)
                call seed_patch_dofs(parent,offsets(i),patch%cmp_dof_map())
            end select
            if (.not.allocated(active)) cycle
            active(offsets(i)+1:offsets(i+1)) = .true.
            nc = patch_control_counts(patches(i))
            do d=1,size(nc)
                call restrict_active_axis(active(offsets(i)+1:offsets(i+1)),nc,d,patch_axis_knots(patches(i),d))
            end do
        end do
        do i=1,nconnection
            associate(conn => connections(i))
                if (conn%continuity<0 .or. allocated(conn%prolongation)) cycle
                pa = conn%patch_a
                pb = conn%patch_b
                a = interpolated_trace_ids(patch_control_counts(patches(pa)), &
                    patch_axis_knots(patches(pa),(conn%side_a+1)/2),conn,.false.)
                if (size(a)==0) cycle
                b = interpolated_trace_ids(patch_control_counts(patches(pb)), &
                    patch_axis_knots(patches(pb),(conn%side_b+1)/2),conn,.true.)
                if (size(a)/=size(b)) cycle
                if (size(a)>1 .and. (patch_is_rational(patches(pa)) .or. patch_is_rational(patches(pb)))) then
                    call merge_trace_dofs(parent,offsets,conn,a,b,patch_weights(patches(pa),a),patch_weights(patches(pb),b),merged)
                else
                    call merge_trace_dofs(parent,offsets,conn,a,b,merged=merged)
                end if
                if (present(shareable)) shareable(i)=merged
            end associate
        end do
        call disjoint_set_map(parent,map,active,include_aliases)
    end subroutine

    !> Assemble a scalar objective gradient in independent control positions and log weights.
    !! Local partials use patch storage order. A nested trace couples fine positions to coarse
    !! weights, so supply its coordinate partial even when requesting only the weight gradient.
    !! Positive map entries identify direct copies. Blended fine controls and inactive controls
    !! without an active representative have ID zero.
    pure subroutine assemble_design_gradient(patches,connections,nconnection,offsets,valid, &
            local_dXc,local_dlogWc,dXc,dlogWc,err,map)
        class(*), optional, intent(in) :: patches(:)
        type(multipatch_connection), optional, intent(in) :: connections(:)
        integer, intent(in) :: nconnection,offsets(:)
        logical, intent(in) :: valid
        real(rk), optional, intent(in) :: local_dXc(:,:),local_dlogWc(:)
        real(rk), optional, allocatable, intent(out) :: dXc(:,:),dlogWc(:)
        type(debug), optional, intent(out) :: err
        integer, optional, allocatable, intent(out) :: map(:)
        type(debug) :: status
        integer, allocatable :: ids(:)
        real(rk), allocatable :: gx(:,:),gw(:)
        integer :: i,d,width(size(offsets)-1),nactive,nstorage
        logical :: ok,nested,have_x,have_w,want_x,want_w

        have_x=present(local_dXc)
        have_w=present(local_dlogWc)
        want_x=present(dXc)
        want_w=present(dlogWc)
        nstorage=offsets(size(offsets))
        call status%reset()
        if (.not.valid) then
            call status%set(code=205,severity=1,category="forcad_multipatch",location="assemble_sensitivities", &
                message="Invalid multipatch geometry.",suggestion="Check the geometry and interfaces.")
        else if ((want_x.and..not.have_x) .or. (have_w.and..not.want_w) .or. &
                (want_w.and..not.(have_x.or.have_w)) .or. (have_x.and..not.(want_x.or.want_w)) .or. &
                .not.(want_x.or.want_w.or.present(map))) then
            call status%set(code=208,severity=1,category="forcad_multipatch",location="assemble_sensitivities", &
                message="Inconsistent local partials and requested design gradients.", &
                suggestion="Supply local partials for the requested gradient, or request the map alone.")
        else
            ok=.true.
            if (present(local_dlogWc)) ok=size(local_dlogWc)==nstorage.and.all(ieee_is_finite(local_dlogWc))
            if (present(local_dXc)) then
                width=[(patch_ncoord(patches(i)),i=1,size(width))]
                ok=ok.and.all(shape(local_dXc)==[nstorage,maxval(width)]).and.all(ieee_is_finite(local_dXc))
                if (ok) then
                    do i=1,size(offsets)-1
                        d=width(i)
                        if (any(local_dXc(offsets(i)+1:offsets(i+1),d+1:)/=0.0_rk)) ok=.false.
                    end do
                end if
            end if
            if (.not.ok) then
                call status%set(code=208,severity=1,category="forcad_multipatch",location="assemble_sensitivities", &
                    message="Invalid local derivative shape, padding or value.", &
                    suggestion="Use finite storage-ordered arrays and zero padding for absent coordinates.")
            else
                call multipatch_design_map(patches,connections,nconnection,offsets,ids,status, &
                    want_w.or..not.want_x,nested)
            end if
        end if
        if (status%ok) then
            nactive=max(0,maxval(ids))
            ok=.true.
            if (have_x.and.(want_x.or.nested)) then
                allocate(gx(nactive,size(local_dXc,2)))
                do d=1,size(gx,2)
                    call sum_design_column(ids,local_dXc(:,d),gx(:,d),ok)
                    if (.not.ok) exit
                end do
            end if
            if (want_w.and.ok) then
                allocate(gw(nactive),source=0.0_rk)
                if (present(local_dlogWc)) call sum_design_column(ids,local_dlogWc,gw,ok)
            end if
            if (.not.ok) call status%set(code=208,severity=1,category="forcad_multipatch",location="assemble_sensitivities", &
                message="An assembled geometry derivative exceeds the working real range.",suggestion="Rescale the objective.")
            if (ok.and.nested) call nested_design_pullback(patches,connections,nconnection,offsets,ids,gx,gw,status)
        end if
        if (status%ok) then
            if (present(dXc)) call move_alloc(gx,dXc)
            if (present(dlogWc)) call move_alloc(gw,dlogWc)
            if (present(map)) call move_alloc(ids,map)
        else
            if (present(dXc)) allocate(dXc(0,0))
            if (present(dlogWc)) dlogWc=[real(rk)::]
            if (present(map)) map=[integer::]
        end if
        if (present(err)) err=status
    end subroutine

    !> Compensated scatter for one control coordinate or the log-weight column.
    pure subroutine sum_design_column(ids,local,global,ok)
        integer, intent(in) :: ids(:)
        real(rk), intent(in) :: local(:)
        real(rk), intent(out) :: global(:)
        logical, intent(out) :: ok
        real(rk) :: correction(size(global)),total,tail,remainder
        integer :: i,id

        ! The absolute sum stays within range on the fast path, in any storage order.
        if (maxval(abs(local))>=huge(total)/(2.0_rk*max(1,size(local)))) then
            call sum_scaled_column(ids,local,global,ok)
            return
        end if
        global=0.0_rk
        correction=0.0_rk
        do i=1,size(ids)
            id=ids(i)
            if (id==0) cycle
            call two_sum(global(id),local(i),total,tail)
            global(id)=total
            if (tail==0.0_rk) cycle
            call two_sum(correction(id),tail,total,remainder)
            if (remainder/=0.0_rk) then
                call sum_scaled_column(ids,local,global,ok)
                return
            end if
            correction(id)=total
        end do
        global=global+correction
        ok=all(ieee_is_finite(global))
    end subroutine

    !> For extreme partials keep separate radix exponents until a complete class has been summed.
    pure subroutine sum_scaled_column(ids,local,global,ok)
        integer, intent(in) :: ids(:)
        real(rk), intent(in) :: local(:)
        real(rk), intent(out) :: global(:)
        logical, intent(out) :: ok
        integer :: head(size(global)),next(size(ids)),count(size(global)),id,i,n
        integer, allocatable :: powers(:)
        real(rk), allocatable :: parts(:)

        head=0
        count=0
        do i=1,size(ids)
            id=ids(i)
            if (id==0.or.local(i)==0.0_rk) cycle
            next(i)=head(id)
            head(id)=i
            count(id)=count(id)+1
        end do
        allocate(parts(maxval(count)),powers(maxval(count)))
        global=0.0_rk
        ok=.false.
        do id=1,size(global)
            n=0
            i=head(id)
            do while (i/=0)
                call add_scaled_term(local(i),0,parts,powers,n)
                i=next(i)
            end do
            call scaled_sum_value(parts,powers,n,1.0_rk,global(id),ok)
            if (.not.ok) return
        end do
        ok=.true.
    end subroutine

    !> At a directly shared C0 trace, add the same displacement and log-weight increment
    !! at each paired control. This preserves their positions and the boundary weight ratio.
    !! Inactive copies of participating periodic controls must follow the same variation.
    pure subroutine multipatch_design_map(patches,connections,nconnection,offsets,map,err,with_weights,nested)
        class(*), optional, intent(in) :: patches(:)
        type(multipatch_connection), optional, intent(in) :: connections(:)
        integer, intent(in) :: nconnection,offsets(:)
        integer, allocatable, intent(out) :: map(:)
        type(debug), intent(out) :: err
        logical, intent(in) :: with_weights
        logical, intent(out) :: nested
        logical, allocatable :: shareable(:)
        integer :: i

        call err%reset()
        nested=.false.
        call multipatch_dof_map(patches,connections,nconnection,offsets,.true.,map,shareable,include_aliases=.true.)
        do i=1,nconnection
            associate(conn=>connections(i))
                if (conn%continuity<0) cycle
                if (allocated(conn%prolongation)) then
                    nested=.true.
                    cycle
                end if
                if (conn%continuity>0 .or. .not.shareable(i)) then
                    call err%set(code=204,severity=1,category="forcad_multipatch",location="assemble_sensitivities", &
                        message="Independent geometry variables require matching or nested C0 interfaces.", &
                        suggestion="Use the linearized interface equations for constrained geometry variations.")
                else if (with_weights .and. conn%has_projective_jet()) then
                    if (patch_is_rational(patches(conn%patch_a)) .and. patch_is_rational(patches(conn%patch_b))) cycle
                    call err%set(code=204,severity=1,category="forcad_multipatch",location="assemble_sensitivities", &
                        message="Varying polynomial weights can change the gauge of a prescribed projective jet.", &
                        suggestion="Request coordinate sensitivities alone, or use a join with an implicit weight scale.")
                end if
                if (.not.err%ok) then
                    map=[integer ::]
                    return
                end if
            end associate
        end do
    end subroutine

    !> Collapse single-parent rows, then apply the transpose of each remaining rational refinement row.
    !! In reverse topological order every child's complete adjoint is available before its parents.
    pure subroutine nested_design_pullback(patches,connections,nconnection,offsets,ids,gx,gw,err)
        class(*), intent(in) :: patches(:)
        type(multipatch_connection), intent(in) :: connections(:)
        integer, intent(in) :: nconnection,offsets(:)
        integer, allocatable, intent(inout) :: ids(:)
        real(rk), allocatable, intent(inout) :: gx(:,:),gw(:)
        type(debug), intent(inout) :: err
        type design_row
            integer :: target,coarse,interface,first,needed
            integer, allocatable :: parent(:)
            real(rk), allocatable :: alpha(:)
        end type
        type(design_row), allocatable :: rows(:)
        integer, allocatable :: group(:),forest(:),owner(:),consumer(:),order(:),ptr(:),cursor(:),ones(:),roots(:)
        real(rk), allocatable :: weights(:),adj(:),x(:,:),terms(:),paths(:,:)
        integer :: i,j,k,r,h,c,f,n,nc,nrow,lo,hi,source,target,head,tail,d,ndim,first,last
        real(rk) :: denom,value(1)
        integer :: power,maxrow
        logical :: changed,ok,weight

        nc=maxval(ids)
        nrow=0
        maxrow=0
        do i=1,nconnection
            if (allocated(connections(i)%prolongation)) nrow=nrow+size(connections(i)%prolongation%fine)
        end do
        allocate(rows(nrow))
        r=0
        do i=1,nconnection
            if (.not.allocated(connections(i)%prolongation)) cycle
            associate(conn=>connections(i),t=>connections(i)%prolongation)
                c=merge(conn%patch_a,conn%patch_b,t%coarse_a)
                f=merge(conn%patch_b,conn%patch_a,t%coarse_a)
                do j=1,size(t%fine)
                    r=r+1
                    lo=t%ptr(j)
                    hi=t%ptr(j+1)-1
                    associate(row=>rows(r))
                        row%coarse=c
                        row%interface=i
                        row%first=lo
                        row%target=ids(offsets(f)+t%fine(j))
                        row%parent=ids(offsets(c)+t%coarse(lo:hi))
                        if (lo>hi.or.row%target<=0.or.any(row%parent<=0)) goto 204
                        if (hi>lo) then
                            weights=patch_weights(patches(c),t%coarse(lo:hi))
                            call weighted_trace_sum(t%val(lo:hi),weights,denom,power)
                            row%alpha=scale(fraction(t%val(lo:hi))*fraction(weights)/denom, &
                                exponent(t%val(lo:hi))+exponent(weights)-power)
                        end if
                    end associate
                    maxrow=max(maxrow,hi-lo+1)
                end do
            end associate
        end do

        group=[(i,i=1,nc)]
        forest=group
        do
            changed=.false.
            do r=1,nrow
                associate(row=>rows(r))
                    target=group(row%target)
                    source=group(row%parent(1))
                    if (target==source .or. .not.all(group(row%parent)==source)) cycle
                    call disjoint_set_union(forest,row%target,row%parent(1))
                    changed=.true.
                end associate
            end do
            if (.not.changed) exit
            call disjoint_set_map(forest,group)
        end do

        n=maxval(group)
        allocate(owner(n),order(n),ptr(n+1),source=0)
        ptr(1)=1
        do r=1,nrow
            associate(row=>rows(r))
                row%target=group(row%target)
                row%parent=group(row%parent)
                if (all(row%parent==row%target)) cycle
                if (any(row%alpha==0.0_rk)) goto 208
                owner(row%target)=-1
                row%needed=size(row%parent)
                do j=1,size(row%parent)
                    source=row%parent(j)
                    ptr(source+1)=ptr(source+1)+1
                end do
            end associate
        end do
        ! Each class stores its local partial followed by the contributions from rows that use it.
        do i=1,n
            if (ptr(i+1)>=huge(i)-ptr(i)) goto 208
            ptr(i+1)=ptr(i)+1+ptr(i+1)
        end do
        allocate(consumer(ptr(n+1)-1))
        cursor=ptr(:n)
        do r=1,nrow
            associate(row=>rows(r))
                if (all(row%parent==row%target)) cycle
                do j=1,size(row%parent)
                    source=row%parent(j)
                    cursor(source)=cursor(source)+1
                    consumer(cursor(source))=r
                end do
            end associate
        end do

        ! Resolve a junction from the first stencil whose parents are all known.
        roots=pack([(i,i=1,n)],owner==0)
        tail=size(roots)
        order(:tail)=roots
        head=0
        do while (head<tail)
            head=head+1
            i=order(head)
            do j=ptr(i)+1,ptr(i+1)-1
                r=consumer(j)
                target=rows(r)%target
                if (owner(target)>0) cycle
                rows(r)%needed=rows(r)%needed-1
                if (rows(r)%needed>0) cycle
                owner(target)=r
                tail=tail+1
                order(tail)=target
            end do
        end do
        if (tail/=n) goto 204
        deallocate(consumer)
        order=order(n:1:-1)

        ! At a patch junction two different insertion paths may define the same control.
        ! Compare their action on independent controls only when the immediate stencils differ.
        do r=1,nrow
            target=rows(r)%target
            h=owner(target)
            if (h==r.or.all(rows(r)%parent==target)) cycle
            if (same_design_row(rows(h)%parent,rows(h)%alpha,rows(r)%parent,rows(r)%alpha)) cycle
            if (.not.allocated(paths)) allocate(paths(n,2))
            paths=0.0_rk
            paths(target,1)=1.0_rk
            do j=1,size(rows(r)%parent)
                source=rows(r)%parent(j)
                paths(source,2)=paths(source,2)+rows(r)%alpha(j)
            end do
            do k=1,n
                i=order(k)
                h=owner(i)
                if (h==0.or.all(paths(i,:)==0.0_rk)) cycle
                do j=1,size(rows(h)%parent)
                    source=rows(h)%parent(j)
                    paths(source,:)=paths(source,:)+rows(h)%alpha(j)*paths(i,:)
                end do
            end do
            if (any(abs(paths(roots,1)-paths(roots,2))>128*epsilon(denom)*max(paths(roots,1),paths(roots,2)))) goto 204
        end do
        ndim=0
        if (allocated(gx)) ndim=size(gx,2)
        if (ndim>0.or.allocated(gw)) then
            allocate(x(n,ndim),adj(n),terms(ptr(n+1)-1),ones(maxval(ptr(2:)-ptr(:n))))
            ones=1
            if (allocated(weights)) deallocate(weights)
            allocate(weights(maxrow))
        end if
        do d=1,ndim+merge(1,0,allocated(gw))
            weight=d>ndim
            if (weight) then
                call sum_design_column(group,gw,adj,ok)
            else
                call sum_design_column(group,gx(:,d),adj,ok)
            end if
            if (.not.ok) goto 208
            cursor=ptr(:n)
            terms(cursor)=adj
            do k=1,n
                target=order(k)
                first=ptr(target)
                last=cursor(target)
                if (last>first) then
                    call sum_design_column(ones(:last-first+1),terms(first:last),value,ok)
                    if (.not.ok) goto 208
                    adj(target)=value(1)
                end if
                r=owner(target)
                if (r==0) cycle
                associate(row=>rows(r))
                    if (weight.and.ndim>0) then
                        associate(t=>connections(row%interface)%prolongation)
                            lo=row%first
                            hi=lo+size(row%parent)-1
                            call nested_row_weight_terms(patches(row%coarse),t%coarse(lo:hi),row%alpha,x(target,:),adj(target), &
                                weights(:size(row%parent)),ok)
                        end associate
                        if (.not.ok) goto 208
                    else
                        weights(:size(row%parent))=row%alpha*adj(target)
                    end if
                    do j=1,size(row%parent)
                        source=row%parent(j)
                        cursor(source)=cursor(source)+1
                        terms(cursor(source))=weights(j)
                    end do
                end associate
            end do
            if (.not.weight) x(:,d)=adj
        end do
        if (allocated(gx)) gx=x(roots,:)
        if (allocated(gw)) gw=adj(roots)
        owner=0
        owner(roots)=[(i,i=1,size(roots))]
        do i=1,size(ids)
            if (ids(i)>0) ids(i)=owner(group(ids(i)))
        end do
        return
204     call err%set(code=204,severity=1,category="forcad_multipatch",location="assemble_sensitivities", &
            message="Nested control dependencies have a cycle or incompatible definitions.", &
            suggestion="Supply nested interfaces with a consistent coarse-to-fine design hierarchy.")
        return
208     call err%set(code=208,severity=1,category="forcad_multipatch",location="assemble_sensitivities", &
            message="A nested design coefficient or adjoint is outside the working real range.", &
            suggestion="Rescale the weights or objective.")
    end subroutine

    !> Different orientations or periodic aliases may give the same row in class numbering.
    pure logical function same_design_row(a,wa,b,wb) result(ok)
        integer, intent(in) :: a(:),b(:)
        real(rk), intent(in) :: wa(:),wb(:)
        integer :: i
        real(rk) :: lhs,rhs
        ok=.false.
        do i=1,size(b)
            if (.not.any(a==b(i))) return
        end do
        do i=1,size(a)
            lhs=sum(wa,mask=a==a(i))
            rhs=sum(wb,mask=b==a(i))
            if (abs(lhs-rhs)>128*epsilon(lhs)*max(lhs,rhs)) return
        end do
        ok=.true.
    end function

    !> A degree-p basis participates when its support (U_i,U_{i+p+1}) overlaps the active domain.
    pure subroutine restrict_active_axis(active,nc,dir,knot)
        logical, intent(inout) :: active(:)
        integer, intent(in) :: nc(:), dir
        real(rk), intent(in) :: knot(:)
        integer :: n, p, lower, upper, stride, base
        n = nc(dir)
        p = size(knot)-n-1
        lower = count(knot(p+2:n+p+1)<=knot(p+1))
        upper = count(knot(:n)<knot(n+1))
        if (lower==0 .and. upper==n) return
        stride = product(nc(:dir-1))
        do base=0,size(active)-1,n*stride
            active(base+1:base+lower*stride) = .false.
            active(base+upper*stride+1:base+n*stride) = .false.
        end do
    end subroutine

    !> Seed a patch's disjoint sets from its independent, compact local numbering.
    pure subroutine seed_patch_dofs(parent, offset, local)
        integer, intent(inout) :: parent(:)
        integer, intent(in) :: offset, local(:)
        integer, allocatable :: first(:)
        integer :: i, id, n

        n = maxval(local)
        if (n<=0 .or. n==size(local)) return
        allocate(first(n),source=0)
        do i = 1, size(local)
            id = local(i)
            if (first(id)==0) first(id) = offset+i
            parent(offset+i) = first(id)
        end do
    end subroutine

    !> Oriented control IDs only when the endpoint interpolates a single normal layer.
    pure function interpolated_trace_ids(nc, knot, conn, second, native) result(ids)
        integer, contiguous, intent(in) :: nc(:)
        real(rk), contiguous, intent(in) :: knot(:)
        type(multipatch_connection), intent(in) :: conn
        logical, intent(in) :: second
        logical, optional, intent(in) :: native !! Keep the patch's storage order for a separate orientation map.
        integer, allocatable :: ids(:)
        integer :: side, n, layer
        logical :: orient

        orient = second.and..not.requested(native)
        side = merge(conn%side_b,conn%side_a,second)
        n = nc((side+1)/2)
        layer = interpolatory_boundary_layer(knot,n,size(knot)-n-1,side)
        if (layer<0) then
            allocate(ids(0))
        else
            select case (size(nc))
            case (1)
                call boundary_layer_index(nc(1),side,layer,.false.,ids)
            case (2)
                call boundary_layer_index(nc,side,layer,orient.and.conn%reverse,ids)
            case (3)
                call boundary_layer_index(nc,side,layer,orient.and.conn%reverse,orient.and.conn%swap,orient.and.conn%flip,ids)
            case default
                allocate(ids(0))
            end select
        end if
    end function

    !> Merge paired interpolatory controls with matching polynomial or rational trace bases.
    pure subroutine merge_trace_dofs(parent, offsets, conn, a, b, Wa, Wb, merged)
        integer, contiguous, intent(inout) :: parent(:)
        integer, intent(in) :: offsets(:), a(:), b(:)
        type(multipatch_connection), intent(in) :: conn
        real(rk), contiguous, optional, intent(in) :: Wa(:), Wb(:)
        logical, optional, intent(out) :: merged
        real(rk) :: ratio
        integer :: i, power
        logical :: compatible

        if (present(merged)) merged=.false.
        if (size(a)/=size(b) .or. (present(Wa).neqv.present(Wb))) return
        if (present(Wa)) then
            call multipatch_projective_weight_scale(Wa,Wb,128.0_rk*epsilon(1.0_rk),ratio,compatible,power)
            if (.not.compatible) return
        end if
        if (present(merged)) merged=.true.
        do i = 1, size(a)
            call disjoint_set_union(parent,offsets(conn%patch_a)+a(i),offsets(conn%patch_b)+b(i))
        end do
    end subroutine

    !> Compare the oriented denominator coefficients W_A = factor * W_B.
    !! Polynomial patches have unit denominator, including stored constant weights.
    !! Only active normal layers enter the endpoint contraction on each rational side.
    pure subroutine connection_weight_scale(patches, conn, tol, weight_scale, ok, scale_power)
        class(*), intent(in) :: patches(:)
        type(multipatch_connection), intent(in) :: conn
        real(rk), intent(in) :: tol
        real(rk), intent(out) :: weight_scale
        logical, intent(out) :: ok
        integer, intent(out) :: scale_power
        integer :: k, id, dir, j, n, lo, hi
        integer, allocatable :: plane(:), offsets(:), powers(:,:)
        real(rk), allocatable :: basis(:,:), trace(:,:)

        weight_scale = 1.0_rk
        scale_power = 0
        ok = .true.
        do k = 1, 2
            id = merge(conn%patch_a,conn%patch_b,k==1)
            dir = (merge(conn%side_a,conn%side_b,k==1)+1)/2
            if (.not.patch_is_rational(patches(id))) cycle
            if (allocated(conn%prolongation)) then
                n=size(conn%prolongation%fine)
            else
                call trace_weight_stencil(patch_control_counts(patches(id)),patch_axis_knots(patches(id),dir), &
                    conn,k==2,plane,offsets,basis)
                n=size(plane)
            end if
            if (.not.allocated(trace)) then
                allocate(trace(n,2),source=1.0_rk)
                allocate(powers(n,2),source=0)
            end if
            if (n/=size(trace,1)) then
                ok = .false.
                return
            end if
            if (allocated(conn%prolongation)) then
                associate(transfer=>conn%prolongation)
                    if (transfer%coarse_a.eqv.(k==1)) then
                        ! Compare T*w_coarse with w_fine in the same fine basis.
                        do j=1,n
                            lo=transfer%ptr(j)
                            hi=transfer%ptr(j+1)-1
                            call weighted_trace_sum(transfer%val(lo:hi),patch_weights(patches(id),transfer%coarse(lo:hi)), &
                                trace(j,k),powers(j,k))
                        end do
                    else
                        trace(:,k)=patch_weights(patches(id),transfer%fine)
                    end if
                end associate
            else
                do j=1,n
                    call weighted_trace_sum(basis(0,:),patch_weights(patches(id),plane(j)+offsets),trace(j,k),powers(j,k))
                end do
            end if
        end do
        if (allocated(trace)) call multipatch_projective_weight_scale(trace(:,1),trace(:,2),tol,weight_scale,ok,scale_power,powers)
    end subroutine

    !> Select the oriented normal basis and control offsets for one boundary trace.
    pure subroutine trace_weight_stencil(nc, knot, conn, second, plane, offsets, basis)
        integer, contiguous, intent(in) :: nc(:)
        real(rk), contiguous, intent(in) :: knot(:)
        type(multipatch_connection), intent(in) :: conn
        logical, intent(in) :: second
        integer, allocatable, intent(out) :: plane(:), offsets(:)
        real(rk), allocatable, intent(out) :: basis(:,:)
        real(rk) :: length
        integer :: j, stride

        call endpoint_stencil(nc, knot, merge(conn%side_b,conn%side_a,second), 0, &
            second.and.conn%reverse, second.and.conn%swap, second.and.conn%flip, plane, stride, length, basis)
        offsets = [(j*stride,j=0,size(basis,2)-1)]
    end subroutine

    !> Apply a factor before returning the leading components of a signed expansion to rk.
    pure subroutine scaled_sum_value(parts,powers,n,factor,value,ok)
        real(rk), intent(in) :: parts(:),factor
        integer, intent(in) :: powers(:),n
        real(rk), intent(out) :: value
        logical, intent(out) :: ok
        integer :: power
        ok=.true.
        value=0.0_rk
        if (n==0.or.factor==0.0_rk) return
        value=parts(n)
        if (n>1) then
            if (powers(n)-powers(n-1)<=digits(value)+1) value=value+scale(parts(n-1),powers(n-1)-powers(n))
        end if
        value=value*fraction(factor)
        power=powers(n)+exponent(factor)
        ok=power+exponent(value)<=maxexponent(value)
        if (ok) value=scale(value,power)
    end subroutine

    !> alpha_j [bar_eta + bar_x . (P_j - sum_k alpha_k P_k)] for one refined control.
    pure subroutine nested_row_weight_terms(patch,local,alpha,bar_x,bar_eta,terms,ok)
        class(*), intent(in) :: patch
        integer, intent(in) :: local(:)
        real(rk), intent(in) :: alpha(:),bar_x(:),bar_eta
        real(rk), intent(out) :: terms(:)
        logical, intent(out) :: ok
        real(rk) :: points(size(local)),delta(size(local),size(bar_x))
        real(rk) :: center(size(delta)+1),parts(size(center)+size(bar_x))
        integer :: shifts(size(local),size(bar_x)),center_power(size(center)),powers(size(parts))
        integer :: i,d,n,nd,ncenter,anchor

        ok=.true.
        if (size(local)==1.or.all(bar_x==0.0_rk)) then
            terms=alpha*bar_eta
            return
        end if
        ok=.false.
        nd=min(size(bar_x),patch_ncoord(patch))
        anchor=maxloc(alpha,dim=1)
        ncenter=0
        call add_scaled_term(bar_eta,0,center,center_power,ncenter)
        ! Project displacements from a common origin before subtracting their weighted mean.
        do d=1,nd
            select type(p=>patch)
            type is(nurbs_surface)
                points=p%get_Xc(local,d)
            type is(nurbs_volume)
                points=p%get_Xc(local,d)
            class default
                return
            end select
            call scaled_difference(points(anchor),points,delta(:,d),shifts(:,d))
            shifts(:,d)=shifts(:,d)+exponent(delta(:,d))+exponent(bar_x(d))
            delta(:,d)=fraction(delta(:,d))*fraction(bar_x(d))
            do i=1,size(local)
                call add_scaled_term(-fraction(alpha(i))*delta(i,d),exponent(alpha(i))+shifts(i,d), &
                    center,center_power,ncenter)
            end do
        end do
        do i=1,size(local)
            n=ncenter
            parts(:n)=center(:n)
            powers(:n)=center_power(:n)
            do d=1,nd
                call add_scaled_term(delta(i,d),shifts(i,d),parts,powers,n)
            end do
            call scaled_sum_value(parts,powers,n,alpha(i),terms(i),ok)
            if (.not.ok) return
        end do
    end subroutine

    !> Positive endpoint denominator sum, returned in a power-of-radix gauge.
    pure subroutine weighted_trace_sum(basis, weights, value, power)
        real(rk), intent(in) :: basis(:), weights(:)
        real(rk), intent(out) :: value
        integer, intent(out) :: power
        real(rk) :: terms(size(basis))
        integer :: powers(size(basis))

        value = 0.0_rk
        power = 0
        if (.not.all(ieee_is_finite(basis)) .or. any(basis<0.0_rk) .or. all(basis==0.0_rk)) return
        terms = fraction(basis)*fraction(weights)
        powers = exponent(basis)+exponent(weights)
        power = maxval(powers,mask=basis>0.0_rk)
        value = sum(scale(terms,powers-power))
    end subroutine

    !> Compare positive denominator coefficients in a common exponent for each ratio.
    pure subroutine multipatch_projective_weight_scale(Wa, Wb, tol, weight_scale, ok, scale_power, powers)
        real(rk), contiguous, intent(in) :: Wa(:), Wb(:)
        real(rk), intent(in) :: tol
        real(rk), intent(out) :: weight_scale
        logical, intent(out) :: ok
        integer, optional, intent(out) :: scale_power
            !! The factor is weight_scale * radix**scale_power. Without this output a scalar factor must fit.
        integer, optional, intent(in) :: powers(:,:)
            !! Binary exponents of the two oriented trace coefficient arrays when kept in scaled form.
        integer :: i, power, ref_power, common
        real(rk) :: f, ref, a, b, rtol

        weight_scale = 1.0_rk
        if (present(scale_power)) scale_power = 0
        ok = .false.
        if (size(Wa)==0 .or. size(Wa)/=size(Wb) .or. .not.ieee_is_finite(tol) .or. tol<0.0_rk) return
        if (present(powers)) then
            if (any(shape(powers)/=[size(Wa),2])) return
        end if
        rtol = max(tol,128.0_rk*epsilon(1.0_rk))
        do i = 1, size(Wa)
            if (.not.ieee_is_finite(Wa(i)) .or. .not.ieee_is_finite(Wb(i))) return
            if (Wa(i)<=0.0_rk .or. Wb(i)<=0.0_rk) return
            f = fraction(Wa(i))/fraction(Wb(i))
            power = exponent(Wa(i))-exponent(Wb(i))+exponent(f)
            if (present(powers)) power = power+powers(i,1)-powers(i,2)
            f = fraction(f)
            if (i==1) then
                ref = f
                ref_power = power
            else
                common = max(power,ref_power)
                a = scale(f,power-common)
                b = scale(ref,ref_power-common)
                if (abs(a-b)>rtol*max(a,b)) return
            end if
        end do
        if (ref_power>=minexponent(ref) .and. ref_power<=maxexponent(ref)) then
            weight_scale = scale(ref,ref_power)
        else if (present(scale_power)) then
            weight_scale = ref
            scale_power = ref_power
        else
            if (ref_power<minexponent(ref)-digits(ref)+1 .or. ref_power>maxexponent(ref)) return
            a = scale(ref,ref_power)
            if (abs(scale(a,-ref_power)-ref)>rtol*ref) return
            weight_scale = a
        end if
        ok = .true.
    end subroutine
    !===============================================================================

    !> Store physical controls, padding missing coordinates with the common origin.
    pure subroutine store_constraint_points(points,offset,Xc)
        real(rk), intent(inout) :: points(:,:)
        integer, intent(in) :: offset
        real(rk), intent(in) :: Xc(:,:)
        integer :: d, last
        last = offset+size(Xc,1)
        points(offset+1:last,:size(Xc,2)) = Xc
        do d = size(Xc,2)+1,size(points,2)
            points(offset+1:last,d) = points(1,d)
        end do
    end subroutine

    !> Test each centered homogeneous CSR residual in a common radix scale.
    !! Reuse the validation copy of each coefficient and exponent to store
    !! its product with the control weight. The final component checks the
    !! denominator whenever weights are present. For each component the
    !! signed sum, its absolute-term bound and `atol` share one exponent.
    pure subroutine check_constraint_residuals(rowptr,col,val,powers,points,weights,rtol,atol,ok)
        logical, intent(out) :: ok
        integer, intent(in) :: rowptr(:), col(:)
        integer, intent(inout) :: powers(:)
        real(rk), intent(inout) :: val(:)
        real(rk), intent(in) :: points(:,:), weights(:), rtol, atol
        real(rk) :: w, delta, term, residual, magnitude, roundoff, relative
        integer :: j, row, c, power, delta_power, term_power, ncoord
        logical :: homogeneous

        ok = .false.
        if (.not. all(ieee_is_finite(val))) return
        homogeneous = size(weights)>0
        ncoord = size(points,2)
        do j = 1, size(val)
            w = 1.0_rk
            if (homogeneous) w = weights(col(j))
            powers(j) = powers(j) + exponent(val(j)) + exponent(w)
            val(j) = fraction(val(j))*fraction(w)
        end do
        do row = 1, size(rowptr)-1
            roundoff = real(rowptr(row+1)-rowptr(row),rk)*epsilon(1.0_rk)
            if (roundoff<1.0_rk) roundoff = roundoff/(1.0_rk-roundoff)
            relative = min(1.0_rk,max(rtol,roundoff))
            do c = 1, ncoord+merge(1,0,homogeneous)
                residual = 0.0_rk
                magnitude = 0.0_rk
                power = -huge(0)
                if (atol>0.0_rk) power = exponent(atol)
                do j = rowptr(row), rowptr(row+1)-1
                    delta = 1.0_rk
                    delta_power = 0
                    if (c<=ncoord) call scaled_difference(points(1,c),points(col(j),c),delta,delta_power)
                    if (delta==0.0_rk .or. val(j)==0.0_rk) cycle
                    term_power = powers(j) + delta_power + exponent(delta)
                    if (term_power>power) then
                        if (magnitude>0.0_rk) then
                            residual = scale(residual,power-term_power)
                            magnitude = scale(magnitude,power-term_power)
                        end if
                        power = term_power
                    end if
                    term = scale(val(j)*fraction(delta),term_power-power)
                    residual = residual + term
                    magnitude = magnitude + abs(term)
                end do
                if (magnitude==0.0_rk) cycle
                if (abs(residual)>scale(atol,-power)+relative*magnitude) return
            end do
        end do
        ok = .true.
    end subroutine
end module forcad_multipatch
