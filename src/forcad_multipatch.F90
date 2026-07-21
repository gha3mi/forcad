!> author: Seyed Ali Ghasemi
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
!! orientation. Geometric connections compose the patch-B normal coordinate
!! with a scalar transition map \(\phi\). The stored transition jet is
!!
!! \[
!! \boldsymbol\phi_\Gamma=
!! [\phi'(0),\phi''(0),\ldots,\phi^{(n)}(0)].
!! \]
!!
!! For corresponding mapped tangential coordinates
!! \(\boldsymbol\tau\), the order-\(r\) interface equation is
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
!! parameter-transition model. For rational curves, surfaces, and volumes, the
!! algebraic rows act on a selected homogeneous lift and permit one constant
!! projective scale between boundary weight vectors. This is sufficient for the
!! corresponding Euclidean rational join, but it is not necessary in the most
!! general projective formulation, where the homogeneous lift can be rescaled
!! by a nonconstant function whose derivative jet also enters the equations.
!!
!! For surfaces and volumes, the implemented parameter map has the restricted
!! separable form
!! \(\Phi(\boldsymbol\tau,\eta_A)
!!   =(T\boldsymbol\tau,\phi(\eta_A))\), where \(T\) is the fixed orientation
!! map. It excludes tangentially varying normal jets and normal-dependent
!! tangential coordinates. Consequently, surface and volume `geometric`
!! connections enforce \(G^n\) within this normal-only model; they do not cover
!! every multivariate \(G^n\) reparameterization allowed by the general theory.
!!
!! Tangential orientation is represented by reversal, swap, and per-axis
!! flips. Rank-specific `connect` methods validate metadata and compatible trace
!! spaces; rank-specific `is_valid` methods additionally check physical traces,
!! rational projective weights, and assembled derivative residuals. A
!! continuity value of `-1` records topology without merging degrees of freedom;
!! values `0:n` request continuity through that derivative order.
!!
!! **General geometric-continuity reference:** P. Kiciak, "Conditions for
!! geometric continuity between polynomial and rational surface patches,"
!! *Computer Aided Geometric Design* 13 (1996), 709--741. The surface and
!! volume implementation here uses the restricted normal-only transition
!! described above.
!! [doi:10.1016/0167-8396(96)00006-4](https://doi.org/10.1016/0167-8396(96)00006-4).
module forcad_multipatch

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use forcad_kinds, only: rk

    implicit none

    private
    public multipatch_connection, multipatch_compatible_trace_space, multipatch_projective_weight_scale, &
        multipatch_growth_capacity, multipatch_chain_rule_coefficients, multipatch_valid_reparameterization

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
    !! normalized local increment at the interface.
    !!
    !! @warning [[multipatch_connection:set]] is a low-level metadata setter. It
    !! does not have access to patches and cannot validate topology or spline
    !! compatibility. Prefer the rank-specific `connect` binding.
    !! @endwarning
    type multipatch_connection
        integer, private :: patch_a = 0 !! One-based first patch identifier.
        integer, private :: side_a = 0 !! First patch boundary identifier.
        integer, private :: patch_b = 0 !! One-based second patch identifier.
        integer, private :: side_b = 0 !! Second patch boundary identifier.
        integer, private :: continuity = 0 !! Requested order, or `-1` for no continuity constraint.
        logical, private :: geometric = .false. !! False for normalized \(C^n\), true for the scalar normal-transition model.
        real(rk), allocatable, private :: reparameterization(:) !! Normal transition derivatives 1 through \(n\).
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
        procedure :: is_geometric => get_connection_geometric !!> Report whether a scalar geometric transition is selected.
        procedure :: is_reversed => get_connection_reverse !!> Report tangential reversal.
        procedure :: is_swapped => get_connection_swap !!> Report tangential-axis exchange.
        procedure :: is_flipped => get_connection_flip !!> Report a selected tangential-axis flip.
    end type
    !===============================================================================

contains

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
        geometric, reparameterization)
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
        integer :: n

        this%patch_a = patch_a
        this%side_a = side_a
        this%patch_b = patch_b
        this%side_b = side_b
        this%continuity = 0
        this%geometric = .false.
        this%reverse = .false.
        this%swap = .false.
        this%flip = .false.
        if (allocated(this%reparameterization)) deallocate(this%reparameterization)
        if (present(continuity)) this%continuity = continuity
        if (present(geometric)) this%geometric = geometric
        if (present(reparameterization)) then
            allocate(this%reparameterization(size(reparameterization)), source=reparameterization)
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
    !> Report whether the connection uses the scalar geometric transition model.
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
    !> Validate the transition jet for a parametric or geometric connection.
    !! For \(G^n\), `reparameterization(k)` is \(\phi^{(k)}(0)\) in the scalar
    !! normal-only model documented by this module. For `n>0`, the array must
    !! have exactly `n` finite entries and its first derivative must satisfy
    !! `normal_sign*reparameterization(1)>sqrt(tiny(1.0_rk))`. For `n=0`, the array
    !! must be absent or empty. Higher entries are otherwise unconstrained; any
    !! finite jet with an accepted first derivative defines a local scalar
    !! diffeomorphism jet, but this routine does not validate a more general
    !! multivariate reparameterization.
    pure logical function multipatch_valid_reparameterization(connection, normal_sign) result(ok)
        type(multipatch_connection), intent(in) :: connection
            !! Connection metadata to validate.
        integer, intent(in) :: normal_sign
            !! Expected sign of the first normalized derivative.
        integer :: continuity

        continuity = connection%continuity
        if (.not. connection%geometric) then
            ok = .not. allocated(connection%reparameterization)
            return
        end if
        if (continuity < 0) then
            ok = .false.
            return
        end if
        if (continuity == 0) then
            ok = .not. allocated(connection%reparameterization) .or. size(connection%reparameterization) == 0
            return
        end if
        if (.not. allocated(connection%reparameterization)) then
            ok = .false.
            return
        end if
        ok = size(connection%reparameterization) == continuity
        if (.not. ok) return
        ok = all(ieee_is_finite(connection%reparameterization))
        if (.not. ok) return
        ok = real(normal_sign, rk)*connection%reparameterization(1) > sqrt(tiny(1.0_rk))
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
    !! by the univariate Faa di Bruno formula. For surfaces and volumes these
    !! coefficients apply only to the scalar normal composition described in
    !! the module documentation.
    !! See the module-level Kiciak reference for the general geometric-
    !! continuity theory and the documented restriction of this model.
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
    !! second trace in opposite orientation. Multiplicities are therefore
    !! preserved by the entrywise comparison.
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
            !! Absolute tolerance in normalized parameter coordinates.
        real(rk) :: start_a, end_a, start_b, end_b, ta, tb, eps
        integer :: i, j

        ok = .false.
        if (degree_a /= degree_b .or. nc_a /= nc_b) return
        if (degree_a < 0 .or. nc_a < 1) return
        if (size(knot_a) /= nc_a + degree_a + 1) return
        if (size(knot_b) /= nc_b + degree_b + 1) return

        start_a = knot_a(degree_a + 1)
        end_a = knot_a(nc_a + 1)
        start_b = knot_b(degree_b + 1)
        end_b = knot_b(nc_b + 1)
        if (end_a <= start_a .or. end_b <= start_b) return

        eps = max(tol, 128.0_rk*epsilon(1.0_rk))
        do i = 1, size(knot_a)
            ta = (knot_a(i) - start_a)/(end_a - start_a)
            if (reverse) then
                j = size(knot_b) - i + 1
                tb = 1.0_rk - (knot_b(j) - start_b)/(end_b - start_b)
            else
                tb = (knot_b(i) - start_b)/(end_b - start_b)
            end if
            if (abs(ta - tb) > eps) return
        end do
        ok = .true.
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Test projective equivalence of two rational boundary weight vectors.
    !!
    !! Homogeneous control data represent the same Euclidean geometry when all
    !! weights and weighted coordinates differ by one nonzero scalar. This
    !! routine checks only the weight part: it determines the candidate from the
    !! first weights and verifies \(W_a(i)=sW_b(i)\) componentwise. Coordinate
    !! compatibility must be established separately.
    !!
    pure subroutine multipatch_projective_weight_scale(Wa, Wb, tol, weight_scale, ok)
        real(rk), intent(in), contiguous :: Wa(:)
            !! First boundary weight vector.
        real(rk), intent(in), contiguous :: Wb(:)
            !! Second boundary weight vector of the same size.
        real(rk), intent(in) :: tol
            !! Relative comparison tolerance, augmented by roundoff scaling.
        real(rk), intent(out) :: weight_scale
            !! Projective scale \(s\), or one when validation fails.
        logical, intent(out) :: ok
            !! True only when one scale relates every corresponding weight.
        real(rk) :: ref, wtol
        integer :: i

        weight_scale = 1.0_rk
        ok = .false.
        if (size(Wa) /= size(Wb) .or. size(Wa) < 1) return
        if (abs(Wb(1)) <= tiny(1.0_rk)) return

        weight_scale = Wa(1)/Wb(1)
        ok = .true.
        do i = 1, size(Wa)
            ref = max(abs(Wa(i)), abs(weight_scale*Wb(i)))
            wtol = max(tol*ref, 128.0_rk*epsilon(1.0_rk)*ref)
            if (abs(Wa(i) - weight_scale*Wb(i)) > wtol) then
                weight_scale = 1.0_rk
                ok = .false.
                return
            end if
        end do
    end subroutine
    !===============================================================================
end module forcad_multipatch
