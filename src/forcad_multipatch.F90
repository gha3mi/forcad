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
    use forcad_kinds, only: rk

    implicit none

    private
    public multipatch_connection, multipatch_compatible_trace_space, multipatch_projective_weight_scale, &
        multipatch_growth_capacity, multipatch_chain_rule_coefficients, multipatch_composition_coefficients, &
        multipatch_valid_reparameterization

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
        if (present(reparameterization)) then
            allocate(this%reparameterization(size(reparameterization)), source=reparameterization)
        end if
        if (present(transition_jet)) then
            if (size(transition_jet) > 0) then
                allocate(this%transition_jet(&
                    size(transition_jet,1),size(transition_jet,2),size(transition_jet,3),size(transition_jet,4)),&
                    source=transition_jet)
            end if
        end if
        if (present(projective_jet)) then
            if (size(projective_jet) > 0) then
                allocate(this%projective_jet(&
                    size(projective_jet,1),size(projective_jet,2),size(projective_jet,3)),&
                    source=projective_jet)
            end if
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
        real(rk) :: choose
        integer :: a1, a2, a3, j, r, continuity, dimension

        dimension = size(transition,1)
        continuity = size(transition,2)
        allocate(coefficient(0:continuity,0:continuity,0:continuity,0:continuity), source=0.0_rk)
        coefficient(0,0,0,0) = 1.0_rk

        do r = 1, continuity
            do a3 = 0, r
                do a2 = 0, r - a3
                    do a1 = 0, r - a2 - a3
                        if (a1 + a2 + a3 == 0) cycle
                        if (a1 > 0 .and. dimension >= 1) then
                            choose = 1.0_rk
                            do j = 1, r
                                coefficient(r,a1,a2,a3) = coefficient(r,a1,a2,a3) + &
                                    choose*transition(1,j)*coefficient(r-j,a1-1,a2,a3)
                                choose = choose*real(r-j, rk)/real(j, rk)
                            end do
                        end if
                        if (a2 > 0 .and. dimension >= 2) then
                            choose = 1.0_rk
                            do j = 1, r
                                coefficient(r,a1,a2,a3) = coefficient(r,a1,a2,a3) + &
                                    choose*transition(2,j)*coefficient(r-j,a1,a2-1,a3)
                                choose = choose*real(r-j, rk)/real(j, rk)
                            end do
                        end if
                        if (a3 > 0 .and. dimension >= 3) then
                            choose = 1.0_rk
                            do j = 1, r
                                coefficient(r,a1,a2,a3) = coefficient(r,a1,a2,a3) + &
                                    choose*transition(3,j)*coefficient(r-j,a1,a2,a3-1)
                                choose = choose*real(r-j, rk)/real(j, rk)
                            end do
                        end if
                    end do
                end do
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
            allocate(transition_jet(&
                size(this%transition_jet,1),size(this%transition_jet,2),&
                size(this%transition_jet,3),size(this%transition_jet,4)),&
                source=this%transition_jet)
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
            allocate(projective_jet(&
                size(this%projective_jet,1),size(this%projective_jet,2),size(this%projective_jet,3)),&
                source=this%projective_jet)
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
            ok = .not. full_transition .and. &
                (.not. allocated(connection%reparameterization) .or. size(connection%reparameterization) == 0)
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
    !! transition jet; [[multipatch_composition_coefficients]] handles the
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
