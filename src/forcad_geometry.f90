!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Tensor-product construction and interpolation operations for NURBS geometry.
!!
!! The public generics construct a higher-rank NURBS object without sampling
!! the source geometry. Polynomial source objects remain polynomial except for
!! revolution, which necessarily introduces the rational quadratic weights of
!! exact circular arcs.
!!
!! **Constructions**
!!
!! - `extrude` forms
!!   \(\mathbf{S}(u,v)=\mathbf{C}(u)+v\mathbf d\) or
!!   \(\mathbf{V}(u,v,w)=\mathbf{S}(u,v)+w\mathbf d\), with degree one and
!!   active interval \([0,1]\) in the new direction.
!! - `revolve` constructs exact rational quadratic-arc control layers about a
!!   three-dimensional axis. Arc endpoints are rotated source controls; each
!!   intermediate control lies at the intersection of endpoint tangents. For
!!   unit axis \(\widehat{\mathbf a}\), axis point \(\mathbf a_0\),
!!   \(\mathbf r=(\mathbf I-\widehat{\mathbf a}\widehat{\mathbf a}^{T})
!!   (\mathbf P-\mathbf a_0)\), and axial projection
!!   \(\mathbf h=\widehat{\mathbf a}\widehat{\mathbf a}^{T}
!!   (\mathbf P-\mathbf a_0)\), Rodrigues' formula gives
!!
!!   \[
!!   \mathcal R_\theta(\mathbf P)=\mathbf a_0+\mathbf h+
!!   \cos\theta\,\mathbf r+
!!   \sin\theta\,(\widehat{\mathbf a}\mathbin{\times}\mathbf r).
!!   \]
!!
!!   The signed angle is split into exact rational quadratic arcs no larger
!!   than \(\pi/2\).
!! - `sweep` forms the translational tensor product
!!   \(\mathbf S(u,v)=\mathbf P(u)+\mathbf Q(v)-\mathbf o\), with
!!   tensor-product weights \(w^P_iw^Q_j\). It intentionally keeps the profile
!!   frame fixed; it is not a Frenet-, rotation-minimizing-, or guide-rail
!!   sweep.
!! - `loft` globally interpolates compatible sections in homogeneous
!!   coordinates. Section degree, control-net shape, coordinate dimension, and
!!   inherited wrapping mode must agree. Knot arrays must have equal shape and
!!   agree within the precision-scaled tolerance documented by
!!   [[same_knot_vector]]; the first section's knots are retained.
!!
!! For a loft through homogeneous section controls \(\mathbf H_i^{(\ell)}\)
!! at parameters \(t_\ell\), the new controls satisfy
!!
!! \[
!! \sum_{j=1}^{n_s}N_{j,q}(t_\ell)\mathbf Q_{ij}
!!   =\mathbf H_i^{(\ell)},\qquad \ell=1,\ldots,n_s.
!! \]
!!
!! Thus the supplied homogeneous control data are reproduced to the accepted
!! linear-system residual; no sampled fitting approximation is introduced. If
!! all section knots are identical, this reproduces each section geometry at
!! its supplied loft parameter. If knots differ only within the accepted
!! tolerance, the output still uses the first section's knots and exact equality
!! to every independently parameterized input section is not asserted. A result
!! is accepted only when all interpolated homogeneous weights remain finite and
!! strictly positive.
!!
!! Invalid input produces an invalid result whose `err` member describes the
!! failure. No partially constructed geometry is returned.
!!
!! **Reference**
!!
!! L. Piegl and W. Tiller, *The NURBS Book*, 2nd ed., Springer,
!! 1997, Chapters 8, 9, and 11.
!! [doi:10.1007/978-3-642-59223-2](https://doi.org/10.1007/978-3-642-59223-2).
module forcad_geometry

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use forcad_kinds, only: rk
    use forcad_nurbs_curve, only: nurbs_curve
    use forcad_nurbs_surface, only: nurbs_surface
    use forcad_nurbs_volume, only: nurbs_volume
    use forcad_utils, only: basis_bspline, solve

    implicit none

    private
    public extrude, revolve, sweep, loft

    real(rk), parameter :: linear_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]

    !> Construct a degree-one extrusion of a curve or surface.
    !!
    !! Curve input returns a surface; surface input returns a volume. The source
    !! directions are copied exactly and the new direction has knot vector
    !! `[0,0,1,1]`, degree one, and no parameter wrapping.
    interface extrude
        module procedure extrude_curve
        module procedure extrude_surface
    end interface

    !> Construct an exact rational revolution of a curve or surface.
    !!
    !! Curve input returns a surface and surface input returns a volume. The
    !! output adds a quadratic rational direction whose control count and knot
    !! multiplicities depend only on the number of at-most-quarter-circle arcs.
    interface revolve
        module procedure revolve_curve
        module procedure revolve_surface
    end interface

    !> Construct a fixed-frame translational sweep along a NURBS curve.
    !!
    !! The first argument is a curve or surface profile and the second argument
    !! is the spine curve. Source knot vectors, degrees, and wrapping policies
    !! are retained in their corresponding output directions.
    interface sweep
        module procedure sweep_curve
        module procedure sweep_surface
    end interface

    !> Interpolate compatible curve or surface sections in homogeneous space.
    !!
    !! Curve sections return a surface and surface sections return a volume.
    !! The output interpolation direction uses an averaged clamped knot vector;
    !! its degree defaults to `min(3,size(sections)-1)`.
    interface loft
        module procedure loft_curves
        module procedure loft_surfaces
    end interface

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Extrude a NURBS curve exactly into a degree-one tensor-product surface.
    !!
    !! The source knot vector, degree, control points, weights, and parameter
    !! wrapping are preserved in direction 1. Direction 2 uses the clamped
    !! linear knot vector `[0,0,1,1]`; the two control layers are separated by
    !! `vector` and carry identical weights.
    !!
    !! In homogeneous coordinates
    !! \(\mathbf H_i=(w_i\mathbf P_i,w_i)\), the output layers are
    !! \(\mathbf H_{i0}=\mathbf H_i\) and
    !! \(\mathbf H_{i1}=(w_i(\mathbf P_i+\mathbf d),w_i)\).
    !!
    pure function extrude_curve(curve, vector) result(surface)
        type(nurbs_curve), intent(in) :: curve
            !! Valid source curve.
        real(rk), intent(in), contiguous :: vector(:)
            !! Finite, nonzero extrusion vector. Its size determines the output coordinate dimension and must be between
            !! the source dimension and three.
        type(nurbs_surface) :: surface
        real(rk), allocatable :: Xc(:,:), Wc(:), knot(:), Xc_out(:,:), Wc_out(:)
        integer :: degree, nc, ndim, i, j, d, n

        if (.not. curve%err%ok) then
            call surface%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "Cannot extrude an invalid NURBS curve.",&
                location   = "extrude_curve",&
                suggestion = "Construct a valid curve before calling extrude.")
            return
        end if

        Xc = curve%get_Xc()
        if (size(vector) < size(Xc,2) .or. size(vector) > 3 .or. &
            .not. all(ieee_is_finite(vector)) .or. all(vector == 0.0_rk)) then
            call surface%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "The extrusion vector is invalid or has insufficient dimension.",&
                location   = "extrude_curve",&
                suggestion = "Pass a finite nonzero vector with at least the curve coordinate dimension.")
            return
        end if

        Wc = curve%get_Wc()
        if (size(Wc) == 0) then
            deallocate(Wc)
            allocate(Wc(size(Xc,1)), source=1.0_rk)
        end if
        knot = curve%get_knot()
        degree = curve%get_degree()
        nc = size(Xc,1)
        ndim = size(vector)
        allocate(Xc_out(2*nc,ndim), Wc_out(2*nc))

        do concurrent (j = 1:2, i = 1:nc, d = 1:ndim) local(n)
            n = i + (j-1)*nc
            if (d <= size(Xc,2)) then
                Xc_out(n,d) = Xc(i,d)
            else
                Xc_out(n,d) = 0.0_rk
            end if
            if (j == 2) Xc_out(n,d) = Xc_out(n,d) + vector(d)
        end do
        do concurrent (j = 1:2, i = 1:nc)
            Wc_out(i+(j-1)*nc) = Wc(i)
        end do

        call surface%set(&
            knot1    = knot,&
            knot2    = linear_knot,&
            Xc       = Xc_out,&
            Wc       = Wc_out,&
            degree   = [degree,1],&
            wrap_parameters = [curve%get_parameter_wrapping(),.false.])
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Extrude a NURBS surface exactly into a degree-one tensor-product volume.
    !!
    !! Source directions 1 and 2 are copied unchanged. Direction 3 uses
    !! `[0,0,1,1]`, and the second control layer is translated by `vector`.
    !! The construction is exact because both layers retain the source weights.
    !!
    pure function extrude_surface(surface, vector) result(volume)
        type(nurbs_surface), intent(in) :: surface
            !! Valid source surface.
        real(rk), intent(in), contiguous :: vector(:)
            !! Finite, nonzero extrusion vector. Its size determines the output coordinate dimension and must be between
            !! the source dimension and three.
        type(nurbs_volume) :: volume
        real(rk), allocatable :: Xc(:,:), Wc(:), knot1(:), knot2(:), Xc_out(:,:), Wc_out(:)
        integer :: degree(2), nc, ndim, i, j, d, n

        if (.not. surface%err%ok) then
            call volume%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "Cannot extrude an invalid NURBS surface.",&
                location   = "extrude_surface",&
                suggestion = "Construct a valid surface before calling extrude.")
            return
        end if

        Xc = surface%get_Xc()
        if (size(vector) < size(Xc,2) .or. size(vector) > 3 .or. &
            .not. all(ieee_is_finite(vector)) .or. all(vector == 0.0_rk)) then
            call volume%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "The extrusion vector is invalid or has insufficient dimension.",&
                location   = "extrude_surface",&
                suggestion = "Pass a finite nonzero vector with at least the surface coordinate dimension.")
            return
        end if

        Wc = surface%get_Wc()
        if (size(Wc) == 0) then
            deallocate(Wc)
            allocate(Wc(size(Xc,1)), source=1.0_rk)
        end if
        knot1 = surface%get_knot(1)
        knot2 = surface%get_knot(2)
        degree = surface%get_degree()
        nc = size(Xc,1)
        ndim = size(vector)
        allocate(Xc_out(2*nc,ndim), Wc_out(2*nc))

        do concurrent (j = 1:2, i = 1:nc, d = 1:ndim) local(n)
            n = i + (j-1)*nc
            if (d <= size(Xc,2)) then
                Xc_out(n,d) = Xc(i,d)
            else
                Xc_out(n,d) = 0.0_rk
            end if
            if (j == 2) Xc_out(n,d) = Xc_out(n,d) + vector(d)
        end do
        do concurrent (j = 1:2, i = 1:nc)
            Wc_out(i+(j-1)*nc) = Wc(i)
        end do

        call volume%set(&
            knot1    = knot1,&
            knot2    = knot2,&
            knot3    = linear_knot,&
            Xc       = Xc_out,&
            Wc       = Wc_out,&
            degree   = [degree,1],&
            wrap_parameters = [surface%get_parameter_wrapping(1),surface%get_parameter_wrapping(2),.false.])
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Revolve a NURBS curve exactly into a rational quadratic-arc surface.
    !! The angle is signed and measured in radians. The axis is defined by a
    !! three-dimensional point and any finite nonzero direction vector.
    !! Every arc spans at most \(\pi/2\); its middle-layer weight is multiplied
    !! by \(\cos(\Delta\theta/2)\). The source parameter is direction 1 and the
    !! revolution parameter is direction 2.
    !!
    !! For an arc of signed angle \(\Delta\theta\), endpoint layers use the
    !! original section weight \(w_i\), while the tangent-intersection layer
    !! uses \(w_i\cos(\Delta\theta/2)\). This is the exact projective conic
    !! construction, not a polynomial approximation. One- and two-coordinate
    !! source controls are embedded in three dimensions by appending zeros.
    !!
    !! See Piegl and Tiller, *The NURBS Book*, 2nd ed., Section 8.4.
    pure function revolve_curve(curve, axis_point, axis_direction, angle) result(surface)
        type(nurbs_curve), intent(in) :: curve
            !! Valid source curve with at most three coordinate components.
        real(rk), intent(in), contiguous :: axis_point(:)
            !! Point on the axis, shape `[3]`.
        real(rk), intent(in), contiguous :: axis_direction(:)
            !! Finite nonzero axis direction, shape `[3]`; it need not be normalized.
        real(rk), intent(in) :: angle
            !! Nonzero signed revolution angle in radians.
        type(nurbs_surface) :: surface
        real(rk), allocatable :: Xc(:,:), Wc(:), knot(:), Xc_out(:,:), Wc_out(:), knot_revolve(:)
        integer :: degree

        if (.not. curve%err%ok) then
            call surface%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "Cannot revolve an invalid NURBS curve.",&
                location   = "revolve_curve",&
                suggestion = "Construct a valid curve before calling revolve.")
            return
        end if

        Xc = curve%get_Xc()
        if (.not. valid_revolution(Xc, axis_point, axis_direction, angle)) then
            call surface%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "The revolution axis, angle, or source coordinate dimension is invalid.",&
                location   = "revolve_curve",&
                suggestion = "Pass finite three-component axis data, a nonzero axis and angle, and at most 3D coordinates.")
            return
        end if

        Wc = curve%get_Wc()
        if (size(Wc) == 0) then
            deallocate(Wc)
            allocate(Wc(size(Xc,1)), source=1.0_rk)
        end if
        knot = curve%get_knot()
        degree = curve%get_degree()
        call revolution_control_net(&
            Xc             = Xc,&
            Wc             = Wc,&
            axis_point     = axis_point,&
            axis_direction = axis_direction,&
            angle          = angle,&
            Xc_out         = Xc_out,&
            Wc_out         = Wc_out,&
            knot           = knot_revolve)
        call surface%set(&
            knot1    = knot,&
            knot2    = knot_revolve,&
            Xc       = Xc_out,&
            Wc       = Wc_out,&
            degree   = [degree,2],&
            wrap_parameters = [curve%get_parameter_wrapping(),.false.])
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Revolve a NURBS surface exactly into a rational quadratic-arc volume.
    !! The angle is signed and measured in radians. The axis is defined by a
    !! three-dimensional point and any finite nonzero direction vector.
    !! Source directions 1 and 2 are retained; the revolution parameter becomes
    !! direction 3. Every rational arc spans at most \(\pi/2\). One- and
    !! two-coordinate source controls are embedded in three dimensions by
    !! appending zeros.
    !!
    !! See Piegl and Tiller, *The NURBS Book*, 2nd ed., Section 8.4.
    pure function revolve_surface(surface, axis_point, axis_direction, angle) result(volume)
        type(nurbs_surface), intent(in) :: surface
            !! Valid source surface with at most three coordinate components.
        real(rk), intent(in), contiguous :: axis_point(:)
            !! Point on the axis, shape `[3]`.
        real(rk), intent(in), contiguous :: axis_direction(:)
            !! Finite nonzero axis direction, shape `[3]`; it need not be normalized.
        real(rk), intent(in) :: angle
            !! Nonzero signed revolution angle in radians.
        type(nurbs_volume) :: volume
        real(rk), allocatable :: Xc(:,:), Wc(:), knot1(:), knot2(:), Xc_out(:,:), Wc_out(:), knot_revolve(:)
        integer :: degree(2)

        if (.not. surface%err%ok) then
            call volume%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "Cannot revolve an invalid NURBS surface.",&
                location   = "revolve_surface",&
                suggestion = "Construct a valid surface before calling revolve.")
            return
        end if

        Xc = surface%get_Xc()
        if (.not. valid_revolution(Xc, axis_point, axis_direction, angle)) then
            call volume%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "The revolution axis, angle, or source coordinate dimension is invalid.",&
                location   = "revolve_surface",&
                suggestion = "Pass finite three-component axis data, a nonzero axis and angle, and at most 3D coordinates.")
            return
        end if

        Wc = surface%get_Wc()
        if (size(Wc) == 0) then
            deallocate(Wc)
            allocate(Wc(size(Xc,1)), source=1.0_rk)
        end if
        knot1 = surface%get_knot(1)
        knot2 = surface%get_knot(2)
        degree = surface%get_degree()
        call revolution_control_net(&
            Xc             = Xc,&
            Wc             = Wc,&
            axis_point     = axis_point,&
            axis_direction = axis_direction,&
            angle          = angle,&
            Xc_out         = Xc_out,&
            Wc_out         = Wc_out,&
            knot           = knot_revolve)
        call volume%set(&
            knot1    = knot1,&
            knot2    = knot2,&
            knot3    = knot_revolve,&
            Xc       = Xc_out,&
            Wc       = Wc_out,&
            degree   = [degree,2],&
            wrap_parameters = [surface%get_parameter_wrapping(1),surface%get_parameter_wrapping(2),.false.])
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Check the common revolution input invariants.
    pure logical function valid_revolution(Xc, axis_point, axis_direction, angle) result(ok)
        real(rk), intent(in), contiguous :: Xc(:,:), axis_point(:), axis_direction(:)
        real(rk), intent(in) :: angle

        ok = .false.
        if (size(Xc,2) < 1 .or. size(Xc,2) > 3 .or. &
            size(axis_point) /= 3 .or. size(axis_direction) /= 3) return
        if (.not. all(ieee_is_finite(axis_point)) .or. &
            .not. all(ieee_is_finite(axis_direction)) .or. .not. ieee_is_finite(angle)) return

        ok = maxval(abs(axis_direction)) > 0.0_rk .and. angle /= 0.0_rk
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct exact quadratic rational-arc control layers for a revolution.
    pure subroutine revolution_control_net(Xc, Wc, axis_point, axis_direction, angle, Xc_out, Wc_out, knot)
        real(rk), intent(in), contiguous :: Xc(:,:), Wc(:), axis_point(:), axis_direction(:)
        real(rk), intent(in) :: angle
        real(rk), allocatable, intent(out) :: Xc_out(:,:), Wc_out(:), knot(:)
        real(rk) :: axis(3), point(3), projection(3), radial(3), tangent(3), rotated(3)
        real(rk) :: pi, step, theta, arc_weight, factor, axis_scale
        integer :: narc, nlayer, nc, i, j, k, n

        pi = acos(-1.0_rk)
        narc = ceiling(abs(angle)/(0.5_rk*pi))
        nlayer = 2*narc + 1
        nc = size(Xc,1)
        step = angle/real(narc,rk)
        arc_weight = cos(0.5_rk*step)
        axis_scale = maxval(abs(axis_direction))
        axis = axis_direction/axis_scale
        axis = axis/sqrt(dot_product(axis,axis))
        allocate(Xc_out(nc*nlayer,3), Wc_out(nc*nlayer), knot(nlayer+3))

        knot(1:3) = 0.0_rk
        do k = 1, narc - 1
            knot(2*k+2:2*k+3) = real(k,rk)/real(narc,rk)
        end do
        knot(size(knot)-2:) = 1.0_rk

        do concurrent (j = 1:nlayer, i = 1:nc) &
            local(point, projection, radial, tangent, rotated, theta, factor, n)
            point = 0.0_rk
            point(1:size(Xc,2)) = Xc(i,:)
            projection = axis_point + dot_product(point-axis_point,axis)*axis
            radial = point - projection
            tangent = [axis(2)*radial(3)-axis(3)*radial(2), &
                axis(3)*radial(1)-axis(1)*radial(3), &
                axis(1)*radial(2)-axis(2)*radial(1)]
            if (mod(j,2) == 1) then
                theta = real((j-1)/2,rk)*step
                factor = 1.0_rk
            else
                theta = (real(j/2,rk)-0.5_rk)*step
                factor = 1.0_rk/arc_weight
            end if
            rotated = projection + factor*(cos(theta)*radial+sin(theta)*tangent)
            n = i + (j-1)*nc
            Xc_out(n,:) = rotated
            if (mod(j,2) == 1) then
                Wc_out(n) = Wc(i)
            else
                Wc_out(n) = arc_weight*Wc(i)
            end if
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Sweep a NURBS curve by exact translation along a NURBS spine.
    !! The profile orientation remains fixed. The optional origin is subtracted
    !! from the spine before translation. The output control net is the tensor
    !! sum of both control polygons and its weights are the pairwise products
    !! \(w^{P}_i w^{Q}_j\).
    !!
    !! The resulting homogeneous tensor product evaluates exactly to
    !! \(\mathbf P(u)+\mathbf Q(v)-\mathbf o\), provided both source weight
    !! denominators are nonzero.
    !!
    pure function sweep_curve(profile, spine, origin) result(surface)
        type(nurbs_curve), intent(in) :: profile
            !! Valid profile curve; output direction 1.
        type(nurbs_curve), intent(in) :: spine
            !! Valid translation curve; output direction 2.
        real(rk), intent(in), contiguous, optional :: origin(:)
            !! Optional finite translation origin. If present, its size must equal the output coordinate dimension.
        type(nurbs_surface) :: surface
        real(rk), allocatable :: profile_Xc(:,:), spine_Xc(:,:), profile_Wc(:), spine_Wc(:)
        real(rk), allocatable :: profile_knot(:), spine_knot(:), Xc_out(:,:), Wc_out(:)
        real(rk) :: coordinate
        integer :: profile_degree, spine_degree, np, ns, ndim, i, j, d, n
        logical :: invalid_origin

        if (.not. profile%err%ok .or. .not. spine%err%ok) then
            call surface%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "Cannot sweep an invalid NURBS profile or spine.",&
                location   = "sweep_curve",&
                suggestion = "Construct valid profile and spine curves before calling sweep.")
            return
        end if

        profile_Xc = profile%get_Xc()
        spine_Xc = spine%get_Xc()
        ndim = max(size(profile_Xc,2),size(spine_Xc,2))
        invalid_origin = .false.
        if (present(origin)) invalid_origin = size(origin) /= ndim .or. .not. all(ieee_is_finite(origin))
        if (ndim > 3 .or. invalid_origin) then
            call surface%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "The sweep coordinate dimension or spine origin is invalid.",&
                location   = "sweep_curve",&
                suggestion = "Use at most 3D geometry and pass a finite origin matching the output dimension.")
            return
        end if

        profile_Wc = profile%get_Wc()
        if (size(profile_Wc) == 0) then
            deallocate(profile_Wc)
            allocate(profile_Wc(size(profile_Xc,1)), source=1.0_rk)
        end if
        spine_Wc = spine%get_Wc()
        if (size(spine_Wc) == 0) then
            deallocate(spine_Wc)
            allocate(spine_Wc(size(spine_Xc,1)), source=1.0_rk)
        end if
        profile_knot = profile%get_knot()
        spine_knot = spine%get_knot()
        profile_degree = profile%get_degree()
        spine_degree = spine%get_degree()
        np = size(profile_Xc,1)
        ns = size(spine_Xc,1)
        allocate(Xc_out(np*ns,ndim), Wc_out(np*ns))

        do concurrent (j = 1:ns, i = 1:np, d = 1:ndim) local(coordinate, n)
            coordinate = 0.0_rk
            if (present(origin)) coordinate = -origin(d)
            if (d <= size(profile_Xc,2)) coordinate = coordinate + profile_Xc(i,d)
            if (d <= size(spine_Xc,2)) coordinate = coordinate + spine_Xc(j,d)
            n = i + (j-1)*np
            Xc_out(n,d) = coordinate
        end do
        do concurrent (j = 1:ns, i = 1:np)
            Wc_out(i+(j-1)*np) = profile_Wc(i)*spine_Wc(j)
        end do

        call surface%set(&
            knot1    = profile_knot,&
            knot2    = spine_knot,&
            Xc       = Xc_out,&
            Wc       = Wc_out,&
            degree   = [profile_degree,spine_degree],&
            wrap_parameters = [profile%get_parameter_wrapping(),spine%get_parameter_wrapping()])
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Sweep a NURBS surface by exact translation along a NURBS spine.
    !! The profile orientation remains fixed. The optional origin is subtracted
    !! from the spine before translation. The surface directions become output
    !! directions 1 and 2; the spine becomes direction 3.
    !!
    pure function sweep_surface(profile, spine, origin) result(volume)
        type(nurbs_surface), intent(in) :: profile
            !! Valid profile surface.
        type(nurbs_curve), intent(in) :: spine
            !! Valid translation curve.
        real(rk), intent(in), contiguous, optional :: origin(:)
            !! Optional finite translation origin. If present, its size must equal the output coordinate dimension.
        type(nurbs_volume) :: volume
        real(rk), allocatable :: profile_Xc(:,:), spine_Xc(:,:), profile_Wc(:), spine_Wc(:)
        real(rk), allocatable :: profile_knot1(:), profile_knot2(:), spine_knot(:)
        real(rk), allocatable :: Xc_out(:,:), Wc_out(:)
        real(rk) :: coordinate
        integer :: profile_degree(2), spine_degree, np, ns, ndim, i, j, d, n
        logical :: invalid_origin

        if (.not. profile%err%ok .or. .not. spine%err%ok) then
            call volume%err%set(&
                code       = 102,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "Cannot sweep an invalid NURBS profile or spine.",&
                location   = "sweep_surface",&
                suggestion = "Construct a valid profile surface and spine curve before calling sweep.")
            return
        end if

        profile_Xc = profile%get_Xc()
        spine_Xc = spine%get_Xc()
        ndim = max(size(profile_Xc,2),size(spine_Xc,2))
        invalid_origin = .false.
        if (present(origin)) invalid_origin = size(origin) /= ndim .or. .not. all(ieee_is_finite(origin))
        if (ndim > 3 .or. invalid_origin) then
            call volume%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "The sweep coordinate dimension or spine origin is invalid.",&
                location   = "sweep_surface",&
                suggestion = "Use at most 3D geometry and pass a finite origin matching the output dimension.")
            return
        end if

        profile_Wc = profile%get_Wc()
        if (size(profile_Wc) == 0) then
            deallocate(profile_Wc)
            allocate(profile_Wc(size(profile_Xc,1)), source=1.0_rk)
        end if
        spine_Wc = spine%get_Wc()
        if (size(spine_Wc) == 0) then
            deallocate(spine_Wc)
            allocate(spine_Wc(size(spine_Xc,1)), source=1.0_rk)
        end if
        profile_knot1 = profile%get_knot(1)
        profile_knot2 = profile%get_knot(2)
        spine_knot = spine%get_knot()
        profile_degree = profile%get_degree()
        spine_degree = spine%get_degree()
        np = size(profile_Xc,1)
        ns = size(spine_Xc,1)
        allocate(Xc_out(np*ns,ndim), Wc_out(np*ns))

        do concurrent (j = 1:ns, i = 1:np, d = 1:ndim) local(coordinate, n)
            coordinate = 0.0_rk
            if (present(origin)) coordinate = -origin(d)
            if (d <= size(profile_Xc,2)) coordinate = coordinate + profile_Xc(i,d)
            if (d <= size(spine_Xc,2)) coordinate = coordinate + spine_Xc(j,d)
            n = i + (j-1)*np
            Xc_out(n,d) = coordinate
        end do
        do concurrent (j = 1:ns, i = 1:np)
            Wc_out(i+(j-1)*np) = profile_Wc(i)*spine_Wc(j)
        end do

        call volume%set(&
            knot1    = profile_knot1,&
            knot2    = profile_knot2,&
            knot3    = spine_knot,&
            Xc       = Xc_out,&
            Wc       = Wc_out,&
            degree   = [profile_degree,spine_degree],&
            wrap_parameters = [profile%get_parameter_wrapping(1),profile%get_parameter_wrapping(2),spine%get_parameter_wrapping()])
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Interpolate compatible NURBS curve sections by a NURBS surface.
    !! Section knot vectors, degrees, control counts and coordinate dimensions
    !! must match, with knot values compared by [[same_knot_vector]]. Parameters
    !! must be finite and strictly increasing. The default loft degree is
    !! `min(3,number_of_sections-1)`. Interpolation is performed on homogeneous
    !! coordinates. With identical section knots, rational section geometry is
    !! reproduced at each supplied section parameter to the accepted linear
    !! residual. When knots differ within tolerance, the first knot vector is
    !! retained and only homogeneous control-data reproduction is asserted.
    !! For each source control index \(i\), the routine solves
    !! \(A\mathbf Q_i=\mathbf H_i\), where
    !! \(A_{\ell j}=N_{j,q}(t_\ell)\).
    !!
    !! Uses the averaged open knot vector and global interpolation construction
    !! of Piegl and Tiller, *The NURBS Book*, 2nd ed., Section 9.2.1. The result
    !! is rejected if the solve residual exceeds the criterion documented by
    !! [[interpolate_loft]] or any interpolated weight is nonpositive/nonfinite.
    pure function loft_curves(sections, parameters, degree) result(surface)
        type(nurbs_curve), intent(in), contiguous :: sections(:)
            !! Two or more valid and mutually compatible curves.
        real(rk), intent(in), contiguous, optional :: parameters(:)
            !! Optional strictly increasing section parameters, one per section. Uniform values on `[0,1]` are used when
            !! absent.
        integer, intent(in), optional :: degree
            !! Optional loft-direction degree in `[1,size(sections)-1]`.
        type(nurbs_surface) :: surface
        real(rk), allocatable :: Xc(:,:), knot(:), section_Xc(:,:), section_Wc(:), section_knot(:)
        real(rk), allocatable :: parameters_(:), data(:,:), control(:,:), loft_knot(:), Xc_out(:,:), Wc_out(:)
        integer :: nsection, section_degree, loft_degree, nc, ndim, section, i, d, n
        logical :: ok

        nsection = size(sections)
        if (nsection < 2) then
            call surface%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "A curve loft requires at least two sections.",&
                location   = "loft_curves",&
                suggestion = "Pass an array containing two or more compatible NURBS curves.")
            return
        end if
        do section = 1, nsection
            if (.not. sections(section)%err%ok) then
                call surface%err%set(&
                    code       = 102,&
                    severity   = 1,&
                    category   = "forcad_geometry",&
                    message    = "Cannot loft an invalid NURBS curve section.",&
                    location   = "loft_curves",&
                    suggestion = "Construct every section successfully before calling loft.")
                return
            end if
        end do

        Xc = sections(1)%get_Xc()
        knot = sections(1)%get_knot()
        section_degree = sections(1)%get_degree()
        nc = size(Xc,1)
        ndim = size(Xc,2)
        if (ndim > 3) then
            call surface%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "Curve lofting supports geometry coordinates through dimension three.",&
                location   = "loft_curves",&
                suggestion = "Use one-, two-, or three-dimensional curve control points.")
            return
        end if

        allocate(parameters_(nsection))
        if (present(parameters)) then
            if (size(parameters) /= nsection .or. .not. all(ieee_is_finite(parameters))) then
                call surface%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_geometry",&
                    message    = "Loft parameters must provide one finite value per curve section.",&
                    location   = "loft_curves",&
                    suggestion = "Pass a finite parameter array with size equal to the number of sections.")
                return
            end if
            parameters_ = parameters
        else
            do concurrent (section = 1:nsection)
                parameters_(section) = real(section-1,rk)/real(nsection-1,rk)
            end do
        end if
        if (any(parameters_(2:) <= parameters_(:nsection-1))) then
            call surface%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "Loft parameters must be strictly increasing.",&
                location   = "loft_curves",&
                suggestion = "Remove duplicate parameters and order the sections along the loft direction.")
            return
        end if

        loft_degree = min(3,nsection-1)
        if (present(degree)) loft_degree = degree
        if (loft_degree < 1 .or. loft_degree >= nsection) then
            call surface%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "The loft degree must be between one and number_of_sections-1.",&
                location   = "loft_curves",&
                suggestion = "Reduce the degree or provide more curve sections.")
            return
        end if

        allocate(data(nsection,nc*(ndim+1)))
        do section = 1, nsection
            section_Xc = sections(section)%get_Xc()
            section_knot = sections(section)%get_knot()
            if (any(shape(section_Xc) /= shape(Xc)) .or. &
                sections(section)%get_degree() /= section_degree .or. &
                .not. same_knot_vector(section_knot,knot) .or. &
                sections(section)%get_parameter_wrapping() .neqv. sections(1)%get_parameter_wrapping()) then
                call surface%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_geometry",&
                    message    = "Curve loft sections do not share one compatible NURBS space.",&
                    location   = "loft_curves",&
                    suggestion = "Match section degrees, knot vectors, control counts, and coordinate dimensions.")
                return
            end if
            section_Wc = sections(section)%get_Wc()
            if (size(section_Wc) == 0) then
                deallocate(section_Wc)
                allocate(section_Wc(nc), source=1.0_rk)
            end if
            do concurrent (i = 1:nc, d = 1:ndim)
                data(section,(i-1)*(ndim+1)+d) = section_Xc(i,d)*section_Wc(i)
            end do
            do concurrent (i = 1:nc)
                data(section,i*(ndim+1)) = section_Wc(i)
            end do
        end do

        call interpolate_loft(parameters_, loft_degree, data, loft_knot, control, ok)
        if (.not. ok) then
            call surface%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "The homogeneous curve-loft interpolation system could not be solved accurately.",&
                location   = "loft_curves",&
                suggestion = "Improve parameter spacing, reduce the loft degree, or remove duplicate sections.")
            return
        end if

        allocate(Wc_out(nc*nsection))
        do concurrent (section = 1:nsection, i = 1:nc) local(n)
            n = i + (section-1)*nc
            Wc_out(n) = control(section,i*(ndim+1))
        end do
        if (.not. all(ieee_is_finite(Wc_out) .and. Wc_out > 0.0_rk)) then
            call surface%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "Curve-loft interpolation produced a nonpositive homogeneous weight.",&
                location   = "loft_curves",&
                suggestion = "Use degree one or provide section weights with a smoother variation.")
            return
        end if
        allocate(Xc_out(nc*nsection,ndim))
        do concurrent (section = 1:nsection, i = 1:nc, d = 1:ndim) local(n)
            n = i + (section-1)*nc
            Xc_out(n,d) = control(section,(i-1)*(ndim+1)+d)/Wc_out(n)
        end do

        call surface%set(&
            knot1    = knot,&
            knot2    = loft_knot,&
            Xc       = Xc_out,&
            Wc       = Wc_out,&
            degree   = [section_degree,loft_degree],&
            wrap_parameters = [sections(1)%get_parameter_wrapping(),.false.])
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Interpolate compatible NURBS surface sections by a NURBS volume.
    !! Section knot vectors, degrees, control counts and coordinate dimensions
    !! must match, with directional knot values compared by
    !! [[same_knot_vector]]. Parameters must be finite and strictly increasing.
    !! The default loft degree is `min(3,number_of_sections-1)`.
    !! Interpolation is performed independently for every homogeneous
    !! control-net coordinate, including the weight. Identical section knots
    !! yield section-geometry interpolation to the accepted solve residual. If
    !! knots differ within tolerance, the first surface's knots are retained and
    !! only homogeneous control-data reproduction is asserted.
    !! For each flattened source control index \(i\), the routine solves
    !! \(A\mathbf Q_i=\mathbf H_i\), where
    !! \(A_{\ell j}=N_{j,q}(t_\ell)\).
    !!
    !! Uses the averaged open knot vector and global interpolation construction
    !! of Piegl and Tiller, *The NURBS Book*, 2nd ed., Section 9.2.1. The result
    !! is rejected if the solve residual exceeds the criterion documented by
    !! [[interpolate_loft]] or any interpolated weight is nonpositive/nonfinite.
    pure function loft_surfaces(sections, parameters, degree) result(volume)
        type(nurbs_surface), intent(in), contiguous :: sections(:)
            !! Two or more valid and mutually compatible surfaces.
        real(rk), intent(in), contiguous, optional :: parameters(:)
            !! Optional strictly increasing section parameters, one per section. Uniform values on `[0,1]` are used when
            !! absent.
        integer, intent(in), optional :: degree
            !! Optional loft-direction degree in `[1,size(sections)-1]`.
        type(nurbs_volume) :: volume
        real(rk), allocatable :: Xc(:,:), knot1(:), knot2(:)
        real(rk), allocatable :: section_Xc(:,:), section_Wc(:), section_knot1(:), section_knot2(:)
        real(rk), allocatable :: parameters_(:), data(:,:), control(:,:), loft_knot(:), Xc_out(:,:), Wc_out(:)
        integer :: nsection, section_degree(2), section_nc(2), loft_degree, nc, ndim, section, i, d, n
        logical :: ok

        nsection = size(sections)
        if (nsection < 2) then
            call volume%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "A surface loft requires at least two sections.",&
                location   = "loft_surfaces",&
                suggestion = "Pass an array containing two or more compatible NURBS surfaces.")
            return
        end if
        do section = 1, nsection
            if (.not. sections(section)%err%ok) then
                call volume%err%set(&
                    code       = 102,&
                    severity   = 1,&
                    category   = "forcad_geometry",&
                    message    = "Cannot loft an invalid NURBS surface section.",&
                    location   = "loft_surfaces",&
                    suggestion = "Construct every section successfully before calling loft.")
                return
            end if
        end do

        Xc = sections(1)%get_Xc()
        knot1 = sections(1)%get_knot(1)
        knot2 = sections(1)%get_knot(2)
        section_degree = sections(1)%get_degree()
        section_nc = sections(1)%get_nc()
        nc = size(Xc,1)
        ndim = size(Xc,2)
        if (ndim > 3) then
            call volume%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "Surface lofting supports geometry coordinates through dimension three.",&
                location   = "loft_surfaces",&
                suggestion = "Use one-, two-, or three-dimensional surface control points.")
            return
        end if

        allocate(parameters_(nsection))
        if (present(parameters)) then
            if (size(parameters) /= nsection .or. .not. all(ieee_is_finite(parameters))) then
                call volume%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_geometry",&
                    message    = "Loft parameters must provide one finite value per surface section.",&
                    location   = "loft_surfaces",&
                    suggestion = "Pass a finite parameter array with size equal to the number of sections.")
                return
            end if
            parameters_ = parameters
        else
            do concurrent (section = 1:nsection)
                parameters_(section) = real(section-1,rk)/real(nsection-1,rk)
            end do
        end if
        if (any(parameters_(2:) <= parameters_(:nsection-1))) then
            call volume%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "Loft parameters must be strictly increasing.",&
                location   = "loft_surfaces",&
                suggestion = "Remove duplicate parameters and order the sections along the loft direction.")
            return
        end if

        loft_degree = min(3,nsection-1)
        if (present(degree)) loft_degree = degree
        if (loft_degree < 1 .or. loft_degree >= nsection) then
            call volume%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "The loft degree must be between one and number_of_sections-1.",&
                location   = "loft_surfaces",&
                suggestion = "Reduce the degree or provide more surface sections.")
            return
        end if

        allocate(data(nsection,nc*(ndim+1)))
        do section = 1, nsection
            section_Xc = sections(section)%get_Xc()
            section_knot1 = sections(section)%get_knot(1)
            section_knot2 = sections(section)%get_knot(2)
            if (any(shape(section_Xc) /= shape(Xc)) .or. &
                any(sections(section)%get_degree() /= section_degree) .or. &
                any(sections(section)%get_nc() /= section_nc) .or. &
                .not. same_knot_vector(section_knot1,knot1) .or. &
                .not. same_knot_vector(section_knot2,knot2) .or. &
                sections(section)%get_parameter_wrapping(1) .neqv. sections(1)%get_parameter_wrapping(1) .or. &
                sections(section)%get_parameter_wrapping(2) .neqv. sections(1)%get_parameter_wrapping(2)) then
                call volume%err%set(&
                    code       = 101,&
                    severity   = 1,&
                    category   = "forcad_geometry",&
                    message    = "Surface loft sections do not share one compatible NURBS space.",&
                    location   = "loft_surfaces",&
                    suggestion = "Match section degrees, knot vectors, control counts, and coordinate dimensions.")
                return
            end if
            section_Wc = sections(section)%get_Wc()
            if (size(section_Wc) == 0) then
                deallocate(section_Wc)
                allocate(section_Wc(nc), source=1.0_rk)
            end if
            do concurrent (i = 1:nc, d = 1:ndim)
                data(section,(i-1)*(ndim+1)+d) = section_Xc(i,d)*section_Wc(i)
            end do
            do concurrent (i = 1:nc)
                data(section,i*(ndim+1)) = section_Wc(i)
            end do
        end do

        call interpolate_loft(parameters_, loft_degree, data, loft_knot, control, ok)
        if (.not. ok) then
            call volume%err%set(&
                code       = 108,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "The homogeneous surface-loft interpolation system could not be solved accurately.",&
                location   = "loft_surfaces",&
                suggestion = "Improve parameter spacing, reduce the loft degree, or remove duplicate sections.")
            return
        end if

        allocate(Wc_out(nc*nsection))
        do concurrent (section = 1:nsection, i = 1:nc) local(n)
            n = i + (section-1)*nc
            Wc_out(n) = control(section,i*(ndim+1))
        end do
        if (.not. all(ieee_is_finite(Wc_out) .and. Wc_out > 0.0_rk)) then
            call volume%err%set(&
                code       = 100,&
                severity   = 1,&
                category   = "forcad_geometry",&
                message    = "Surface-loft interpolation produced a nonpositive homogeneous weight.",&
                location   = "loft_surfaces",&
                suggestion = "Use degree one or provide section weights with a smoother variation.")
            return
        end if
        allocate(Xc_out(nc*nsection,ndim))
        do concurrent (section = 1:nsection, i = 1:nc, d = 1:ndim) local(n)
            n = i + (section-1)*nc
            Xc_out(n,d) = control(section,(i-1)*(ndim+1)+d)/Wc_out(n)
        end do

        call volume%set(&
            knot1    = knot1,&
            knot2    = knot2,&
            knot3    = loft_knot,&
            Xc       = Xc_out,&
            Wc       = Wc_out,&
            degree   = [section_degree,loft_degree],&
            wrap_parameters = [sections(1)%get_parameter_wrapping(1),sections(1)%get_parameter_wrapping(2),.false.])
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compare same-shaped knot vectors with a precision-scaled absolute tolerance.
    !! The vectors are accepted exactly when
    !!
    !! \[
    !! \max_i|U_i-V_i|\le64\,\epsilon_{rk}
    !! \max\!\left(1,\max_i|U_i|,\max_i|V_i|\right).
    !! \]
    !!
    !! This is a numerical compatibility test, not exact knot identity.
    pure logical function same_knot_vector(knot_a, knot_b) result(same)
        real(rk), intent(in), contiguous :: knot_a(:), knot_b(:)
        real(rk) :: knot_scale

        same = .false.
        if (size(knot_a) /= size(knot_b)) return
        if (size(knot_a) == 0) then
            same = .true.
            return
        end if
        knot_scale = max(1.0_rk,maxval(abs(knot_a)),maxval(abs(knot_b)))
        same = maxval(abs(knot_a-knot_b)) <= 64.0_rk*epsilon(1.0_rk)*knot_scale
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Interpolate homogeneous section data with an averaged open knot vector.
    !! Uses global B-spline interpolation from Piegl and Tiller, The NURBS Book,
    !! 2nd edition, Section 9.2.1. For zero-based section parameters
    !! \(t_0<\cdots<t_n\) and degree \(q\), each endpoint is repeated \(q+1\)
    !! times and the internal knots are
    !!
    !! \[
    !! U_{j+q}=\frac1q\sum_{i=j}^{j+q-1}t_i,
    !! \qquad j=1,\ldots,n-q.
    !! \]
    !!
    !! The square collocation system is solved by [[solve]]. `ok` is true only
    !! for a finite solution whose maximum componentwise interpolation residual
    !! satisfies
    !!
    !! \[
    !! r_{\max}\le4096\,\epsilon_{rk}(n+1)
    !! \max\!\left(1,\max_{\ell,c}|H_{\ell c}|\right).
    !! \]
    pure subroutine interpolate_loft(parameters, degree, data, knot, control, ok)
        real(rk), intent(in), contiguous :: parameters(:), data(:,:)
        integer, intent(in) :: degree
        real(rk), allocatable, intent(out) :: knot(:), control(:,:)
        logical, intent(out) :: ok
        real(rk), allocatable :: collocation(:,:)
        real(rk) :: residual, data_scale, fitted_value
        integer :: nsection, section, j, column

        ok = .false.
        nsection = size(parameters)
        if (nsection < 2 .or. degree < 1 .or. degree >= nsection .or. &
            size(data,1) /= nsection) then
            allocate(knot(0), control(0,0))
            return
        end if

        allocate(knot(nsection+degree+1))
        knot(1:degree+1) = parameters(1)
        do j = 1, nsection - degree - 1
            knot(degree+1+j) = sum(parameters(j+1:j+degree))/real(degree,rk)
        end do
        knot(nsection+1:) = parameters(nsection)

        allocate(collocation(nsection,nsection))
        do section = 1, nsection
            collocation(section,:) = basis_bspline(parameters(section),knot,nsection,degree)
        end do
        control = solve(collocation,data)
        if (size(control,1) /= nsection .or. size(control,2) /= size(data,2)) return
        if (.not. all(ieee_is_finite(control))) return

        residual = 0.0_rk
        do column = 1, size(data,2)
            do section = 1, nsection
                fitted_value = dot_product(collocation(section,:),control(:,column))
                residual = max(residual,abs(fitted_value-data(section,column)))
            end do
        end do
        data_scale = max(1.0_rk,maxval(abs(data)))
        ok = residual <= 4096.0_rk*epsilon(1.0_rk)*real(nsection,rk)*data_scale
    end subroutine
    !===============================================================================
end module forcad_geometry
