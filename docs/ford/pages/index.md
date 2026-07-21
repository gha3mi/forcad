---
title: NURBS Theory and ForCAD Usage
author: Seyed Ali Ghasemi
---

This chapter connects the mathematical definition of B-splines and NURBS to
the public ForCAD interface. It is intended both as a compact theoretical
reference and as a guide for constructing, evaluating, refining, analyzing,
and exporting geometry correctly.

The mathematical statements below use the conventions implemented by ForCAD.
Where the implementation intentionally provides an approximation, a restricted
model, or a local numerical candidate rather than a mathematical guarantee,
that distinction is stated explicitly.

## 1. Public interface and conventions

Applications should normally import the public facade rather than implementation
modules:

```fortran
use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume
```

Additional public facilities may be imported as needed:

```fortran
use forcad, only: extrude, revolve, sweep, loft
use forcad, only: nurbs_multipatch_curve, nurbs_multipatch_surface, &
    nurbs_multipatch_volume
use forcad, only: ndgrid, hexahedron_Xc, solve
```

The principal geometry types are [[nurbs_curve]], [[nurbs_surface]], and
[[nurbs_volume]]. They represent parametric mappings of rank one, two, and
three. Rank is not the same as physical dimension: a curve may be embedded in
two or three dimensions, and a surface used by the physical IGA routines is
embedded in three dimensions.

ForCAD uses the following conventions throughout its public API:

| Quantity | Convention |
| --- | --- |
| Real kind | All public floating-point geometry uses `rk`. |
| Degree | `degree=p` means polynomial degree \(p\); spline order is \(p+1\). |
| Index base | Public control, element, patch, and quadrature indices are one-based. |
| Coordinates | `Xc(A,d)` stores coordinate \(d\) of flattened control point \(A\). |
| Tensor order | The first parametric direction varies fastest. |
| Weights | Omitted weights mean \(w_A=1\); supplied weights must be finite and strictly positive. |
| Diagnostics | Recoverable failures are reported through the public `err` component. |

For a surface control index \((i,j)\) and a volume control index \((i,j,k)\),
the flattened one-based row is

\[
A_{2D}=i+n_1(j-1),
\qquad
A_{3D}=i+n_1\bigl((j-1)+n_2(k-1)\bigr).
\]

This ordering applies to control points, weights, tensor-product basis columns,
and element connectivity.

### 1.1 Minimal curve workflow

The following complete program constructs a nonuniform quadratic NURBS curve,
evaluates one point, and samples the active parameter interval:

```fortran
program first_nurbs_curve

    use forcad, only: rk, nurbs_curve

    implicit none

    type(nurbs_curve) :: curve
    real(rk), parameter :: knot(8) = [&
        0.0_rk,0.0_rk,0.0_rk,0.35_rk,0.80_rk,1.0_rk,1.0_rk,1.0_rk]
    real(rk), parameter :: Wc(5) = [1.0_rk,0.8_rk,1.3_rk,0.9_rk,1.0_rk]
    real(rk) :: Xc(5,3)
    real(rk), allocatable :: point(:)

    Xc(1,:) = [0.0_rk,0.0_rk,0.0_rk]
    Xc(2,:) = [0.8_rk,1.1_rk,0.2_rk]
    Xc(3,:) = [1.7_rk,0.4_rk,0.8_rk]
    Xc(4,:) = [2.6_rk,-0.5_rk,0.3_rk]
    Xc(5,:) = [3.4_rk,0.1_rk,0.0_rk]

    call curve%set(&
        knot   = knot,&
        Xc     = Xc,&
        Wc     = Wc,&
        degree = 2)

    if (curve%err%ok) then
        point = curve%cmp_Xg(0.42_rk)
        call curve%create(res = 101)
        write(*,"(a,3(1x,es14.6))") "C(0.42):", point
    end if

    call curve%finalize()

end program first_nurbs_curve
```

The relation `size(knot)=size(Xc,1)+degree+1` holds here:
\(8=5+2+1\). Passing `degree` explicitly is recommended for unclamped or
otherwise ambiguous knot vectors.

## 2. B-spline spaces and knot vectors

Let

\[
U=(u_0,u_1,\ldots,u_m),
\qquad
u_i\leq u_{i+1},
\]

be a finite nondecreasing knot vector. With \(n+1\) control variables and
degree \(p\geq0\), the fundamental size relation is

\[
m=n+p+1,
\qquad
|U|=(n+1)+p+1.
\]

ForCAD stores complete knot vectors, including repeated endpoint knots and any
extension knots outside the active interval. The active parameter interval is

\[
\Omega_\xi=[u_p,u_{n+1}].
\]

Only nonzero spans \([u_i,u_{i+1}]\) inside this interval define elements.
Repeated knots produce zero-length intervals, which are not elements.

### 2.1 Cox-de Boor basis

Using zero-based mathematical indices, the degree-zero basis is

\[
N_{i,0}(u)=
\begin{cases}
1, & u_i\leq u<u_{i+1},\\
0, & \text{otherwise},
\end{cases}
\]

with the active right endpoint evaluated by the conventional limiting value.
For \(p>0\),

\[
N_{i,p}(u)=
\frac{u-u_i}{u_{i+p}-u_i}N_{i,p-1}(u)
+
\frac{u_{i+p+1}-u}{u_{i+p+1}-u_{i+1}}N_{i+1,p-1}(u).
\]

A term with zero denominator is defined as zero. On a valid active interval,
the basis has the central properties

\[
N_{i,p}(u)\geq0,
\qquad
\sum_iN_{i,p}(u)=1,
\qquad
\operatorname{supp}N_{i,p}\subseteq[u_i,u_{i+p+1}].
\]

At a parameter that is not a knot, at most \(p+1\) univariate basis functions
are nonzero. This local support is the main reason spline evaluation and IGA
assembly can be performed without dense global basis work.

### 2.2 Multiplicity and continuity

If an interior knot has multiplicity \(s\), a degree-\(p\) B-spline basis is
generally \(C^{p-s}\) there:

| Multiplicity | Generic continuity |
| ---: | ---: |
| \(s=1\) | \(C^{p-1}\) |
| \(1<s<p\) | \(C^{p-s}\) |
| \(s=p\) | \(C^0\) |
| \(s=p+1\) | Discontinuous across the breakpoint |

This is a property of the spline space. A particular geometry can be smoother
than the space when its control data satisfy additional relations. Conversely,
derivatives above \(C^{p-s}\) at the knot do not possess a unique classical
two-sided value; an evaluator necessarily returns a side selected by its span
convention.

ForCAD reports stored-knot multiplicities and the corresponding \(p-s\)
metadata through `get_multiplicity` and `get_continuity`.

### 2.3 Knot-vector families

The descriptions below are independent properties, not mutually exclusive
enumeration values:

| Property | Mathematical meaning | ForCAD use |
| --- | --- | --- |
| Uniform | Consecutive distinct knots have constant spacing. | Supply the complete vector or use an open-uniform `set` overload. |
| Nonuniform | Distinct spans have unequal lengths. | Supply the complete vector explicitly. |
| Open/clamped | Each active endpoint has multiplicity \(p+1\). | End control points are interpolated for a valid clamped curve. |
| Unclamped | One or both active endpoints have lower multiplicity and exterior knots may be present. | Supply `degree` explicitly when inference is ambiguous. |
| Repeated interior | An interior value occurs \(s>1\) times. | Controls continuity through \(C^{p-s}\). |
| Bezier | No interior breakpoints and both endpoints have multiplicity \(p+1\). | `set(Xc=...,Wc=...)` constructs this form. |
| Periodic-form | Exterior knots and repeated control data define a periodic spline space and seam. | The caller supplies that structure; `wrap_parameters` only enables modulo evaluation. |

Parameter wrapping is not a periodic-geometry constructor. With active interval
\([a,b]\), wrapping maps an evaluation parameter to the same residue class
modulo \(b-a\). The caller remains responsible for compatible exterior knots,
control-point closure, weights, and derivative matching at the seam.

```fortran
call curve%set(&
    knot            = periodic_knot,&
    Xc              = periodic_Xc,&
    Wc              = periodic_Wc,&
    degree          = 2,&
    wrap_parameters = .true.)
```

See [[curve_knot_clamped]], [[curve_knot_unclamped]],
[[curve_knot_repeated]], and [[curve_knot_periodic]] for focused constructions.

### 2.4 Constructing a space from breakpoints

ForCAD can generate a knot vector from strictly increasing breakpoints and
requested continuity values. At breakpoint \(i\),

\[
s_i=p-c_i.
\]

Values \(-1\leq c_i\leq p-1\) retain the breakpoint with multiplicity
\(p-c_i\); `c_i=-1` therefore creates multiplicity \(p+1\). The value
`c_i=p` omits that breakpoint.

```fortran
call curve%set(&
    Xth_dir   = [0.0_rk,0.25_rk,0.60_rk,1.0_rk],&
    degree    = 3,&
    continuity = [-1,2,1,-1])
```

This call defines the spline space. If control points or weights are supplied,
their sizes must match the control count implied by the generated knots.

## 3. Rational basis functions and NURBS geometry

Let \(w_i>0\) be control weights. The univariate rational basis is

\[
R_i(u)=\frac{w_iN_{i,p}(u)}{W(u)},
\qquad
W(u)=\sum_jw_jN_{j,p}(u).
\]

The corresponding NURBS curve is

\[
\mathbf C(u)=\sum_iR_i(u)\mathbf P_i.
\]

Positive weights and partition of unity imply

\[
R_i(u)\geq0,
\qquad
\sum_iR_i(u)=1.
\]

Equal weights cancel from numerator and denominator, so the geometry reduces
to a polynomial B-spline. `is_rational()` reports whether nonuniform stored
weights affect the basis within a roundoff-scaled comparison; it does not mean
merely that a weight array was supplied.

### 3.1 Homogeneous coordinates

Define projective control points

\[
\mathbf H_i=(w_i\mathbf P_i,w_i).
\]

The polynomial B-spline in homogeneous space is

\[
\mathbf H(u)=\sum_iN_{i,p}(u)\mathbf H_i
             =(\mathbf X_w(u),W(u)),
\]

and dehomogenization gives

\[
\mathbf C(u)=\frac{\mathbf X_w(u)}{W(u)}.
\]

This formulation is fundamental. Exact conics, rational refinement, rational
lofting, and rational multipatch constraints must transform weighted
coordinates and weights consistently; applying a polynomial control-point
formula only to Euclidean coordinates is generally incorrect.

### 3.2 Rational derivatives

For \(A_i=w_iN_i\), univariate derivatives can be written recursively as

\[
R_i^{(k)}=
\frac{1}{W}
\left[
A_i^{(k)}-
\sum_{j=1}^{k}{k\choose j}W^{(j)}R_i^{(k-j)}
\right].
\]

Tensor-product mixed derivatives follow the corresponding multi-index
Leibniz recurrence. Unlike polynomial B-spline derivatives, rational basis
derivatives do not in general vanish merely because the total derivative order
exceeds the polynomial degree. Where the derivatives exist,

\[
\sum_A R_A=1,
\qquad
\sum_A \partial^{\boldsymbol\alpha}R_A=0
\quad\text{for }|\boldsymbol\alpha|>0.
\]

These identities are useful numerical checks, but their residual should be
interpreted relative to precision, derivative scale, degree, and knot spacing.

## 4. Curves, surfaces, and volumes

For tensor-product geometry, let

\[
B_{ij}(u,v)=N_{i,p}(u)M_{j,q}(v),
\qquad
B_{ijk}(u,v,w)=N_{i,p}(u)M_{j,q}(v)L_{k,r}(w).
\]

The rational tensor-product basis is

\[
R_A(\boldsymbol\xi)=
\frac{w_AB_A(\boldsymbol\xi)}
     {\sum_Bw_BB_B(\boldsymbol\xi)}.
\]

ForCAD represents

\[
\begin{aligned}
\mathbf C(u) &= \sum_iR_i(u)\mathbf P_i,\\
\mathbf S(u,v) &= \sum_{i,j}R_{ij}(u,v)\mathbf P_{ij},\\
\mathbf V(u,v,w) &= \sum_{i,j,k}R_{ijk}(u,v,w)\mathbf P_{ijk}.
\end{aligned}
\]

The essential public shapes are:

| Geometry | Directional controls | Flattened `Xc` | Parameter for scalar evaluation |
| --- | --- | --- | --- |
| Curve | \(n_1\) | `[n1,ncoord]` | scalar \(u\) |
| Surface | \((n_1,n_2)\) | `[n1*n2,ncoord]` | `[u,v]` |
| Volume | \((n_1,n_2,n_3)\) | `[n1*n2*n3,ncoord]` | `[u,v,w]` |

A Bezier or rational Bezier tensor product is constructed when only control
counts, controls, and optional weights are supplied:

```fortran
call surface%set(&
    nc = [3,3],&
    Xc = surface_Xc,&
    Wc = surface_Wc)

call volume%set(&
    nc = [2,2,2],&
    Xc = volume_Xc)
```

Explicit NURBS surfaces and volumes receive one complete knot vector per
direction:

```fortran
call surface%set(&
    knot1  = knot_u,&
    knot2  = knot_v,&
    Xc     = surface_Xc,&
    Wc     = surface_Wc,&
    degree = [p,q])

call volume%set(&
    knot1  = knot_u,&
    knot2  = knot_v,&
    knot3  = knot_w,&
    Xc     = volume_Xc,&
    Wc     = volume_Wc,&
    degree = [p,q,r])
```

Passing explicit degrees avoids ambiguity for non-open knot vectors. The row
count must equal the product of the directional control counts inferred from
the knot-size relations.

## 5. Evaluation, basis values, and derivatives

### 5.1 Geometry evaluation

`cmp_Xg` evaluates geometry directly and does not require a sampling cache:

```fortran
curve_point   = curve%cmp_Xg(0.35_rk)
surface_point = surface%cmp_Xg([0.35_rk,0.60_rk])
volume_point  = volume%cmp_Xg([0.35_rk,0.60_rk,0.20_rk])
```

Each scalar evaluation visits only the tensor product of locally active
univariate controls, nominally \(\prod_d(p_d+1)\). Outside the active interval,
an unwrapped parameter has zero basis support. With wrapping enabled, it is
first mapped modulo the active period.

### 5.2 Dense and active basis interfaces

The scalar `basis` overload returns one value per control variable:

```fortran
call surface%basis(&
    Xt  = [0.35_rk,0.60_rk],&
    Tgc = basis)
```

`derivative_order_active` returns only the local tensor-product basis values or
derivatives and the first active index in each direction:

```fortran
call surface%derivative_order_active(&
    Xt           = [0.35_rk,0.60_rk],&
    order        = [1,1],&
    first_active = first_active,&
    Dgc          = mixed_derivative)
```

For a surface, `mixed_derivative` has
\((p+1)(q+1)\) entries. `first_active=[i_0,j_0]` identifies the first active
control in each direction; local direction 1 varies fastest. Active interfaces
are preferred in element kernels because they avoid allocating and zeroing a
vector of total control-point size.

### 5.3 Derivative interfaces

| Binding | Result |
| --- | --- |
| `derivative` | All first parametric derivatives, optionally with basis values. |
| `derivative2` | The complete parametric Hessian, optionally with first derivatives and basis values. |
| `derivative_order` | One arbitrary derivative order for a curve or one mixed multi-index for a tensor product. |
| `derivative_order_active` | The same requested derivative restricted to local support. |

For a surface scalar evaluation, `derivative` returns `[ncp,2]`. For a volume
it returns `[ncp,3]`. The complete scalar Hessian is packed as

\[
\mathtt{d2Tgc((a-1)n_{cp}+A,b)}
=\frac{\partial^2R_A}{\partial\xi_a\partial\xi_b},
\]

with shape `[2*ncp,2]` for a surface and `[3*ncp,3]` for a volume. Both mixed
entries are populated. When only one mixed derivative is required,
`derivative_order` is clearer and uses less memory:

```fortran
call volume%derivative_order(&
    Xt    = [0.35_rk,0.60_rk,0.20_rk],&
    order = [1,0,2],&
    Dgc   = derivative_uww)
```

The curve equivalents use a scalar order. See [[curve_derivatives]],
[[surface_derivatives]], and [[volume_derivatives]].

### 5.4 Sampling and cache semantics

`create` samples a Cartesian parameter grid and caches parameters, geometry
points, and visualization connectivity:

```fortran
call curve%create(res = 121)
call surface%create(&
    res1 = 61,&
    res2 = 41)
call volume%create(&
    res1 = 31,&
    res2 = 25,&
    res3 = 17)
```

Directional arrays can be supplied instead of uniform resolutions. Cached data
are queried with `get_Xt`, `get_Xg`, and `get_ng`. `create` is required before
operations that consume cached geometry, including `export_Xg` and the
sample-based `nearest_point`.

Control modification, weight modification, fitting, transformation, and
refinement do not all rebuild every derived cache. After changing geometry,
call `create` again before using cached samples. After changing the spline
space or control count, also recompute any connectivity retained by the
application.

### 5.5 Mapping externally ordered points

`put_to_nurbs` is a cache-mapping operation, not fitting or closest-point
projection. For a curve it maps column 1 of an input array affinely from its
data range to the active parameter interval. A surface uses columns 1 and 2,
and a volume uses columns 1 through 3 independently:

\[
\xi_{id}=\xi_{a,d}+
\frac{X_{id}-X_{\min,d}}{X_{\max,d}-X_{\min,d}}
(\xi_{b,d}-\xi_{a,d}).
\]

It then evaluates one geometry point per input row and replaces the sampled
geometry cache. Each mapped input column must have a nonzero range. Remaining
columns affect neither parameters nor geometry evaluation. Use a fitting or
nearest-point method instead when input coordinates must be approximated or
projected geometrically.

## 6. Elements and isogeometric ansatz data

### 6.1 Nonzero-span elements

Every nonzero active knot span is a univariate element. Tensor products of
these spans form surface and volume elements. If direction \(d\) has \(n_{e,d}\)
nonzero spans, then

\[
n_e=\prod_dn_{e,d},
\qquad
n_{loc}=\prod_d(p_d+1).
\]

`cmp_elem()` returns `[nelem,nlocal]` one-based control connectivity with the
first parametric direction varying fastest:

```fortran
elem = volume%cmp_elem()
```

This IGA connectivity is distinct from control-net, sampled-geometry, and
parameter-grid visualization connectivity. The distinction is demonstrated by
[[curve_connectivity]], [[surface_connectivity]], and [[volume_connectivity]].

### 6.2 Mapping derivatives to physical space

For a volume mapping \(\mathbf x=\mathbf V(\boldsymbol\xi)\), define

\[
\mathbf J=\frac{\partial\mathbf x}{\partial\boldsymbol\xi}.
\]

The physical gradient and weighted volume contribution are

\[
\nabla_xR_A=\mathbf J^{-T}\nabla_\xi R_A,
\qquad
dV=\omega_g\,|\det\mathbf J|\,|\det\mathbf J_{e,\hat\xi}|.
\]

For a surface in three dimensions, with
\(\mathbf J=[\mathbf S_{,u}\ \mathbf S_{,v}]\) and metric
\(\mathbf G=\mathbf J^T\mathbf J\),

\[
dA=\omega_g\sqrt{\det\mathbf G}\,|\det\mathbf J_{e,\hat\xi}|,
\qquad
\nabla_\Gamma R_A=\mathbf J\mathbf G^{-1}\nabla_\xi R_A.
\]

For a curve,

\[
dL=\omega_g\|\mathbf C_{,u}\|\,|du/d\hat\xi|,
\qquad
\nabla_\Gamma R_A=
\frac{R_{A,u}}{\mathbf C_{,u}\cdot\mathbf C_{,u}}\mathbf C_{,u}.
\]

`ansatz` evaluates local basis values, physical gradients, an optional physical
Hessian, and the complete quadrature-weighted differential measure:

```fortran
call surface%ansatz(&
    ie         = element_id,&
    ig         = quadrature_id,&
    Tgc        = basis,&
    dTgc_dXg   = gradient,&
    d2Tgc_dXg2 = hessian,&
    dA         = dA,&
    ngauss     = [4,4])
```

Default quadrature counts are `degree+1` in each direction. `ig` is a flattened
one-based tensor quadrature index. Singular tangents, metrics, or Jacobians are
reported through the geometry diagnostic rather than silently inverted.

Typical assembly uses

\[
\mathbf K_e=\sum_g\mathbf B_g^T\mathbf D\mathbf B_g\,d\Omega_g,
\qquad
\mathbf f_e=\sum_g\mathbf R_g^T f(\mathbf x_g)\,d\Omega_g,
\]

and scatters these local contributions through `cmp_elem()`. ForCAD provides
the geometry, basis, derivatives, connectivity, and measure; the physical weak
form, boundary conditions, global matrix format, and scalable linear solver
remain application responsibilities. See [[curve_ansatz]], [[surface_ansatz]],
[[volume_ansatz]], [[surface_iga_poisson]], and [[volume_iga_poisson]].

## 7. Geometry-preserving refinement

Refinement changes the spline representation while seeking to preserve the
same geometry on the active domain.

### 7.1 Knot insertion

Inserting a knot enlarges the spline space without changing degree. For one
insertion \(\bar u\in[u_k,u_{k+1})\), the new controls have the local affine
form

\[
\mathbf Q_i=\alpha_i\mathbf P_i+(1-\alpha_i)\mathbf P_{i-1},
\]

with coefficients determined by knots and degree. Rational data are refined in
homogeneous coordinates. The requested final multiplicity must remain valid.

```fortran
call curve%insert_knots(&
    Xth = [0.25_rk,0.50_rk,0.75_rk],&
    r   = [1,1,1])

call surface%insert_knots(&
    dir = 1,&
    Xth = [0.30_rk,0.65_rk],&
    r   = [1,1])
```

### 7.2 Degree elevation

Degree elevation increases \(p\) while preserving the active-domain mapping.
Knot multiplicities are elevated consistently, and rational geometry is
transformed homogeneously:

```fortran
call curve%elevate_degree(t = 1)
call volume%elevate_degree(&
    dir = 2,&
    t   = 2)
```

For positive elevation, a nonclamped representation is converted to an
equivalent open-clamped representation on its active interval. Exterior
extension knots and periodic seam constraints are not retained. Therefore,
applications that depend on a particular periodic representation must rebuild
and verify that representation after elevation.

### 7.3 Knot removal

Knot removal attempts to reduce the representation while remaining within an
internal roundoff-scaled control test. Requests are processed in order, and the
longest accepted prefix of each request is removed. For rational geometry the
test acts on homogeneous controls; it is not a certified Euclidean Hausdorff
error bound.

```fortran
call curve%remove_knots(&
    Xth = [0.50_rk],&
    r   = [1])
```

Removal requires exact equality with an existing stored interior knot and does
not remove active endpoint knots.

### 7.4 Refinement operators and verification

The optional scalar operator `Bs` maps old Euclidean control coordinates to
new coordinates after the new weights are fixed. `B` is its coordinate-block
expansion. For rational data these are not the homogeneous refinement
operators and depend on old and new weights.

The correct regression test for refinement is geometric:

\[
e_\infty=
\max_{u\in\mathcal U}
\|\mathbf C_{before}(u)-\mathbf C_{after}(u)\|_\infty.
\]

Use a parameter set that includes span interiors, repeated knots approached
from valid sides, and active endpoints. Scale the acceptance tolerance with
coordinate magnitude, machine precision, degree, and the number and
conditioning of refinement operations. [[curve_refine]], [[surface_refine]],
and [[volume_refine]] demonstrate the rank-specific API.

## 8. Geometric constructions and transformations

The public generics [extrude](../interface/extrude.html),
[revolve](../interface/revolve.html), [sweep](../interface/sweep.html), and
[loft](../interface/loft.html) construct a higher-rank object.

| Construction | Mathematical model | Exactness and restriction |
| --- | --- | --- |
| `extrude` | \(\mathbf S(u,v)=\mathbf C(u)+v\mathbf d\) or \(\mathbf V(u,v,w)=\mathbf S(u,v)+w\mathbf d\) | Exact degree-one tensor direction. |
| `revolve` | Rotation about a 3D axis using rational quadratic arcs | Exact circular revolution; the angle is split into arcs no larger than \(\pi/2\). |
| `sweep` | \(\mathbf S(u,v)=\mathbf P(u)+\mathbf Q(v)-\mathbf o\) | Exact translational tensor product with a fixed profile frame; not a Frenet, rotation-minimizing, or guide-rail sweep. |
| `loft` | Homogeneous interpolation of compatible sections | Interpolates supplied section controls at supplied parameters; sections must have compatible degree, knots, dimensions, and wrapping policy. |

### 8.1 Built-in primitive geometries

The rank-specific types also provide compact constructors for common validated
control nets:

| Geometry rank | Constructors | Representation |
| --- | --- | --- |
| Curve | `set_circle`, `set_half_circle`, `set_C` | Exact planar rational conic templates. |
| Surface | `set_tetragon`, `set_ring`, `set_half_ring`, `set_C` | Bilinear rectangular or exact rational annular templates. |
| Volume | `set_hexahedron`, `set_ring`, `set_half_ring`, `set_C` | Tensor-product blocks or rational annular extrusions. |

For `set_circle`, the represented radius is `abs(radius)`. The legacy
`set_half_circle` template instead represents radius `abs(radius)/2`; this
distinction is intentional in the current API and must be accounted for when
matching a prescribed semicircle. Signed radius-like arguments can change
template orientation, and zero scales or zero extrusion lengths create
degenerate geometry rather than a regular domain.

```fortran
call curve%set_circle(&
    center = [0.0_rk,0.0_rk,0.0_rk],&
    radius = 2.0_rk)
call surface%set_ring(&
    center  = [0.0_rk,0.0_rk,0.0_rk],&
    radius1 = 1.0_rk,&
    radius2 = 2.0_rk)
call volume%set_hexahedron(&
    L  = [2.0_rk,1.0_rk,0.5_rk],&
    nc = [4,3,2])
```

### 8.2 Rank-raising constructor usage

Example construction:

```fortran
call profile%set_circle(&
    center = [0.0_rk,0.0_rk,0.0_rk],&
    radius = 1.0_rk)

cylinder = extrude(&
    curve  = profile,&
    vector = [0.0_rk,0.0_rk,3.0_rk])
```

For a revolution, Rodrigues' formula is used about unit axis
\(\widehat{\mathbf a}\):

\[
\mathcal R_\theta(\mathbf P)=\mathbf a_0+\mathbf h+
\cos\theta\,\mathbf r+
\sin\theta(\widehat{\mathbf a}\times\mathbf r),
\]

where \(\mathbf r\) is the radial component and \(\mathbf h\) the axial
component of \(\mathbf P-\mathbf a_0\).

Rotation and translation bindings act separately on controls (`rotate_Xc`,
`translate_Xc`) and cached samples (`rotate_Xg`, `translate_Xg`). Transforming
controls does not guarantee that an existing sample cache has been updated;
resample with `create` when the cache must represent the transformed geometry.

Focused examples are [[surface_extrude]], [[surface_revolve]],
[[surface_sweep]], [[surface_loft]], [[volume_extrude]], [[volume_revolve]],
[[volume_sweep]], and [[volume_loft]].

## 9. Fitting and nearest-point projection

### 9.1 Polynomial least-squares fitting

Given parameters \(\boldsymbol\xi_g\) and data
\(\mathbf y_g\), polynomial fitting minimizes

\[
\min_{\mathbf P}\|\mathbf A\mathbf P-\mathbf Y\|_F^2,
\qquad
A_{gA}=N_A(\boldsymbol\xi_g).
\]

`lsq_fit_bspline` uses the spline space already stored in the object and solves
the structured normal equations. This is memory-efficient for local support,
but normal equations square the condition number. Poor parameterization,
excessive degree, clustered knots, or rank-deficient data can therefore reduce
accuracy.

```fortran
call curve%lsq_fit_bspline(&
    Xt    = parameters,&
    Xdata = data_points,&
    ndata = size(parameters))
```

### 9.2 Nonlinear NURBS fitting

NURBS fitting also optimizes positive weights. ForCAD parameterizes weights in
logarithmic form, alternates control and weight updates, and uses damping and
regularization. The problem is nonlinear and generally nonconvex, so the result
is a converged numerical candidate, not a guaranteed global optimum. Inspect
`err`, residuals, weight range, and sensitivity to initial space and
regularization. See [[curve_fit_bspline]], [[curve_fit_nurbs]], and the
corresponding surface and volume examples.

### 9.3 Nearest points

The sample-based method

```fortran
call curve%nearest_point(&
    point_Xg   = query,&
    nearest_Xg = nearest_point,&
    nearest_Xt = nearest_parameter,&
    id         = sample_id)
```

searches only cached samples, so its error is limited by the sampling grid.
Call `create` first.

The continuous method minimizes

\[
f(\boldsymbol\xi)=
\frac12\|\mathbf x(\boldsymbol\xi)-\mathbf x_q\|^2
\]

on the active parameter box using streamed coarse seeds, bounded Newton steps,
and backtracking. It returns a bounded local candidate. It is not a certified
global closest point, and convergence can fail near singular Jacobians,
multiple equally near regions, or poorly scaled parameterizations.

```fortran
call curve%nearest_point2(&
    point_Xg   = query,&
    tol        = 1.0e-11_rk,&
    maxit      = 40,&
    nearest_Xt = nearest_parameter,&
    nearest_Xg = nearest_point)
```

Use multiple independent initial regions or a problem-specific spatial search
when a global certificate is required. See [[curve_nearest]],
[[surface_nearest]], and [[volume_nearest]].

## 10. Geometric measures and quadrature

ForCAD integrates each nonzero active span or span cell independently with
Gauss-Legendre quadrature:

\[
L_h=\sum_{e,g}\omega_g\|\mathbf C'(u_{eg})\|\,|J_e|,
\]

\[
A_h=\sum_{e,g}\omega_g
\|\mathbf S_{,u}\times\mathbf S_{,v}\|\,|J_e|,
\]

\[
V_h=\sum_{e,g}\omega_g
|\det\mathbf V_{,(u,v,w)}|\,|J_e|.
\]

```fortran
call curve%cmp_length(&
    length = length,&
    ngauss = 8)
call surface%cmp_area(&
    area    = area,&
    ngauss  = [6,6])
call volume%cmp_volume(&
    volume  = volume_measure,&
    ngauss  = [5,5,5])
```

The default is `degree+1` points per direction. Although Gaussian quadrature is
exact for sufficiently low-degree polynomial integrands, arc length, area, and
the absolute Jacobian determinant are generally nonpolynomial, even for
polynomial geometry. These routines therefore return numerical approximations,
not universal exact measures. Increase quadrature order and verify convergence.

The volume routine uses an absolute determinant, so inverted or folded cells
can contribute positive measure. A positive result is not proof that a mapping
is globally injective. See [[curve_cmp_length]], [[surface_cmp_area]], and
[[volume_cmp_volume]].

## 11. Multipatch topology and continuity

The multipatch types own copies of complete single-patch objects and record
oriented interfaces:

| Rank | Type | Boundary labels |
| --- | --- | --- |
| Curve | [[nurbs_multipatch_curve]] | `left`, `right` |
| Surface | [[nurbs_multipatch_surface]] | `u_min`, `u_max`, `v_min`, `v_max` |
| Volume | [[nurbs_multipatch_volume]] | Surface labels plus `w_min`, `w_max` |

`add_patch` returns a one-based patch identifier. `connect` validates metadata
and trace-space compatibility. `is_valid` must then be called before analysis
to verify geometry traces, orientation, rational projective compatibility, and
the requested derivative residuals.

### 11.1 Parametric continuity

Let \(\eta_A\) and \(\eta_B\) be normalized local increments along the two
stored normal parameter axes. ForCAD's parametric \(C^n\) convention imposes

\[
\left.\frac{\partial^r f_A}{\partial\eta_A^r}\right|_\Gamma
=
\sigma^r
\left.\frac{\partial^r f_B}{\partial\eta_B^r}\right|_\Gamma,
\qquad r=0,\ldots,n,
\]

after applying tangential orientation. The sign \(\sigma\) accounts for whether
the two stored side axes point consistently across the interface. Normalizing
by each active parameter interval includes the required interval-scale factors.

### 11.2 Geometric continuity

For geometric continuity, patch B is composed with a scalar normal transition
\(\eta_B=\phi(\eta_A)\). The order-\(r\) chain rule is

\[
\frac{d^r}{d\eta_A^r}f_B(\phi(\eta_A))
=\sum_{k=0}^{r}B_{r,k}
\bigl(\phi'(0),\ldots,\phi^{(r-k+1)}(0)\bigr)
\frac{\partial^kf_B}{\partial\eta_B^k},
\]

where \(B_{r,k}\) is a partial exponential Bell polynomial. The user supplies
the transition jet

\[
[\phi'(0),\phi''(0),\ldots,\phi^{(n)}(0)]
\]

through `reparameterization` and sets `geometric=.true.`.

For curves this is the complete local univariate reparameterization model. For
surfaces and volumes, ForCAD implements the restricted separable map

\[
\Phi(\boldsymbol\tau,\eta)=
(T\boldsymbol\tau,\phi(\eta)).
\]

It does not represent all multivariate \(G^n\) joins: tangentially varying
normal jets and normal-dependent tangential coordinates are outside the current
model. Rational constraints use a homogeneous lift with one constant
projective interface scale; the most general nonconstant projective rescaling
is not implemented.

### 11.3 Degrees of freedom and sparse constraints

`cmp_dof_map` merges only corresponding \(C^0\) interface controls into compact
global identifiers. Higher-order continuity remains an explicit homogeneous
constraint system

\[
\mathbf A\mathbf q=\mathbf0.
\]

`cmp_dof_constraint` returns one-based CSR arrays `rowptr`, `col`, and `val`.
The caller may eliminate constraints, use Lagrange multipliers, or pass them to
an external constrained solver. For rational geometry, raw rows apply to
homogeneous controls. A rational scalar field uses weighted numerator
coefficients \(w_iq_i\), together with a compatible denominator.

```fortran
call patches%add_patch(&
    patch = left,&
    id    = left_id)
call patches%add_patch(&
    patch = right,&
    id    = right_id)
call patches%connect(&
    patch_a    = left_id,&
    side_a     = "u_max",&
    patch_b    = right_id,&
    side_b     = "u_min",&
    continuity = 1)

if (patches%is_valid()) then
    dof_map = patches%cmp_dof_map()
    call patches%cmp_dof_constraint(&
        rowptr = rowptr,&
        col    = col,&
        val    = val)
end if
```

Surface edges can be reversed. Volume faces can additionally exchange their
two tangential axes and flip either mapped axis. The connection must encode
these orientations correctly before controls, derivatives, or element degrees
of freedom are compared.

See [[curve_multipatch]], [[surface_multipatch]], [[volume_multipatch]],
[[curve_iga_beam_c1]], [[surface_multipatch_gn]], and
[[volume_multipatch_gn]].

## 12. Export and visualization

Single-patch geometry can export four distinct representations:

| Method | Data represented |
| --- | --- |
| `export_Xc` | Control polygon, net, or lattice. |
| `export_Xg` | Geometry samples cached by `create`. |
| `export_Xth` | Parameter-space span mesh. |
| `export_Xth_in_Xg` | Parameter-span lines or cells evaluated in physical space. |

VTK point fields are supplied as one column per field with matching names:

```fortran
call surface%create(&
    res1 = 81,&
    res2 = 41)
call surface%export_Xg(&
    filename    = "vtk/surface_temperature.vtk",&
    point_data  = temperature,&
    field_names = ["temperature"])
```

`show` visualizes files that have already been exported; it does not perform
the export. Scalar fields in the sampled `Xg` file are discovered by the
PyVista interface, where the user can select fields and colormaps.

```fortran
call surface%export_Xc("vtk/surface_Xc.vtk")
call surface%export_Xth_in_Xg("vtk/surface_Xth.vtk")
call surface%show(&
    vtkfile_Xc        = "vtk/surface_Xc.vtk",&
    vtkfile_Xg        = "vtk/surface_temperature.vtk",&
    vtkfile_Xth_in_Xg = "vtk/surface_Xth.vtk")
```

Curves and surfaces export IGES entities 126 and 128. Trivariate NURBS volume
export is not represented by the current IGES backend; `export_iges` on a
volume reports this as a recoverable diagnostic. Multipatch VTK methods write
one file per patch with a common stem, preserving patch ownership for
post-processing.

## 13. Diagnostics and numerical contracts

Every geometry and multipatch object owns a structured `err` state. The normal
workflow is:

1. Construct or mutate the object.
2. Check `object%err%ok` before consuming outputs.
3. Treat zero-sized allocatable outputs or zero measures after failure as
   failure results, not valid geometry.
4. Resolve or reset the diagnostic according to the `fordebug` contract before
   continuing with methods that require a clear state.

Important numerical distinctions are:

| Operation | What is guaranteed or returned |
| --- | --- |
| Basis evaluation | A value under the stored span and wrapping conventions. |
| Knot insertion | Geometry preservation to floating-point roundoff for valid input. |
| Degree elevation | Active-domain geometry preservation; nonclamped representation may change. |
| Knot removal | Acceptance by an internal homogeneous/control test, not a certified Euclidean bound. |
| Least-squares fit | A numerical least-squares solution subject to conditioning and solver checks. |
| NURBS fit | A damped nonlinear candidate, not a global optimum. |
| `nearest_point` | Nearest cached sample. |
| `nearest_point2` | Bounded local numerical candidate. |
| Length/area/volume | Elementwise quadrature approximation. |
| Multipatch \(G^n\) | Full scalar transition for curves; restricted normal-only model for surfaces and volumes. |

Tolerance selection should be dimensionally meaningful. For a geometry scale
\(L\), a basic roundoff scale is

\[
\tau=C\,\epsilon_{rk}\max(1,L),
\]

where \(C\) must also account for degree, derivative order, knot-span size,
conditioning, and the number of accumulated operations. A fixed absolute
tolerance cannot be expected to work equally for micrometer and kilometer
models.

## 14. Accuracy, memory, and performance guidance

The mathematically fastest interface depends on the requested output:

| Task | Preferred approach |
| --- | --- |
| One geometry point | `cmp_Xg`; no sample cache is needed. |
| Element assembly | `cmp_elem` plus `ansatz` or active derivatives. |
| One mixed derivative | `derivative_order_active` when global zeros are unnecessary. |
| Visualization grid | One `create`, then reuse cached `Xg` and connectivity. |
| Repeated integration | Reuse geometry and choose the lowest quadrature order that passes a convergence study. |
| Large multipatch system | Use CSR constraints and a sparse external PDE solver. |

For one point in rank \(r\), the active tensor basis has

\[
n_{loc}=\prod_{d=1}^{r}(p_d+1)
\]

entries. Direct point evaluation and active derivatives scale with local
support, while a dense basis result necessarily stores all global control
columns. For very large control nets, repeatedly requesting dense basis or
Hessian arrays can dominate both memory traffic and allocation cost.

Sampling stores \(O(n_g n_{coord})\) geometry data in addition to parameter and
connectivity arrays. Do not call `create` solely to evaluate isolated points.
Likewise, optional refinement matrices can be much larger than the refined
geometry itself; request them only when the application needs an explicit
control-variable transfer operator.

The library uses standard-conforming `pure`, `elemental`, and `do concurrent`
kernels where their independence requirements hold. Parallel compiler flags
change execution strategy, not mathematical semantics. Validate numerical
agreement between debug and optimized builds, and do not use unsafe
floating-point reassociation when reproducibility or robust diagnostics matter.

### 14.1 Selected public utilities

Three lower-level facilities are intentionally re-exported through the facade:

| Utility | Contract and appropriate scale |
| --- | --- |
| `ndgrid` | Forms all pairs or triples from directional vectors with direction 1 varying fastest. The result can be passed to scalar-point APIs or used to reproduce tensor ordering. |
| `hexahedron_Xc` | Generates an axis-aligned, uniformly spaced control lattice on signed extents. Its input-shape and count preconditions are low-level; `volume%set_hexahedron` provides object diagnostics. |
| `solve` | Solves a small or medium dense system. It attempts unpivoted Cholesky for a numerically symmetric matrix and otherwise falls back to partial-pivoting LU. |

`solve` allocates dense work, costs \(O(n^3+n^2n_{rhs})\), does not estimate a
condition number, and performs no iterative refinement. A shape or accepted
pivot failure returns a zero-sized result. It is useful for compact examples
and local dense systems, but it is not a sparse, distributed, or production
large-scale IGA solver.

```fortran
call ndgrid(&
    X_dir1 = parameters_u,&
    X_dir2 = parameters_v,&
    Xt     = parameter_pairs)

solution = solve(&
    A = matrix,&
    B = right_hand_sides)
```

## 15. Choosing the relevant example

The examples are intentionally narrower than this chapter. They are executable
usage references for individual public features:

| Topic | Curve | Surface | Volume |
| --- | --- | --- | --- |
| Basic construction | [[curve_basic]] | [[surface_basic]] | [[volume_basic]] |
| Basis functions | [[curve_basis]] | [[surface_basis]] | [[volume_basis]] |
| Derivatives | [[curve_derivatives]] | [[surface_derivatives]] | [[volume_derivatives]] |
| Element connectivity | [[curve_connectivity]] | [[surface_connectivity]] | [[volume_connectivity]] |
| IGA ansatz | [[curve_ansatz]] | [[surface_ansatz]] | [[volume_ansatz]] |
| Refinement | [[curve_refine]] | [[surface_refine]] | [[volume_refine]] |
| Nearest point | [[curve_nearest]] | [[surface_nearest]] | [[volume_nearest]] |
| Least-squares fit | [[curve_fit_bspline]] | [[surface_fit_bspline]] | [[volume_fit_bspline]] |
| NURBS fit | [[curve_fit_nurbs]] | [[surface_fit_nurbs]] | [[volume_fit_nurbs]] |
| Geometric measure | [[curve_cmp_length]] | [[surface_cmp_area]] | [[volume_cmp_volume]] |
| Multipatch | [[curve_multipatch]] | [[surface_multipatch]] | [[volume_multipatch]] |
| \(G^n\) multipatch | [[curve_multipatch_gn]] | [[surface_multipatch_gn]] | [[volume_multipatch_gn]] |

For a new scientific application, a reliable progression is:

1. Verify the knot-size, degree, control-count, and positive-weight invariants.
2. Check partition of unity and derivative sums at representative parameters.
3. Compare `cmp_Xg` against known points or exact geometry.
4. Validate element connectivity and Jacobian signs or metric rank.
5. Perform refinement-invariance and quadrature-convergence studies.
6. For multipatches, call `is_valid` and measure CSR constraint residuals.
7. Only then benchmark optimized and parallel builds on representative sizes.

## 16. References

The implementation documentation identifies the algorithm used by each
nontrivial routine. The following sources provide the principal theoretical
foundation for this chapter:

1. L. Piegl and W. Tiller, *The NURBS Book*, 2nd ed., Springer, 1997.
   [doi:10.1007/978-3-642-59223-2](https://doi.org/10.1007/978-3-642-59223-2).
2. M. G. Cox, "The Numerical Evaluation of B-Splines," *Journal of the
   Institute of Mathematics and Its Applications* 10 (1972), 134-149.
   [doi:10.1093/imamat/10.2.134](https://doi.org/10.1093/imamat/10.2.134).
3. C. de Boor, "On Calculating with B-Splines," *Journal of Approximation
   Theory* 6 (1972), 50-62.
   [doi:10.1016/0021-9045(72)90080-9](https://doi.org/10.1016/0021-9045(72)90080-9).
4. W. Boehm, "Inserting New Knots into B-Spline Curves,"
   *Computer-Aided Design* 12 (1980), 199-201.
   [doi:10.1016/0010-4485(80)90154-2](https://doi.org/10.1016/0010-4485(80)90154-2).
5. H. Prautzsch, "Degree Elevation of B-Spline Curves,"
   *Computer Aided Geometric Design* 1 (1984), 193-198.
   [doi:10.1016/0167-8396(84)90031-1](https://doi.org/10.1016/0167-8396(84)90031-1).
6. W. Tiller, "Knot-removal Algorithms for NURBS Curves and Surfaces,"
   *Computer-Aided Design* 24 (1992), 445-453.
   [doi:10.1016/0010-4485(92)90012-Y](https://doi.org/10.1016/0010-4485(92)90012-Y).
7. T. J. R. Hughes, J. A. Cottrell, and Y. Bazilevs, "Isogeometric Analysis:
   CAD, Finite Elements, NURBS, Exact Geometry and Mesh Refinement,"
   *Computer Methods in Applied Mechanics and Engineering* 194 (2005),
   4135-4195.
   [doi:10.1016/j.cma.2004.10.008](https://doi.org/10.1016/j.cma.2004.10.008).
8. P. Kiciak, "Conditions for Geometric Continuity between Polynomial and
   Rational Surface Patches," *Computer Aided Geometric Design* 13 (1996),
   709-741.
   [doi:10.1016/0167-8396(96)00006-4](https://doi.org/10.1016/0167-8396(96)00006-4).
9. L. Armijo, "Minimization of Functions Having Lipschitz Continuous First
   Partial Derivatives," *Pacific Journal of Mathematics* 16 (1966), 1-3.
   [doi:10.2140/pjm.1966.16.1](https://doi.org/10.2140/pjm.1966.16.1).
