title: Public Fortran API

Import the names your program uses from the [`forcad` module](../../module/forcad.html). Use `real(rk)` for model data and the `_rk` suffix for real literals.

```fortran
use forcad, only: rk, nurbs_curve
```

[TOC]

## Geometry and construction

| Task | Public interface |
| --- | --- |
| Define, evaluate, edit, or measure a patch | [`nurbs_curve`](../../type/nurbs_curve.html), [`nurbs_surface`](../../type/nurbs_surface.html), [`nurbs_volume`](../../type/nurbs_volume.html) |
| Construct geometry from existing profiles | [`extrude`](../../interface/extrude.html), [`revolve`](../../interface/revolve.html), [`sweep`](../../interface/sweep.html), [`loft`](../../interface/loft.html) |
| Group patches and register their interfaces | [`nurbs_multipatch_curve`](../../type/nurbs_multipatch_curve.html), [`nurbs_multipatch_surface`](../../type/nurbs_multipatch_surface.html), [`nurbs_multipatch_volume`](../../type/nurbs_multipatch_volume.html), [`multipatch_connection`](../../type/multipatch_connection.html) |
| Integrate a physical quantity | [Curve length](../../type/nurbs_curve.html#boundprocedure-cmp_length), [surface area](../../type/nurbs_surface.html#boundprocedure-cmp_area), [volume](../../type/nurbs_volume.html#boundprocedure-cmp_volume) |

Call a generic such as `extrude` by its public name. Its reference page lists the argument types for each supported overload.

## Grids and numerical utilities

| Name | Result |
| --- | --- |
| [`ndgrid`](../../interface/ndgrid.html) | Tensor grid coordinates from two or three coordinate vectors |
| [`hexahedron_Xc`](../../proc/hexahedron_xc.html) | Control coordinates for an axis-aligned rectangular volume |
| [`rotation`](../../proc/rotation.html) | A three-dimensional rotation matrix from three angles in degrees |
| [`solve`](../../proc/solve.html) | Solution of a dense square system \(AX=B\) with one or more right-hand sides |
| [`rk`](../../module/forcad_kinds.html#variable-rk) | Working real kind for coordinates, weights, and parameters |

## Patch data and evaluation

Let \(N\) be the total number of control points and \(d\) the number of physical coordinates per point. For a surface with \(n_u\) controls along U and \(n_v\) along V, \(N=n_u n_v\). For a volume, \(N=n_u n_v n_w\).

| Input to `set` | Shape and convention |
| --- | --- |
| Control points, `Xc` | A `real(rk)` array of shape `(N, d)`. Each row is a point. U varies fastest, followed by V, then W. |
| Weights, `Wc` | An optional `real(rk)` vector of length \(N\) in the same order as `Xc`. Each weight must be finite and positive. Omit `Wc` for unit weights. |
| Explicit knots | For \(n\) control points and degree \(p\) in one direction, supply \(n+p+1\) nondecreasing knots. With one-based indexing, the active interval is \([\xi_{p+1},\,\xi_{n+1}]\). |

For a volume control with one-based indices \(i_u,i_v,i_w\), the corresponding row of `Xc` and `Wc` is

\[
i=i_u+n_u(i_v-1)+n_u n_v(i_w-1).
\]

For a surface, omit the last term. For example, the volume control \((2,2,2)\) is row 8 when each direction has two controls. [Move this control in Fortran](../fortran/move-volume-control.html).

A polynomial patch can store its unit weights implicitly. For a valid patch with no stored weights, `get_Wc()` returns an empty array. Studio shows their effective value of one in the control editor.

Use `cmp_Xg` to evaluate positions at specified parameters. Use `create` to retain a sampled grid for `export_Xg`. After `modify_Xc` or `modify_Wc`, call `create()` to refresh the positions on a previously sampled grid.

See [Evaluate a surface point and its tangents](../fortran/evaluate-a-surface.html) for physical derivatives and [Evaluate the rational basis](../fortran/evaluate-rational-basis.html) for the contribution of each control point.

## Point queries, fitting, and analysis

| Public call | What it returns | Worked example |
| --- | --- | --- |
| `nearest_point` / `nearest_point2` | A cached sample index and position / a location on the continuous patch with its parameters | [Compare the two surface queries](../fortran/project-a-point.html), [locate a point in a volume](../fortran/locate-a-point-in-a-volume.html) |
| `lsq_fit_bspline` | Updated controls fitted to physical observations at supplied parameters | [Fit a profile](../fortran/fit-a-curve.html), [fit panel heights on a grid](../fortran/fit-a-surface.html) |
| `ansatz_element` | Local basis values, physical gradients and quadrature weights, with optional geometry derivatives | [Integrate a field](../fortran/integrate-a-surface-field.html), [differentiate arc length](../fortran/curve-length-sensitivity.html) |
| `assemble_sensitivities` | Patch-local derivatives combined for independent multipatch controls and log weights | [Length gradient at a shared control](../fortran/shared-control-sensitivity.html), [volume gradients at a shared face](../fortran/shared-face-volume-sensitivity.html) |
| `cmp_elem(shared=.true.)` / `interface_quadrature` | Global element connectivity / paired interface parameters and integration weights | [Solve across two volumes and verify face flux](../iga/solve-two-volumes-poisson.html) |

## Diagnostics

Check `patch%err%ok` after a construction or edit. If it is false, inspect the diagnostic, correct the input, and call `patch%err%reset()` before retrying. See [Evaluate your first NURBS curve](../tutorials/first-curve.html) for a complete program and [Studio and Fortran equivalents](studio-fortran.html) for the matching Studio commands.
