title: Studio and Fortran equivalents

Find a Studio command and a worked example using the public ForCAD API. A [curve](../../type/nurbs_curve.html) has one parameter direction. On a [surface](../../type/nurbs_surface.html), U and V correspond to `dir=1` and `dir=2`. A [volume](../../type/nurbs_volume.html) adds W at `dir=3`.

[TOC]

## Evaluate a patch

| ForCAD Studio | Public Fortran API |
| --- | --- |
| **Inspector → Selection and probe** · surface position and tangents | [Evaluate `cmp_Xg` at a parameter pair](../fortran/evaluate-a-surface.html). The optional derivative orders select the U or V tangent. |
| **Patch → Parameter panels** · basis plot | [Evaluate a curve’s rational basis](../fortran/evaluate-rational-basis.html) with `basis`. |
| **Snap → Snap to closest point** · while dragging a control point | [Project a physical point onto a surface](../fortran/project-a-point.html) with `nearest_point2`. For a solid, Studio snaps to an exterior face. A [Fortran volume query](../fortran/locate-a-point-in-a-volume.html) also searches its interior. |

## Edit a patch

| ForCAD Studio | Public Fortran API |
| --- | --- |
| **Patch → Points + weights** · [Move a control point](../how-to/change-control-point.html) | `modify_Xc` replaces one coordinate by control index and physical direction, or accepts the whole control array. [Move the Z coordinate of a volume corner](../fortran/move-volume-control.html). |
| **Patch → Points + weights** · [Change a weight](../how-to/change-weight.html) | `call curve%modify_Wc(W=2.0_rk, num=3)`. [Run the same edit in Fortran](../fortran/change-curve-weight.html). |
| **Patch → Insert knot** · [Set value and copies](../how-to/insert-knot.html) | `call curve%insert_knots(Xth=[0.5_rk], r=[1])`. [Check the midpoint after insertion](../fortran/insert-a-knot.html). |
| [**Patch → Remove knot**](../how-to/remove-knot.html) | [Read the accepted count](../fortran/remove-a-knot.html) from `remove_knots` when only one copy can be removed. |
| [**Patch → Elevate degree**](../how-to/elevate-degree.html) | [Raise the U degree of an annulus](../fortran/elevate-degree.html) with `elevate_degree`. |
| [**Patch → Reduce degree**](../how-to/reduce-degree.html) | [Read the accepted reduction and error bound](../fortran/reduce-degree.html) from `reduce_degree`. |
| [**Patch → Refine patch**](../how-to/refine-patch.html) | [Subdivide each knot span](../fortran/refine-a-patch.html) with `insert_knots`. For a degree increase in the same edit, call `elevate_degree` first. |
| [**Patch → Set domain**](../how-to/set-domain.html) | [Map a volume's U interval to new bounds](../fortran/set-domain.html) with `reparameterize`. |
| **Patch → Split patch** · [Divide a patch](../how-to/split-a-patch.html) | `call curve%split(u, lower, upper)` produces two patches. Surface and volume splits also take a parameter direction. [Check both halves of a unit cube](../fortran/split-a-volume.html). |
| **Model → Measure** · [Measure the selection](../how-to/measure-geometry.html) | `cmp_length`, `cmp_area`, or `cmp_volume` integrates a patch of the corresponding rank. [Measurement examples](../fortran/measure-a-patch.html). |

## Construct geometry

| ForCAD Studio | Public Fortran API |
| --- | --- |
| [**Create → Extrude**](../how-to/extrude-a-profile.html) | [`extrude`](../../interface/extrude.html) constructs a surface from a curve or a volume from a surface. [Run a volume extrusion](../fortran/extrude-a-surface.html). |
| [**Create → Revolve**](../how-to/revolve-a-profile.html) | [`revolve`](../../interface/revolve.html) revolves a profile about an axis. [Build a quarter-cylinder wall](../fortran/revolve-a-profile.html). |
| [**Create → Sweep**](../how-to/sweep-a-profile.html) | [`sweep`](../../interface/sweep.html) translates a profile along a path. [Sweep a tilted line along a circle](../fortran/sweep-a-profile.html). |
| [**Create → Loft**](../how-to/loft-sections.html) | [`loft`](../../interface/loft.html) builds geometry through ordered sections. [Loft between two circles](../fortran/loft-sections.html). |

## Join and mesh patches

| ForCAD Studio | Public Fortran API |
| --- | --- |
| [**Multipatch → New multipatch**](../how-to/join-volume-patches.html) | For a volume assembly, use [`nurbs_multipatch_volume`](../../type/nurbs_multipatch_volume.html). [Register two members and connect their U faces](../fortran/join-volume-patches.html). |
| **Multipatch → Check interfaces** | The assembly's `is_valid` checks the registered connections. |
| **Mesh → Generate mesh** · [one volume](../how-to/mesh-a-volume.html) | `create` samples the patch. `get_Xg` returns nodes and `get_elem_Xg_vis` returns cell connectivity. [Generate 24 hex cells in Fortran](../fortran/sample-volume-mesh.html). |
| **Mesh → Generate mesh** · [conforming volume assembly](../how-to/mesh-a-multipatch.html) | Assembly `create` accepts per-patch sampling counts or a target size and can synchronize joined directions. Studio then assigns common node IDs on shared faces. [Export the complete reconstruction](../how-to/export-fortran.html). |

## Exchange and project operations

A dash indicates a Studio operation with no direct public ForCAD call.

| ForCAD Studio | Public Fortran API |
| --- | --- |
| **Project → Export IGES** · [curve or surface](../how-to/export-iges.html) | [Write a circle with `export_iges`](../fortran/export-curve-iges.html). |
| **Project → Export VTK** · [sampled patch geometry](../how-to/export-geometry-vtk.html) | [Write a single sampled patch with `create` and `export_Xg`](../fortran/export-volume-vtk.html). |
| **Mesh → Export mesh (VTK)** · [generated mesh](../how-to/export-mesh-vtk.html) | — |
| **Project → Save** · [editable FCAD project](../how-to/save-project.html) | — |
| **Project → Export Fortran** · [generate a standalone program](../how-to/export-fortran.html) | — |
| **Undo** a project edit | — |
