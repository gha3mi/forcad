---
name: patch-editing
description: Inspect and edit native control coordinates, weights, knots, degrees, ranges, and exact patch sections.
---

# Edit native patches

1. Read the degrees, domains, knot vectors, and control data with `inspect_patch`. Use `inspect_assembly` for shared controls. Use the IDs returned for the intended target and keep patch-local control IDs distinct from assembly-wide control IDs.
2. Use `edit_patch` with the appropriate operation and active parameter direction. Control positions use model units. Knot and domain values use the patch’s native parameter coordinates. Use `set_probe` to evaluate a parameter position when a physical point needs confirmation.
3. Run `preview_controls` with the planned control or weight edits and review dependent changes. Batch related control edits when they should share an Undo step. Use `preview_refinement` before assembly knot insertion or degree elevation, using the same starting member and operation as the edit.
4. Reinspect changed quantities and check diagnostics. For a split, bounded extraction, face, or section, verify which native objects the operation produced.
