---
name: multipatch
description: Design, join, inspect, and refine native assemblies while preserving interface orientation and continuity.
---

# Assemble and validate patches

1. Identify member patches of the same dimension and the boundaries intended to meet. Use `inspect_patch` for native domains and `inspect_assembly` for existing members, joins, and shared controls. Read interface diagnostics before editing an existing assembly.
2. Use `assemble` with explicit members and interfaces. Verify the paired boundaries and tangential orientation from the actual geometry. Specify only the continuity the design supports. For a geometric continuity constraint, use the required normal transition derivatives and verify the result in ForCAD Studio.
3. Inspect the updated definition and check reported validation. Report unjoined members or unconstrained interfaces as such. When changing a connected assembly, preview control or refinement adjustments before applying them and verify the affected member IDs afterward. Check the linked mesh state when the edit affects its source.
