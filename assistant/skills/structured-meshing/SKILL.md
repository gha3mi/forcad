---
name: structured-meshing
description: Generate or replace structured meshes using native patch counts or physical target size, and verify conforming volume meshes.
---

# Size and verify a structured mesh

1. Resolve the patch, assembly, or existing mesh ID. Check source geometry and assembly interfaces. For a conforming hex mesh, verify that the source is a volume assembly whose interfaces can share nodes.
2. Use `preview_mesh` for estimated counts from a target size in model units. For exact counts, choose uniform sizing or supply per-patch counts in native parameter order. For conforming output, match counts along corresponding interface directions, including swapped or reversed directions. Use synchronization where supported.
3. Apply `mesh` to the exact source or mesh being replaced. Inspect the result for the requested mesh kind, sizing, and diagnostics. Use `set_mesh_view` to inspect the volume interior and `pick_mesh_cell` to verify a visible cell as needed. Report conforming output only when the native result supports that claim.
