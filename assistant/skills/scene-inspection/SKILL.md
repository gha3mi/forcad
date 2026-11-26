---
name: scene-inspection
description: Read the live design, identify exact patch or mesh features, and inspect or adjust the ForCAD Studio viewport and selection.
---

# Inspect the scene

1. Use `inspect` for objects, selection, units, mesh summaries, and diagnostics. Follow `next_offset` when the needed collections have more results. Use `inspect_patch` for exact patch data, `inspect_assembly` for members and shared controls, and `measure` for physical measurements in the display units.
2. Read `inspect_view` for the actual camera, drawable bounds, panels, and display settings. Use `select` for one object or `set_selection` for exact collections and subfeatures. Use `select_region` only when selecting by viewport position. Its points are canvas-relative, and box drag direction changes the selection rule.
3. Use `set_camera`, `set_workspace`, and `set_parameter_view` for view changes. Discover current dialog-free commands with `list_commands` and pass a returned command name to `execute`. Inspect again to verify the resulting view or selection. Keep viewport selection and camera changes separate from geometry edits.
