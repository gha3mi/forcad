---
name: project-exchange
description: Work with ForCAD Studio projects, tabs, and project-relative CAD imports, exports, and workspace captures.
---

# Exchange files and manage project tabs

1. Use `list_project_tabs` to identify the active project and unsaved edits. Use its current tab ID and revision for `switch_project_tab` or `close_project_tab`. Use `new_project` or `open_project` with a separate tab when the current scene must remain open. Inspect again before using object IDs in the active project.
2. Use `import_file` to merge an FCAD project or import IGES into the current design. Confirm the IGES source units or the intent to preserve the numeric coordinates. Verify the imported objects and their scale.
3. Use `save_project` to save the working project. Use `export_file` with a new project-relative path for a supported selection export, project snapshot, or viewport capture. Use `.forcad/screenshots/<name>.png` for viewport captures unless the user chooses a location. Choose the format for the actual content. Confirm the returned path and report it to the user.
