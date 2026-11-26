module studio_command_presentation
    use studio_commands, only: &
        find_command, interaction_guide_command, about_command, toggle_assistant_command, &
        diagnose_with_ai_command, fullscreen_command, canvas_command, restore_layout_command, &
        reset_controls_command, reset_presentation_command, home_command, new_project_command, &
        open_project_command, save_project_command, save_as_command, close_project_tab_command, &
        project_folder_command, units_command, undo_patch_command, redo_patch_command, screenshot_command, &
        import_project_command, import_iges_command, export_vtk_command, export_iges_command, &
        export_fortran_command, copy_geometry_command, paste_geometry_command, fit_scene_command, &
        fit_around_panels_command, select_tool_command, orbit_tool_command, clear_selection_command, &
        box_selection_command, lasso_selection_command, project_isolate_command, project_show_all_command, &
        project_hide_selection_command, scene_exit_isolation_command, scene_reveal_command, &
        scene_select_mesh_source_command, inspect_patch_command, show_inspector_command, &
        toggle_inspector_command, toggle_parameter_panels_command, toggle_patch_inspection_command, &
        show_basis_command, snap_off_command, snap_vertex_command, snap_closest_command, transform_command, &
        translate_geometry_command, rotate_geometry_command, scale_geometry_command, affine_geometry_command, &
        measure_geometry_command, color_command, translucent_command, edit_controls_command, &
        edit_interval_command, insert_knot_command, remove_knot_command, edit_knot_values_command, &
        elevate_degree_command, reduce_degree_command, refine_command, split_patch_command, &
        reverse_patch_command, extract_patch_command, extract_section_command, extract_face_command, &
        new_multipatch_command, edit_multipatch_command, validate_multipatch_command, &
        remove_multipatch_command, open_mesh_dialog_command, mesh_edges_command, mesh_pick_command, &
        mesh_only_command, mesh_shaded_command, export_mesh_command, physical_knots_command, &
        world_grid_command, clean_display_command, control_structure_command, interactive_quality_command, &
        shaded_display_command, shaded_edges_display_command, wireframe_display_command, &
        fast_shaded_display_command, viewport_layer_base, parameter_lines_layer, toggle_panel_command_base, &
        viewport_style_names, viewport_shaded, viewport_edges, viewport_wireframe, viewport_fast_shaded
    use studio_icons, only: &
        icon_bounding_vertex, icon_box_selection, icon_close, icon_copy, icon_degree_elevate, &
        icon_degree_reduce, icon_export, icon_extract_face, icon_extract_section, icon_face_grid, &
        icon_fit_panels, icon_fullscreen, icon_generate_mesh, icon_isolate, icon_lasso_selection, &
        icon_maximize, icon_measure, icon_mesh, icon_mesh_pick, icon_mesh_shaded, icon_paste, icon_refresh, &
        icon_restore, icon_scene_reference, icon_screenshot, icon_settings, viewport_style_icons
    implicit none(type,external)
    private
    public :: command_visual

    type :: command_presentation_t
        character(:),allocatable :: caption
        integer :: glyph
    end type command_presentation_t

contains

    pure function command_visual(command,fallback_label,fallback_icon) result(visual)
        integer,intent(in) :: command,fallback_icon
        character(*),intent(in) :: fallback_label
        type(command_presentation_t) :: visual

        select case(command)
        case(find_command)
            visual=command_presentation_t('Find command',27)
        case(interaction_guide_command)
            visual=command_presentation_t('Help',16)
        case(about_command)
            visual=command_presentation_t('About ForCAD Studio',17)
        case(toggle_assistant_command)
            visual=command_presentation_t('Assistant',29)
        case(diagnose_with_ai_command)
            visual=command_presentation_t('Ask Assistant…',29)
        case(fullscreen_command)
            visual=command_presentation_t('Fullscreen',icon_fullscreen)
        case(canvas_command)
            visual=command_presentation_t('Canvas Only',icon_maximize)
        case(restore_layout_command)
            visual=command_presentation_t('Restore panels',icon_refresh)
        case(reset_controls_command)
            visual=command_presentation_t('Reset view and tools',icon_refresh)
        case(reset_presentation_command)
            visual=command_presentation_t('Restore workspace defaults',icon_settings)
        case(home_command)
            visual=command_presentation_t('Home',86)
        case(new_project_command)
            visual=command_presentation_t('New project',81)
        case(open_project_command)
            visual=command_presentation_t('Open project',82)
        case(save_project_command)
            visual=command_presentation_t('Save project',3)
        case(save_as_command)
            visual=command_presentation_t('Save project as',3)
        case(close_project_tab_command)
            visual=command_presentation_t('Close tab',icon_close)
        case(project_folder_command)
            visual=command_presentation_t('Working folder',82)
        case(units_command)
            visual=command_presentation_t('Project units',84)
        case(undo_patch_command)
            visual=command_presentation_t('Undo',4)
        case(redo_patch_command)
            visual=command_presentation_t('Redo',5)
        case(screenshot_command)
            visual=command_presentation_t('Screenshot',icon_screenshot)
        case(import_project_command)
            visual=command_presentation_t('Import FCAD',83)
        case(import_iges_command)
            visual=command_presentation_t('Import IGES',83)
        case(export_vtk_command)
            visual=command_presentation_t('Export VTK',icon_export)
        case(export_iges_command)
            visual=command_presentation_t('Export IGES',icon_export)
        case(export_fortran_command)
            visual=command_presentation_t('Export Fortran',icon_export)
        case(copy_geometry_command)
            visual=command_presentation_t('Copy',icon_copy)
        case(paste_geometry_command)
            visual=command_presentation_t('Paste',icon_paste)
        case(fit_scene_command)
            visual=command_presentation_t('Fit all',19)
        case(fit_around_panels_command)
            visual=command_presentation_t('Fit area',icon_fit_panels)
        case(select_tool_command)
            visual=command_presentation_t('Select',18)
        case(orbit_tool_command)
            visual=command_presentation_t('Orbit',28)
        case(clear_selection_command)
            visual=command_presentation_t('Clear selection',icon_close)
        case(box_selection_command)
            visual=command_presentation_t('Box selection',icon_box_selection)
        case(lasso_selection_command)
            visual=command_presentation_t('Lasso selection',icon_lasso_selection)
        case(project_isolate_command)
            visual=command_presentation_t('Isolate selection',icon_isolate)
        case(project_show_all_command)
            visual=command_presentation_t('Show all',31)
        case(project_hide_selection_command)
            visual=command_presentation_t('Hide selected',32)
        case(scene_exit_isolation_command)
            visual=command_presentation_t('Exit isolation',icon_restore)
        case(scene_reveal_command)
            visual=command_presentation_t('Reveal in Scene',2)
        case(scene_select_mesh_source_command)
            visual=command_presentation_t('Select mesh source',icon_scene_reference)
        case(inspect_patch_command)
            visual=command_presentation_t('Inspect parameters',33)
        case(show_inspector_command)
            visual=command_presentation_t('Inspector',40)
        case(toggle_inspector_command)
            visual=command_presentation_t('Inspector',40)
        case(toggle_parameter_panels_command)
            visual=command_presentation_t('Parameter panels',33)
        case(toggle_patch_inspection_command)
            visual=command_presentation_t('Parameter panels',33)
        case(show_basis_command)
            visual=command_presentation_t('Basis Plot',23)
        case(snap_off_command)
            visual=command_presentation_t('Snap off',0)
        case(snap_vertex_command)
            visual=command_presentation_t('Snap to vertices',icon_bounding_vertex)
        case(snap_closest_command)
            visual=command_presentation_t('Snap to closest point',93)
        case(transform_command)
            visual=command_presentation_t('Exact transform…',38)
        case(translate_geometry_command)
            visual=command_presentation_t('Move',9)
        case(rotate_geometry_command)
            visual=command_presentation_t('Rotate',60)
        case(scale_geometry_command)
            visual=command_presentation_t('Scale',61)
        case(affine_geometry_command)
            visual=command_presentation_t('Affine matrix…',62)
        case(measure_geometry_command)
            visual=command_presentation_t('Measure',icon_measure)
        case(color_command)
            visual=command_presentation_t('Object color…',37)
        case(translucent_command)
            visual=command_presentation_t('Translucent',21)
        case(edit_controls_command)
            visual=command_presentation_t('Points + weights',22)
        case(edit_interval_command)
            visual=command_presentation_t('Set domain',63)
        case(insert_knot_command)
            visual=command_presentation_t('Insert knot',65)
        case(remove_knot_command)
            visual=command_presentation_t('Remove knot',66)
        case(edit_knot_values_command)
            visual=command_presentation_t('Knot values',68)
        case(elevate_degree_command)
            visual=command_presentation_t('Elevate degree',icon_degree_elevate)
        case(reduce_degree_command)
            visual=command_presentation_t('Reduce degree',icon_degree_reduce)
        case(refine_command)
            visual=command_presentation_t('Refine patch',64)
        case(split_patch_command)
            visual=command_presentation_t('Split patch',90)
        case(reverse_patch_command)
            visual=command_presentation_t('Reverse direction',92)
        case(extract_patch_command)
            visual=command_presentation_t('Extract range',91)
        case(extract_section_command)
            visual=command_presentation_t('Extract section',icon_extract_section)
        case(extract_face_command)
            visual=command_presentation_t('Extract face',icon_extract_face)
        case(new_multipatch_command)
            visual=command_presentation_t('New multipatch…',87)
        case(edit_multipatch_command)
            visual=command_presentation_t('Edit multipatch…',88)
        case(validate_multipatch_command)
            visual=command_presentation_t('Check interfaces',89)
        case(remove_multipatch_command)
            visual=command_presentation_t('Dissolve multipatch',36)
        case(open_mesh_dialog_command)
            visual=command_presentation_t('Generate mesh…',icon_generate_mesh)
        case(mesh_edges_command)
            visual=command_presentation_t('Show meshes',31)
        case(mesh_pick_command)
            visual=command_presentation_t('Pick mesh cells',icon_mesh_pick)
        case(mesh_only_command)
            visual=command_presentation_t('Meshes only',icon_mesh)
        case(mesh_shaded_command)
            visual=command_presentation_t('Filled mesh faces',icon_mesh_shaded)
        case(export_mesh_command)
            visual=command_presentation_t('Export mesh (VTK)',icon_export)
        case(physical_knots_command)
            visual=command_presentation_t('Knot lines',63)
        case(world_grid_command)
            visual=command_presentation_t('Ground grid',12)
        case(clean_display_command)
            visual=command_presentation_t('Clean preview',8)
        case(control_structure_command)
            visual=command_presentation_t('Control structure',22)
        case(interactive_quality_command)
            visual=command_presentation_t('Interactive Quality',icon_settings)
        case(viewport_layer_base)
            visual=command_presentation_t('Fill selected faces',8)
        case(viewport_layer_base+1)
            visual=command_presentation_t('Face grid',icon_face_grid)
        case(viewport_layer_base+2)
            visual=command_presentation_t('Control net',22)
        case(viewport_layer_base+3)
            visual=command_presentation_t('Control points',icon_bounding_vertex)
        case(viewport_layer_base+parameter_lines_layer-1)
            visual=command_presentation_t('Parameter lines',24)
        case(toggle_panel_command_base+1)
            visual=command_presentation_t('Scene',2)
        case(toggle_panel_command_base+5)
            visual=command_presentation_t('Parameter Space',63)
        case(toggle_panel_command_base+6)
            visual=command_presentation_t('Basis Plot',23)
        case(shaded_display_command)
            visual=command_presentation_t(trim(viewport_style_names(viewport_shaded)),viewport_style_icons(viewport_shaded))
        case(shaded_edges_display_command)
            visual=command_presentation_t(trim(viewport_style_names(viewport_edges)),viewport_style_icons(viewport_edges))
        case(wireframe_display_command)
            visual=command_presentation_t(trim(viewport_style_names(viewport_wireframe)),viewport_style_icons(viewport_wireframe))
        case(fast_shaded_display_command)
            visual=command_presentation_t(trim(viewport_style_names(viewport_fast_shaded)),viewport_style_icons(viewport_fast_shaded))
        case default
            visual=command_presentation_t(fallback_label,fallback_icon)
        end select
    end function command_visual

end module studio_command_presentation
