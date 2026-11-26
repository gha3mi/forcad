! Native Qt widget composition for the ForCAD Studio workspace.
module studio_ui
    use studio_inspector, only: build_inspector_widgets, update_inspector_section
    use studio_project_widgets, only: build_project_home, edit_project_units
    use studio_theme, only: studio_style_sheet, studio_tooltip, studio_command_tooltip
    use studio_status_hints, only: set_status_hint,watch_status_hint
    use studio_ribbon, only: build_ribbon, &
        parameter_directions,camera_commands,camera_names
    use, intrinsic :: iso_c_binding, only: c_associated
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_forcad, only: affine_transform_t, rotation_matrix
    use studio_geometry_operations, only: read_geometry_target, commit_document_transform, commit_patch_range, &
        rename_document_patch, commit_document_visibility
    use studio_edit_history, only: document_undo_stamp, undo_document_edit, document_can_redo, redo_document_edit, document_undo_result_t, &
        undo_mesh, undo_material, undo_assembly, undo_patch, undo_project_merge, undo_properties, undo_deletion, prune_document_redo
    use qt, only: qt_string, QAbstractItemView_NoEditTriggers, QAbstractItemView_SelectRows, &
        QAbstractItemView_SingleSelection, QAction, QCloseEvent, QColor, &
        QColorDialog_getColor, QComboBox, QContextMenuEvent, QContextMenuEvent_Keyboard, QDialog, &
        QDialogButtonBox, QDialogButtonBox_AcceptRole, QDialogButtonBox_Cancel, QDialog_Accepted, &
        QDoubleSpinBox, QFileDialog, QFileDialog_AcceptSave, QFileDialog_ExistingFile, &
        QFileDialog_getOpenFileName, QFileDialog_getSaveFileName, QFormLayout, QGridLayout, QHBoxLayout, &
        QIcon, QInputDialog_getText, QKeyEvent, QLabel, QLineEdit, QMainWindow, QMenu, &
        QMessageBox_Cancel, QMessageBox_Discard, QMessageBox_No, QMessageBox_Save, QMessageBox_Yes, &
        QMessageBox_information, QMessageBox_question, QMouseEvent, QObject, &
        QPixmap, QPoint, QPushButton, QRect, QRegularExpressionValidator, QResizeEvent, &
        QSizePolicy_Ignored, QSizePolicy_Preferred, QSpinBox, QSplitter, QStackedWidget, QTabWidget, &
        QTableWidget, QTableWidgetItem, QToolButton, QToolButton_InstantPopup, QTreeWidgetItem, &
        QVBoxLayout, QWidget, Qt_AlignCenter, Qt_AlignRight, Qt_AlignVCenter, Qt_ClickFocus, &
        Qt_Horizontal, Qt_ItemIsEnabled, Qt_ItemIsSelectable, Qt_Key_Delete, Qt_Key_Escape, Qt_Key_F2, Qt_NoModifier, Qt_LeftButton, &
        Qt_PlainText, Qt_RightButton, Qt_ScrollBarAlwaysOff, Qt_StrongFocus, Qt_ToolButtonIconOnly, &
        Qt_ToolButtonTextBesideIcon, Qt_Vertical, QtfEventFilter, QtfWidget
    use studio_selection, only: selection_target_t, selection_replace, &
        selection_patch, selection_assembly, selection_control, selection_mesh, &
        pick_multipatch,pick_control_points,pick_automatic
    use studio_workspace_selection, only: primary_workspace_target, patch_in_selection, mesh_in_selection, pointer_selection_mode
    use studio_state, only: studio_t, studio_action_binding_t, dp, knot_select_tool, knot_insert_tool, &
        knot_remove_tool, mesh_exterior, mesh_layer, studio_mesh_count, &
        initialize_state, adjust_mesh_selection_after_removal, clear_mesh_selection, repaint_all, &
        studio_mesh_cell_label, studio_mesh_valid, area_click, area_box, area_lasso, default_assistant_width
    use studio_icons, only: line_icon,icon_close,icon_minimize,icon_maximize,icon_restore,icon_more, &
        icon_fit_selection,icon_fit_panels,icon_export,icon_mesh,icon_refresh, &
        icon_delete,icon_extract_section,icon_extract_face,icon_generate_mesh,icon_mesh_size,icon_mesh_pick, &
        icon_bounding_face,icon_bounding_edge,icon_bounding_vertex
    use studio_render, only: connect_canvases, domain_knot_fraction, domain_knot_point, reset_domain_view, zoom_domain_view
    use studio_canvas_input, only: cancel_area_drag
    use studio_preferences, only: clear_preferences
    use studio_workspace_preferences, only: save_preferences
    use studio_units, only: unitless, unit_mm, unit_symbol, length_caption, model_length, display_length
    use studio_project_file, only: reset_studio_project,reset_project_document,open_project_document,save_project_document
    use studio_document, only: studio_document_t, document_has_changes,document_display_name,mark_document_saved,document_has_patch, &
        document_patch_in_assembly,record_document_properties,clear_document_history
    use studio_document_deletion, only: can_delete_document_patch, delete_document_patch, delete_document_assembly, &
        delete_document_mesh, deleted_mesh_t, document_deletion_error, delete_document_selection
    use studio_project_folder_dialog, only: choose_project_folder
    use studio_presentation_text, only: entity_name, integer_text, studio_product_name
    use studio_project_tree, only: create_project_tree_widgets, build_project_rows, update_project_mesh_rows, &
        mesh_project_id, sync_scene_visibility_icons, sync_scene_visibility_batch, select_project_patch_row, &
        update_project_summary, project_item_identity, apply_project_filter, reveal_project_row, &
        project_row_visibility, project_target_row, scene_visibility_column
    use studio_command_finder, only: find_studio_command
    use studio_measure_dialog, only: show_geometry_measure
    use studio_commands, only: measure_geometry_command, extract_section_command, extract_face_command, &
        viewport_edges,viewport_wireframe,viewport_style_commands,viewport_style_names, &
        shaded_display_command,shaded_edges_display_command,fast_shaded_display_command,control_structure_command
    use studio_commands, only: parameter_lines_layer, viewport_layer_count, default_viewport_layers
    use studio_commands, only: close_panel_command_base, toggle_panel_command_base, workspace_panels, shaded_layer_command, wireframe_layer_command, control_net_layer_command, &
        close_window_command, minimize_window_command, maximize_window_command, fullscreen_command, reset_presentation_command, &
        inspector_section_command_base, interaction_guide_command, basis_u_command, basis_du_command, basis_dw_command, &
        knots_u_command, knots_w_command, select_curve_command, select_volume_command, clean_display_command, wireframe_display_command
    use studio_commands, only: &
        show_assistant_command, toggle_assistant_command, toggle_inspector_command, project_folder_command, resume_project_command, &
        sample_project_command, diagnose_with_ai_command, edit_interval_axis_base, fit_scene_command, &
        fit_selection_command, inspect_patch_command, toggle_patch_inspection_command, viewport_layer_base, basis_legend_command, &
        basis_probe_command, world_grid_command, fit_around_panels_command, domain_details_command, &
        domain_reset_command, domain_probe_command, basis_functions_command, basis_derivative_command, &
        domain_zoom_out_command, domain_zoom_in_command, apply_interval_command, refine_command, interactive_quality_command, &
        domain_select_command, domain_insert_command, domain_remove_command, translate_geometry_command, &
        affine_geometry_command, transform_exact_base, transform_apply_command, transform_done_command, &
        primitive_command, construction_command, project_duplicate_command, project_rename_command, &
        project_remove_command, project_isolate_command, project_show_all_command, project_expand_command, &
        project_collapse_command, project_show_branch, project_hide_branch, project_hide_selection_command, split_patch_command, &
        scene_exit_isolation_command, scene_reveal_command, scene_select_mesh_source_command, &
        extract_patch_command, reverse_patch_command, snap_off_command, snap_vertex_command, &
        snap_closest_command, generate_mesh_command, open_mesh_dialog_command, mesh_patch_counts_command, mesh_target_size_command, &
        mesh_edges_command, export_mesh_command, &
        mesh_pick_command, mesh_only_command, physical_knots_command, mesh_shaded_command, &
        home_command, new_project_command, open_project_command, &
        save_project_command, save_as_command, units_command, edit_controls_command, &
        import_project_command, copy_geometry_command, paste_geometry_command, close_project_tab_command, &
        undo_patch_command, redo_patch_command, insert_knot_command, remove_knot_command, elevate_degree_command, &
        reduce_degree_command, add_patch_command, edit_knot_values_command, edit_interval_command, &
        new_multipatch_command, edit_multipatch_command, validate_multipatch_command, remove_multipatch_command, &
        insert_assembly_knot_command,elevate_assembly_degree_command, &
        export_iges_command, color_command, translucent_command, transform_command, &
        snap_command, canvas_command, pick_mode_base, pick_controls_command, pick_multipatch_command, &
        clear_selection_command, export_vtk_command, export_fortran_command, &
        back_view_command, bottom_view_command, left_view_command, right_view_command, &
        toggle_patch_visibility_command, reset_controls_command, screenshot_command, import_iges_command, &
        isometric_view_command, zoom_in_command, zoom_out_command, restore_layout_command, &
        front_view_command, top_view_command, select_tool_command, orbit_tool_command, box_selection_command, lasso_selection_command, &
        reset_probe_command, parameter_grid_command, basis_support_command, show_knots_command, &
        show_basis_command, show_inspector_command, toggle_parameter_panels_command, find_command
    use studio_assistant, only: show_studio_assistant, sync_studio_assistant, &
        studio_assistant_is_docked, studio_assistant_is_visible, hide_studio_assistant, &
        ask_studio_diagnostic
    use studio_diagnostics, only: studio_diagnostic_summary, record_studio_diagnostic
    use studio_assistant_actions, only: assistant_action_context

    use studio_single_patch, only: read_single_patch, read_single_patch_scene, read_single_patch_plot, &
        read_single_patch_slice, set_active_patch, active_patch_id, &
        single_patch_count, single_patch_name, single_patch_controls, read_single_patch_controls, commit_single_patch_controls, &
        single_patch_revision, initialize_studio_patch_catalog, &
        commit_single_patch_interval, commit_single_patch_refinement, single_patch_refinement_limits, single_patch_available, &
        add_single_patch, &
        commit_single_patch_reverse,studio_native_patch,publish_native_patches
    use studio_patch_range_dialog, only: edit_studio_patch_range
    use studio_control_editor, only: edit_single_patch_controls, edit_multipatch_controls, control_index_text=>index_text
    use studio_assembly_refinement_dialog, only: edit_assembly_refinement
    use studio_multipatch_controls, only: multipatch_controls_t,read_document_multipatch_controls, &
        commit_document_multipatch_controls,joined_patch_assembly,assembly_control_id,set_multipatch_control_displacement
    use studio_knot_editor, only: edit_single_patch_knots
    use studio_degree_editor, only: edit_single_patch_degree
    use studio_create_patch, only: create_studio_patch
    use studio_primitive_dialog, only: create_studio_primitive
    use studio_forcad, only: primitive_names
    use studio_construction, only: construction_names
    use studio_construction_dialog, only: create_studio_construction
    use studio_parameter_editor, only: edit_single_patch_knot_values, parse_parameter_number, number_text, parameter_number_pattern
    use studio_viewport_camera, only: fit_workspace_selection,set_interactive_quality,set_viewport_orientation
    use studio_knot_interaction, only: connect_knot_interactions
    use studio_live_viewport, only: zoom_live_viewport, frame_workspace_scene, &
        fit_around_panels, fit_framed_scene_in_clear_area, workspace_view_bounds, workspace_view_hit, &
        live_viewport_point, pick_scene_patch, &
        transform_gizmo_available, transform_gizmo_pivot, fit_mesh_selection, cancel_transform_drag
    use studio_scene, only: commit_studio_material
    use studio_scene, only: refresh_studio_scene, studio_scene_is_visible, studio_scene_set_visible, studio_scene_count, &
        frame_studio_scene, get_studio_material, set_studio_material, studio_scene_set_visibility
    use studio_multipatch, only: multipatch_draft,multipatch_join
    use studio_multipatch_collection, only: count_multipatches, save_multipatch, refresh_multipatch, &
        multipatch_undo_stamp
    use studio_multipatch_dialog, only: edit_studio_multipatch
    use studio_transform_dialog, only: edit_studio_transform
    use studio_mesh, only: mesh_axis_cell_limit
    use studio_mesh_dialog, only: edit_mesh_patch_counts
    use studio_workspace_widgets, only: mesh_sizing_t
    use studio_mesh_operations, only: mesh_source_t, mesh_publication_t, read_document_mesh_source, find_document_mesh, &
        generate_document_mesh, mesh_source_missing, mesh_build_failed
    use studio_exchange, only: studio_exchange_suffix,export_studio_mesh
    use studio_document_exchange, only: import_document_iges,export_document_geometry
    implicit none(type, external)
    private
    integer,parameter,public :: studio_min_width=1180, studio_min_height=740
    integer,parameter,public :: studio_default_width=1536, studio_default_height=900
    character(*),parameter,public :: studio_min_size_message='Minimum supported size is 1180 x 740'
    character(12),parameter :: mesh_presentation_names(mesh_exterior:mesh_layer)=[character(12):: &
        'Exterior','Full lattice','Cell layer']

    public :: set_workspace_units,capture_workspace
    public :: build_ui, refresh_inspector, refresh_knot_table, sync_controls, invoke_command, show_ready, set_panel_visibility
    public :: select_entity, open_studio_scene, show_scene_multipatch, rebuild_project_tree
    public :: show_home, on_studio_close, mark_project_saved, open_project_path, new_studio_project
    public :: load_project_workspace, open_project_tab, import_studio_project, copy_project_selection, paste_project_selection
    public :: activate_project_tab, close_project_tab, start_project_tab
    public :: select_patch,select_project_mesh,select_multipatch_group,refresh_rotated_scene, &
        on_assistant_model_created,on_assistant_geometry_changed,available_assistant_tools

    public :: generate_workspace_mesh
    public :: set_workspace_material,select_workspace_mesh_cell,set_workspace_mesh_view,set_workspace_mesh_section,has_visible_meshes
    public :: set_patch_visibility,set_mesh_visibility,rename_project_patch,set_parameter_probe
    public :: import_document_geometry,export_geometry_selection
    public :: apply_patch_range, apply_patch_section, delete_project_patch, delete_project_assembly, delete_project_mesh
    public :: delete_project_selection
    public :: apply_geometry_transform, apply_workspace_selection, apply_workspace_selection_batch, clear_workspace_selection

    character(*),parameter::studio_project_filter=studio_product_name//' project (*.fcad)'


    character(len=6),parameter :: transform_captions(3)=[character(len=6)::'Move','Rotate','Scale']
    integer,parameter :: wide_workspace_with_chat=1610
    logical,parameter :: clean_layers(viewport_layer_count)=[.true.,.false.,.false.,.false.,.false.]
    logical,parameter :: display_layer_presets(viewport_layer_count,2)= &
        reshape([clean_layers,default_viewport_layers],[viewport_layer_count,2])
    interface
        recursive module subroutine invoke_command(s,command)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: command
        end subroutine invoke_command
    end interface

    interface
        module subroutine build_ui(s,empty_document)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            logical, intent(in), optional :: empty_document
        end subroutine build_ui

        module subroutine sync_window_controls(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine sync_window_controls

        module subroutine refresh_workspace_status(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine refresh_workspace_status

        module logical function assistant_shares_inspector_space(s) result(shares_space)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
        end function assistant_shares_inspector_space

        recursive module subroutine sync_assistant_layout(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine sync_assistant_layout

        module subroutine show_ready(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine show_ready

        module subroutine show_scene_context(s,target,at,item)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            type(selection_target_t), intent(in) :: target
            type(QPoint), intent(in) :: at
            type(QTreeWidgetItem), intent(in), optional :: item
        end subroutine show_scene_context

        module subroutine sync_controls(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine sync_controls

        module subroutine restore_workspace_layout(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine restore_workspace_layout

        module subroutine leave_viewport_focus(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine leave_viewport_focus

        module subroutine toggle_canvas(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine toggle_canvas

        module subroutine set_panel_visibility(s,index,visible)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            integer, intent(in) :: index
            logical, intent(in) :: visible
        end subroutine set_panel_visibility

        module subroutine sync_project_column(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine sync_project_column
    end interface

    interface
        module subroutine show_home(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine show_home

        module subroutine enter_workspace(s,start_assistant)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            logical, intent(in) :: start_assistant
        end subroutine enter_workspace

        module logical function document_patch(s,kind,variant) result(present_in_project)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
            integer, intent(in) :: kind, variant
        end function document_patch

        module subroutine new_studio_project(s,folder)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            character(*), intent(in) :: folder
        end subroutine new_studio_project

        module subroutine build_document_tabs(s,row)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
            type(QHBoxLayout),intent(inout) :: row
        end subroutine build_document_tabs

        module subroutine ensure_project_tab(s)
            implicit none(type, external)
            type(studio_t),intent(inout) :: s
        end subroutine ensure_project_tab

        module subroutine start_project_tab(s,assistant_folder)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
            character(*),intent(in),optional :: assistant_folder
        end subroutine start_project_tab

        module subroutine update_project_tabs(s)
            implicit none(type, external)
            type(studio_t),intent(in) :: s
        end subroutine update_project_tabs

        module subroutine refresh_project_tab_actions(s)
            implicit none(type, external)
            type(studio_t),intent(in) :: s
        end subroutine refresh_project_tab_actions

        module subroutine close_active_project_tab(s)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
        end subroutine close_active_project_tab

        module logical function confirm_all_project_tabs(s) result(allowed)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
        end function confirm_all_project_tabs

        module logical function confirm_stored_document(s,index) result(allowed)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
            integer,intent(in) :: index
        end function confirm_stored_document

        module integer function find_project_tab(s,path) result(tab)
            implicit none(type, external)
            type(studio_t),intent(in) :: s
            character(*),intent(in) :: path
        end function find_project_tab

        module subroutine open_project_tab(s,path,error,remembered_path,project_directory,sample_title)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
            character(*),intent(in) :: path
            character(:),allocatable,intent(out) :: error
            character(*),intent(in),optional :: remembered_path,project_directory,sample_title
        end subroutine open_project_tab

        module subroutine activate_project_tab(s,index,keep_assistant)
            implicit none(type,external)
            type(studio_t),intent(inout),target :: s
            integer,intent(in) :: index
            logical,intent(in),optional :: keep_assistant
        end subroutine activate_project_tab

        module subroutine close_project_tab(s,index)
            implicit none(type,external)
            type(studio_t),intent(inout),target :: s
            integer,intent(in) :: index
        end subroutine close_project_tab

        module subroutine import_studio_project(s,path,error,source_name,imported)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
            character(*),intent(in) :: path
            character(:),allocatable,intent(out) :: error
            character(*),intent(in),optional :: source_name
            integer,intent(out),optional :: imported
        end subroutine import_studio_project

        module subroutine run_project_import(s)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
        end subroutine run_project_import

        module subroutine copy_project_selection(s,failure)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
            character(:),allocatable,intent(out),optional :: failure
        end subroutine copy_project_selection

        module subroutine paste_project_selection(s,failure,same_folder)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
            character(:),allocatable,intent(out),optional :: failure
            logical,intent(in),optional :: same_folder
        end subroutine paste_project_selection

        module logical function project_clipboard_in_scope(s) result(current)
            implicit none(type,external)
            type(studio_t),intent(in) :: s
        end function project_clipboard_in_scope

        module subroutine prepare_document_view(s)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
        end subroutine prepare_document_view

        module subroutine pick_document_selection(s)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
        end subroutine pick_document_selection

        module subroutine refresh_project_merge_view(s)
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
        end subroutine refresh_project_merge_view

        module subroutine reset_document(s,templates)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            logical, intent(in) :: templates
        end subroutine reset_document

        module subroutine sync_project_title(s)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
        end subroutine sync_project_title

        module subroutine mark_project_saved(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine mark_project_saved

        module logical function confirm_project_change(s) result(allowed)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end function confirm_project_change

        module logical function write_project(s,save_as) result(saved)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            logical, intent(in) :: save_as
        end function write_project

        module subroutine choose_project(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine choose_project

        module subroutine open_project_path(s,path)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            character(*), intent(in) :: path
        end subroutine open_project_path

        module subroutine load_project_workspace(s,path,error,remembered_path,project_directory)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            character(*), intent(in) :: path
            character(:), allocatable, intent(out) :: error
            character(*), intent(in), optional :: remembered_path,project_directory
        end subroutine load_project_workspace

        module subroutine on_studio_close(event,data)
            implicit none(type, external)
            type(QCloseEvent), intent(in) :: event
            class(*), intent(inout), optional :: data
        end subroutine on_studio_close

        module subroutine choose_units(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine choose_units

        module subroutine set_workspace_units(s,model_unit,display_unit,error)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            integer, intent(in) :: model_unit,display_unit
            character(:), allocatable, intent(out), optional :: error
        end subroutine set_workspace_units

        module logical function capture_workspace(s,path) result(saved)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            character(*), intent(in) :: path
        end function capture_workspace

        module integer function active_length_unit(s) result(id)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
        end function active_length_unit

        module function project_unit_status(s) result(label)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
            character(:), allocatable :: label
        end function project_unit_status

        module subroutine rebuild_project_tree(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine rebuild_project_tree

        module subroutine sync_project_meshes(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine sync_project_meshes

        module subroutine select_project_mesh(s,id)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: id
        end subroutine select_project_mesh

        module subroutine refresh_multipatch_card(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine refresh_multipatch_card

        module subroutine show_multipatch_inspector_summary(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine show_multipatch_inspector_summary

        module subroutine select_multipatch_group(s,id,reveal)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: id
            logical, intent(in), optional :: reveal
        end subroutine select_multipatch_group

        module subroutine select_project_patch(s,kind)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            integer, intent(in) :: kind
        end subroutine select_project_patch

        module subroutine sync_patch_visibility_button(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine sync_patch_visibility_button

        module subroutine exit_scene_isolation(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine exit_scene_isolation

        module subroutine toggle_scene_selection_visibility(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine toggle_scene_selection_visibility

        module subroutine hide_scene_selection(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine hide_scene_selection

        module subroutine reveal_scene_selection(s,target)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            type(selection_target_t), intent(in), optional :: target
        end subroutine reveal_scene_selection

        module subroutine select_scene_mesh_source(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine select_scene_mesh_source

        module logical function multipatch_members_shown(s) result(visible)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
        end function multipatch_members_shown

        module subroutine forget_workspace_selection(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine forget_workspace_selection

        module subroutine apply_workspace_selection(s,target,mode)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            type(selection_target_t), intent(in) :: target
            integer, intent(in) :: mode
        end subroutine apply_workspace_selection

        module subroutine apply_workspace_selection_batch(s,targets,mode)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            type(selection_target_t), intent(in) :: targets(:)
            integer, intent(in) :: mode
        end subroutine apply_workspace_selection_batch

        module subroutine clear_workspace_selection(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine clear_workspace_selection

        module function tree_selection_target(s,item) result(target)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
            type(QTreeWidgetItem), intent(in) :: item
            type(selection_target_t) :: target
        end function tree_selection_target

        module subroutine sync_selection_rows(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine sync_selection_rows

        module subroutine on_tree_click(item,column,data)
            implicit none(type, external)
            type(QTreeWidgetItem), intent(in) :: item
            integer, intent(in) :: column
            class(*), intent(inout), optional :: data
        end subroutine on_tree_click

        module subroutine on_tree_selection(item,data)
            implicit none(type, external)
            type(QTreeWidgetItem), intent(in) :: item
            class(*), intent(inout), optional :: data
        end subroutine on_tree_selection

        module subroutine on_project_key(event,data)
            implicit none(type, external)
            type(QKeyEvent), intent(in) :: event
            class(*), intent(inout), optional :: data
        end subroutine on_project_key

        module subroutine on_project_mouse_press(event,data)
            implicit none(type, external)
            type(QMouseEvent), intent(in) :: event
            class(*), intent(inout), optional :: data
        end subroutine on_project_mouse_press

        module subroutine on_project_context(event,data)
            implicit none(type, external)
            type(QContextMenuEvent), intent(in) :: event
            class(*), intent(inout), optional :: data
        end subroutine on_project_context

        module subroutine run_project_delete(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine run_project_delete

        module subroutine delete_project_selection(s,targets,removed,error)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            type(selection_target_t), intent(in) :: targets(:)
            integer, intent(out) :: removed
            character(:), allocatable, intent(out) :: error
        end subroutine delete_project_selection

        module subroutine delete_project_patch(s,kind,variant,error)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: kind, variant
            character(:), allocatable, intent(out) :: error
        end subroutine delete_project_patch

        module subroutine delete_project_assembly(s,id,error)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: id
            character(:), allocatable, intent(out) :: error
        end subroutine delete_project_assembly

        module subroutine delete_project_mesh(s,id,error)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: id
            character(:), allocatable, intent(out) :: error
        end subroutine delete_project_mesh

        module subroutine undo_project_edit(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine undo_project_edit

        module subroutine redo_project_edit(s)
            implicit none(type, external)
            type(studio_t),intent(inout),target::s
        end subroutine redo_project_edit

        module subroutine run_project_edit(s,command)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: command
        end subroutine run_project_edit

        module subroutine rename_project_patch(s,kind,variant,name,error,changed)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            integer, intent(in) :: kind, variant
            character(*), intent(in) :: name
            character(:), allocatable, intent(out) :: error
            logical, optional, intent(out) :: changed
        end subroutine rename_project_patch

        module subroutine set_patch_visibility(s,kind,variants,visible,changed,error)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            integer, intent(in) :: kind, variants(:)
            logical, intent(in) :: visible
            logical, intent(out) :: changed
            character(:), allocatable, intent(out), optional :: error
        end subroutine set_patch_visibility

        module subroutine update_project_visibility(s,show_all)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            logical, intent(in) :: show_all
        end subroutine update_project_visibility

        module subroutine on_tree_open(item,column,data)
            implicit none(type, external)
            type(QTreeWidgetItem), intent(in) :: item
            integer, intent(in) :: column
            class(*), intent(inout), optional :: data
        end subroutine on_tree_open

        module subroutine open_studio_scene(s,keep_tool)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            logical, intent(in), optional :: keep_tool
        end subroutine open_studio_scene

        module subroutine show_scene_multipatch(s,id)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: id
        end subroutine show_scene_multipatch

        module subroutine run_geometry_import(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine run_geometry_import

        module subroutine import_document_geometry(s,path,kind,variant,count,skipped,error,source_unit,needs_source_unit)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            character(*), intent(in) :: path
            integer, intent(out) :: kind, variant, count, skipped
            character(:), allocatable, intent(out) :: error
            integer, intent(in), optional :: source_unit
            logical, intent(out), optional :: needs_source_unit
        end subroutine import_document_geometry

        module function choose_export_path(s,title,base,filter,suffix,alternate) result(path)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
            character(*), intent(in) :: title, base, filter, suffix
            character(*), intent(in), optional :: alternate
            character(:), allocatable :: path
        end function choose_export_path

        module subroutine run_geometry_export(s,vtk)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            logical, intent(in) :: vtk
        end subroutine run_geometry_export

        module subroutine run_fortran_export(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine run_fortran_export

        module subroutine export_geometry_selection(s,path,error)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            character(*), intent(in) :: path
            character(:), allocatable, intent(out) :: error
        end subroutine export_geometry_selection
    end interface

    interface
        module function compact_scalar(value) result(text)
            implicit none(type, external)
            real(dp), intent(in) :: value
            character(len=:), allocatable :: text
        end function compact_scalar

        module subroutine select_entity(s,kind,variant)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: kind
            integer, intent(in), optional :: variant
        end subroutine select_entity

        module subroutine select_patch(s,kind,variant,reveal)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: kind
            integer, intent(in), optional :: variant
            logical, intent(in) :: reveal
        end subroutine select_patch

        module subroutine refresh_inspector(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine refresh_inspector

        module subroutine on_control_selection(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine on_control_selection

        module subroutine refresh_edit_actions(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine refresh_edit_actions

        module logical function can_reduce_degree(s) result(available)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
        end function can_reduce_degree

        module logical function can_remove_knot(s,any_direction) result(available)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
            logical, intent(in), optional :: any_direction
        end function can_remove_knot

        module subroutine refresh_studio_diagnostics(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine refresh_studio_diagnostics

        module subroutine refresh_inline_interval(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine refresh_inline_interval

        module subroutine refresh_domain_readout(s)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
        end subroutine refresh_domain_readout

        module subroutine refresh_knot_table(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine refresh_knot_table

        module subroutine on_knot_double_click(row,column,data)
            implicit none(type, external)
            integer, intent(in) :: row, column
            class(*), intent(inout), optional :: data
        end subroutine on_knot_double_click

        module subroutine on_knot_table_selection(row,column,previous_row,previous_column,data)
            implicit none(type, external)
            integer, intent(in) :: row, column, previous_row, previous_column
            class(*), intent(inout), optional :: data
        end subroutine on_knot_table_selection

        module subroutine on_domain(checked,data)
            implicit none(type, external)
            logical, intent(in) :: checked
            class(*), intent(inout), optional :: data
        end subroutine on_domain

        module subroutine on_probe(value,data)
            implicit none(type, external)
            real(dp), intent(in) :: value
            class(*), intent(inout), optional :: data
        end subroutine on_probe

        module subroutine set_parameter_probe(s,normalized)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            real(dp), intent(in) :: normalized(3)
        end subroutine set_parameter_probe
    end interface

    interface
        module subroutine refresh_shared_controls(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine refresh_shared_controls

        module subroutine edit_workspace_controls(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine edit_workspace_controls

        module subroutine on_control_open(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine on_control_open

        module subroutine on_control_move(s,point,position)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: point
            real(dp), intent(in) :: position(3)
        end subroutine on_control_move

        module subroutine on_interval_apply(data)
            implicit none(type, external)
            class(*), intent(inout), optional :: data
        end subroutine on_interval_apply

        module subroutine run_knot_editor(s,removing)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            logical, intent(in) :: removing
        end subroutine run_knot_editor

        module subroutine run_patch_tool(s,command)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: command
        end subroutine run_patch_tool

        module subroutine apply_patch_range(s,axis,lower,upper,splitting,created,error)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: axis
            real(dp), intent(in) :: lower, upper
            logical, intent(in) :: splitting
            integer, intent(out) :: created
            character(:), allocatable, intent(out) :: error
        end subroutine apply_patch_range

        module subroutine run_patch_section(s,face)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            logical, intent(in) :: face
        end subroutine run_patch_section

        module subroutine apply_patch_section(s,axis,parameter,created,error)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: axis
            real(dp), intent(in) :: parameter
            integer, intent(out) :: created
            character(:), allocatable, intent(out) :: error
        end subroutine apply_patch_section

        module subroutine run_degree_editor(s,reducing)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            logical, intent(in) :: reducing
        end subroutine run_degree_editor

        module subroutine run_patch_interval(s,axis_hint)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: axis_hint
        end subroutine run_patch_interval

        module subroutine apply_inline_interval(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine apply_inline_interval

        module subroutine run_refinement(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine run_refinement

        module subroutine run_knot_values_editor(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine run_knot_values_editor

        module subroutine on_knot_interaction(s,value,operation,error,changed)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            real(dp), intent(in) :: value
            character(*), intent(in) :: operation, error
            logical, intent(in) :: changed
        end subroutine on_knot_interaction

        module logical function can_restructure_patch(s) result(allowed)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
        end function can_restructure_patch
    end interface

    interface
        module subroutine update_assembly_ribbon(s,has_group)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
            logical, intent(in) :: has_group
        end subroutine update_assembly_ribbon

        module subroutine update_mesh_ribbon(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine update_mesh_ribbon

        module subroutine volume_mesh_style_changed(index,data)
            implicit none(type, external)
            integer, intent(in) :: index
            class(*), intent(inout), optional :: data
        end subroutine volume_mesh_style_changed

        module subroutine volume_mesh_axis_changed(index,data)
            implicit none(type, external)
            integer, intent(in) :: index
            class(*), intent(inout), optional :: data
        end subroutine volume_mesh_axis_changed

        module subroutine volume_mesh_layer_changed(index,data)
            implicit none(type, external)
            integer, intent(in) :: index
            class(*), intent(inout), optional :: data
        end subroutine volume_mesh_layer_changed

        module subroutine mesh_section_axis_changed(index,data)
            implicit none(type, external)
            integer, intent(in) :: index
            class(*), intent(inout), optional :: data
        end subroutine mesh_section_axis_changed

        module subroutine mesh_section_offset_changed(index,data)
            implicit none(type, external)
            integer, intent(in) :: index
            class(*), intent(inout), optional :: data
        end subroutine mesh_section_offset_changed

        module subroutine set_workspace_mesh_view(s,mode,axis,layer,error)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            integer, intent(in) :: mode
            integer, intent(in), optional :: axis,layer
            character(:), allocatable, intent(out) :: error
        end subroutine set_workspace_mesh_view

        module subroutine set_workspace_mesh_section(s,axis,percent,error)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            integer, intent(in) :: axis,percent
            character(:), allocatable, intent(out) :: error
        end subroutine set_workspace_mesh_section

        module logical function has_visible_meshes(s,minimum_rank) result(shown)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
            integer, intent(in), optional :: minimum_rank
        end function has_visible_meshes

        module subroutine sync_selection_filter(s)
            implicit none(type, external)
            type(studio_t),intent(inout)::s
        end subroutine sync_selection_filter

        module subroutine sync_snap_controls(s)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
        end subroutine sync_snap_controls

        module subroutine sync_transform_controls(s)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
        end subroutine sync_transform_controls

        module subroutine reset_transform_values(s)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
        end subroutine reset_transform_values

        module subroutine activate_transform(s,mode)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: mode
        end subroutine activate_transform

        module subroutine on_transform_exit(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine on_transform_exit

        module subroutine apply_transform_values(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine apply_transform_values

        module subroutine on_transform_drag(s,linear,offset,pivot)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            real(dp), intent(in) :: linear(3,3), offset(3), pivot(3)
        end subroutine on_transform_drag

        module subroutine on_assistant_geometry_changed(s,kind,variant)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: kind, variant
        end subroutine on_assistant_geometry_changed

        module subroutine on_assistant_model_created(s,group,kind,variant,mesh_id)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: group, kind, variant, mesh_id
        end subroutine on_assistant_model_created

        module function available_assistant_tools(s) result(context)
            implicit none(type, external)
            type(studio_t), intent(in) :: s
            type(assistant_action_context) :: context
        end function available_assistant_tools

        module subroutine refresh_rotated_scene(s)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine refresh_rotated_scene

        module subroutine run_multipatch_editor(s,id)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: id
        end subroutine run_multipatch_editor

        module subroutine run_geometry_transform(s,operation)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in), optional :: operation
        end subroutine run_geometry_transform

        module subroutine apply_geometry_transform(s,transform,fit_after,changed_out,error_out,reveal)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            type(affine_transform_t), intent(in) :: transform
            logical, intent(in) :: fit_after
            logical, intent(out), optional :: changed_out
            character(:), allocatable, intent(out), optional :: error_out
            logical, intent(in), optional :: reveal
        end subroutine apply_geometry_transform

        module subroutine set_workspace_material(s,changed,error,rgb,opacity)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            logical, intent(out) :: changed
            character(:), allocatable, intent(out) :: error
            integer, intent(in), optional :: rgb(3), opacity
        end subroutine set_workspace_material

        module subroutine select_workspace_mesh_cell(s,id,number,part,cell,error)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: id, number
            integer, intent(out) :: part, cell
            character(:), allocatable, intent(out) :: error
        end subroutine select_workspace_mesh_cell

        module subroutine change_material(s,change_color)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            logical, intent(in) :: change_color
        end subroutine change_material

        module subroutine generate_workspace_mesh(s,requested_cells,mesh_id,failure,conforming,patch_cells,source_override,target_size,synchronize,prepared_sizing)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            integer, intent(in), optional :: requested_cells(:)
            integer, intent(out), optional :: mesh_id
            character(:), allocatable, intent(out), optional :: failure
            logical, intent(in), optional :: conforming
            integer, intent(in), optional :: patch_cells(:,:)
            type(mesh_source_t), intent(in), optional :: source_override
            real(dp), intent(in), optional :: target_size
            logical, intent(in), optional :: synchronize
            type(mesh_sizing_t), intent(in), optional :: prepared_sizing
        end subroutine generate_workspace_mesh

        module subroutine run_patch_mesh_dialog(s,by_target_size,generate)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            logical, intent(in) :: by_target_size
            logical, intent(in), optional :: generate
        end subroutine run_patch_mesh_dialog

        module subroutine clamp_mesh_layer(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine clamp_mesh_layer

        module subroutine save_selected_mesh(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine save_selected_mesh

        module subroutine toggle_selected_mesh_visibility(s)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
        end subroutine toggle_selected_mesh_visibility

        module subroutine set_mesh_visibility(s,id,visible,changed,error)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            integer, intent(in) :: id
            logical, intent(in) :: visible
            logical, intent(out) :: changed
            character(:), allocatable, intent(out), optional :: error
        end subroutine set_mesh_visibility

        module subroutine on_mesh_cell_selected(s,id,part,cell,face)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: id, part, cell, face
        end subroutine on_mesh_cell_selected

        module subroutine show_created_patch(s,kind)
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: kind
        end subroutine show_created_patch

        module subroutine change_snap_mode(s,command)
            implicit none(type, external)
            type(studio_t), intent(inout) :: s
            integer, intent(in) :: command
        end subroutine change_snap_mode
    end interface
contains

    subroutine selected_geometry(s,kind,patch_ids,name)
        type(studio_t), intent(in) :: s
        integer, intent(out) :: kind
        integer, allocatable, intent(out) :: patch_ids(:)
        character(:), allocatable, intent(out) :: name

        if(s%selected_mesh>0)then
            kind=s%entity_kind
            allocate(patch_ids(0))
            name=''
            return
        end if
        call read_geometry_target(s%document,s%entity_kind,selected_patch_id(s),s%selected_multipatch,kind,patch_ids,name)
    end subroutine selected_geometry

    function mesh_source_selection(s) result(source)
        type(studio_t), intent(in) :: s
        type(mesh_source_t) :: source

        source=read_document_mesh_source(s%document,s%entity_kind,selected_patch_id(s),s%selected_multipatch,s%selected_mesh)
    end function mesh_source_selection

    integer function mesh_context_record(s) result(id)
        type(studio_t), intent(in) :: s

        id=s%selected_mesh
        if(id>0)return
        id=find_document_mesh(s%document,s%entity_kind,selected_patch_id(s),s%selected_multipatch)
    end function mesh_context_record

    integer function selected_patch_id(s) result(id)
        type(studio_t), intent(in) :: s

        id=0
        if(s%live_scene%valid)id=active_patch_id(s%entity_kind,catalog=s%document%patches)
    end function selected_patch_id

    pure integer function pick_action_index(command) result(index)
        integer, intent(in) :: command

        index=command-pick_mode_base+1
        if(command==pick_controls_command)index=pick_control_points+1
        if(command==pick_multipatch_command)index=pick_multipatch+1
    end function pick_action_index

end module studio_ui
