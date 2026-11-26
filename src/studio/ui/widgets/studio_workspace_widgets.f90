module studio_workspace_widgets
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: iso_c_binding, only: c_double
    use qt, only: QWidget, QToolButton, QRegularExpressionValidator, QLineEdit, QPushButton, QAction, QLabel, &
        QComboBox, QSlider, QSpinBox, QDoubleSpinBox, QTabWidget, QTreeWidget, QTreeWidgetItem, QIcon
    use studio_commands, only: viewport_style_count
    implicit none(type, external)
    private

    type, public :: studio_ribbon_controls
        type(QAction) :: assistant_toggle_action,alternate_redo_action
        type(QAction) :: assembly_connect_action,assembly_fit_action,assembly_visibility_action
        type(QAction) :: assembly_refine_actions(2)
        type(QAction) :: mesh_generate_action,mesh_edges_action,mesh_export_action,mesh_pick_action
        type(QAction) :: mesh_only_action,mesh_shaded_action,physical_knot_action
        type(QToolButton) :: volume_mesh_modes(3)
        type(QComboBox) :: volume_mesh_axis
        type(QSpinBox) :: volume_mesh_layer
        type(QWidget) :: volume_mesh_controls,volume_mesh_layer_controls
        type(QAction) :: transform_tools(4)
        type(QAction) :: patch_tool_actions(5)
        type(QAction) :: refinement_action
        type(QAction) :: parameters_ribbon
        type(QAction) :: profile_actions(4),measure_action,view_fit_action,view_visibility_action,interactive_quality_action
    end type studio_ribbon_controls

    type, public :: project_tree_row
        type(QTreeWidgetItem) :: item
        character(:), allocatable :: key
        integer :: category=0, rank=0, id=0, branch=0, assembly=0
        integer :: mesh_cells=-1
        logical :: stale=.false.
    end type project_tree_row

    type, public :: project_tree_state
        character(:), allocatable :: key
        logical :: expanded=.false., selected=.false., current=.false.
    end type project_tree_state

    type, public :: project_tree_view
        type(QTreeWidget) :: tree
        type(QWidget) :: viewport
        type(QTreeWidgetItem) :: geometry_branch, multipatch_branch, mesh_branch, entity_branches(3)
        type(project_tree_row), allocatable :: rows(:)
        type(project_tree_state), allocatable :: expanded_before_filter(:)
        integer(int64) :: stamp=0_int64, document_revision=-1_int64
        integer :: row_count=0, filter_mode=0, applied_filter=0, scroll_before_filter=0
        character(:), allocatable :: search_text, applied_search, reveal_key
        logical :: changing_filter=.false., visibility_dirty=.true.
        type(QLineEdit) :: search
        type(QAction) :: filters(9)
        type(QToolButton) :: filter_button, exit_isolation
        type(QWidget) :: empty_state
        type(QLabel) :: empty_title
        type(QPushButton) :: empty_clear
        type(QIcon) :: visible_icon, hidden_icon, mixed_icon, stale_mesh_icon, reference_icon, assembly_icon, folder_icon
        type(QIcon) :: mesh_row_icon, geometry_icons(3)
    end type project_tree_view

    type, public :: mesh_sizing_t
        integer :: rank=0,assembly_id=0
        integer(int64) :: document_revision=-1_int64
        integer, allocatable :: members(:),cells(:,:)
        real(c_double) :: target_size=0._c_double
        logical :: per_patch=.false.,conforming=.false.,synchronize=.true.
        integer(int64), allocatable :: preview_revisions(:)
    end type mesh_sizing_t

    type :: canvas_mesh_controls
        type(QWidget) :: panel, layer_controls, section_position
        type(QComboBox) :: presentation, layer_axis, section_axis
        type(QSlider) :: layer_index, section_offset
        type(QLabel) :: layer_value, section_value
    end type canvas_mesh_controls

    type, public :: workspace_widgets_t
        type(QPushButton) :: maximize_button
        type(QAction) :: fullscreen_action
        logical :: restore_maximized=.false., window_restore_icon=.false.
        type(QPushButton) :: resume_button
        type(QPushButton) :: leave_canvas_button
        type(QLabel) :: project_caption
        type(QToolButton) :: find_command_button
        type(QTabWidget) :: ribbon_tabs
        type(studio_ribbon_controls) :: ribbon
        type(project_tree_view) :: project_browser
        type(QAction) :: multipatch_actions(3)
        logical :: syncing_volume_mesh_controls=.false.
        type(mesh_sizing_t), allocatable :: mesh_drafts(:)
        integer :: mesh_divisions(3,3)=reshape([32,1,1,12,12,1,8,8,8],[3,3])
        type(QAction) :: material_color_action, translucent_action, transform_action
        type(QWidget) :: transform_strip
        type(QLabel) :: transform_caption
        type(QDoubleSpinBox) :: transform_inputs(3)
        type(QPushButton) :: transform_snap_button
        type(QAction) :: export_iges_action, export_vtk_action, canvas_action, pick_actions(6), snap_mode_actions(3)
        type(QToolButton) :: viewport_snap,viewport_pick,viewport_fit,viewport_orientation,viewport_visibility
        type(QAction) :: viewport_snap_action,viewport_click_action,camera_view_actions(7),free_orbit_action
        type(QAction) :: viewport_hide_action,viewport_isolate_action,viewport_exit_isolation_action
        type(QAction) :: viewport_style_actions(viewport_style_count)
        integer :: displayed_viewport_style=0,displayed_viewport_orientation=-1
        type(canvas_mesh_controls) :: canvas_mesh
        type(QAction) :: area_selection_actions(2)
        type(QAction) :: entity_menu_actions(3)
        type(QAction) :: edit_controls_action,undo_patch_action,insert_knot_action,remove_knot_action
        type(QAction) :: elevate_degree_action, reduce_degree_action
        type(QAction) :: edit_knot_values_action, edit_interval_action
        type(QPushButton) :: domain_probe_button,domain_knots_button,domain_groups_button
        type(QPushButton) :: domain_knot_buttons(3),domain_values_button,domain_refine_button
        type(QPushButton) :: interval_apply
        type(QLineEdit) :: interval_lower,interval_upper
        type(QLabel) :: interval_label
        type(QRegularExpressionValidator) :: interval_validator
        integer :: interval_context(3)=0
        integer(int64) :: interval_revision=-1_int64
        type(QWidget) :: domain_knot_tools,domain_groups
        type(QAction) :: scene_selection_visibility_action
        type(QAction) :: project_inspect_action,project_fit_action
        type(QAction) :: project_duplicate_action,project_rename_action,project_remove_action
        type(QPushButton) :: insert_knot_button,remove_knot_button
        logical :: inspector_covered_by_chat=.false.
        logical :: inspector_requested=.true.
        logical :: updating_assistant_layout=.false.
    end type workspace_widgets_t
end module studio_workspace_widgets
