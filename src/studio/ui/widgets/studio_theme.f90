module studio_theme
    use qt, only: QColor,QPalette,QPalette_Window,QPushButton,QVBoxLayout,QWidget
    use studio_icons, only: stylesheet_icon, icon_chevron, icon_chevron_up, icon_chevron_right, icon_check, &
        icon_scene_branch_closed, icon_scene_branch_open
    use studio_commands, only: canvas_command, restore_layout_command, toggle_patch_inspection_command, find_command, &
        home_command, new_project_command, open_project_command, save_project_command, save_as_command, project_folder_command, &
        import_project_command, import_iges_command, export_iges_command, export_vtk_command, export_fortran_command, export_mesh_command, &
        close_project_tab_command, screenshot_command, units_command, clean_display_command, control_structure_command, &
        primitive_command, add_patch_command, extrude_command, revolve_command, sweep_command, loft_command, &
        translate_geometry_command, rotate_geometry_command, scale_geometry_command, affine_geometry_command, transform_command, &
        edit_controls_command, edit_interval_command, insert_knot_command, remove_knot_command, edit_knot_values_command, &
        show_knots_command, domain_details_command, show_basis_command, elevate_degree_command, reduce_degree_command, refine_command, &
        split_patch_command, extract_patch_command, reverse_patch_command, new_multipatch_command, edit_multipatch_command, &
        validate_multipatch_command, remove_multipatch_command, insert_assembly_knot_command,elevate_assembly_degree_command, &
        generate_mesh_command, open_mesh_dialog_command, mesh_pick_command, &
        extract_section_command, extract_face_command, &
        mesh_only_command, mesh_shaded_command, mesh_edges_command, fit_scene_command, fit_selection_command, fit_around_panels_command, &
        select_tool_command, orbit_tool_command, color_command, translucent_command, toggle_patch_visibility_command, interactive_quality_command, &
        physical_knots_command, project_duplicate_command, project_remove_command, project_rename_command, project_isolate_command, &
        project_show_all_command, scene_exit_isolation_command, clear_selection_command, measure_geometry_command, &
        undo_patch_command, redo_patch_command, toggle_assistant_command, toggle_inspector_command, reset_controls_command, reset_presentation_command, &
        viewport_layer_base, parameter_lines_layer, fast_shaded_display_command, wireframe_display_command, &
        shaded_edges_display_command, shaded_display_command
    use studio_forcad_primitives, only: point_shape, line_shape, circle_shape, arc_shape, rectangle_shape, box_shape, cylinder_shape, tube_shape
    implicit none(type, external)
    private

    public :: studio_style_sheet, studio_tooltip, studio_command_tooltip, studio_style_primary, build_opaque_panel

    character(*), parameter, public :: studio_frame_edge='#3b4d5a', studio_divider='#2e3e4b'
    character(*), parameter, public :: studio_heading='#202c37', studio_secondary_text='#b4c8d7'
    character(*), parameter, public :: studio_selected_fill='#244653', studio_active_edge='#76d8d4'
    character(*), parameter, public :: studio_panel_fill='#19232d', studio_ribbon_fill='#19242e', studio_messages_fill='#1a2732'
    character(*), parameter, public :: studio_viewport_controls_fill='#182630'
    character(*), parameter :: popup_fill='#1b2832'

contains

    subroutine build_opaque_panel(panel,frame,fill,frame_name,parent)
        type(QWidget), intent(out) :: panel,frame
        character(*), intent(in) :: fill,frame_name
        class(QWidget), intent(in), optional :: parent
        type(QPalette) :: palette
        type(QVBoxLayout) :: layout

        panel=QWidget(parent)
        palette=panel%palette()
        call palette%setColor(QPalette_Window,QColor(fill))
        call panel%setPalette(palette)
        call palette%delete()
        call panel%setAutoFillBackground(.true.)
        frame=QWidget(panel)
        call frame%setObjectName(frame_name)
        layout=QVBoxLayout(panel)
        call layout%setContentsMargins(0,0,0,0)
        call layout%addWidget(frame)
    end subroutine build_opaque_panel

    subroutine studio_style_primary(button)
        class(QPushButton), intent(in) :: button

        call button%setProperty('primary',.true.)
        call button%repolish()
    end subroutine studio_style_primary

    function studio_style_sheet() result(css)
        character(:), allocatable :: css,down_arrow,up_arrow,closed_branch,open_branch,check_mark
        character(:), allocatable :: closed_scene_branch,open_scene_branch

        down_arrow=stylesheet_icon(icon_chevron,8,5,'image')
        up_arrow=stylesheet_icon(icon_chevron_up,8,5,'image')
        closed_branch=stylesheet_icon(icon_chevron_right,12,12,'image')
        open_branch=stylesheet_icon(icon_chevron,12,12,'image')
        closed_scene_branch=stylesheet_icon(icon_scene_branch_closed,9,9,'image')
        open_scene_branch=stylesheet_icon(icon_scene_branch_open,9,9,'image')
        check_mark=stylesheet_icon(icon_check,14,14,'image')
        css = &
            'QWidget { color:#edf3f7; font-family:"Inter","Noto Sans","Adwaita Sans","Sans Serif"; font-size:14px; }'// &
            'QMainWindow,QWidget#workspace,QWidget#root { background:#10161c; }'// &
            'QWidget#studioHome { background:qradialgradient(cx:.56,cy:.25,radius:.88,fx:.56,fy:.25,'// &
            'stop:0 #223845,stop:.56 #15222c,stop:1 #111820); }'// &
            'QLabel#studioHomeTitle { color:#f3f8fc; font-size:44px; font-weight:700; }'// &
            'QLabel#studioHomeSubtitle { color:#b4c8d6; font-size:16px; }'// &
            'QLabel#studioHomeSection { font-size:17px; font-weight:700; }'// &
            'QWidget#studioHomeCard { background:#192630; border:1px solid #344a59; }'// &
            'QPushButton[welcomeAction="true"] { text-align:left; padding:10px 14px; '// &
            'background:#202f3a; border-color:#3b505e; color:#e7f2f9; font-size:14px; }'// &
            'QPushButton[welcomeAction="true"]:hover { background:#2b414e; border-color:#6a94aa; }'// &
            'QPushButton[primary="true"]:enabled { background:'//studio_selected_fill//'; color:#f0ffff;'// &
            'border-color:'//studio_active_edge//'; }'// &
            'QPushButton[primary="true"]:enabled:hover { background:#2d5561; border-color:#a1e8e4; }'// &
            'QPushButton[welcomeAction="true"]:focus,QPushButton[primary="true"]:enabled:focus { border-color:#9cf5e9; }'// &
            'QDialog { background:#19232d; border:1px solid '//studio_frame_edge//'; }'// &
            'QDialogButtonBox,QMessageBox { dialogbuttonbox-buttons-have-icons:0; }'// &
            'QDialogButtonBox QPushButton,QMessageBox QPushButton { min-width:68px; min-height:24px; }'// &
            'QDialog QPushButton:default:enabled { border-color:#76c4c8; }'// &
            'QDialog QPushButton:focus:enabled { border-color:#b0f1ec; }'// &
            'QDialog QTableWidget:focus,QDialog QListWidget:focus { border-color:'//studio_active_edge//'; }'// &
            'QDialog QHeaderView::section { color:#dbe9f3; font-weight:600; }'// &
            'QDialog QGroupBox { color:#dfebf4; font-weight:600; margin-top:10px; padding-top:5px;'// &
            'border:1px solid '//studio_divider//'; background:#1b2933; }'// &
            'QDialog QGroupBox::title { subcontrol-origin:margin; subcontrol-position:top left; left:8px; padding:0 4px; }'// &
            'QDialog QTabWidget::pane { border:1px solid '//studio_divider//'; background:#18232c; }'// &
            'QDialog QTabBar::tab:hover:!selected { background:#2a3b47; }'// &
            'QDialog QListWidget,QDialog QTableWidget {'// &
            'selection-background-color:'//studio_selected_fill//'; selection-color:#f5fbff; }'// &
            'QDialog QLabel#studioDialogTitle,QDialog QLabel[dialogRole="title"] {'// &
            'color:#f4f9fd; font-size:19px; font-weight:700; padding:5px 9px;'// &
            'background:qlineargradient(x1:0,y1:0,x2:1,y2:0,stop:0 #263945,stop:1 #202c37);'// &
            'border:0; border-left:2px solid '//studio_active_edge//'; border-bottom:1px solid #354957; }'// &
            'QMessageBox QLabel#qt_msgbox_label { color:#f1f8fd; font-weight:600; }'// &
            'QLabel[dialogRole="detail"] { color:'//studio_secondary_text//'; }'// &
            'QLabel[dialogRole="context"] { color:#daedf6; background:#1f333e;'// &
            'border-left:2px solid '//studio_active_edge//'; padding:5px 8px; }'// &
            'QLabel[dialogRole="error"] { color:#ffd5ce; background:#382831; border:1px solid #69505a;'// &
            'border-left:2px solid #e49b9d; padding:6px 8px; }'// &
            'QWidget#brand { background:#171f28; border-bottom:1px solid '//studio_divider//'; }'// &
            'QLabel#brandTitle { font-size:17px; font-weight:650; background:transparent; }'// &
            'QLabel#studioDocumentTitle { color:'//studio_secondary_text//'; padding-left:8px; font-size:13px; font-weight:500; }'// &
            'QWidget#panel { border:1px solid '//studio_frame_edge//'; background:'//studio_panel_fill//'; }'// &
            'QWidget#panelHeader { font-weight:700; font-size:14px;'// &
            'border:0; border-bottom:1px solid '//studio_divider//'; background:'//studio_heading//'; }'// &
            'QLabel#panelTitle { color:#f2f7fb; font-size:14px; font-weight:700; border:0; background:transparent; }'// &
            'QWidget#transformStrip { background:#1b2934; border:1px solid #435b6b; }'// &
            'QLabel#transformCaption { font-weight:650; color:#daedf9; }'// &
            'QLabel#sectionHeader { font-size:14px; font-weight:600; padding:8px 12px;'// &
            'background:#22303b; border:0; }'// &
            'QWidget#inspectorSection { background:#1b2630; border:1px solid '//studio_divider//'; }'// &
            'QWidget#entityContext { background:#1c2a34; border:1px solid '//studio_divider//';'// &
            'border-left:3px solid '//studio_active_edge//'; }'// &
            'QLabel#contextIdentity { font-size:16px; font-weight:700; color:#edf3f7; }'// &
            'QLabel#contextDimensions { font-size:12px; color:#bfd0dc; }'// &
            'QLabel#studioNormalizedProbe { font-size:12px; color:#bfd0dc; }'// &
            'QLabel[probeAxis] { font-size:12px; font-weight:700; }'// &
            'QLabel[probeAxis="U"],QPushButton[parameterAxis="U"] { color:#8ce0db; }'// &
            'QLabel[probeAxis="V"],QPushButton[parameterAxis="V"] { color:#a2c7fb; }'// &
            'QLabel[probeAxis="W"],QPushButton[parameterAxis="W"] { color:#eda9da; }'// &
            'QPushButton#sectionToggle { text-align:left; font-size:13px; font-weight:600;'// &
            'padding:3px 9px; min-height:20px; border:1px solid transparent; background:#22303b; }'// &
            'QPushButton#sectionToggle:hover { background:#2a3944; }'// &
            'QPushButton#sectionToggle:checked { border:1px solid transparent; padding:3px 9px; background:'//studio_heading//'; }'// &
            'QPushButton#sectionToggle:focus { border-color:#8de5e1; }'// &
            'QPushButton#probeReset { background:#22313d; padding:2px 6px; }'// &
            'QPushButton#probeReset:hover { background:#314756; }'// &
            'QPushButton#studioRepairDiagnostic { margin-bottom:7px; }'// &
            'QLabel { background:transparent; border:0; }'// &
            'QLabel#muted { color:'//studio_secondary_text//'; }'// &
            'QWidget#studioSelectedControl { background:#1c2f39; border-top:1px solid #345463; }'// &
            'QLabel#muted:disabled { color:#758391; }'// &
            'QMenu { color:#edf3f7; background:'//popup_fill//'; border:1px solid '//studio_frame_edge//'; padding:4px; }'// &
            'QMenu::item { padding:5px 27px 5px 10px; border-left:2px solid transparent; }'// &
            'QMenu::item:selected:enabled { color:#f5fbff; background:'//studio_selected_fill//';'// &
            'border-left-color:'//studio_active_edge//'; }'// &
            'QMenu#studioViewportContext::item:disabled, QMenu#studioProjectContextMenu::item:disabled {'// &
            'color:'//studio_secondary_text//'; }'// &
            'QMenu[viewportMenu="true"]::item:checked { color:#e5ffff; background:#21444d;'// &
            'border-left-color:'//studio_active_edge//'; }'// &
            'QMenu[viewportMenu="true"]::indicator:checked { '//check_mark//' width:13px; height:13px; }'// &
            'QMenu::separator { background:'//studio_divider//'; height:1px; margin:3px 7px; }'// &
            'QTabWidget#studioRibbonTabs::pane { border:0; background:transparent; }'// &
            'QWidget#studioRibbonPage { background:transparent; border:0; }'// &
            'QScrollArea#studioRibbonStrip { border-top:1px solid '//studio_divider//'; }'// &
            'QScrollArea#studioRibbonStrip QScrollBar:horizontal { height:7px; }'// &
            'QScrollArea#studioRibbonStrip QScrollBar::handle:horizontal:hover { background:#6f8e9e; }'// &
            'QTabWidget#studioRibbonTabs QTabBar { background:transparent; border:0; }'// &
            'QTabWidget#studioRibbonTabs QTabBar::tab { color:#bdcddd; background:transparent;'// &
            'padding:5px 11px; margin-right:3px; min-height:19px; '// &
            'border:1px solid transparent; border-bottom-width:2px; font-size:13px; font-weight:600; }'// &
            'QTabWidget#studioRibbonTabs QTabBar::tab:hover { background:#223846; color:#edf8ff; }'// &
            'QTabWidget#studioRibbonTabs QTabBar::tab:selected { color:#edf8ff; background:#243642;'// &
            'border-bottom-color:'//studio_active_edge//'; }'// &
            'QTabWidget#studioRibbonTabs QTabBar::tab:selected:focus { border-color:'//studio_active_edge//'; }'// &
            'QWidget#studioViewportControlsFrame { background:'//studio_viewport_controls_fill// &
            '; border:1px solid '//studio_frame_edge//'; }'// &
            'QWidget#studioViewerDivider { background:'//studio_divider//'; }'// &
            'QWidget#studioCanvasMeshControls { border-top:1px solid '//studio_divider//'; }'// &
            'QLabel#studioCanvasMeshValue { color:#c5dbeb; font-weight:600; }'// &
            'QSlider[meshScrubber="true"] { background:transparent; border:0; }'// &
            'QSlider[meshScrubber="true"]::groove:horizontal { height:4px; background:#385061; }'// &
            'QSlider[meshScrubber="true"]::add-page:horizontal { background:#385061; }'// &
            'QSlider[meshScrubber="true"]::sub-page:horizontal { background:#629ca6; }'// &
            'QSlider[meshScrubber="true"]::handle:horizontal { width:9px; margin:-5px 0; background:#a6e9ea; }'// &
            'QSlider[meshScrubber="true"]::handle:horizontal:hover,QSlider[meshScrubber="true"]:focus::handle:horizontal {'// &
            'background:#defcff; }'// &
            'QWidget#studioPlotControls { border-bottom:1px solid '//studio_divider//'; }'// &
            'QToolButton { border:1px solid transparent; padding:3px;'// &
            'background:transparent; } QToolButton:hover { background:#314252; }'// &
            'QToolButton:checked { background:'//studio_selected_fill//'; border-color:'//studio_active_edge//'; }'// &
            'QToolButton#studioDomainKnotMenu { padding:0; }'// &
            'QToolButton#studioDomainKnotMenu::menu-indicator { image:none; width:0; }'// &
            'QToolButton#studioProjectActionsButton,QToolButton#studioSceneFilterButton { padding:0; }'// &
            'QToolButton#studioProjectActionsButton::menu-indicator,'// &
            'QToolButton#studioSceneFilterButton::menu-indicator { image:none; width:0; }'// &
            'QToolButton#studioSceneFilterButton { background:#202f3a; border-color:#3a4f60; }'// &
            'QToolButton#studioSceneFilterButton:hover { background:#304352; border-color:#5f7f90; }'// &
            'QToolButton#studioSceneFilterButton:checked { background:#214853; border-color:#62b6ba; }'// &
            'QToolButton#studioSceneFilterButton:focus { border-color:#8de5e1; }'// &
            'QToolButton#studioSceneExitIsolation { color:#f7d8a4; background:#352e25;'// &
            'border-color:#6c5943; margin:0 8px 4px; }'// &
            'QToolButton#studioSceneExitIsolation:hover { background:#493b2b; border-color:#ac895b; }'// &
            'QLabel#studioSceneEmptyTitle { color:#c9d9e4; font-weight:600; padding-bottom:3px; }'// &
            'QToolButton[viewerSettings="true"] { padding:0; background:transparent; }'// &
            'QToolButton[viewerSettings="true"]:hover,QToolButton[viewerSettings="true"]:pressed { background:#314756; }'// &
            'QToolButton[viewerSettings="true"]::menu-indicator { image:none; width:0; }'// &
            'QToolButton[viewportMenu="true"] { padding:0 11px 0 1px; }'// &
            'QToolButton[viewportMenu="true"]:focus { border-color:'//studio_active_edge//'; }'// &
            'QToolButton[viewportMenu="true"]::menu-indicator { '//down_arrow// &
            'subcontrol-origin:border; subcontrol-position:right center; right:3px; width:7px; height:5px; }'// &
            'QToolButton[viewportSplit="true"] { padding:0 12px 0 0; }'// &
            'QToolButton[viewportSplit="true"]::menu-button { width:11px; border-left:1px solid '//studio_divider//'; }'// &
            'QToolButton[viewportSplit="true"]::menu-button:hover { background:#314756; }'// &
            'QToolButton[viewportSplit="true"]::menu-arrow { '//down_arrow//' width:7px; height:5px; }'// &
            'QToolButton#studioViewportSnapButton:checked,QToolButton#studioViewportVisibilityButton:checked {'// &
            'background:'//studio_selected_fill//'; border-color:'//studio_active_edge//'; }'// &
            'QToolButton#studioViewportPickButton { background:#202c37;'// &
            'border-color:'//studio_divider//'; }'// &
            'QToolButton#studioViewportPickButton:hover { background:#314756; border-color:'//studio_frame_edge//'; }'// &
            'QToolButton#studioViewportPickButton:focus { border-color:'//studio_active_edge//'; }'// &
            'QPushButton#studioExitCanvasOnly { padding:0; }'// &
            'QWidget#studioDomainKnotTools { background:#19232d; border-top:1px solid '//studio_divider//'; }'// &
            'QToolButton:focus,QPushButton:focus { border-color:#8de5e1; }'// &
            'QPushButton { border:1px solid #405463; min-height:22px;'// &
            'padding:2px 9px; background:#25323d; }'// &
            'QPushButton:hover { background:#30414d; border-color:#6d8b9c; }'// &
            'QPushButton:pressed { background:#1d4b53; }'// &
            'QPushButton:checked { background:#173e46; border:1px solid #58cfc9; }'// &
            'QWidget#studioRibbonFrame { background:'//studio_ribbon_fill//'; border:1px solid '//studio_frame_edge//'; }'// &
            'QWidget#ribbonGroup { background:transparent; border:0; border-left:1px solid '//studio_divider//'; }'// &
            'QWidget#ribbonGroup[firstRibbonGroup="true"] { border-left:0; }'// &
            'QWidget#studioRibbonQuickDivider { background:'//studio_divider//'; }'// &
            'QToolButton[ribbonQuickAction="true"] { padding:3px 7px; font-size:13px;'// &
            ' color:#dce9f1; background:transparent; border:1px solid transparent; }'// &
            'QToolButton[ribbonQuickAction="true"]:hover { background:#2b3d4a; }'// &
            'QToolButton[ribbonQuickAction="true"]:checked { background:'//studio_selected_fill//'; border-color:'//studio_active_edge//'; }'// &
            'QToolButton[ribbonQuickAction="true"]:focus { border-color:#8de5e1; }'// &
            'QToolButton[ribbonQuickAction="true"]:disabled { color:#8296a7; }'// &
            'QLabel#ribbonCaption { color:#c4d5e0; font-size:12px; font-weight:600;'// &
            'padding:0 4px; background:transparent; }'// &
            'QToolButton[ribbonTile="true"] { padding:2px 4px; '// &
            'font-size:13px; font-weight:600; color:#e9f2fa; }'// &
            'QToolButton[ribbonStacked="true"] { padding:0 4px; }'// &
            'QScrollArea[compactCommandRibbon="true"] QToolButton[ribbonStacked="true"] { padding:0; }'// &
            'QToolButton[ribbonTile="true"]:hover { background:#32495a; border-color:#537388; }'// &
            'QToolButton[ribbonTile="true"]:pressed { background:#203743; border-color:#76a5b7; }'// &
            'QToolButton[ribbonTile="true"]:checked { background:'//studio_selected_fill//'; border-color:'//studio_active_edge//'; }'// &
            'QToolButton[ribbonTile="true"]:focus { border-color:#9ee7eb; }'// &
            'QToolButton[ribbonTile="true"]:disabled { color:#8296a7; background:transparent; border-color:transparent; }'// &
            'QToolButton[ribbonPrimary="true"] { background:#244b53; border-color:#447a80; }'// &
            'QToolButton[ribbonPrimary="true"]:hover { background:#2c5d65; border-color:#6fa5aa; }'// &
            'QToolButton[ribbonTile="true"]::menu-indicator { '//down_arrow// &
            'width:8px; height:5px; subcontrol-origin:padding; subcontrol-position:top right; top:4px; right:2px; }'// &
            'QToolButton[ribbonTile="true"]::menu-button { width:13px; border:0; border-left:1px solid #3a5366; }'// &
            'QToolButton[ribbonTile="true"]::menu-button:hover { background:#32495a; }'// &
            'QToolButton[ribbonTile="true"]::menu-arrow { '//down_arrow//' width:8px; height:5px; }'// &
            'QPushButton#axisTab { min-height:26px; padding:1px 8px; font-size:13px;'// &
            'background:#202c37; border-color:'//studio_frame_edge//'; }'// &
            'QPushButton#axisTab:hover { background:#314756; }'// &
            'QPushButton#axisTab:checked { background:'//studio_selected_fill//'; border-color:'//studio_active_edge//'; }'// &
            'QPushButton#axisTab:focus { border-color:#9ee7eb; }'// &
            'QPushButton#entity:checked { border-bottom:3px solid #58cfc9; padding-bottom:0; font-weight:600; }'// &
            'QPushButton#titleControl,QPushButton#panelClose { padding:0; border:1px solid transparent;'// &
            'background:transparent; }'// &
            'QPushButton#titleControl:hover,QPushButton#panelClose:hover { background:#2c4050; }'// &
            'QPushButton#titleControl[windowClose="true"]:hover { background:#703744; }'// &
            'QPushButton#titleControl:focus,QPushButton#panelClose:focus { border-color:#8de5e1; }'// &
            'QPushButton#viewportTool,QPushButton[viewerTool="true"] {'// &
            'min-height:22px; padding:0; background:transparent; border:1px solid transparent; }'// &
            'QPushButton#viewportTool:hover,QPushButton[viewerTool="true"]:hover { background:#314756; }'// &
            'QPushButton#viewportTool:checked,QPushButton[viewerTool="true"]:checked {'// &
            'background:'//studio_selected_fill//'; border-color:'//studio_active_edge//'; }'// &
            'QPushButton#viewportTool:focus,QPushButton[viewerTool="true"]:focus { border-color:#9ee7eb; }'// &
            'QWidget:disabled,QPushButton:disabled { color:#758391; }'// &
            'QDoubleSpinBox,QSpinBox,QLineEdit,QTextEdit,QPlainTextEdit { background:#121c25; border:1px solid #465c6b;'// &
            'min-height:20px; padding:2px 7px; selection-background-color:#256176; }'// &
            'QDoubleSpinBox:hover,QSpinBox:hover,QLineEdit:hover { border-color:#6a8799; }'// &
            'QDoubleSpinBox:focus,QSpinBox:focus,QLineEdit:focus,QTextEdit:focus,QPlainTextEdit:focus { border-color:#8de5e1; }'// &
            'QLineEdit#studioSceneSearch { background:#151f29; border-color:#465c6e; padding:2px 8px; }'// &
            'QLineEdit#studioSceneSearch:hover { border-color:#607f93; }'// &
            'QLineEdit#studioSceneSearch:focus { border-color:#83d8d8; }'// &
            'QComboBox { background:#1a2933; border:1px solid #465c6b; '// &
            'padding:3px 24px 3px 8px; min-height:22px; }'// &
            'QComboBox:hover,QComboBox:focus { border-color:#8de5e1; }'// &
            'QComboBox::drop-down { width:22px; border:0; }'// &
            'QComboBox::down-arrow { '//down_arrow//' width:8px; height:5px; }'// &
            'QComboBox QAbstractItemView { background:'//popup_fill//'; border:1px solid '//studio_frame_edge//';'// &
            'selection-background-color:'//studio_selected_fill//'; selection-color:#edf3f7; }'// &
            'QDoubleSpinBox::up-button,QSpinBox::up-button { width:14px; border:0; background:#1c2832; }'// &
            'QDoubleSpinBox::down-button,QSpinBox::down-button { width:14px; border:0; background:#1c2832; }'// &
            'QDoubleSpinBox::up-arrow,QSpinBox::up-arrow { '//up_arrow// &
            'width:7px; height:4px; }'// &
            'QDoubleSpinBox::down-arrow,QSpinBox::down-arrow { '//down_arrow// &
            'width:7px; height:4px; }'// &
            'QCheckBox { spacing:7px; background:transparent; min-height:19px; }'// &
            'QCheckBox::indicator { width:13px; height:13px; border:1px solid #bac8d4;'// &
            'background:#17212a; }'// &
            'QCheckBox::indicator:checked {'// &
            'background:#187f8b; border:1px solid #58cfc9;'// &
            check_mark//' }'// &
            'QCheckBox:focus::indicator { border-color:#8de5e1; }'// &
            'QTreeWidget { background:transparent; border:0; padding:4px 2px; outline:0;'// &
            'show-decoration-selected:1; selection-background-color:#173e46; }'// &
            'QTreeWidget::item { height:23px; border:1px solid transparent; }'// &
            'QTreeWidget::item:selected { border:1px solid #58cfc9; }'// &
            'QTreeWidget::branch { background:transparent; }'// &
            'QTreeWidget::branch:has-children:closed { '//closed_branch//' }'// &
            'QTreeWidget::branch:has-children:open { '//open_branch//' }'// &
            'QTreeWidget#studioProjectTree { border:1px solid transparent; padding:0 3px 2px;'// &
            'alternate-background-color:#1b2630; selection-background-color:#234a53;'// &
            'qproperty-indentation:12; qproperty-textElideMode:ElideMiddle; qproperty-verticalScrollMode:ScrollPerPixel; }'// &
            'QTreeWidget#studioProjectTree:focus { border-color:'//studio_active_edge//'; }'// &
            'QTreeWidget#studioProjectTree::branch:has-children:closed { '//closed_scene_branch//' }'// &
            'QTreeWidget#studioProjectTree::branch:has-children:open { '//open_scene_branch//' }'// &
            'QTreeWidget#studioProjectTree::item { height:24px; border:0; padding:0 1px; }'// &
            'QTreeWidget#studioProjectTree::item:hover:!selected { background:#273946; }'// &
            'QTreeWidget#studioProjectTree::item:selected,QTreeWidget#studioProjectTree::branch:selected {'// &
            'background:'//studio_selected_fill//'; color:#f5fbff; border:0; }'// &
            'QTableWidget { background:#17222c; alternate-background-color:#1e2b36; gridline-color:#354653; border:1px solid #3d4f5c;'// &
            'selection-background-color:#265d88; }'// &
            'QTableWidget::item { padding:1px 3px; }'// &
            'QTableWidget#knotTable { color:#ecf4fa; font-size:13px; background:#17222c; alternate-background-color:#20303a;'// &
            'gridline-color:'//studio_divider//'; border-color:'//studio_frame_edge//'; qproperty-verticalScrollMode:ScrollPerItem;'// &
            'selection-background-color:'//studio_selected_fill//'; selection-color:#f5fbff; }'// &
            'QTableWidget#knotTable:focus { border-color:'//studio_active_edge//'; }'// &
            'QTableWidget#knotTable QHeaderView { font-size:12px; font-weight:600; }'// &
            'QTableWidget#knotTable QHeaderView::section { color:#c6d8e5; background:'//studio_heading//'; padding:0 4px; }'// &
            'QTableWidget#knotTable::item { padding:0 4px; }'// &
            'QTableWidget#knotTable::item:hover:!selected { background:#283d4a; }'// &
            'QTableCornerButton::section,QAbstractScrollArea::corner { background:#26323d; border:0; }'// &
            'QHeaderView { background:#18232c; }'// &
            'QHeaderView::section { background:#26323d; border:0; border-right:1px solid #3a4855;'// &
            'border-bottom:1px solid #3a4855; padding:3px; }'// &
            'QTabWidget::pane { border:0; background:transparent; }'// &
            'QTabBar::tab { background:#222e39; border:1px solid #35424e; '// &
            'padding:6px 7px; margin:0; }'// &
            'QTabBar::tab:selected { background:#173e46; border-color:#58cfc9; }'// &
            'QWidget#studioPanelResizeGuide { background:'//studio_active_edge//'; }'// &
            'QSplitterHandle { qproperty-autoFillBackground:true; }'// &
            'QSplitter::handle { background:#182630; }'// &
            'QSplitter::handle:horizontal { border-left:3px solid #182630; border-right:3px solid #182630;'// &
            'background:'//studio_divider//'; }'// &
            'QWidget#studioStatusContent { background:#151e27; border-top:1px solid '//studio_frame_edge//' }'// &
            'QLabel#studioStatusWarning { color:#ffcd89; }'// &
            'QLabel#studioStatusGuidance { color:#d4e4ef; }'// &
            'QWidget#studioStatusContent QToolButton { color:#c1d5e1; background:transparent; border:0; '// &
            'border-left:1px solid #334653; padding:3px 8px; min-height:19px; }'// &
            'QWidget#studioStatusContent QToolButton:hover,QWidget#studioStatusContent QToolButton:focus { '// &
            'background:#273a49; border-left-color:#65b6bd; color:#f1f7fb; }'// &
            'QToolButton#studioStatusAssistant { color:#b9e2e3; }'// &
            'QLabel#studioAssistantElapsed { color:#a9c4d6; padding:0 3px; }'// &
            'QWidget#studioMessagesFrame { background:'//studio_messages_fill//'; border:1px solid #4a6476; }'// &
            'QLabel#studioMessagesTitle { color:#e5eef5; font-weight:700; }'// &
            'QPlainTextEdit#studioMessagesList { background:#14212b; border:1px solid #364b5b; padding:6px; }'// &
            'QScrollBar:vertical { width:8px; background:#17212a; }'// &
            'QScrollBar::handle:vertical { min-height:24px; background:#4c6474; }'// &
            'QScrollBar::add-line:vertical,QScrollBar::sub-line:vertical { height:0; }'// &
            'QScrollBar:horizontal { height:8px; background:#17212a; }'// &
            'QScrollBar::handle:horizontal { min-width:24px; background:#4c6474; }'// &
            'QScrollBar::add-line:horizontal,QScrollBar::sub-line:horizontal { width:0; }'// &
            'QScrollBar::add-page:horizontal,QScrollBar::sub-page:horizontal { background:#17212a; }'// &
            'QScrollArea { background:transparent; border:0; }'// &
            'QScrollArea>QWidget>QWidget { background:transparent; }'
    end function studio_style_sheet

    function studio_tooltip(title, description, shortcut) result(tip)
        character(*), intent(in) :: title
        character(*), intent(in), optional :: description, shortcut
        character(:), allocatable :: tip

        tip=trim(title)
        if(len(tip)==0)return
        if(present(description))then
            if(len_trim(description)>0)tip=tip//' · '//trim(description)
        end if
        if(present(shortcut))then
            if(len_trim(shortcut)>0)tip=tip//' · '//trim(shortcut)
        end if
    end function studio_tooltip

    function studio_command_tooltip(command, title, shortcut) result(tip)
        integer, intent(in) :: command
        character(*), intent(in) :: title
        character(*), intent(in), optional :: shortcut
        character(:), allocatable :: tip, description

        select case(command)
        case(0); description='This tool is planned for a later ForCAD Studio release.'
        case(home_command); description='Return to the start page. The current project remains available to resume.'
        case(new_project_command); description='Start an empty project and choose its working folder.'
        case(open_project_command); description='Open an existing ForCAD Studio project.'
        case(save_project_command); description='Save the current geometry, meshes, and project settings.'
        case(save_as_command); description='Save the current project under a new name or location.'
        case(close_project_tab_command); description='Close the active project tab. ForCAD Studio checks for unsaved changes.'
        case(project_folder_command); description='Choose where the project stores its inputs and generated files.'
        case(import_project_command); description='Add geometry and meshes from another ForCAD Studio project to this one.'
        case(import_iges_command); description='Add curves and surfaces from an IGES file to the current scene.'
        case(export_iges_command); description='Write the selected curve or surface to an IGES file.'
        case(export_vtk_command); description='Write the selected geometry or multipatch to a VTK file.'
        case(export_fortran_command); description='Export the project geometry, multipatches, and meshes as a standalone Fortran program.'
        case(export_mesh_command); description='Export the generated mesh nodes and cell connectivity as VTK.'
        case(screenshot_command); description='Save a workspace image in the project’s screenshot folder.'
        case(units_command); description='Choose the length unit used for dimensions and measurements.'
        case(primitive_command); description='Choose a basic shape and enter its dimensions and placement.'
        case(primitive_command+point_shape); description='Place a point at an exact XYZ position.'
        case(primitive_command+line_shape); description='Create a straight segment between two XYZ endpoints.'
        case(primitive_command+circle_shape); description='Create an XY-plane circle from a center and radius.'
        case(primitive_command+arc_shape); description='Create a circular segment from a center, radius, start angle, and sweep.'
        case(primitive_command+rectangle_shape); description='Create a rectangular surface parallel to XY from a corner, width, and height.'
        case(primitive_command+box_shape); description='Create a box volume from a corner and its three side lengths.'
        case(primitive_command+cylinder_shape); description='Create a cylindrical surface along Z from a base center, radius, and height.'
        case(primitive_command+tube_shape); description='Create a hollow volume along Z from two radii and a height.'
        case(add_patch_command); description='Create a new patch from an existing patch of the same kind.'
        case(extrude_command); description='Extend a curve or surface along a specified direction and distance.'
        case(revolve_command); description='Rotate a curve or surface profile around an axis through a specified angle.'
        case(sweep_command); description='Move a curve or surface profile along a guide curve.'
        case(loft_command); description='Build geometry through an ordered series of curves or surfaces.'
        case(translate_geometry_command); description='Move the selection using viewport handles or exact coordinates.'
        case(rotate_geometry_command); description='Rotate the selection with viewport handles or exact angles.'
        case(scale_geometry_command); description='Resize the selection using viewport handles or axis scale factors.'
        case(affine_geometry_command); description='Enter an exact affine matrix, offset, and optional pivot.'
        case(transform_command); description='Choose Move, Rotate, Scale, or Affine and enter exact values for the selection.'
        case(edit_controls_command); description='Edit control-point coordinates and weights for a patch or multipatch.'
        case(edit_interval_command); description='Choose a direction and set its parameter bounds.'
        case(insert_knot_command); description='Choose a direction and insert a knot while preserving the patch shape.'
        case(remove_knot_command); description='Choose a direction and remove a knot or reduce its multiplicity within a tolerance.'
        case(edit_knot_values_command); description='Choose a direction and edit its interior knot values.'
        case(show_knots_command); description='Open Parameter Space with knot controls on the active axis.'
        case(domain_details_command); description='Show exact knot values and multiplicities for the active axis.'
        case(show_basis_command); description='Inspect basis functions and first derivatives along the active axis.'
        case(elevate_degree_command); description='Choose a direction and increase its degree while preserving the patch shape.'
        case(reduce_degree_command); description='Choose a direction and lower its degree within a specified tolerance.'
        case(refine_command); description='Choose one or all directions to split knot spans and optionally elevate degree.'
        case(split_patch_command); description='Choose a direction and divide the patch at an exact parameter value.'
        case(extract_patch_command); description='Create a patch from a range in a direction you choose.'
        case(reverse_patch_command); description='Choose a direction to reverse its parameter traversal.'
        case(new_multipatch_command); description='Connect existing patches by explicitly pairing their boundaries.'
        case(edit_multipatch_command); description='Change the selected multipatch members or interface definitions.'
        case(validate_multipatch_command); description='Check geometry matching and continuity across multipatch interfaces.'
        case(insert_assembly_knot_command); description='Insert knots in a starting direction and propagate to connected directions.'
        case(elevate_assembly_degree_command); description='Elevate a starting degree and propagate to connected directions.'
        case(remove_multipatch_command); description='Remove the assembly and keep its member patches.'
        case(generate_mesh_command); description='Generate or replace the selected geometry’s structured mesh.'
        case(open_mesh_dialog_command); description='Choose cell counts or a target size, then generate or replace the structured mesh.'
        case(extract_section_command); description='Create an exact native section at a literal U, V, or W value.'
        case(extract_face_command); description='Create an editable surface from a selected volume boundary.'
        case(mesh_pick_command); description='Select individual cells in visible meshes. Press Esc to clear the selection.'
        case(mesh_only_command); description='Hide geometry temporarily and show only visible meshes.'
        case(mesh_shaded_command); description='Switch between filled mesh faces and wireframe.'
        case(mesh_edges_command); description='Show or hide meshes that belong to the selection.'
        case(interactive_quality_command); description='Use a lighter preview while scrolling and restore full quality when zoom stops.'
        case(clean_display_command); description='Show the shaded scene with CAD control guides turned off.'
        case(control_structure_command); description='Show the control net, control points, and parameter lines for the selection.'
        case(fast_shaded_display_command); description='Keep filled faces with lighter rendering, even when the camera stops.'
        case(wireframe_display_command); description='Draw geometry and mesh edges without filled faces.'
        case(shaded_edges_display_command); description='Show filled faces and geometry edges at full quality when the camera stops.'
        case(shaded_display_command); description='Show filled faces at full quality when the camera stops.'
        case(fit_scene_command); description='Frame all visible geometry and meshes in the viewport.'
        case(fit_selection_command); description='Frame the current selection without changing its geometry.'
        case(fit_around_panels_command); description='Frame the scene within the viewport area clear of panels.'
        case(select_tool_command); description='Pick objects or elements in the viewport.'
        case(orbit_tool_command); description='Drag to rotate the camera around the displayed scene.'
        case(color_command); description='Change the display color of the selected patch or multipatch.'
        case(translucent_command); description='Toggle transparency to see geometry behind the selection.'
        case(toggle_patch_visibility_command); description='Show or hide the selected geometry in the scene.'
        case(physical_knots_command); description='Show knot locations on selected curves and knot-span boundaries on selected surfaces and volumes.'
        case(viewport_layer_base); description='Show filled faces on selected surfaces and volumes.'
        case(viewport_layer_base+1); description='Show a sampled grid across visible surface and volume faces.'
        case(viewport_layer_base+2); description='Connect the control points of the selected geometry.'
        case(viewport_layer_base+3); description='Show control points on the selected geometry.'
        case(viewport_layer_base+parameter_lines_layer-1)
            description='Show evenly sampled lines on the selected surface or volume, or sample marks and a direction arrow on the selected curve.'
        case(project_duplicate_command); description='Create a separate copy of the selected patch.'
        case(project_rename_command); description='Change the selected scene object’s name.'
        case(project_remove_command); description='Remove the selected objects from the scene. Undo restores them.'
        case(project_isolate_command); description='Show only the selection until you exit isolation.'
        case(project_show_all_command); description='Show all scene objects and exit isolation.'
        case(scene_exit_isolation_command); description='Restore the visibility settings from before isolation.'
        case(clear_selection_command); description='Deselect every item in the scene and viewport.'
        case(measure_geometry_command); description='Measure the selected geometry, edge, or face in project units.'
        case(undo_patch_command); description='Reverse the last ForCAD Studio edit.'
        case(redo_patch_command); description='Reapply the last undone edit.'
        case(find_command); description='Search ForCAD Studio commands and run them from the keyboard.'
        case(toggle_assistant_command); description='Show or hide the Assistant.'
        case(toggle_inspector_command); description='Show or hide geometry, mesh, and parameter details.'
        case(toggle_patch_inspection_command); description='Show or hide Parameter Space and the Basis Plot together.'
        case(canvas_command); description='Hide the ribbon and panels to focus on the viewport. Repeat to restore them.'
        case(restore_layout_command); description='Bring the workspace panels back to the default arrangement.'
        case(reset_controls_command); description='Reset the camera, display settings, and active tools.'
        case(reset_presentation_command); description='Reset the view and tools, restore the panel layout, and clear saved workspace preferences.'
        case default; description=''
        end select
        tip=studio_tooltip(title,description,shortcut)
    end function studio_command_tooltip
end module studio_theme
