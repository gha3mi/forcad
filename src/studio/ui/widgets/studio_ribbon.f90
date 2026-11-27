module studio_ribbon
    use studio_status_hints, only: set_status_hint, set_tab_status_hint, watch_status_hint, register_status_menu
    use qt, only: QObject,QAction,QWidget,QToolButton,QMenu,QLabel,QHBoxLayout,QVBoxLayout,QGridLayout,QTabWidget, &
        QSizePolicy_Fixed,QSizePolicy_Ignored,QScrollArea,Qt_ScrollBarAsNeeded,Qt_ScrollBarAlwaysOff, &
        QToolButton_InstantPopup,Qt_AlignCenter,Qt_AlignTop,Qt_AlignRight,Qt_AlignVCenter, &
        Qt_ToolButtonIconOnly,Qt_ToolButtonTextBesideIcon, &
        Qt_ToolButtonTextUnderIcon
    use studio_state, only: studio_t
    use studio_theme, only: studio_tooltip,studio_ribbon_fill,build_opaque_panel
    use studio_workspace_widgets, only: studio_ribbon_controls
    use studio_icons, only: line_icon,icon_close,icon_fit_selection,icon_fit_panels,icon_export,icon_mesh, &
        icon_measure,icon_delete,icon_box_selection,icon_lasso_selection,icon_mesh_pick,icon_mesh_shaded, &
        icon_bounding_vertex,viewport_style_icons, &
        icon_degree_elevate,icon_degree_reduce,icon_extract_section,icon_extract_face,icon_generate_mesh, &
        icon_camera_front,icon_camera_top,icon_refresh,icon_maximize
    use studio_viewport_camera, only: camera_names,camera_icons
    use studio_forcad, only: primitive_names,point_shape,line_shape,circle_shape,arc_shape,rectangle_shape,cylinder_shape,box_shape,tube_shape
    use studio_commands, only: measure_geometry_command, toggle_panel_command_base, new_project_command,undo_patch_command,redo_patch_command,edit_controls_command,translate_geometry_command, &
        transform_command,toggle_assistant_command,home_command,open_project_command,save_project_command,save_as_command, &
        project_folder_command,export_iges_command,export_vtk_command,export_fortran_command,units_command, &
        import_project_command,copy_geometry_command,paste_geometry_command,close_project_tab_command, &
        toggle_patch_inspection_command, &
        insert_knot_command,remove_knot_command,edit_knot_values_command,edit_interval_command,elevate_degree_command, &
        reduce_degree_command,open_mesh_dialog_command,mesh_edges_command, &
        mesh_pick_command,mesh_shaded_command,mesh_only_command,export_mesh_command,fit_scene_command,fit_selection_command, &
        fit_around_panels_command,color_command,translucent_command,toggle_patch_visibility_command,restore_layout_command,clean_display_command, &
        viewport_layer_base,viewport_layer_count,parameter_lines_layer,world_grid_command,interactive_quality_command, &
        viewport_style_count,viewport_style_commands,viewport_style_names,viewport_shaded,viewport_edges,viewport_wireframe, &
        viewport_fast_shaded, &
        physical_knots_command,canvas_command,refine_command,split_patch_command,extract_patch_command,reverse_patch_command, &
        extract_section_command,extract_face_command, &
        new_multipatch_command,edit_multipatch_command,validate_multipatch_command,remove_multipatch_command, &
        insert_assembly_knot_command,elevate_assembly_degree_command, &
        extrude_command,revolve_command,sweep_command,loft_command,primitive_command,add_patch_command, &
        clear_selection_command,box_selection_command,lasso_selection_command,pick_multipatch_command,pick_controls_command, &
        project_duplicate_command,project_rename_command,project_remove_command, &
        back_view_command,left_view_command,right_view_command,bottom_view_command,about_command,interaction_guide_command,model_ribbon_tab, &
        reset_controls_command,reset_presentation_command,close_window_command,select_curve_command,select_surface_command,select_volume_command
    use studio_commands, only: screenshot_command, import_iges_command, isometric_view_command, front_view_command, &
        top_view_command, zoom_in_command, zoom_out_command, select_tool_command, orbit_tool_command, shaded_layer_command, &
        wireframe_layer_command, control_net_layer_command, parameter_grid_command, basis_support_command
    implicit none(type, external)
    private

    public :: studio_ribbon_controls, build_ribbon
    public :: parameter_directions, camera_commands, camera_names

    integer, parameter :: ribbon_tile_height=52, ribbon_stacked_height=24
    character(len=3),parameter :: parameter_directions='UVW'
    integer,parameter :: camera_commands(7)=[isometric_view_command,front_view_command,back_view_command,left_view_command, &
        right_view_command,top_view_command,bottom_view_command]

    abstract interface
        subroutine ribbon_action_factory(parent,text,shortcut,s,command,icon_kind,action)
            import :: QObject, studio_t, QAction
            implicit none(type, external)
            class(QObject), intent(in) :: parent
            character(*), intent(in) :: text,shortcut
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: command,icon_kind
            type(QAction), intent(out) :: action
        end subroutine ribbon_action_factory

        subroutine ribbon_menu_factory(menu,text,shortcut,s,command,icon_kind,menu_item)
            import :: QMenu, studio_t, QAction
            implicit none(type, external)
            type(QMenu), intent(inout) :: menu
            character(*), intent(in) :: text,shortcut
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: command,icon_kind
            type(QAction), intent(out), optional :: menu_item
        end subroutine ribbon_menu_factory

    end interface

    integer, parameter :: compact_group_spacing=2

contains

    subroutine build_ribbon(s,make_command,menu_action,fullscreen_action)
        type(studio_t),intent(inout),target::s
        type(studio_ribbon_controls),pointer :: controls
        procedure(ribbon_action_factory)::make_command
        procedure(ribbon_menu_factory)::menu_action
        type(QAction),intent(in)::fullscreen_action
        type(QTabWidget)::tabs
        type(QAction)::action,save_action,edit_controls_action
        type(QGridLayout)::ribbon_layout
        type(QWidget)::frame,quick
        type(QHBoxLayout)::quick_row
        integer::i
        integer,parameter::codes(17)=[new_project_command,import_iges_command,screenshot_command,undo_patch_command,redo_patch_command, &
            select_tool_command,orbit_tool_command,zoom_out_command,zoom_in_command,isometric_view_command,shaded_layer_command, &
            wireframe_layer_command,front_view_command,top_view_command,parameter_grid_command,control_net_layer_command,basis_support_command]
        integer,parameter::icons(17)=[81,83,26,4,5,18,28,6,7,25,8,21,icon_camera_front,icon_camera_top,10,22,12]
        character(len=18),parameter::labels(17)=[character(len=18)::'New project…','Import IGES…','Screenshot', &
            'Undo modeling edit','Redo','Select','Orbit','Zoom out','Zoom in','Isometric view','Shaded','Wireframe', &
            'Front view','Top view','Parameter grid','Control net','Support']
        character(len=12),parameter::keys(17)=[character(len=12)::'Ctrl+N','Ctrl+I','Ctrl+Alt+S','Ctrl+Z','Ctrl+Shift+Z','', &
            '','Ctrl+-','Ctrl++','R','','','','','','','']
        character(len=26),parameter::transform_buttons(4)=[character(len=26)::'studioTranslateButton', &
            'studioRotateButton','studioScaleButton','studioAffineButton']
        character(len=6),parameter::transform_captions(4)=[character(len=6)::'Move','Rotate','Scale','Affine']
        character(len=20),parameter::transform_labels(4)=[character(len=20)::'Move','Rotate','Scale','Affine matrix…']
        integer,parameter::transform_icons(4)=[9,60,61,62]
        character(64),parameter :: tab_hints(7)=[character(64):: &
            'Open and save projects, exchange files, and set project units.', &
            'Create basic shapes or build geometry from profiles.', &
            'Select, transform, and measure scene objects.', &
            'Edit a patch’s parameters, control points, knots, and degree.', &
            'Connect patches, check joins, and edit shared controls.', &
            'Generate, inspect, display, and export structured meshes.', &
            'Set the camera, appearance, and workspace panels.']

        controls=>s%widgets%ribbon
        call build_opaque_panel(s%command_ribbon,frame,studio_ribbon_fill,'studioRibbonFrame',s%surface)
        call s%command_ribbon%setObjectName('studioRibbon')
        ribbon_layout=QGridLayout(frame)
        call ribbon_layout%setContentsMargins(1,1,1,1)
        call ribbon_layout%setSpacing(0)
        tabs=QTabWidget(frame)
        s%widgets%ribbon_tabs=tabs
        call ribbon_layout%addWidget(tabs,0,0)
        call tabs%setObjectName('studioRibbonTabs')
        call tabs%setAccessibleName('Ribbon categories')
        do i=1,size(s%toolbar_actions)
            call make_command(s%surface,trim(labels(i)),trim(keys(i)),s,codes(i),icons(i),action)
            s%toolbar_actions(i)=action
            if(i==6.or.i==7.or.i==11.or.i==12.or.i>=15)call s%toolbar_actions(i)%setCheckable(.true.)
        end do

        call make_command(s%surface,'Redo','Ctrl+Y',s,redo_patch_command,0,controls%alternate_redo_action)
        call make_command(s%surface,'Save project','Ctrl+S',s,save_project_command,3,save_action)
        call make_command(s%surface,'Quit','Ctrl+Q',s,close_window_command,0,action)
        call make_command(s%surface,'Parameter panels','Shift+F11',s,toggle_patch_inspection_command,33,controls%parameters_ribbon)
        call controls%parameters_ribbon%setCheckable(.true.)

        call build_project_ribbon(tabs)
        call build_create_ribbon(tabs)
        call build_model_ribbon(tabs)
        call build_patch_ribbon(tabs)
        call build_multipatch_ribbon(tabs)
        call build_mesh_ribbon(tabs)
        call build_view_ribbon(tabs)
        do i=1,size(tab_hints)
            call set_tab_status_hint(tabs,i-1,trim(tab_hints(i)))
        end do
        call tabs%setCurrentIndex(model_ribbon_tab)

        quick=QWidget(frame)
        call quick%setObjectName('studioRibbonQuickActions')
        quick_row=QHBoxLayout(quick)
        call quick_row%setContentsMargins(2,0,6,0)
        call quick_row%setSpacing(4)
        call build_quick_actions(quick,quick_row)
        call ribbon_layout%addWidget(quick,0,0,alignment=ior(Qt_AlignTop,Qt_AlignRight))
        call quick%raise()

    contains

        subroutine build_quick_actions(parent,row)
            type(QWidget),intent(in)::parent
            type(QHBoxLayout),intent(inout)::row
            type(QAction)::help,about

            call quick_action(parent,row,save_action,'studioQuickSaveButton')
            call quick_action(parent,row,s%toolbar_actions(4),'studioUndoControlsButton')
            call quick_action(parent,row,s%toolbar_actions(5),'studioRedoControlsButton')
            call quick_divider(parent,row)
            call quick_action(parent,row,fullscreen_action,'studioFullscreenButton')
            call s%panel_actions(8)%setIconText('Inspector')
            call quick_action(parent,row,s%panel_actions(8),'studioQuickInspectorButton',labeled=.true.)
            call make_command(s%surface,'Assistant','Ctrl+Shift+A',s,toggle_assistant_command,29,controls%assistant_toggle_action)
            call controls%assistant_toggle_action%setCheckable(.true.)
            call quick_action(parent,row,controls%assistant_toggle_action,'studioAssistantButton',labeled=.true.)
            call quick_divider(parent,row)
            call make_command(s%surface,'Help','F1',s,interaction_guide_command,16,help)
            call set_status_hint(help,studio_tooltip('Help','Open the interaction guide and keyboard shortcuts.','F1'))
            call quick_action(parent,row,help,'studioHelpButton')
            call make_command(s%surface,'About ForCAD Studio','',s,about_command,17,about)
            call about%setObjectName('studioAboutAction')
            call quick_action(parent,row,about,'studioAboutButton')
        end subroutine build_quick_actions

        subroutine build_create_ribbon(tabs)
            type(QTabWidget),intent(inout)::tabs
            type(QHBoxLayout)::groups,row
            type(QWidget)::group
            type(QAction)::copy_action

            call ribbon_page(tabs,'Create',30,groups)
            call ribbon_group(groups,'Points & curves',group,row)
            call primitive_tile(row,group,point_shape,'Point','studioCreatePointButton')
            call primitive_tile(row,group,line_shape,'Line','studioCreateLineButton')
            call primitive_tile(row,group,circle_shape,'Circle','studioCreateCircleButton')
            call primitive_tile(row,group,arc_shape,'Arc','studioCreateArcButton')
            call ribbon_group(groups,'Surfaces',group,row)
            call primitive_tile(row,group,rectangle_shape,'Rectangle','studioCreateRectangleButton')
            call primitive_tile(row,group,cylinder_shape,'Cylinder','studioCreateCylinderButton')
            call ribbon_group(groups,'Volumes',group,row)
            call primitive_tile(row,group,box_shape,'Box','studioCreateBoxButton')
            call primitive_tile(row,group,tube_shape,'Tube','studioCreateTubeButton')
            call ribbon_group(groups,'From profiles',group,row)
            call ribbon_command(row,group,'Extrude…','',extrude_command,77,'Extrude','studioExtrudeButton',controls%profile_actions(1))
            call ribbon_command(row,group,'Revolve…','',revolve_command,78,'Revolve','studioRevolveButton',controls%profile_actions(2))
            call ribbon_command(row,group,'Sweep…','',sweep_command,79,'Sweep','studioSweepButton',controls%profile_actions(3))
            call ribbon_command(row,group,'Loft…','',loft_command,80,'Loft','studioLoftButton',controls%profile_actions(4))
            call ribbon_group(groups,'Reuse',group,row)
            call ribbon_command(row,group,'Copy a starting patch…','',add_patch_command,11,'Copy patch','studioCopyPatchButton',copy_action)
            call copy_action%setObjectName('studioCopyFromPatchAction')
            call groups%addStretch()
        end subroutine build_create_ribbon

        subroutine primitive_tile(row,group,shape,caption,name)
            type(QHBoxLayout),intent(inout)::row
            type(QWidget),intent(in)::group
            integer,intent(in)::shape
            character(*),intent(in)::caption,name

            call ribbon_command(row,group,trim(primitive_names(shape))//'…','',primitive_command+shape,68+shape,caption,name)
        end subroutine primitive_tile

        subroutine build_model_ribbon(tabs)
            type(QTabWidget),intent(inout)::tabs
            type(QHBoxLayout)::groups,row
            type(QVBoxLayout)::column
            type(QWidget)::group
            integer::index

            call ribbon_page(tabs,'Model',24,groups)
            call ribbon_group(groups,'Selection',group,row)
            call ribbon_action(row,group,s%toolbar_actions(6),'studioModelSelectButton')
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Box selection','',box_selection_command,icon_box_selection,'Box', &
                'studioBoxSelectionButton',s%widgets%area_selection_actions(1),column=column)
            call ribbon_command(row,group,'Lasso selection','',lasso_selection_command,icon_lasso_selection,'Lasso', &
                'studioLassoSelectionButton',s%widgets%area_selection_actions(2),column=column)
            call s%widgets%area_selection_actions(1)%setCheckable(.true.)
            call s%widgets%area_selection_actions(2)%setCheckable(.true.)
            call set_status_hint(s%widgets%area_selection_actions(1),studio_tooltip('Box selection', &
                'Drag left to right for enclosed items or right to left for crossed items. Shift adds, Ctrl toggles, Esc cancels.'))
            call set_status_hint(s%widgets%area_selection_actions(2),studio_tooltip('Lasso selection', &
                'Draw around or across items in the viewport. Shift adds, Ctrl toggles, Esc cancels.'))
            call ribbon_command(row,group,'Clear selection','',clear_selection_command,icon_close,'Clear','studioModelClearSelectionButton')
            call ribbon_group(groups,'Recall',group,row)
            call ribbon_command(row,group,'Recall last selected curve','Ctrl+1',select_curve_command,23, &
                'Curve','studioRecallCurveButton')
            call ribbon_command(row,group,'Recall last selected surface','Ctrl+2',select_surface_command,24, &
                'Surface','studioRecallSurfaceButton')
            call ribbon_command(row,group,'Recall last selected volume','Ctrl+3',select_volume_command,8, &
                'Volume','studioRecallVolumeButton')
            call ribbon_group(groups,'Transform',group,row)
            do index=1,size(controls%transform_tools)
                call make_command(s%surface,trim(transform_labels(index)),'',s,translate_geometry_command+index-1, &
                    transform_icons(index),controls%transform_tools(index))
                if(index<4)call controls%transform_tools(index)%setCheckable(.true.)
                call ribbon_action(row,group,controls%transform_tools(index),trim(transform_buttons(index)), &
                    trim(transform_captions(index)))
            end do
            call ribbon_command(row,group,'Exact transform…','Ctrl+T',transform_command,38,'Exact transform…','studioNumericTransformButton')
            call ribbon_group(groups,'Objects',group,row)
            call ribbon_command(row,group,'Duplicate patch','',project_duplicate_command,11,'Duplicate', &
                'studioModelDuplicateButton',s%widgets%project_duplicate_action)
            call ribbon_command(row,group,'Rename patch…','',project_rename_command,68,'Rename', &
                'studioModelRenameButton',s%widgets%project_rename_action)
            call ribbon_command(row,group,'Delete selection','',project_remove_command,icon_delete,'Delete', &
                'studioModelDeleteButton',s%widgets%project_remove_action)
            call s%widgets%project_remove_action%setObjectName('studioProjectDeleteAction')
            call ribbon_group(groups,'Analysis',group,row)
            call ribbon_command(row,group,'Measure geometry…','',measure_geometry_command,icon_measure,'Measure', &
                'studioMeasureButton',controls%measure_action)
            call groups%addStretch()
        end subroutine build_model_ribbon

        subroutine ribbon_command(row,group,text,shortcut,command,icon_kind,display_text,object_name,action_out,primary,column,minimum_width)
            type(QHBoxLayout),intent(inout)::row
            type(QWidget),intent(in)::group
            character(*),intent(in)::text,shortcut
            integer,intent(in)::command,icon_kind
            character(*),intent(in),optional::display_text,object_name
            type(QAction),intent(out),optional::action_out
            logical,intent(in),optional::primary
            integer,intent(in),optional::minimum_width
            type(QVBoxLayout),intent(inout),optional::column
            type(QAction)::action
            type(QToolButton)::tile

            call make_command(s%surface,text,shortcut,s,command,icon_kind,action)
            call ribbon_action(row,group,action,object_name,display_text,tile,column)
            if(present(primary))call tile%setProperty('ribbonPrimary',primary)
            if(present(minimum_width))call tile%setMinimumWidth(minimum_width)
            if(present(action_out))action_out=action
        end subroutine ribbon_command

        subroutine build_project_ribbon(tabs)
            type(QTabWidget),intent(inout)::tabs
            type(QHBoxLayout)::groups,row
            type(QVBoxLayout)::column
            type(QWidget)::group
            type(QMenu)::menu

            call ribbon_page(tabs,'Project',82,groups)
            call ribbon_group(groups,'Files',group,row)
            call ribbon_command(row,group,'Home','',home_command,86,object_name='studioHomeButton')
            call ribbon_column(row,column)
            call ribbon_action(row,group,s%toolbar_actions(1),'studioNewProjectButton','New',column=column)
            call ribbon_command(row,group,'Open project…','Ctrl+O',open_project_command,82,'Open','studioOpenProjectButton',column=column)
            call ribbon_column(row,column)
            call ribbon_action(row,group,save_action,'studioSaveProjectButton','Save',column=column)
            call ribbon_command(row,group,'Save project as…','Ctrl+Shift+S',save_as_command,3,'Save as','studioSaveProjectAsButton',column=column)
            call ribbon_command(row,group,'Close project tab','Ctrl+W',close_project_tab_command,icon_close,'Close tab', &
                'studioCloseProjectTabButton',s%project_tabs%close)
            call s%project_tabs%close%setEnabled(.false.)
            call ribbon_group(groups,'Exchange',group,row)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Import FCAD project into this model…','',import_project_command,83,'Import FCAD', &
                'studioImportProjectButton',column=column)
            call ribbon_action(row,group,s%toolbar_actions(2),'studioImportButton','Import IGES',column=column)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Export geometry as VTK…','',export_vtk_command,icon_export,'Export VTK', &
                'studioExportVtkButton',column=column)
            call ribbon_command(row,group,'Export geometry as IGES…','',export_iges_command,icon_export,'Export IGES', &
                'studioExportIgesButton',column=column)
            call ribbon_command(row,group,'Export project as Fortran…','',export_fortran_command,icon_export,'Export Fortran', &
                'studioExportFortranButton')
            call ribbon_action(row,group,s%toolbar_actions(3),'studioCaptureButton','Screenshot')
            call ribbon_group(groups,'Clipboard',group,row)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Copy selected objects','Ctrl+C',copy_geometry_command,11,'Copy','studioCopyGeometryButton', &
                s%project_tabs%copy,column=column)
            call ribbon_command(row,group,'Paste objects into this project','Ctrl+V',paste_geometry_command,83,'Paste','studioPasteGeometryButton', &
                s%project_tabs%paste,column=column)
            call set_status_hint(s%project_tabs%copy,studio_tooltip('Copy selected objects', &
                'Copy patches, multipatches or current meshes between project tabs.','Ctrl+C'))
            call set_status_hint(s%project_tabs%paste,studio_tooltip('Paste objects', &
                'Add the copied objects to the current project.','Ctrl+V'))
            call s%project_tabs%copy%setEnabled(.false.)
            call s%project_tabs%paste%setEnabled(.false.)
            call ribbon_group(groups,'Workspace',group,row)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Working folder…','',project_folder_command,82,'Folder','studioProjectFolderButton',column=column)
            call ribbon_command(row,group,'Project units…','',units_command,84,'Units','studioProjectUnitsButton',column=column)
            call ribbon_menu(row,group,'Options',33,menu,'studioApplicationMenuButton', &
                description='Reset the view and tools or restore the workspace defaults.')
            call menu_action(menu,'Reset view and tools','Ctrl+R',s,reset_controls_command,0)
            call menu_action(menu,'Restore workspace defaults','',s,reset_presentation_command,0)
            call groups%addStretch()
        end subroutine build_project_ribbon

        subroutine build_patch_ribbon(tabs)
            type(QTabWidget),intent(inout)::tabs
            type(QHBoxLayout)::groups,row
            type(QVBoxLayout)::column
            type(QWidget)::group,page

            call ribbon_page(tabs,'Patch',33,groups)
            call groups%setSpacing(compact_group_spacing)
            page=tabs%widget(tabs%count()-1)
            call page%setProperty('compactCommandRibbon',.true.)
            call ribbon_group(groups,'Control',group,row)
            call ribbon_command(row,group,'Control points and weights…','Ctrl+E',edit_controls_command,22, &
                'Points + weights','studioEditControlsButton',edit_controls_action)

            call ribbon_group(groups,'Parameters',group,row)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Set patch parameter range…','',edit_interval_command,63,'Set domain', &
                'studioPatchRangeButton',column=column)
            call ribbon_action(row,group,controls%parameters_ribbon,'studioParametersButton','Parameter panels',column=column)

            call ribbon_group(groups,'Knots',group,row)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Insert knot…','',insert_knot_command,65,'Insert knot','studioInsertKnotButton',column=column)
            call ribbon_command(row,group,'Remove knot…','',remove_knot_command,66,'Remove knot','studioRemoveKnotButton',column=column)
            call ribbon_command(row,group,'Edit knot values…','',edit_knot_values_command,68,'Knot values','studioEditKnotValuesButton')

            call ribbon_group(groups,'Degree & spans',group,row)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Elevate degree…','',elevate_degree_command,icon_degree_elevate,'Elevate degree','studioElevateDegreeButton',column=column)
            call ribbon_command(row,group,'Reduce degree…','',reduce_degree_command,icon_degree_reduce,'Reduce degree','studioReduceDegreeButton',column=column)
            call ribbon_command(row,group,'Refine patch…','',refine_command,64,'Refine patch', &
                'studioRefinePatchButton',controls%refinement_action)

            call ribbon_group(groups,'Topology',group,row)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Split patch…','',split_patch_command,90,'Split patch', &
                'studioSplitPatchButton',controls%patch_tool_actions(1),column=column)
            call ribbon_command(row,group,'Reverse direction…','',reverse_patch_command,92,'Reverse direction', &
                'studioReversePatchButton',controls%patch_tool_actions(3),column=column)
            call controls%patch_tool_actions(1)%setObjectName('studioSplitPatchAction')
            call controls%patch_tool_actions(3)%setObjectName('studioReversePatchAction')
            call ribbon_group(groups,'Extract',group,row)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Extract parameter range…','',extract_patch_command,91,'Extract range', &
                'studioExtractPatchButton',controls%patch_tool_actions(2),column=column)
            call ribbon_command(row,group,'Extract section…','',extract_section_command,icon_extract_section,'Extract section', &
                'studioExtractSectionButton',controls%patch_tool_actions(4),column=column)
            call ribbon_command(row,group,'Extract volume face…','',extract_face_command,icon_extract_face,'Extract face', &
                'studioExtractFaceButton',controls%patch_tool_actions(5))
            call controls%patch_tool_actions(2)%setObjectName('studioExtractPatchAction')
            call controls%patch_tool_actions(4)%setObjectName('studioExtractSectionAction')
            call controls%patch_tool_actions(5)%setObjectName('studioExtractFaceAction')

            call groups%addStretch()
        end subroutine build_patch_ribbon

        subroutine build_mesh_ribbon(tabs)
            type(QTabWidget),intent(inout)::tabs
            type(QHBoxLayout)::groups,row
            type(QVBoxLayout)::column
            type(QWidget)::group

            call ribbon_page(tabs,'Mesh',icon_mesh,groups)

            call ribbon_group(groups,'Build',group,row)
            call ribbon_command(row,group,'Generate mesh…','',open_mesh_dialog_command,icon_generate_mesh, &
                'Generate…','studioGenerateMeshButton',controls%mesh_generate_action,primary=.true.)

            call ribbon_group(groups,'Inspect',group,row)
            call ribbon_command(row,group,'Pick mesh cells','',mesh_pick_command,icon_mesh_pick, &
                'Pick mesh cells','studioPickMeshCellButton',controls%mesh_pick_action)
            call controls%mesh_pick_action%setCheckable(.true.)

            call ribbon_group(groups,'Display',group,row)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Show meshes only','',mesh_only_command,icon_mesh, &
                'Meshes only','studioMeshOnlyButton',controls%mesh_only_action,column=column)
            call ribbon_command(row,group,'Show or hide meshes for the selection','',mesh_edges_command,31,'Show meshes', &
                'studioMeshVisibilityButton',controls%mesh_edges_action,column=column)
            call ribbon_command(row,group,'Shade mesh faces','',mesh_shaded_command,icon_mesh_shaded, &
                'Filled mesh faces','studioMeshShadedButton',controls%mesh_shaded_action)
            call controls%mesh_only_action%setCheckable(.true.)
            call controls%mesh_shaded_action%setCheckable(.true.)
            call controls%mesh_edges_action%setCheckable(.true.)
            call controls%mesh_edges_action%setObjectName('studioMeshVisibilityAction')

            call ribbon_group(groups,'Volume mesh',group,row)
            controls%volume_mesh_controls=QWidget(group)
            call row%addWidget(controls%volume_mesh_controls)

            call ribbon_group(groups,'Export',group,row)
            call ribbon_command(row,group,'Export mesh (VTK)…','',export_mesh_command,icon_export, &
                'Export mesh (VTK)','studioExportMeshButton',controls%mesh_export_action)
            call groups%addStretch()
        end subroutine build_mesh_ribbon

        subroutine build_view_ribbon(tabs)
            type(QTabWidget),intent(inout)::tabs
            type(QHBoxLayout)::groups,row
            type(QVBoxLayout)::column
            type(QWidget)::group,page
            type(QAction)::action
            type(QToolButton)::shaded_button,guide_button
            integer::i,panel
            integer,parameter::panel_ids(4)=[1,8,5,6],panel_icons(4)=[2,40,63,23]
            character(len=20),parameter::panel_names(4)=[character(len=20)::'Scene','Inspector','Parameter Space','Basis Plot']
            character(len=32),parameter::panel_buttons(4)=[character(len=32)::'studioViewSceneButton','studioViewInspectorButton', &
                'studioViewParameterSpaceButton','studioViewBasisPlotButton']
            integer,parameter::camera_toolbar(7)=[10,13,0,0,0,14,0]
            character(len=26),parameter::camera_buttons(7)=[character(len=26)::'studioCameraIsoButton','studioCameraFrontButton', &
                'studioCameraBackButton','studioCameraLeftButton','studioCameraRightButton','studioCameraTopButton','studioCameraBottomButton']
            integer,parameter::layer_icons(viewport_layer_count)=[8,21,22,icon_bounding_vertex,24]
            character(len=20),parameter::layer_labels(viewport_layer_count)=[character(len=20)::'Shaded geometry','Wire overlay', &
                'Control net','Control points','Parameter lines']
            character(len=20),parameter::presets(*)=[character(len=20)::'Clean preview','Control structure']
            character(len=32),parameter::preset_buttons(*)=[character(len=32):: &
                'studioViewCleanPresetButton','studioViewStructurePresetButton']
            integer,parameter::preset_icons(*)=[8,22]
            character(len=28),parameter::style_action_names(viewport_style_count)=[character(len=28):: &
                'studioViewportShaded','studioViewportShadedEdges','studioViewportWireframe','studioViewportFastShaded']

            do i=1,size(s%viewport_layers)
                call make_command(s%surface,trim(layer_labels(i)),'',s,viewport_layer_base+i-1,layer_icons(i),s%viewport_layers(i))
                call s%viewport_layers(i)%setCheckable(.true.)
            end do

            call ribbon_page(tabs,'View',25,groups)
            call groups%setSpacing(compact_group_spacing)
            page=tabs%widget(tabs%count()-1)
            call page%setProperty('compactCommandRibbon',.true.)
            call ribbon_group(groups,'Camera',group,row)
            call ribbon_column(row,column)
            call ribbon_action(row,group,s%toolbar_actions(7),'studioOrbitToolButton',column=column)
            call ribbon_command(row,group,'Fit scene','F',fit_scene_command,19,'Fit all','studioViewFitSceneButton',column=column)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Fit selection','Shift+F',fit_selection_command,icon_fit_selection,'Fit selection', &
                'studioViewFitSelectionButton',controls%view_fit_action,column=column,minimum_width=102)
            call ribbon_command(row,group,'Fit around panels','',fit_around_panels_command,icon_fit_panels,'Fit area', &
                'studioFitClearAreaButton',column=column)
            call ribbon_column(row,column)
            call compact_icon_action(row,group,s%toolbar_actions(9),'studioViewZoomInButton',column)
            call compact_icon_action(row,group,s%toolbar_actions(8),'studioViewZoomOutButton',column)

            do i=1,size(camera_commands)
                if(i==1.or.mod(i,2)==0)call ribbon_column(row,column)
                if(camera_toolbar(i)>0)then
                    action=s%toolbar_actions(camera_toolbar(i))
                else
                    call make_command(s%surface,trim(camera_names(i)),'',s,camera_commands(i),camera_icons(i),action)
                end if
                call action%setCheckable(.true.)
                s%widgets%camera_view_actions(i)=action
                call compact_icon_action(row,group,action,trim(camera_buttons(i)),column)
                if(i==1)then
                    call make_command(s%surface,'Interactive Quality','',s,interactive_quality_command,33,controls%interactive_quality_action)
                    call controls%interactive_quality_action%setCheckable(.true.)
                    call controls%interactive_quality_action%setChecked(s%interactive_quality)
                    call compact_icon_action(row,group,controls%interactive_quality_action,'studioInteractiveQualityButton',column)
                end if
            end do

            call ribbon_group(groups,'Display',group,row)
            associate(styles=>s%widgets%viewport_style_actions)
                do i=viewport_shaded,viewport_fast_shaded
                    call make_command(s%surface,trim(viewport_style_names(i)),'',s,viewport_style_commands(i),viewport_style_icons(i),styles(i))
                    call styles(i)%setCheckable(.true.)
                    call styles(i)%setObjectName(trim(style_action_names(i)))
                end do
                call ribbon_column(row,column)
                call ribbon_action(row,group,styles(viewport_shaded),'studioShadedFacesButton',tile=shaded_button,column=column)
                call shaded_button%setMinimumWidth(82)
                call ribbon_action(row,group,styles(viewport_edges),'studioShadedButton',tile=shaded_button,column=column)
                call shaded_button%setMinimumWidth(82)
                call ribbon_column(row,column)
                call ribbon_action(row,group,styles(viewport_fast_shaded),'studioFastShadedButton',column=column)
                call ribbon_action(row,group,styles(viewport_wireframe),'studioWireframeButton',column=column)
            end associate
            call ribbon_column(row,column)
            do i=1,size(presets)
                call make_command(s%surface,trim(presets(i)),'',s,clean_display_command+i-1,preset_icons(i),s%viewport_presets(i))
                call s%viewport_presets(i)%setCheckable(.true.)
                call compact_icon_action(row,group,s%viewport_presets(i),trim(preset_buttons(i)),column)
            end do

            call ribbon_group(groups,'Guides',group,row)
            call ribbon_column(row,column)
            call compact_icon_action(row,group,s%viewport_layers(3),'studioControlNetButton',column)
            call compact_icon_action(row,group,s%viewport_layers(4),'studioControlPointsButton',column)
            call ribbon_column(row,column)
            call ribbon_action(row,group,s%viewport_layers(parameter_lines_layer),'studioParameterLinesButton',tile=guide_button,column=column)
            call guide_button%setMinimumWidth(126)
            call ribbon_command(row,group,'Knot boundaries','',physical_knots_command,63, &
                'Knot lines','studioPhysicalKnotsRibbonButton',controls%physical_knot_action,column=column)
            call controls%physical_knot_action%setCheckable(.true.)
            call ribbon_column(row,column)
            call compact_icon_action(row,group,s%viewport_layers(2),'studioWireOverlayButton',column)
            call make_command(s%surface,'Ground grid','',s,world_grid_command,12,s%viewport_grid)
            call s%viewport_grid%setCheckable(.true.)
            call compact_icon_action(row,group,s%viewport_grid,'studioGroundGridButton',column)

            call ribbon_group(groups,'Appearance',group,row)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Object color…','',color_command,37,'Color','studioObjectColorButton',column=column)
            call ribbon_command(row,group,'Translucent','',translucent_command,21,'Translucent','studioTranslucentButton',column=column)
            call ribbon_column(row,column)
            call make_command(s%surface,'Show or hide selection','Ctrl+Shift+H',s,toggle_patch_visibility_command,31,controls%view_visibility_action)
            call ribbon_action(row,group,controls%view_visibility_action,'studioSelectionVisibilityButton','Show',column=column)

            call ribbon_group(groups,'Panels',group,row)
            s%panel_actions(4)=QAction('3D View',s%surface)
            call s%panel_actions(4)%setCheckable(.true.)
            call s%panel_actions(4)%setChecked(.true.)
            do i=1,size(panel_ids)
                if(mod(i,2)==1)call ribbon_column(row,column)
                panel=panel_ids(i)
                call make_command(s%surface,trim(panel_names(i)),'',s,toggle_panel_command_base+panel,panel_icons(i),s%panel_actions(panel))
                call s%panel_actions(panel)%setCheckable(.true.)
                call s%panel_actions(panel)%setChecked(.true.)
                if(panel==5.or.panel==6)then
                    call compact_icon_action(row,group,s%panel_actions(panel),trim(panel_buttons(i)),column)
                else
                    call ribbon_action(row,group,s%panel_actions(panel),trim(panel_buttons(i)),column=column)
                end if
            end do
            call ribbon_column(row,column)
            call make_command(s%surface,'Restore panel layout','',s,restore_layout_command,icon_refresh,action)
            call compact_icon_action(row,group,action,'studioRestoreLayoutButton',column)
            call make_command(s%surface,'Canvas only','Ctrl+Shift+F',s,canvas_command,icon_maximize,action)
            call compact_icon_action(row,group,action,'studioCanvasOnlyButton',column)
            call groups%addStretch()
        end subroutine build_view_ribbon

        subroutine build_multipatch_ribbon(tabs)
            type(QTabWidget),intent(inout)::tabs
            type(QHBoxLayout)::groups,row
            type(QVBoxLayout)::column
            type(QWidget)::group
            type(QAction)::remove_action

            call ribbon_page(tabs,'Multipatch',87,groups)
            call ribbon_group(groups,'Selection',group,row)
            call ribbon_command(row,group,'Multipatches','',pick_multipatch_command,87, &
                'Pick multipatches','studioPickMultipatchesButton')
            call ribbon_group(groups,'Assembly',group,row)
            call ribbon_command(row,group,'New multipatch…','Ctrl+Shift+M',new_multipatch_command,87, &
                'New multipatch…','studioNewMultipatchButton',controls%assembly_connect_action,primary=.true.)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Edit multipatch…','',edit_multipatch_command,88, &
                'Edit multipatch…','studioAssemblyEditButton',column=column)
            call ribbon_command(row,group,'Dissolve multipatch…','',remove_multipatch_command,36, &
                'Dissolve multipatch','studioDissolveMultipatchButton',remove_action,column=column)
            call remove_action%setObjectName('studioDissolveMultipatchAction')

            call ribbon_group(groups,'Controls',group,row)
            call ribbon_command(row,group,'Control points','',pick_controls_command,22, &
                'Pick points','studioPickSharedControlsButton')
            call ribbon_action(row,group,edit_controls_action,'studioAssemblyControlsButton','Points + weights')

            call ribbon_group(groups,'Refine',group,row)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Insert multipatch knot…','',insert_assembly_knot_command,65, &
                'Insert knot','studioAssemblyInsertKnotButton',controls%assembly_refine_actions(1),column=column)
            call ribbon_command(row,group,'Elevate multipatch degree…','',elevate_assembly_degree_command,icon_degree_elevate, &
                'Elevate degree','studioAssemblyElevateButton',controls%assembly_refine_actions(2),column=column)

            call ribbon_group(groups,'Interfaces',group,row)
            call ribbon_command(row,group,'Check interfaces','',validate_multipatch_command,89, &
                'Check interfaces','studioAssemblyValidateButton')

            call ribbon_group(groups,'Members',group,row)
            call ribbon_column(row,column)
            call ribbon_command(row,group,'Fit members','',fit_selection_command,icon_fit_selection, &
                'Fit members','studioFitMultipatchButton',controls%assembly_fit_action,column=column)
            call ribbon_command(row,group,'Show or hide members','',toggle_patch_visibility_command,31, &
                'Show members','studioMultipatchVisibilityButton',controls%assembly_visibility_action,column=column)
            call groups%addStretch()
        end subroutine build_multipatch_ribbon

    end subroutine build_ribbon

    subroutine ribbon_page(tabs,caption,glyph,groups)
        type(QTabWidget),intent(inout)::tabs
        character(*),intent(in)::caption
        integer,intent(in)::glyph
        type(QHBoxLayout),intent(out)::groups
        type(QWidget)::page
        type(QScrollArea)::strip

        page=QWidget()
        call page%setObjectName('studioRibbonPage')
        groups=QHBoxLayout(page)
        call groups%setContentsMargins(4,2,4,0)
        call groups%setSpacing(5)
        strip=QScrollArea(tabs)
        call strip%setObjectName('studioRibbonStrip')
        call strip%setAccessibleName(caption//' commands')
        call strip%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
        call strip%setFixedHeight(82)
        call strip%setHorizontalScrollBarPolicy(Qt_ScrollBarAsNeeded)
        call strip%setVerticalScrollBarPolicy(Qt_ScrollBarAlwaysOff)
        call strip%setWidgetResizable(.true.)
        call strip%setWidget(page)
        call tabs%addTab(strip,caption)
        call tabs%setTabIcon(tabs%count()-1,line_icon(glyph))
    end subroutine ribbon_page

    subroutine ribbon_menu(row,group,label,glyph,menu,object_name,description,column)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::group
        character(*),intent(in)::label
        integer,intent(in)::glyph
        type(QMenu),intent(out)::menu
        character(*),intent(in),optional::object_name,description
        type(QVBoxLayout),intent(inout),optional::column
        type(QToolButton)::button

        button=QToolButton(group)
        if(present(object_name))call button%setObjectName(object_name)
        call button%setText(label)
        call button%setIcon(line_icon(glyph))
        call button%setAccessibleName(label)
        call set_status_hint(button,studio_tooltip(label,description))
        call style_ribbon_tile(button,stacked=present(column))
        menu=QMenu(label,button)
        call register_status_menu(menu)
        call button%setMenu(menu)
        call button%setPopupMode(QToolButton_InstantPopup)
        if(present(column))then
            call column%addWidget(button)
        else
            call row%addWidget(button)
        end if
    end subroutine ribbon_menu

    subroutine quick_divider(parent,row)
        type(QWidget),intent(in)::parent
        type(QHBoxLayout),intent(inout)::row
        type(QWidget)::divider

        divider=QWidget(parent)
        call divider%setObjectName('studioRibbonQuickDivider')
        call divider%setFixedSize(1,18)
        call row%addWidget(divider,alignment=Qt_AlignVCenter)
    end subroutine quick_divider

    subroutine quick_action(parent,row,action,name,labeled)
        type(QWidget),intent(in)::parent
        type(QHBoxLayout),intent(inout)::row
        type(QAction),intent(in)::action
        character(*),intent(in)::name
        logical,intent(in),optional::labeled
        type(QToolButton)::button
        logical::show_label

        show_label=.false.
        if(present(labeled))show_label=labeled
        button=QToolButton(parent)
        call button%setObjectName(name)
        call button%setDefaultAction(action)
        call watch_status_hint(button)
        call button%setAccessibleName(action%text())
        call style_quick_button(button,show_label)
        call row%addWidget(button)
    end subroutine quick_action

    subroutine style_quick_button(button,labeled)
        type(QToolButton),intent(inout)::button
        logical,intent(in)::labeled

        call button%setProperty('ribbonQuickAction',.true.)
        call button%setFixedHeight(32)
        call button%setMinimumWidth(32)
        call button%setIconSize(18,18)
        call button%setToolButtonStyle(merge(Qt_ToolButtonTextBesideIcon,Qt_ToolButtonIconOnly,labeled))
    end subroutine style_quick_button

    subroutine ribbon_group(groups,caption,group,row)
        type(QHBoxLayout),intent(inout)::groups
        character(len=*),intent(in)::caption
        type(QWidget),intent(out)::group
        type(QHBoxLayout),intent(out)::row
        type(QVBoxLayout)::layout
        type(QLabel)::label
        integer::side_margin

        side_margin=merge(2,4,groups%spacing()==compact_group_spacing)
        group=QWidget(); call group%setObjectName('ribbonGroup')
        call group%setProperty('firstRibbonGroup',groups%count()==0)
        call group%setSizePolicy(QSizePolicy_Fixed,QSizePolicy_Fixed)
        layout=QVBoxLayout(group); call layout%setContentsMargins(side_margin,2,side_margin,2); call layout%setSpacing(1)
        row=QHBoxLayout(); call row%setContentsMargins(0,0,0,0); call row%setSpacing(2)
        call layout%addLayout(row)
        label=QLabel(caption); call label%setObjectName('ribbonCaption')
        call label%setAlignment(Qt_AlignCenter)
        call layout%addWidget(label)
        call groups%addWidget(group)
    end subroutine ribbon_group

    subroutine ribbon_column(row,column)
        type(QHBoxLayout),intent(inout)::row
        type(QVBoxLayout),intent(out)::column

        column=QVBoxLayout()
        call column%setContentsMargins(0,0,0,0)
        call column%setSpacing(ribbon_tile_height-2*ribbon_stacked_height)
        call row%addLayout(column)
    end subroutine ribbon_column

    subroutine style_ribbon_tile(button,stacked)
        type(QToolButton),intent(inout)::button
        logical,intent(in),optional::stacked
        logical::compact

        compact=.false.
        if(present(stacked))compact=stacked
        call button%setProperty('ribbonTile',.true.)
        call button%setMinimumWidth(47)
        if(compact)then
            call button%setProperty('ribbonStacked',.true.)
            call button%setFixedHeight(ribbon_stacked_height)
            call button%setIconSize(17,17)
            call button%setToolButtonStyle(Qt_ToolButtonTextBesideIcon)
        else
            call button%setFixedHeight(ribbon_tile_height)
            call button%setIconSize(26,26)
            call button%setToolButtonStyle(Qt_ToolButtonTextUnderIcon)
        end if
    end subroutine style_ribbon_tile

    subroutine compact_icon_action(row,group,action,object_name,column)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::group
        type(QAction),intent(inout)::action
        character(*),intent(in)::object_name
        type(QVBoxLayout),intent(inout)::column
        type(QToolButton)::button

        call ribbon_action(row,group,action,object_name,tile=button,column=column)
        call button%setFixedWidth(27)
        call button%setToolButtonStyle(Qt_ToolButtonIconOnly)
    end subroutine compact_icon_action

    subroutine ribbon_action(row,group,action,object_name,display_text,tile,column)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::group
        type(QAction),intent(inout)::action
        character(len=*),intent(in),optional::object_name,display_text
        type(QToolButton),intent(out),optional::tile
        type(QVBoxLayout),intent(inout),optional::column
        type(QToolButton)::button

        if(present(display_text))call action%setIconText(display_text)
        button=QToolButton(group)
        if(present(object_name))call button%setObjectName(object_name)
        call button%setDefaultAction(action)
        call watch_status_hint(button)
        call style_ribbon_tile(button,stacked=present(column))
        call button%setAccessibleName(action%text())
        if(present(column))then
            call column%addWidget(button)
        else
            call row%addWidget(button)
        end if
        if(present(tile))tile=button
    end subroutine ribbon_action
end module studio_ribbon
