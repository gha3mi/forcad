submodule(studio_ui) studio_ui_layout
    use qt, only: QActionGroup, QButtonGroup, QHeaderView, QHeaderView_Fixed, QHeaderView_Stretch, &
        QPointF, QSlider, QShortcut, Qt_WidgetShortcut, Qt_ToolButtonTextOnly, QToolButton_MenuButtonPopup
    use studio_live_viewport, only: pick_mesh_cell,transform_gizmo_hint
    use studio_viewport_camera, only: sync_viewport_orientation
    use studio_commands, only: patch_ribbon_tab,multipatch_ribbon_tab,mesh_ribbon_tab
    use studio_icons, only: icon_scene_filter,viewport_style_icons,icon_box_selection,icon_lasso_selection,icon_fullscreen
    use studio_viewport_scene, only: pick_selected_control,pick_shared_control,viewport_patch_target
    use studio_selection_status, only: selection_status_text,mesh_counts_text
    use studio_status_hints, only: register_status_hints,register_status_menu
    use studio_theme, only: studio_panel_fill,studio_viewport_controls_fill,build_opaque_panel
    use studio_command_presentation, only: command_visual
    implicit none(type, external)
    integer, parameter :: workspace_side_inset=12, workspace_bottom_inset=10, workspace_gap=8
    integer, parameter :: panel_inset=10
    integer, parameter :: quick_action_label_width=1280, assistant_min_width=318
    character(*), parameter :: center_column_name='studioCenterColumn', inspector_column_name='studioInspectorColumn'
    character(*), parameter :: assistant_space_name='studioAssistantSpace'
    character(*), parameter :: camera_guidance='Orbit: Shift+middle or right drag · Pan: Middle drag · Zoom: Scroll'
    integer, parameter :: viewer_button_width=27, viewport_button_width=34, viewer_button_height=26, viewer_icon_size=18
contains
    module procedure build_ui
        type(QWidget) :: root, brand, left, center, right, lower_row
        type(QWidget) :: p, domain_panel, basis_panel, basis_modes
        type(QVBoxLayout) :: root_layout, workspace_layout, left_layout, p_layout, basis_layout
        type(QHBoxLayout) :: brand_layout, center_layout, lower_layout, tabs
        type(QLabel) :: label
        type(QIcon) :: brand_icon
        type(QPushButton) :: button
        type(QAction) :: find_action, fullscreen_action
        type(QShortcut) :: delete_selection
        type(QtfEventFilter) :: window_events, assistant_columns_events, assistant_space_events
        integer :: i
        character(len=16), parameter :: basis_text(2) = [character(len=16) :: 'Basis functions', 'First derivative']

        call initialize_studio_patch_catalog(s%document%patches)

        if(present(empty_document))then
            if(empty_document)then
                call reset_studio_project(s%document%scene,s%document%multipatches,catalog=s%document%patches)
                s%document%show_templates=.false.
            end if
        end if
        s%window = QMainWindow()
        call s%window%setWindowTitle(studio_product_name)
        call s%window%resize(studio_default_width,studio_default_height)
        call s%window%setMinimumSize(studio_min_width,studio_min_height)
        call s%window%setStyleSheet(studio_style_sheet())
        call s%window%setWindowIcon(line_icon(17))
        call register_status_hints(s%window,on_status_hint,s,s%status)
        root = QWidget(); call root%setObjectName('root')
        root_layout = QVBoxLayout(root)
        call root_layout%setContentsMargins(0,0,0,0); call root_layout%setSpacing(0)
        brand = QWidget(); call brand%setObjectName('brand'); call brand%setFixedHeight(36)
        brand_layout = QHBoxLayout(brand)
        call brand_layout%setContentsMargins(14,0,8,0); call brand_layout%setSpacing(8)
        label=QLabel(''); call label%setFixedSize(23,23)
        call label%setScaledContents(.true.)
        brand_icon=line_icon(17,label); call brand_icon%delete()
        call brand_layout%addWidget(label)
        label = QLabel(studio_product_name)
        call label%setObjectName('brandTitle'); call brand_layout%addWidget(label)
        s%widgets%project_caption=QLabel('')
        call s%widgets%project_caption%setObjectName('studioDocumentTitle')
        call s%widgets%project_caption%setTextFormat(Qt_PlainText)
        call s%widgets%project_caption%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Preferred)
        call brand_layout%addWidget(s%widgets%project_caption,1)
        call build_document_tabs(s,brand_layout)
        call make_command(brand,'Find command','Ctrl+K',s,find_command,27,find_action)
        s%widgets%find_command_button=QToolButton(brand)
        call s%widgets%find_command_button%setDefaultAction(find_action)
        call watch_status_hint(s%widgets%find_command_button)
        call s%widgets%find_command_button%setToolButtonStyle(Qt_ToolButtonTextBesideIcon)
        call s%widgets%find_command_button%setObjectName('studioFindCommandButton')
        call s%widgets%find_command_button%setFixedHeight(26)
        call brand_layout%addWidget(s%widgets%find_command_button)
        call title_button(button,brand_layout,s,minimize_window_command)
        call title_button(s%widgets%maximize_button,brand_layout,s,maximize_window_command)
        call title_button(button,brand_layout,s,close_window_command)
        call root_layout%addWidget(brand)
        s%surface%QtfWidget=QtfWidget()
        s%workspace=s%surface%QtfWidget
        s%panels(4)=s%surface%QWidget
        call s%surface%setMinimumSize(298,230)
        call s%surface%setAccessibleName('Interactive 3D scene')
        call s%surface%setMouseTracking(.true.)
        call s%surface%setObjectName('surfaceCanvas')
        call s%workspace%resizeEvent%connect(on_workspace_resize,s)
        workspace_layout=QVBoxLayout(s%surface)
        call workspace_layout%setContentsMargins(workspace_side_inset,workspace_gap,workspace_side_inset,workspace_bottom_inset)
        call workspace_layout%setSpacing(workspace_gap)
        call make_command(s%window,'Fullscreen','F11',s,fullscreen_command,icon_fullscreen,fullscreen_action)
        call fullscreen_action%setObjectName('studioFullscreenAction')
        call fullscreen_action%setCheckable(.true.)
        call fullscreen_action%setProperty('autoRepeat',.false.)
        call set_status_hint(fullscreen_action,studio_tooltip('Fullscreen','Fill the screen. Repeat to restore the window.','F11'))
        s%widgets%fullscreen_action=fullscreen_action
        call build_ribbon(s,make_command,menu_action,fullscreen_action)
        window_events=QtfEventFilter(s%window)
        call window_events%showEvent%connect(on_window_shown,s)
        call window_events%resizeEvent%connect(on_window_resize,s)
        call build_volume_mesh_controls(s)
        call workspace_layout%addWidget(s%command_ribbon)
        s%columns = QSplitter(Qt_Horizontal)
        call s%columns%setObjectName('studioWorkspaceColumns')
        call s%columns%setChildrenCollapsible(.false.); call s%columns%setHandleWidth(workspace_gap)
        left = QWidget(); s%left_column=left; call left%setMinimumWidth(208)
        left_layout = QVBoxLayout(left)
        call left_layout%setContentsMargins(0,0,0,0); call left_layout%setSpacing(8)
        call create_panel('Scene',p,p_layout,s,1)
        call build_project_tree(p_layout,s)
        call left_layout%addWidget(p,1)
        center = QWidget(); call center%setObjectName(center_column_name)
        center_layout = QHBoxLayout(center)
        call center_layout%setContentsMargins(0,0,0,0); call center_layout%setSpacing(workspace_gap)
        s%center_rows = QSplitter(Qt_Vertical)
        call s%center_rows%setChildrenCollapsible(.false.); call s%center_rows%setHandleWidth(workspace_gap)
        s%viewport_aperture=QtfWidget()
        call s%viewport_aperture%setObjectName('studioViewportAperture')
        call set_status_hint(s%viewport_aperture,camera_guidance)
        call s%viewport_aperture%setMinimumSize(285,230)
        call s%viewport_aperture%resizeEvent%connect(on_viewport_aperture_resize,s)
        call build_viewport_controls(s)
        call s%center_rows%addWidget(s%viewport_aperture)
        lower_row = QWidget(); s%lower_row=lower_row; lower_layout = QHBoxLayout(lower_row)
        call lower_layout%setContentsMargins(0,0,0,0)
        call create_panel('Basis Functions (U)',basis_panel,basis_layout,s,6)
        call basis_panel%setMinimumWidth(260)
        basis_modes=QWidget(); call basis_modes%setObjectName('studioPlotControls')
        tabs=QHBoxLayout(basis_modes)
        call tabs%setContentsMargins(panel_inset,6,panel_inset,4); call tabs%setSpacing(4)
        do i=1,size(s%basis_buttons)
            s%basis_buttons(i)=QPushButton(trim(basis_text(i)))
            call s%basis_buttons(i)%setObjectName('axisTab')
            call s%basis_buttons(i)%setCheckable(.true.)
            call s%basis_buttons(i)%setChecked(i==1)
            call connect_button(s%basis_buttons(i),s,basis_functions_command+i-1)
            call tabs%addWidget(s%basis_buttons(i))
        end do
        call set_status_hint(s%basis_buttons(1),'Basis functions on the active axis')
        call set_status_hint(s%basis_buttons(2),'First derivative on the active axis')
        call tabs%addStretch()
        call basis_layout%addWidget(basis_modes)
        s%basis=QtfWidget(); call s%basis%setMinimumSize(250,155)
        call s%basis%setObjectName('basisCanvas')
        call s%basis%setAccessibleName('Basis functions and derivatives plot')
        call basis_layout%addWidget(s%basis,1)
        call lower_layout%addWidget(basis_panel,51)
        call s%center_rows%addWidget(lower_row)
        call s%center_rows%setSizes([488,293])
        call center_layout%addWidget(s%center_rows,63)
        call create_panel('Domain (U–V)',domain_panel,p_layout,s,5)
        call domain_panel%setMinimumWidth(276)
        call domain_panel%setMaximumWidth(435)
        call build_domain(p_layout,s)
        call center_layout%addWidget(domain_panel,37)
        call create_panel('Inspector',right,p_layout,s,8)
        call right%setObjectName(inspector_column_name)
        call right%setMinimumWidth(330); call right%setMaximumWidth(455)
        call build_inspector(p_layout,s)
        call domain_panel%showEvent%connect(sync_parameter_panels,s)
        call domain_panel%hideEvent%connect(sync_parameter_panels,s)
        call basis_panel%showEvent%connect(sync_parameter_panels,s)
        call basis_panel%hideEvent%connect(sync_parameter_panels,s)
        call right%showEvent%connect(on_inspector_visibility,s)
        call right%hideEvent%connect(on_inspector_visibility,s)
        call s%columns%addWidget(left); call s%columns%addWidget(center); call s%columns%addWidget(right)
        call s%columns%setSizes([264,1041,342]); call s%columns%setStretchFactor(1,1)
        call s%columns%splitterMoved%connect(on_scene_divider_moved,s)
        s%assistant_columns=QSplitter(Qt_Horizontal)
        call s%assistant_columns%setObjectName('studioAssistantColumns')
        call s%assistant_columns%setChildrenCollapsible(.false.)
        call s%assistant_columns%setHandleWidth(workspace_gap)
        s%assistant_space=QWidget()
        call s%assistant_space%setObjectName(assistant_space_name)
        call s%assistant_space%setMinimumWidth(assistant_min_width)
        call s%assistant_columns%addWidget(s%columns)
        call s%assistant_columns%addWidget(s%assistant_space)
        call s%assistant_columns%setStretchFactor(0,1)
        call s%assistant_space%hide()
        call s%assistant_columns%splitterMoved%connect(on_assistant_divider_moved,s)
        assistant_columns_events=QtfEventFilter(s%assistant_columns)
        call assistant_columns_events%resizeEvent%connect(on_assistant_columns_resize,s)
        assistant_space_events=QtfEventFilter(s%assistant_space)
        call assistant_space_events%resizeEvent%connect(on_assistant_space_resize,s)
        call workspace_layout%addWidget(s%assistant_columns,1)
        s%panel_resize_guide=QWidget(s%command_ribbon)
        call s%panel_resize_guide%setObjectName('studioPanelResizeGuide')
        call s%panel_resize_guide%setFixedSize(3,12)
        call s%panel_resize_guide%hide()
        call connect_panel_resize(s,s%columns,[character(24)::center_column_name,inspector_column_name])
        call connect_panel_resize(s,s%assistant_columns,[assistant_space_name])
        s%pages=QStackedWidget()
        call s%pages%addWidget(s%workspace)
        call build_project_home(s,connect_button)
        call root_layout%addWidget(s%pages,1)
        call s%window%setCentralWidget(root)
        call s%status%build(s%window,root)
        call s%status%context%clicked%connect(on_status_selection,s)
        call s%status%unit_button%clicked%connect(on_status_units,s)
        call s%status%assistant%clicked%connect(on_status_assistant,s)
        call set_status_hint(s%status%context,'Open the Inspector for selection details')
        call set_status_hint(s%status%unit_button,'Change project units')
        call set_status_hint(s%status%assistant,'Open the Assistant')
        call s%widgets%ribbon_tabs%currentChanged%connect(on_ribbon_tab_changed,s)
        call show_ready(s)
        call refresh_inspector(s)
        call refresh_knot_table(s)
        call connect_canvases(s,control_move=on_control_move,control_selected=on_control_selection, &
            control_open=on_control_open,selection_changed=apply_workspace_selection,transform_committed=on_transform_drag, &
            viewport_context=show_viewport_context,transform_exited=on_transform_exit,region_selected=apply_workspace_selection_batch)
        delete_selection=QShortcut('Del',s%surface)
        call delete_selection%setContext(Qt_WidgetShortcut)
        call delete_selection%setAutoRepeat(.false.)
        call delete_selection%activated%connect(on_viewport_delete,s)
        call connect_knot_interactions(s,on_knot_interaction,domain_knot_fraction,domain_knot_point)
        call select_entity(s,s%entity_kind)
    end procedure build_ui

    module procedure refresh_workspace_status
        character(:),allocatable :: guidance
        integer :: mesh_id

        if(.not.s%dragging)call s%status%setActivity('')
        if(s%home_active)then
            call s%status%setContext('','','Start a project or open a sample')
            return
        end if
        guidance='Ready'
        if(s%focused_surface.or.s%canvas_only)guidance=camera_guidance
        if(.not.s%canvas_only)then
            select case(s%widgets%ribbon_tabs%currentIndex())
            case(patch_ribbon_tab)
                if(.not.s%widgets%edit_controls_action%isEnabled())guidance='Select one patch to edit its controls or parameters'
            case(multipatch_ribbon_tab)
                if(.not.s%widgets%multipatch_actions(1)%isEnabled())then
                    guidance='Select a multipatch or group patches of the same type'
                    if(.not.s%widgets%ribbon%assembly_connect_action%isEnabled())guidance='Add two patches of the same type to form a multipatch'
                end if
            case(mesh_ribbon_tab)
                guidance='Select one patch, valid multipatch, or mesh with a source'
                if(s%widgets%ribbon%mesh_generate_action%isEnabled())then
                    guidance='Ready to mesh'
                    mesh_id=mesh_context_record(s)
                    if(mesh_id>0)then
                        guidance='Mesh needs update · Remesh to rebuild'
                        if(studio_mesh_valid(s,mesh_id))then
                            guidance='Remesh replaces the existing mesh'
                            if(s%selected_mesh==0)guidance='Mesh · '//mesh_counts_text(s%document%meshes(mesh_id))
                        end if
                    end if
                end if
            case default
                continue
            end select
        end if
        if(s%isolation%count()>0)guidance='Scene isolated · Exit isolation to restore visibility'
        if(s%mesh_only)guidance='Meshes only · Select a visible mesh'
        if(s%mesh_pick)guidance='Select a visible mesh cell'
        select case(s%area_selection_tool)
        case(area_box)
            guidance='Box select · Drag left to right to enclose, right to left to cross'
        case(area_lasso)
            guidance='Lasso select · Draw around or across the objects'
        case default
            continue
        end select
        if(s%transform_mode>0)guidance=transform_gizmo_hint(s)
        call s%status%setContext(selection_status_text(s),project_unit_status(s),guidance, &
            compact=selection_status_text(s,compact=.true.))
    end procedure refresh_workspace_status

    subroutine on_status_selection(data)
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            call invoke_command(s,show_inspector_command)
        end select
    end subroutine on_status_selection

    subroutine on_status_units(data)
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            call invoke_command(s,units_command)
        end select
    end subroutine on_status_units

    subroutine on_status_assistant(data)
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            call invoke_command(s,show_assistant_command)
        end select
    end subroutine on_status_assistant

    subroutine on_ribbon_tab_changed(index,data)
        integer,intent(in) :: index
        class(*),intent(inout),optional :: data

        if(index<0)return
        select type(s=>data)
        type is(studio_t)
            call refresh_workspace_status(s)
        end select
    end subroutine on_ribbon_tab_changed

    subroutine on_status_hint(text,data)
        character(*),intent(in) :: text
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            call s%status%setHint(text)
        end select
    end subroutine on_status_hint

    module procedure sync_window_controls
        type(QToolButton) :: quick_inspector,quick_assistant
        integer :: button_style
        logical :: fullscreen, restore_icon

        button_style=merge(Qt_ToolButtonTextBesideIcon,Qt_ToolButtonIconOnly,s%window%width()>=quick_action_label_width)
        call s%command_ribbon%findChild('studioQuickInspectorButton',quick_inspector)
        call s%command_ribbon%findChild('studioAssistantButton',quick_assistant)
        call quick_inspector%setToolButtonStyle(button_style)
        call quick_assistant%setToolButtonStyle(button_style)
        fullscreen=s%window%boolProperty('fullScreen')
        if(.not.fullscreen.and..not.s%window%isMinimized())s%widgets%restore_maximized=s%window%isMaximized()
        call s%widgets%fullscreen_action%setChecked(fullscreen)
        restore_icon=fullscreen.or.s%window%isMaximized()
        if(restore_icon.eqv.s%widgets%window_restore_icon)return
        s%widgets%window_restore_icon=restore_icon
        call s%widgets%maximize_button%setIcon(line_icon(merge(icon_restore,icon_maximize,restore_icon)))
    end procedure sync_window_controls

    subroutine on_window_shown(data)
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            call sync_window_controls(s)
        end select
    end subroutine on_window_shown

    subroutine on_window_resize(event,data)
        type(QResizeEvent),intent(in) :: event
        class(*),intent(inout),optional :: data

        call on_window_shown(data)
    end subroutine on_window_resize

    recursive subroutine on_workspace_resize(event,data)
        type(QResizeEvent),intent(in)::event
        class(*),intent(inout),optional::data
        select type(s=>data)
        type is(studio_t)
            call sync_assistant_layout(s)
        end select
    end subroutine on_workspace_resize

    subroutine on_viewport_aperture_resize(event,data)
        type(QResizeEvent),intent(in)::event
        class(*),intent(inout),optional::data

        select type(s=>data)
        type is(studio_t)
            call s%workspace%update()
        end select
    end subroutine on_viewport_aperture_resize

    subroutine on_inspector_visibility(data)
        class(*), intent(inout), optional :: data
        logical :: visible

        select type(s=>data)
        type is(studio_t)
            visible=s%panels(8)%isVisible()
            call s%panel_actions(8)%setChecked(visible)
        end select
    end subroutine on_inspector_visibility

    subroutine sync_parameter_panels(data)
        class(*), intent(inout), optional :: data
        logical :: visible

        select type(s=>data)
        type is(studio_t)
            visible=s%panels(5)%isVisible().and.s%panels(6)%isVisible()
            call s%widgets%ribbon%parameters_ribbon%setChecked(visible)
        end select
    end subroutine sync_parameter_panels

    module procedure assistant_shares_inspector_space
        shares_space=studio_assistant_is_docked(s).and.s%window%width()<wide_workspace_with_chat
    end procedure assistant_shares_inspector_space

    subroutine connect_panel_resize(s,splitter,child_names)
        type(studio_t),intent(inout),target :: s
        type(QSplitter),intent(in) :: splitter
        character(*),intent(in) :: child_names(:)
        type(QWidget) :: handles(size(child_names))
        type(QtfEventFilter) :: events
        integer :: i

        do i=1,size(handles)
            call splitter%findChild('qt_splithandle_'//trim(child_names(i)),handles(i))
            if(handles(i)%isNull())return
            if(.not.handles(i)%inherits('QSplitterHandle'))return
        end do
        call splitter%setProperty('opaqueResize',.false.)
        do i=1,size(handles)
            events=QtfEventFilter(handles(i))
            call events%mouseMoveEvent%connect(on_panel_resize_move,s)
            call events%mouseReleaseEvent%connect(on_panel_resize_release,s)
            call events%hideEvent%connect(hide_panel_resize_guide,s)
        end do
    end subroutine connect_panel_resize

    subroutine on_panel_resize_move(event,data)
        type(QMouseEvent),intent(in) :: event
        class(*),intent(inout),optional :: data
        type(QPointF) :: global_position
        type(QPoint) :: position
        integer :: x,y

        select type(s=>data)
        type is(studio_t)
            if(iand(event%buttons(),Qt_LeftButton)==0)then
                call hide_panel_resize_guide(s)
                return
            end if
            global_position=event%globalPosition()
            position=s%command_ribbon%mapFromGlobal(QPoint(nint(global_position%x),nint(global_position%y)))
            x=max(0,min(position%x-1,s%command_ribbon%width()-s%panel_resize_guide%width()))
            y=s%command_ribbon%height()-s%panel_resize_guide%height()
            call s%panel_resize_guide%move(x,y)
            if(.not.s%panel_resize_guide%isVisible())then
                call s%panel_resize_guide%raise()
                call s%panel_resize_guide%show()
            end if
            call event%accept()
        end select
    end subroutine on_panel_resize_move

    subroutine on_panel_resize_release(event,data)
        type(QMouseEvent),intent(in) :: event
        class(*),intent(inout),optional :: data

        if(event%button()==Qt_LeftButton)call hide_panel_resize_guide(data)
    end subroutine on_panel_resize_release

    subroutine hide_panel_resize_guide(data)
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            if(.not.s%panel_resize_guide%isNull())call s%panel_resize_guide%hide()
        end select
    end subroutine hide_panel_resize_guide

    subroutine on_scene_divider_moved(position,index,data)
        integer,intent(in) :: position,index
        class(*),intent(inout),optional :: data
        integer :: width

        select type(s=>data)
        type is(studio_t)
            if(index/=1.or..not.s%focused_surface.or..not.s%left_column%isVisible())return
            width=s%left_column%width()
            s%normal_columns(2)=max(1,s%normal_columns(2)+s%normal_columns(1)-width)
            s%normal_columns(1)=width
        end select
    end subroutine on_scene_divider_moved

    subroutine on_assistant_divider_moved(position,index,data)
        integer,intent(in) :: position,index
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            if(index/=1.or..not.s%assistant_space%isVisible())return
            s%assistant_width=s%assistant_space%width()
            call position_assistant_overlay(s)
        end select
    end subroutine on_assistant_divider_moved

    subroutine on_assistant_columns_resize(event,data)
        type(QResizeEvent),intent(in) :: event
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            call size_assistant_column(s)
        end select
    end subroutine on_assistant_columns_resize

    subroutine on_assistant_space_resize(event,data)
        type(QResizeEvent),intent(in) :: event
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            call position_assistant_overlay(s)
        end select
    end subroutine on_assistant_space_resize

    subroutine size_assistant_column(s)
        type(studio_t),intent(inout) :: s
        integer :: available

        if(.not.s%assistant_space%isVisible())return
        available=s%assistant_columns%width()-workspace_gap
        call s%assistant_columns%setSizes([max(1,available-s%assistant_width),s%assistant_width])
    end subroutine size_assistant_column

    subroutine position_assistant_overlay(s)
        type(studio_t),intent(inout) :: s
        integer :: top,width,card_height

        if(s%assistant_overlay%isNull().or..not.s%assistant_space%isVisible())return
        width=s%assistant_space%width()
        top=s%command_ribbon%y()+s%command_ribbon%height()+workspace_gap
        card_height=max(0,s%surface%height()-top-workspace_bottom_inset)
        call s%assistant_overlay%setGeometry(QRect(s%surface%width()-width-workspace_side_inset,top,width,card_height))
        call s%assistant_overlay%raise()
    end subroutine position_assistant_overlay

    module procedure sync_assistant_layout
        logical :: docked

        if(s%widgets%updating_assistant_layout.or.s%assistant_space%isNull())return
        s%widgets%updating_assistant_layout=.true.
        docked=studio_assistant_is_docked(s)
        if(.not.docked)call s%assistant_space%hide()
        if(assistant_shares_inspector_space(s))then
            s%widgets%inspector_covered_by_chat=.true.
            call s%panels(8)%hide()
        else if(s%widgets%inspector_covered_by_chat)then
            s%widgets%inspector_covered_by_chat=.false.
            call s%panels(8)%setVisible(s%widgets%inspector_requested.and..not.s%focused_surface)
        end if
        if(.not.s%assistant_overlay%isNull())then
            if(docked)then
                if(.not.s%assistant_space%isVisible())then
                    call s%assistant_space%show()
                    call size_assistant_column(s)
                end if
                call position_assistant_overlay(s)
            else
                ! Keep the QWidget owner visible so a detached Qt dock stays on screen.
                call s%assistant_overlay%setGeometry(QRect(-s%assistant_width,-s%surface%height(),0,0))
            end if
        end if
        call s%widgets%ribbon%assistant_toggle_action%setChecked(studio_assistant_is_visible(s))
        s%widgets%updating_assistant_layout=.false.
    end procedure sync_assistant_layout

    module procedure show_ready
        call refresh_workspace_status(s)
    end procedure show_ready

    subroutine title_button(button,layout,s,command)
        type(QPushButton),intent(inout)::button
        type(QHBoxLayout),intent(inout)::layout
        type(studio_t),intent(inout),target::s
        integer,intent(in)::command
        button=QPushButton('')
        call button%setObjectName('titleControl'); call button%setFixedSize(37,27)
        call button%setFocusPolicy(Qt_ClickFocus)
        call button%setIconSize(15,15)
        select case(command)
        case(close_window_command)
            call button%setIcon(line_icon(icon_close))
            call button%setProperty('windowClose',.true.)
            call set_status_hint(button,'Close ForCAD Studio'); call button%setAccessibleName('Close window')
        case(minimize_window_command)
            call button%setIcon(line_icon(icon_minimize))
            call set_status_hint(button,'Minimize ForCAD Studio'); call button%setAccessibleName('Minimize window')
        case(maximize_window_command)
            call button%setIcon(line_icon(icon_maximize))
            call set_status_hint(button,'Maximize or restore window')
            call button%setAccessibleName('Maximize or restore window')
        end select
        call connect_button(button,s,command); call layout%addWidget(button)
    end subroutine title_button

    subroutine create_panel(title,panel,layout,s,index)
        character(len=*),intent(in)::title
        type(QWidget),intent(out)::panel
        type(QVBoxLayout),intent(out)::layout
        type(studio_t),intent(inout),target::s
        integer,intent(in)::index
        type(QWidget)::frame,header
        type(QHBoxLayout)::header_layout
        type(QLabel)::label

        call build_opaque_panel(panel,frame,studio_panel_fill,'panel')
        s%panels(index)=panel
        call panel%setAccessibleName(title)
        if(any(index==[1,5,6]))then
            call panel%showEvent%connect(on_panel_shown,s%panel_actions(index))
            call panel%hideEvent%connect(on_panel_hidden,s%panel_actions(index))
        end if
        layout=QVBoxLayout(frame)
        call layout%setContentsMargins(1,1,1,1); call layout%setSpacing(0)
        header=QWidget(); call header%setObjectName('panelHeader'); call header%setFixedHeight(32)
        header_layout=QHBoxLayout(header)
        call header_layout%setContentsMargins(panel_inset,0,6,0); call header_layout%setSpacing(4)
        label=QLabel(title); call label%setObjectName('panelTitle')
        s%panel_titles(index)=label
        if(index==6)s%basis_title=label
        call header_layout%addWidget(label,1)
        if(index==1)then
            call panel%setAccessibleName('Scene browser')
            call build_project_actions(header_layout,header,s)
        end if
        select case(index)
        case(5)
            call build_domain_camera_tools(header_layout,s)
        case(6)
            call build_basis_display_menu(header_layout,header,s)
        end select
        s%close_buttons(index)=QPushButton('',header)
        call s%close_buttons(index)%setObjectName('panelClose')
        call s%close_buttons(index)%setFixedSize(27,27)
        call s%close_buttons(index)%setIcon(line_icon(icon_close))
        call s%close_buttons(index)%setIconSize(15,15)
        if(index==5.or.index==6)then
            call set_status_hint(s%close_buttons(index),'Hide '//title)
        else
            call set_status_hint(s%close_buttons(index),studio_tooltip('Hide '//title,'Reopen it from the View ribbon.'))
        end if
        call s%close_buttons(index)%setAccessibleName('Close '//title)
        call connect_button(s%close_buttons(index),s,close_panel_command_base+index)
        call header_layout%addWidget(s%close_buttons(index))
        call layout%addWidget(header)
    end subroutine create_panel

    subroutine on_panel_shown(data)
        class(*), intent(inout), optional :: data

        select type(action=>data)
        type is(QAction)
            call action%setChecked(.true.)
        end select
    end subroutine on_panel_shown

    subroutine on_panel_hidden(data)
        class(*), intent(inout), optional :: data

        select type(action=>data)
        type is(QAction)
            call action%setChecked(.false.)
        end select
    end subroutine on_panel_hidden

    subroutine build_project_actions(row,parent,s)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::parent
        type(studio_t),intent(inout),target::s
        type(QToolButton)::button
        type(QMenu)::menu
        type(QAction)::find,reveal

        button=QToolButton(parent)
        call button%setObjectName('studioProjectActionsButton')
        call button%setIcon(line_icon(icon_more))
        call button%setIconSize(17,17)
        call button%setFixedSize(25,25)
        call set_status_hint(button,studio_tooltip('Scene commands','Find, isolate, rename, or organize objects in the scene.'))
        call button%setAccessibleName('Scene commands')
        menu=QMenu('Scene',button)
        call register_status_menu(menu)
        call menu%setObjectName('studioProjectSelectionMenu')
        find=QAction('Find in Scene',menu)
        call find%setShortcut('Ctrl+F')
        call find%setIcon(line_icon(27))
        call s%surface%addAction(find)
        call menu%addAction(find)
        call find%triggered%connect(focus_scene_search,s)
        call menu_action(menu,'Reveal in Scene','',s,scene_reveal_command,2,reveal)
        call reveal%setObjectName('studioSceneRevealAction')
        call menu%addSeparator()
        call menu_action(menu,'Fit selection','',s,fit_selection_command,icon_fit_selection,s%widgets%project_fit_action)
        call menu_action(menu,'Show or hide selection','',s,toggle_patch_visibility_command,31,s%widgets%scene_selection_visibility_action)
        call s%widgets%scene_selection_visibility_action%setObjectName('studioSceneSelectionVisibilityAction')
        call menu_action(menu,'Isolate selection','',s,project_isolate_command,31)
        call menu_action(menu,'Show all objects','',s,project_show_all_command,13)
        call menu%addSeparator()
        call menu_action(menu,'Inspect parameters','',s,inspect_patch_command,33,s%widgets%project_inspect_action)
        call s%widgets%project_fit_action%setObjectName('studioProjectFitAction')
        call s%widgets%project_inspect_action%setObjectName('studioProjectInspectAction')
        call menu%addAction(s%widgets%transform_action)
        call menu%addAction(s%widgets%project_duplicate_action)
        call menu%addAction(s%widgets%project_rename_action)
        call menu%addSeparator()
        call menu_action(menu,'Expand all','',s,project_expand_command,30)
        call menu_action(menu,'Collapse all','',s,project_collapse_command,36)
        call menu%addSeparator()
        call menu%addAction(s%widgets%project_remove_action)
        call button%setMenu(menu)
        call button%setPopupMode(QToolButton_InstantPopup)
        call row%addWidget(button)
    end subroutine build_project_actions

    subroutine focus_scene_search(data)
        class(*),intent(inout),optional::data

        select type(s=>data)
        type is(studio_t)
            if(.not.s%workspace%isVisible())return
            if(s%canvas_only)call invoke_command(s,canvas_command)
            if(.not.s%panels(1)%isVisible())call invoke_command(s,toggle_panel_command_base+1)
            call s%widgets%project_browser%search%setFocus()
        end select
    end subroutine focus_scene_search

    subroutine viewer_menu(row,parent,button,menu,name,label,icon_number,viewport,split)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::parent
        type(QToolButton),intent(out)::button
        type(QMenu),intent(out)::menu
        character(*),intent(in)::name,label
        integer,intent(in)::icon_number
        logical,optional,intent(in)::viewport,split

        button=QToolButton(parent)
        call button%setObjectName(name//'Button')
        call button%setProperty('viewerSettings',.true.)
        if(present(viewport))then
            call button%setProperty('viewportMenu',viewport)
            call button%setFocusPolicy(Qt_StrongFocus)
            call button%setToolButtonStyle(Qt_ToolButtonIconOnly)
        end if
        if(present(split))call button%setProperty('viewportSplit',split)
        call button%setFixedSize(viewer_button_width,viewer_button_height)
        if(present(viewport))call button%setFixedWidth(viewport_button_width)
        call button%setIcon(line_icon(icon_number))
        call button%setIconSize(viewer_icon_size,viewer_icon_size)
        call set_status_hint(button,label)
        call button%setAccessibleName(label)
        menu=QMenu(label,button)
        if(present(viewport))call menu%setProperty('viewportMenu',viewport)
        call register_status_menu(menu)
        call menu%setObjectName(name//'Menu')
        call button%setMenu(menu)
        call button%setPopupMode(QToolButton_InstantPopup)
        if(present(split))then
            if(split)call button%setPopupMode(QToolButton_MenuButtonPopup)
        end if
        call row%addWidget(button)
    end subroutine viewer_menu

    subroutine viewer_menu_heading(menu,title,separator)
        type(QMenu),intent(in) :: menu
        character(*),intent(in) :: title
        logical,optional,intent(in) :: separator
        type(QAction) :: heading

        if(present(separator))then
            if(separator)call menu%addSeparator()
        end if
        heading=QAction(title,menu)
        call heading%setEnabled(.false.)
        call menu%addAction(heading)
    end subroutine viewer_menu_heading

    subroutine build_viewport_snap_menu(row,parent,s)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::parent
        type(studio_t),intent(inout),target :: s
        type(QMenu)::menu

        call viewer_menu(row,parent,s%widgets%viewport_snap,menu,'studioViewportSnap','Snap',39,viewport=.true.,split=.true.)
        call make_command(parent,'Snap','',s,snap_command,39,s%widgets%viewport_snap_action)
        call s%widgets%viewport_snap_action%setCheckable(.true.)
        call s%widgets%viewport_snap%setDefaultAction(s%widgets%viewport_snap_action)
        call menu_action(menu,'Snap off','',s,snap_off_command,0)
        call menu_action(menu,'Snap to vertices','',s,snap_vertex_command,icon_bounding_vertex)
        call menu_action(menu,'Snap to closest point','',s,snap_closest_command,93)
    end subroutine build_viewport_snap_menu

    subroutine build_selection_filter(row,parent,s)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::parent
        type(studio_t),intent(inout),target::s
        type(QMenu)::menu

        call viewer_menu(row,parent,s%widgets%viewport_pick,menu,'studioViewportPick','Selection filter',18,viewport=.true.)
        call viewer_menu_heading(menu,'Method')
        call menu_action(menu,'Click','',s,select_tool_command,18,s%widgets%viewport_click_action)
        call s%widgets%viewport_click_action%setCheckable(.true.)
        call menu%addAction(s%widgets%area_selection_actions(1))
        call menu%addAction(s%widgets%area_selection_actions(2))
        call viewer_menu_heading(menu,'Pick',separator=.true.)
        call menu_action(menu,'Auto','',s,pick_mode_base,18)
        call menu_action(menu,'Faces','',s,pick_mode_base+1,icon_bounding_face)
        call menu_action(menu,'Edges','',s,pick_mode_base+2,icon_bounding_edge)
        call menu_action(menu,'Vertices','',s,pick_mode_base+3,icon_bounding_vertex)
        call menu%addAction(s%widgets%pick_actions(pick_multipatch+1))
        call menu%addAction(s%widgets%pick_actions(pick_control_points+1))
        call menu%addAction(s%widgets%ribbon%mesh_pick_action)
        call menu%addSeparator()
        call menu_action(menu,'Clear selection','',s,clear_selection_command,icon_close)
    end subroutine build_selection_filter

    module procedure sync_selection_filter
        integer,parameter :: pick_icons(0:5)=[18,icon_bounding_face,icon_bounding_edge,icon_bounding_vertex,22,87]
        integer :: i,glyph
        character(:),allocatable :: target,method,hint
        type(QIcon) :: icon
        logical :: picks_mesh,clicking

        if(s%widgets%viewport_pick%isNull())return
        clicking=s%surface_tool==1.and.s%area_selection_tool==area_click.and.s%transform_mode==0
        call s%toolbar_actions(6)%setChecked(clicking)
        call s%widgets%viewport_click_action%setChecked(clicking)
        call s%toolbar_actions(7)%setChecked(s%surface_tool==2)
        do i=1,size(s%widgets%area_selection_actions)
            call s%widgets%area_selection_actions(i)%setChecked(s%area_selection_tool==i)
        end do
        picks_mesh=s%mesh_pick.or.s%mesh_only
        do i=1,size(s%widgets%pick_actions)
            call s%widgets%pick_actions(i)%setEnabled(.not.s%mesh_only)
            call s%widgets%pick_actions(i)%setChecked(.not.picks_mesh.and.i-1==s%feature_pick_mode)
        end do
        call s%widgets%ribbon%mesh_pick_action%setChecked(s%mesh_pick)
        target=s%widgets%pick_actions(s%feature_pick_mode+1)%text()
        if(s%mesh_only)target='Meshes'
        if(s%mesh_pick.and.s%area_selection_tool==area_click)target='Mesh cells'
        if(picks_mesh.and.s%area_selection_tool/=area_click)target='Meshes'
        method='Click'
        if(s%area_selection_tool==area_box)method='Box'
        if(s%area_selection_tool==area_lasso)method='Lasso'
        if(s%surface_tool==2)method='Orbit'
        hint='Selection · '//method//' · '//target
        if(method=='Orbit')hint='Orbit tool active · Choose a selection method to resume picking'
        glyph=pick_icons(s%feature_pick_mode)
        if(picks_mesh)glyph=icon_mesh
        if(s%mesh_pick.and.s%area_selection_tool==area_click)glyph=icon_mesh_pick
        if(s%area_selection_tool==area_box)glyph=icon_box_selection
        if(s%area_selection_tool==area_lasso)glyph=icon_lasso_selection
        if(s%surface_tool==2)glyph=18
        if(s%widgets%viewport_pick%intProperty('studioSelectionGlyph')/=glyph)then
            icon=line_icon(glyph)
            call s%widgets%viewport_pick%setIcon(icon)
            call icon%delete()
            call s%widgets%viewport_pick%setProperty('studioSelectionGlyph',glyph)
        end if
        call set_status_hint(s%widgets%viewport_pick,hint)
        call s%widgets%viewport_pick%setAccessibleName(hint)
    end procedure sync_selection_filter

    subroutine build_viewport_display_menu(row,parent,s)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::parent
        type(studio_t),intent(inout),target::s
        type(QMenu)::menu
        type(QActionGroup) :: styles
        integer::i

        call viewer_menu(row,parent,s%viewport_display,menu,'studioViewportDisplay','3D display',33,viewport=.true.)
        call viewer_menu_heading(menu,'Shading')
        s%widgets%displayed_viewport_style=0
        styles=QActionGroup(menu)
        call styles%setExclusive(.true.)
        do i=1,size(s%widgets%viewport_style_actions)
            call styles%addAction(s%widgets%viewport_style_actions(i))
            call menu%addAction(s%widgets%viewport_style_actions(i))
        end do
        call sync_viewport_style(s)
        call viewer_menu_heading(menu,'Presets',separator=.true.)
        do i=1,size(s%viewport_presets)
            call menu%addAction(s%viewport_presets(i))
        end do
        call viewer_menu_heading(menu,'Guides',separator=.true.)
        do i=1,size(s%viewport_layers)
            call menu%addAction(s%viewport_layers(i))
        end do
        call menu%addAction(s%widgets%ribbon%physical_knot_action)
        call menu%addAction(s%viewport_grid)
    end subroutine build_viewport_display_menu

    subroutine sync_viewport_style(s)
        type(studio_t), intent(inout) :: s
        character(:), allocatable :: hint
        type(QIcon) :: icon
        integer :: i,style

        if(s%viewport_display%isNull())return
        style=s%viewport_style
        if(style==s%widgets%displayed_viewport_style)return
        s%widgets%displayed_viewport_style=style
        do i=1,size(s%widgets%viewport_style_actions)
            call s%widgets%viewport_style_actions(i)%setChecked(i==style)
        end do
        hint='3D display · '//trim(viewport_style_names(style))
        icon=line_icon(viewport_style_icons(style))
        call s%viewport_display%setIcon(icon)
        call icon%delete()
        call s%viewport_display%setAccessibleName(hint)
        call set_status_hint(s%viewport_display,hint)
    end subroutine sync_viewport_style

    subroutine build_domain_camera_tools(row,s)
        type(QHBoxLayout),intent(inout)::row
        type(studio_t),intent(inout),target::s
        type(QPushButton)::button

        call domain_button(row,button,s,domain_zoom_out_command,6,'DomainZoomOut','Zoom out')
        call domain_button(row,button,s,domain_zoom_in_command,7,'DomainZoomIn','Zoom in')
        call domain_button(row,button,s,domain_reset_command,19,'DomainFit','Fit parameter space')
    end subroutine build_domain_camera_tools

    subroutine build_domain_overlays(row,parent,s)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::parent
        type(studio_t),intent(inout),target::s
        type(QToolButton)::button
        integer::i
        integer,parameter::icons(5)=[12,63,24,33,22]
        character(len=22),parameter::options(5)=[character(len=22)::'Parameter grid','Knot divisions', &
            'Elements','Basis support','Selected element']

        do i=1,5
            s%domain_overlays(i)=QAction(trim(options(i)),parent)
            call s%domain_overlays(i)%setIcon(line_icon(icons(i)))
            call s%domain_overlays(i)%setCheckable(.true.)
            call s%domain_overlays(i)%setChecked(s%domain_layers(i))
            call s%domain_overlays(i)%toggled%connect(on_domain,s)
            button=QToolButton(parent)
            call button%setDefaultAction(s%domain_overlays(i))
            call watch_status_hint(button)
            call button%setIconSize(17,17)
            call button%setFixedSize(28,27)
            call button%setAccessibleName('Toggle '//trim(options(i)))
            if(i==1)then
                call button%setObjectName('studioDomainDisplay')
                s%domain_display=button
            end if
            call row%addWidget(button)
        end do
    end subroutine build_domain_overlays

    subroutine build_basis_display_menu(row,parent,s)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::parent
        type(studio_t),intent(inout),target::s
        type(QMenu)::menu
        integer::i

        call viewer_menu(row,parent,s%basis_display,menu,'studioBasisDisplay','Basis Plot display settings',33)
        call menu_action(menu,'Series legend','',s,basis_legend_command,0,s%basis_options(1))
        call menu_action(menu,'Probe guide','',s,basis_probe_command,0,s%basis_options(2))
        do i=1,2
            call s%basis_options(i)%setCheckable(.true.)
        end do
    end subroutine build_basis_display_menu

    subroutine domain_button(row,button,s,command,glyph,name,caption)
        type(QHBoxLayout),intent(inout)::row
        type(QPushButton),intent(out)::button
        type(studio_t),intent(inout),target::s
        integer,intent(in)::command,glyph
        character(len=*),intent(in)::name,caption

        button=QPushButton('')
        call button%setFixedSize(28,27)
        call button%setProperty('viewerTool',.true.)
        call button%setIcon(line_icon(glyph))
        call button%setIconSize(viewer_icon_size,viewer_icon_size)
        call button%setObjectName('studio'//name)
        call set_status_hint(button,caption)
        call button%setAccessibleName(caption)
        call connect_button(button,s,command)
        call row%addWidget(button)
    end subroutine domain_button

    subroutine build_viewport_controls(s)
        type(studio_t),intent(inout),target::s
        type(QWidget)::frame
        type(QVBoxLayout)::controls_layout,stack
        type(QHBoxLayout)::tools_row

        controls_layout=QVBoxLayout(s%viewport_aperture)
        call controls_layout%setContentsMargins(0,0,0,0)
        call build_opaque_panel(s%viewport_controls,frame,studio_viewport_controls_fill, &
            'studioViewportControlsFrame',s%viewport_aperture)
        call s%viewport_controls%setObjectName('viewportControls')
        stack=QVBoxLayout(frame)
        call stack%setContentsMargins(4,4,4,4)
        call stack%setSpacing(6)
        tools_row=QHBoxLayout()
        call tools_row%addStretch()
        call tools_row%setSpacing(2)
        call build_selection_filter(tools_row,frame,s)
        call build_viewport_snap_menu(tools_row,frame,s)
        call add_viewer_divider(tools_row,frame)
        call build_viewport_fit_menu(tools_row,frame,s)
        call build_viewport_orientation_menu(tools_row,frame,s)
        call add_viewer_divider(tools_row,frame)
        call build_viewport_display_menu(tools_row,frame,s)
        call build_viewport_visibility_menu(tools_row,frame,s)
        s%widgets%leave_canvas_button=QPushButton('',frame)
        call s%widgets%leave_canvas_button%setObjectName('studioExitCanvasOnly')
        call s%widgets%leave_canvas_button%setIcon(line_icon(icon_restore))
        call s%widgets%leave_canvas_button%setIconSize(15,15)
        call s%widgets%leave_canvas_button%setFixedSize(viewport_button_width,viewer_button_height)
        call s%widgets%leave_canvas_button%setAccessibleName('Exit Canvas Only')
        call set_status_hint(s%widgets%leave_canvas_button,'Exit Canvas Only · Ctrl+Shift+F')
        call connect_button(s%widgets%leave_canvas_button,s,canvas_command)
        call tools_row%addWidget(s%widgets%leave_canvas_button)
        call s%widgets%leave_canvas_button%hide()
        call stack%addLayout(tools_row)
        call build_canvas_mesh_controls(stack,frame,s)
        call build_transform_strip(stack,frame,s)
        call controls_layout%addWidget(s%viewport_controls,alignment=Qt_AlignRight)
        call controls_layout%addStretch(1)
        call sync_selection_filter(s)
    end subroutine build_viewport_controls

    subroutine build_viewport_fit_menu(row,parent,s)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::parent
        type(studio_t),intent(inout),target::s
        type(QMenu)::menu
        type(QAction)::fit

        call viewer_menu(row,parent,s%widgets%viewport_fit,menu,'studioViewportFit','Fit all',19,viewport=.true.,split=.true.)
        call menu_action(menu,'Fit scene','',s,fit_scene_command,19,fit)
        call set_status_hint(fit,'Fit all · F')
        call s%widgets%viewport_fit%setDefaultAction(fit)
        call set_status_hint(s%widgets%viewport_fit,'Fit all · F')
        call menu%addAction(s%widgets%ribbon%view_fit_action)
        call menu_action(menu,'Fit around panels','',s,fit_around_panels_command,icon_fit_panels)
    end subroutine build_viewport_fit_menu

    subroutine build_viewport_orientation_menu(row,parent,s)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::parent
        type(studio_t),intent(inout),target::s
        type(QMenu)::menu
        type(QActionGroup)::views
        integer::i

        call viewer_menu(row,parent,s%widgets%viewport_orientation,menu,'studioViewportOrientation', &
            'Orientation',25,viewport=.true.)
        views=QActionGroup(menu)
        call views%setExclusive(.true.)
        s%widgets%free_orbit_action=QAction('Free orbit',menu)
        call s%widgets%free_orbit_action%setObjectName('studioViewportFreeOrbit')
        call s%widgets%free_orbit_action%setCheckable(.true.)
        call s%widgets%free_orbit_action%setEnabled(.false.)
        call views%addAction(s%widgets%free_orbit_action)
        call menu%addAction(s%widgets%free_orbit_action)
        call menu%addSeparator()
        do i=1,size(s%widgets%camera_view_actions)
            call views%addAction(s%widgets%camera_view_actions(i))
            call menu%addAction(s%widgets%camera_view_actions(i))
        end do
        call sync_viewport_orientation(s)
    end subroutine build_viewport_orientation_menu

    subroutine build_viewport_visibility_menu(row,parent,s)
        type(QHBoxLayout),intent(inout)::row
        type(QWidget),intent(in)::parent
        type(studio_t),intent(inout),target::s
        type(QMenu)::menu
        type(QAction)::show_all

        call viewer_menu(row,parent,s%widgets%viewport_visibility,menu,'studioViewportVisibility', &
            'Visibility',31,viewport=.true.,split=.true.)
        call s%widgets%viewport_visibility%setCheckable(.true.)
        call s%widgets%viewport_visibility%clicked%connect(on_viewport_visibility,s)
        call menu_action(menu,'Exit isolation','',s,scene_exit_isolation_command,icon_restore,s%widgets%viewport_exit_isolation_action)
        call s%widgets%viewport_exit_isolation_action%setObjectName('studioViewportExitIsolation')
        call s%widgets%viewport_exit_isolation_action%setProperty('visible',.false.)
        call menu_action(menu,'Hide selected','',s,project_hide_selection_command,32,s%widgets%viewport_hide_action)
        call s%widgets%viewport_hide_action%setObjectName('studioViewportHideSelected')
        call menu_action(menu,'Isolate selection','',s,project_isolate_command,13,s%widgets%viewport_isolate_action)
        call s%widgets%viewport_isolate_action%setObjectName('studioViewportIsolateSelection')
        call menu_action(menu,'Show all','',s,project_show_all_command,31,show_all)
        call show_all%setObjectName('studioViewportShowAll')
        call menu%addSeparator()
        call menu%addAction(s%widgets%ribbon%mesh_only_action)
    end subroutine build_viewport_visibility_menu

    subroutine on_viewport_visibility(data)
        class(*),intent(inout),optional::data

        select type(s=>data)
        type is(studio_t)
            if(s%isolation%count()>0)then
                call invoke_command(s,scene_exit_isolation_command)
            else
                call s%widgets%viewport_visibility%setChecked(s%mesh_only)
                call s%widgets%viewport_visibility%showMenu()
            end if
        end select
    end subroutine on_viewport_visibility

    subroutine add_viewer_divider(row,parent)
        type(QHBoxLayout), intent(inout) :: row
        type(QWidget), intent(in) :: parent
        type(QWidget) :: divider

        divider=QWidget(parent)
        call divider%setObjectName('studioViewerDivider')
        call divider%setFixedSize(1,18)
        call row%addWidget(divider,alignment=Qt_AlignVCenter)
    end subroutine add_viewer_divider

    subroutine build_volume_mesh_controls(s)
        type(studio_t),intent(inout),target::s
        type(QVBoxLayout)::fields
        type(QHBoxLayout)::mode_row,layer_row
        type(QButtonGroup)::modes
        integer::axis,mode
        character(31),parameter :: button_names(mesh_exterior:mesh_layer)=[character(31):: &
            'studioVolumeMeshExteriorButton','studioVolumeMeshLatticeButton','studioVolumeMeshLayerButton']
        character(64),parameter :: hints(mesh_exterior:mesh_layer)=[character(64):: &
            'Exterior · Show the outer faces of the volume mesh.', &
            'Full lattice · Show the interior mesh edges as well.', &
            'Cell layer · Inspect one layer of cells along U, V, or W.']

        call s%widgets%ribbon%volume_mesh_controls%setFixedHeight(52)
        call s%widgets%ribbon%volume_mesh_controls%setStyleSheet( &
            'QComboBox,QSpinBox { min-height:20px; padding-top:1px; padding-bottom:1px; }')
        fields=QVBoxLayout(s%widgets%ribbon%volume_mesh_controls)
        call fields%setContentsMargins(1,0,1,0)
        call fields%setSpacing(2)
        modes=QButtonGroup(s%widgets%ribbon%volume_mesh_controls)
        call modes%setExclusive(.true.)
        mode_row=QHBoxLayout()
        call mode_row%setSpacing(2)
        do mode=mesh_exterior,mesh_layer
            associate(button=>s%widgets%ribbon%volume_mesh_modes(mode+1))
                button=QToolButton(s%widgets%ribbon%volume_mesh_controls)
                call button%setObjectName(trim(button_names(mode)))
                call button%setText(trim(mesh_presentation_names(mode)))
                call button%setAccessibleName(trim(mesh_presentation_names(mode))//' volume mesh view')
                call button%setToolButtonStyle(Qt_ToolButtonTextOnly)
                call button%setCheckable(.true.)
                call button%setFixedHeight(25)
                call button%setProperty('ribbonTile',.true.)
                call button%setProperty('ribbonStacked',.true.)
                call button%setFocusPolicy(Qt_StrongFocus)
                call set_status_hint(button,trim(hints(mode)))
                call modes%addButton(button,mode)
                call mode_row%addWidget(button,1)
            end associate
        end do
        call s%widgets%ribbon%volume_mesh_modes(s%mesh_view+1)%setChecked(.true.)
        call modes%idClicked%connect(volume_mesh_style_changed,s)
        call fields%addLayout(mode_row)

        s%widgets%ribbon%volume_mesh_layer_controls=QWidget(s%widgets%ribbon%volume_mesh_controls)
        layer_row=QHBoxLayout(s%widgets%ribbon%volume_mesh_layer_controls)
        call layer_row%setContentsMargins(0,0,0,0)
        call layer_row%setSpacing(4)
        s%widgets%ribbon%volume_mesh_axis=QComboBox(s%widgets%ribbon%volume_mesh_layer_controls)
        call s%widgets%ribbon%volume_mesh_axis%setObjectName('studioVolumeMeshLayerAxis')
        call s%widgets%ribbon%volume_mesh_axis%setAccessibleName('Volume mesh layer direction')
        call set_status_hint(s%widgets%ribbon%volume_mesh_axis,'Direction normal to the visible cell layer')
        call s%widgets%ribbon%volume_mesh_axis%setFixedHeight(25)
        do axis=1,3
            call s%widgets%ribbon%volume_mesh_axis%addItem(parameter_directions(axis:axis))
        end do
        call s%widgets%ribbon%volume_mesh_axis%setCurrentIndex(s%mesh_layer_axis-1)
        call layer_row%addWidget(s%widgets%ribbon%volume_mesh_axis,1)
        s%widgets%ribbon%volume_mesh_layer=QSpinBox(s%widgets%ribbon%volume_mesh_layer_controls)
        call s%widgets%ribbon%volume_mesh_layer%setObjectName('studioVolumeMeshLayer')
        call s%widgets%ribbon%volume_mesh_layer%setAccessibleName('Volume mesh cell layer number')
        call set_status_hint(s%widgets%ribbon%volume_mesh_layer,'Cell layer number')
        call s%widgets%ribbon%volume_mesh_layer%setFixedHeight(25)
        call s%widgets%ribbon%volume_mesh_layer%setRange(1,1)
        call layer_row%addWidget(s%widgets%ribbon%volume_mesh_layer,1)
        call fields%addWidget(s%widgets%ribbon%volume_mesh_layer_controls)
        call fields%addStretch()
        call s%widgets%ribbon%volume_mesh_controls%setObjectName('studioVolumeMeshControls')
        call s%widgets%ribbon%volume_mesh_axis%currentIndexChanged%connect(volume_mesh_axis_changed,s)
        call s%widgets%ribbon%volume_mesh_layer%valueChanged%connect(volume_mesh_layer_changed,s)
    end subroutine build_volume_mesh_controls

    subroutine build_canvas_mesh_controls(stack,parent,s)
        type(QWidget), intent(in) :: parent
        type(QVBoxLayout), intent(inout) :: stack
        type(studio_t), intent(inout), target :: s
        type(QVBoxLayout) :: rows
        type(QHBoxLayout) :: row, layer, section
        integer :: i
        character(*), parameter :: axes='XYZ'

        associate(hud=>s%widgets%canvas_mesh)
            hud%panel=QWidget(parent)
            call hud%panel%setObjectName('studioCanvasMeshControls')
            rows=QVBoxLayout(hud%panel)
            call rows%setContentsMargins(1,5,1,1)
            call rows%setSpacing(4)

            row=QHBoxLayout()
            call row%setSpacing(4)
            hud%presentation=QComboBox(hud%panel)
            call hud%presentation%setObjectName('studioCanvasMeshView')
            call hud%presentation%setAccessibleName('Volume mesh display')
            call set_status_hint(hud%presentation,'Volume mesh · Show the exterior, full lattice, or a cell layer')
            call hud%presentation%setFixedSize(108,26)
            do i=mesh_exterior,mesh_layer
                call hud%presentation%addItem(trim(mesh_presentation_names(i)))
            end do
            call row%addWidget(hud%presentation)

            hud%layer_controls=QWidget(hud%panel)
            layer=QHBoxLayout(hud%layer_controls)
            call layer%setContentsMargins(0,0,0,0)
            call layer%setSpacing(4)
            hud%layer_axis=QComboBox(hud%layer_controls)
            call hud%layer_axis%setObjectName('studioCanvasMeshLayerAxis')
            call hud%layer_axis%setAccessibleName('Volume mesh layer direction')
            call set_status_hint(hud%layer_axis,'Layer direction · U, V or W in the patch parameter space')
            call hud%layer_axis%setFixedSize(47,26)
            do i=1,3
                call hud%layer_axis%addItem(parameter_directions(i:i))
            end do
            call layer%addWidget(hud%layer_axis)
            hud%layer_index=QSlider(Qt_Horizontal,hud%layer_controls)
            call configure_mesh_scrubber(hud%layer_index,'studioCanvasMeshLayer','Cell layer number',1,1)
            call layer%addWidget(hud%layer_index,1)
            hud%layer_value=QLabel('',hud%layer_controls)
            call configure_mesh_value(hud%layer_value)
            call layer%addWidget(hud%layer_value)
            call row%addWidget(hud%layer_controls,1)
            call rows%addLayout(row)

            row=QHBoxLayout()
            call row%setSpacing(4)
            hud%section_axis=QComboBox(hud%panel)
            call hud%section_axis%setObjectName('studioCanvasMeshSectionAxis')
            call hud%section_axis%setAccessibleName('Mesh section plane')
            call set_status_hint(hud%section_axis,'Mesh section · Reveal interior cells along a world X, Y, or Z plane')
            call hud%section_axis%setFixedSize(108,26)
            call hud%section_axis%addItem('No cut')
            do i=1,3
                call hud%section_axis%addItem('Cut '//axes(i:i))
            end do
            call row%addWidget(hud%section_axis)
            hud%section_position=QWidget(hud%panel)
            section=QHBoxLayout(hud%section_position)
            call section%setContentsMargins(0,0,0,0)
            call section%setSpacing(4)
            hud%section_offset=QSlider(Qt_Horizontal,hud%section_position)
            call configure_mesh_scrubber(hud%section_offset,'studioCanvasMeshSectionOffset','Section plane position',0,100)
            call section%addWidget(hud%section_offset,1)
            hud%section_value=QLabel('',hud%section_position)
            call configure_mesh_value(hud%section_value)
            call section%addWidget(hud%section_value)
            call row%addWidget(hud%section_position,1)
            call rows%addLayout(row)

            call hud%presentation%currentIndexChanged%connect(volume_mesh_style_changed,s)
            call hud%layer_axis%currentIndexChanged%connect(volume_mesh_axis_changed,s)
            call hud%layer_index%valueChanged%connect(volume_mesh_layer_changed,s)
            call hud%section_axis%currentIndexChanged%connect(mesh_section_axis_changed,s)
            call hud%section_offset%valueChanged%connect(mesh_section_offset_changed,s)
            call stack%addWidget(hud%panel)
            call hud%panel%hide()
        end associate
    end subroutine build_canvas_mesh_controls

    subroutine configure_mesh_scrubber(slider,name,label,lower,upper)
        type(QSlider), intent(inout) :: slider
        character(*), intent(in) :: name,label
        integer, intent(in) :: lower,upper

        call slider%setObjectName(name)
        call slider%setProperty('meshScrubber',.true.)
        call slider%setAccessibleName(label)
        call set_status_hint(slider,label//' · Drag or use the arrow keys')
        call slider%setMinimumSize(48,24)
        call slider%setRange(lower,upper)
    end subroutine configure_mesh_scrubber

    subroutine configure_mesh_value(label)
        type(QLabel), intent(inout) :: label

        call label%setObjectName('studioCanvasMeshValue')
        call label%setAlignment(Qt_AlignRight+Qt_AlignVCenter)
        call label%setMinimumWidth(32)
    end subroutine configure_mesh_value

    subroutine build_transform_strip(stack,parent,s)
        type(QWidget), intent(in) :: parent
        type(QVBoxLayout),intent(inout)::stack
        type(studio_t),intent(inout),target::s
        type(QVBoxLayout)::layout
        type(QHBoxLayout)::row
        type(QGridLayout)::values
        type(QLabel)::axis_label
        type(QPushButton)::button
        integer::i
        character(3),parameter::axes='XYZ'

        s%widgets%transform_strip=QWidget(parent)
        call s%widgets%transform_strip%setObjectName('transformStrip')
        layout=QVBoxLayout(s%widgets%transform_strip)
        call layout%setContentsMargins(8,7,8,8)
        call layout%setSpacing(6)
        row=QHBoxLayout()
        s%widgets%transform_caption=QLabel('Move')
        call s%widgets%transform_caption%setObjectName('transformCaption')
        call row%addWidget(s%widgets%transform_caption,1)
        s%widgets%transform_snap_button=QPushButton('Snap')
        call s%widgets%transform_snap_button%setObjectName('studioGizmoSnap')
        call s%widgets%transform_snap_button%setCheckable(.true.)
        call s%widgets%transform_snap_button%setIcon(line_icon(39))
        call set_status_hint(s%widgets%transform_snap_button,'Snap drag increments · Hold Shift to snap temporarily')
        call connect_button(s%widgets%transform_snap_button,s,snap_command)
        call row%addWidget(s%widgets%transform_snap_button)
        button=QPushButton('')
        call button%setObjectName('panelClose')
        call button%setFixedSize(26,26)
        call button%setIcon(line_icon(icon_close))
        call button%setIconSize(15,15)
        call set_status_hint(button,'Exit transform · Esc')
        call button%setAccessibleName('Exit transform')
        call connect_button(button,s,transform_done_command)
        call row%addWidget(button)
        call layout%addLayout(row)
        values=QGridLayout()
        call values%setSpacing(4)
        do i=1,3
            axis_label=QLabel(axes(i:i),s%widgets%transform_strip)
            call axis_label%setObjectName('muted')
            call values%addWidget(axis_label,0,i-1)
            s%widgets%transform_inputs(i)=QDoubleSpinBox(s%widgets%transform_strip)
            call s%widgets%transform_inputs(i)%setObjectName('studioGizmoValue'//axes(i:i))
            call s%widgets%transform_inputs(i)%setAccessibleName('Exact transform value along '//axes(i:i))
            call s%widgets%transform_inputs(i)%setDecimals(6)
            call s%widgets%transform_inputs(i)%setFixedSize(78,28)
            call values%addWidget(s%widgets%transform_inputs(i),1,i-1)
        end do
        call layout%addLayout(values)
        row=QHBoxLayout()
        button=QPushButton('Exact…')
        call button%setObjectName('studioGizmoExact')
        call set_status_hint(button,'Open the full numeric editor with pivot controls')
        call connect_button(button,s,transform_exact_base)
        call row%addWidget(button)
        call row%addStretch()
        button=QPushButton('Apply values')
        call button%setObjectName('studioGizmoApply')
        call button%setIcon(line_icon(35))
        call connect_button(button,s,transform_apply_command)
        call row%addWidget(button)
        call layout%addLayout(row)
        call stack%addWidget(s%widgets%transform_strip)
        call s%widgets%transform_strip%hide()
    end subroutine build_transform_strip

    subroutine show_viewport_context(s,x,y)
        type(studio_t), intent(inout), target :: s
        real(dp), intent(in) :: x,y
        type(selection_target_t) :: target
        integer :: kind,id,element,number,mesh,part,cell,face,point

        if(s%mesh_pick.or.s%mesh_only)then
            call pick_mesh_cell(s,x,y,mesh,part,cell,face)
            if(mesh>0)then
                target=selection_target_t(selection_mesh,mesh)
                if(s%mesh_pick)target=selection_target_t(selection_mesh,mesh,part=part,item=cell,feature=face)
            end if
        else
            if(s%feature_pick_mode==pick_automatic.or.s%feature_pick_mode==pick_control_points)then
                point=pick_shared_control(s,x,y)
                if(point>0)then
                    target=selection_target_t(selection_assembly,s%selected_multipatch,item=point)
                else
                    call pick_selected_control(s,x,y,kind,id,point)
                    if(point>0)target=selection_target_t(selection_control,id,kind,item=point)
                end if
            end if
            if(target%category==0)then
                select case(s%feature_pick_mode)
                case(pick_automatic)
                    call pick_scene_patch(s,x,y,kind,id)
                    if(kind>0)then
                        target=selection_target_t(selection_patch,id,kind)
                    else
                        call pick_mesh_cell(s,x,y,mesh,part,cell,face)
                        if(mesh>0)target=selection_target_t(selection_mesh,mesh)
                    end if
                case(pick_multipatch)
                    call pick_scene_patch(s,x,y,kind,id)
                    target=viewport_patch_target(s,kind,id)
                case(pick_control_points)
                    continue
                case default
                    call pick_scene_patch(s,x,y,kind,id,element,number)
                    if(kind>0)target=selection_target_t(selection_patch,id,kind,item=number,feature=element)
                end select
            end if
        end if
        call show_scene_context(s,target,s%surface%mapToGlobal(QPoint(nint(x),nint(y))))
    end subroutine show_viewport_context

    subroutine build_project_tree(layout,s)
        type(QVBoxLayout),intent(inout)::layout
        type(studio_t),intent(inout),target::s
        type(QtfEventFilter)::press_filter,key_filter

        call build_scene_search(layout,s)
        call build_scene_empty(layout,s)
        call create_project_tree_widgets(s%widgets%project_browser,layout)
        press_filter=QtfEventFilter(s%widgets%project_browser%viewport)
        call press_filter%mousePressEvent%connect(on_project_mouse_press,s)
        key_filter=QtfEventFilter(s%widgets%project_browser%tree)
        call key_filter%keyPressEvent%connect(on_project_key,s)
        call key_filter%showEvent%connect(on_scene_tree_shown,s)
        call s%widgets%project_browser%viewport%contextMenuEvent%connect(on_project_context,s)
        call rebuild_project_tree(s)
        call s%widgets%project_browser%tree%currentItemChanged%connect(on_tree_selection,s)
        call s%widgets%project_browser%tree%itemClicked%connect(on_tree_click,s)
        call s%widgets%project_browser%tree%itemDoubleClicked%connect(on_tree_open,s)
        call s%widgets%project_browser%tree%itemExpanded%connect(on_scene_branch_open,s)
    end subroutine build_project_tree

    subroutine on_scene_branch_open(item,data)
        type(QTreeWidgetItem),intent(in) :: item
        class(*),optional,intent(inout) :: data

        select type(s=>data)
        type is(studio_t)
            if(s%updating.or.s%applying_selection.or.item%isNull())return
            if(s%selection%count()>1)call sync_selection_rows(s)
        end select
    end subroutine on_scene_branch_open

    subroutine on_scene_tree_shown(data)
        class(*),optional,intent(inout) :: data

        select type(s=>data)
        type is(studio_t)
            if(s%updating.or.s%applying_selection)return
            if(s%selection%count()>1)call sync_selection_rows(s)
        end select
    end subroutine on_scene_tree_shown

    subroutine build_scene_search(layout,s)
        type(QVBoxLayout),intent(inout)::layout
        type(studio_t),intent(inout),target::s
        type(QWidget)::bar
        type(QHBoxLayout)::row
        type(QMenu)::menu
        type(QActionGroup)::choices
        type(QAction)::exit_action
        integer::i
        character(15),parameter::labels(9)=[character(15):: &
            'All objects','Geometry','Multipatches','Meshes','Hidden','Needs update','Curves','Surfaces','Volumes']
        integer,parameter::menu_order(9)=[1,2,7,8,9,3,4,5,6]

        bar=QWidget(s%panels(1))
        call bar%setObjectName('studioSceneSearchBar')
        row=QHBoxLayout(bar)
        call row%setContentsMargins(panel_inset,6,panel_inset,4)
        call row%setSpacing(5)
        associate(view=>s%widgets%project_browser)
            view%search=QLineEdit(parent=bar)
            call view%search%setObjectName('studioSceneSearch')
            call view%search%setPlaceholderText('Find in scene')
            call view%search%setAccessibleName('Find in scene')
            call set_status_hint(view%search,'Find objects by name · Ctrl+F')
            call view%search%setClearButtonEnabled(.true.)
            call view%search%setFixedHeight(29)
            call view%search%textChanged%connect(on_scene_search,s)
            call row%addWidget(view%search,1)
            view%filter_button=QToolButton(bar)
            call view%filter_button%setObjectName('studioSceneFilterButton')
            call view%filter_button%setIcon(line_icon(icon_scene_filter))
            call view%filter_button%setIconSize(19,19)
            call view%filter_button%setFixedSize(31,29)
            call view%filter_button%setCheckable(.true.)
            call view%filter_button%setAccessibleName('Filter scene objects')
            call set_status_hint(view%filter_button,'Filter scene objects')
            menu=QMenu('Filter scene objects',view%filter_button)
            call register_status_menu(menu)
            choices=QActionGroup(menu)
            call choices%setExclusive(.true.)
            do i=1,size(labels)
                view%filters(i)=QAction(trim(labels(i)),menu)
                call view%filters(i)%setObjectName('studioSceneFilter'//integer_text(i-1))
                call view%filters(i)%setCheckable(.true.)
                call choices%addAction(view%filters(i))
            end do
            do i=1,size(menu_order)
                if(i==2.or.i==8)call menu%addSeparator()
                call menu%addAction(view%filters(menu_order(i)))
            end do
            call view%filters(1)%setChecked(.true.)
            call choices%triggered%connect(on_scene_filter,s)
            call view%filter_button%setMenu(menu)
            call view%filter_button%setPopupMode(QToolButton_InstantPopup)
            call row%addWidget(view%filter_button)
            call layout%addWidget(bar)

            view%exit_isolation=QToolButton(s%panels(1))
            call view%exit_isolation%setObjectName('studioSceneExitIsolation')
            call make_command(s%panels(1),'Exit isolation','',s,scene_exit_isolation_command,icon_restore,exit_action)
            call view%exit_isolation%setDefaultAction(exit_action)
            call view%exit_isolation%setToolButtonStyle(Qt_ToolButtonTextBesideIcon)
            call view%exit_isolation%setAccessibleName('Exit isolation')
            call set_status_hint(view%exit_isolation,'Restore the scene visibility before isolation')
            call view%exit_isolation%setFixedHeight(28)
            call layout%addWidget(view%exit_isolation)
            call view%exit_isolation%hide()
        end associate
    end subroutine build_scene_search

    subroutine build_scene_empty(layout,s)
        type(QVBoxLayout),intent(inout)::layout
        type(studio_t),intent(inout),target::s
        type(QVBoxLayout)::content

        associate(view=>s%widgets%project_browser)
            view%empty_state=QWidget(s%panels(1))
            call view%empty_state%setObjectName('studioSceneEmpty')
            content=QVBoxLayout(view%empty_state)
            call content%setContentsMargins(24,12,24,12)
            call content%setSpacing(7)
            call content%addStretch()
            view%empty_title=QLabel('No geometry yet',view%empty_state)
            call view%empty_title%setObjectName('studioSceneEmptyTitle')
            call view%empty_title%setAlignment(Qt_AlignCenter)
            call content%addWidget(view%empty_title)
            view%empty_clear=QPushButton('Clear filters',view%empty_state)
            call view%empty_clear%setObjectName('studioSceneClearFilters')
            call view%empty_clear%clicked%connect(clear_scene_filters,s)
            call content%addWidget(view%empty_clear)
            call content%addStretch()
            call layout%addWidget(view%empty_state,1)
            call view%empty_state%hide()
        end associate
    end subroutine build_scene_empty

    subroutine on_scene_search(text,data)
        character(*),intent(in)::text
        class(*),intent(inout),optional::data

        select type(s=>data)
        type is(studio_t)
            if(s%widgets%project_browser%changing_filter)return
            s%widgets%project_browser%search_text=text
            call sync_scene_filter_ui(s)
        end select
    end subroutine on_scene_search

    subroutine on_scene_filter(index,data)
        integer,intent(in)::index
        class(*),intent(inout),optional::data

        select type(s=>data)
        type is(studio_t)
            if(s%widgets%project_browser%changing_filter)return
            if(index<0.or.index>=size(s%widgets%project_browser%filters))return
            s%widgets%project_browser%filter_mode=index
            call sync_scene_filter_ui(s)
        end select
    end subroutine on_scene_filter

    subroutine clear_scene_filters(data)
        class(*),intent(inout),optional::data

        select type(s=>data)
        type is(studio_t)
            associate(view=>s%widgets%project_browser)
                view%changing_filter=.true.
                call view%search%clear()
                call view%filters(1)%setChecked(.true.)
                view%filter_mode=0
                view%search_text=''
                view%changing_filter=.false.
            end associate
            call sync_scene_filter_ui(s)
            call s%widgets%project_browser%search%setFocus()
        end select
    end subroutine clear_scene_filters

    subroutine sync_scene_filter_ui(s)
        type(studio_t),intent(inout)::s
        logical::was_updating

        was_updating=s%updating
        s%updating=.true.
        call apply_project_filter(s%widgets%project_browser,s)
        call sync_selection_rows(s)
        s%updating=was_updating
    end subroutine sync_scene_filter_ui

    subroutine build_domain(layout,s)
        type(QVBoxLayout),intent(inout)::layout
        type(studio_t),intent(inout),target::s
        type(QHBoxLayout)::row
        type(QVBoxLayout)::groups_layout
        type(QWidget)::display_row,knot_frame
        type(QHeaderView)::knot_header,knot_rows
        type(QLabel)::dash
        character(len=*),parameter::axes='UVW'
        integer::i

        row=QHBoxLayout()
        call row%setContentsMargins(panel_inset,6,panel_inset,5)
        call row%setSpacing(3)
        s%widgets%domain_probe_button=QPushButton('Probe')
        s%widgets%domain_knots_button=QPushButton('Knots')
        call s%widgets%domain_probe_button%setObjectName('axisTab')
        call s%widgets%domain_knots_button%setObjectName('axisTab')
        call s%widgets%domain_probe_button%setCheckable(.true.)
        call s%widgets%domain_knots_button%setCheckable(.true.)
        call set_status_hint(s%widgets%domain_probe_button,'Move the parameter probe')
        call set_status_hint(s%widgets%domain_knots_button,'Edit knots · Hold Shift temporarily')
        call connect_button(s%widgets%domain_probe_button,s,domain_probe_command)
        call connect_button(s%widgets%domain_knots_button,s,90)
        call row%addWidget(s%widgets%domain_probe_button)
        call row%addWidget(s%widgets%domain_knots_button)
        call row%addStretch()
        do i=1,3
            s%knot_buttons(i)=QPushButton(axes(i:i))
            call s%knot_buttons(i)%setObjectName('axisTab')
            call s%knot_buttons(i)%setProperty('parameterAxis',axes(i:i))
            call s%knot_buttons(i)%setCheckable(.true.)
            call set_status_hint(s%knot_buttons(i),'Activate parameter axis '//axes(i:i))
            call s%knot_buttons(i)%setAccessibleName('Parameter axis '//axes(i:i))
            call connect_button(s%knot_buttons(i),s,39+i)
            call row%addWidget(s%knot_buttons(i))
        end do
        call layout%addLayout(row)

        row=QHBoxLayout()
        call row%setContentsMargins(panel_inset,0,panel_inset,6)
        call row%setSpacing(5)
        s%widgets%interval_label=QLabel('U range')
        call s%widgets%interval_label%setObjectName('muted')
        call row%addWidget(s%widgets%interval_label)
        s%widgets%interval_lower=QLineEdit()
        s%widgets%interval_upper=QLineEdit()
        call s%widgets%interval_lower%setObjectName('studioDomainLower')
        call s%widgets%interval_upper%setObjectName('studioDomainUpper')
        call s%widgets%interval_lower%setAccessibleName('Lower parameter bound')
        call s%widgets%interval_upper%setAccessibleName('Upper parameter bound')
        call s%widgets%interval_lower%setMaxLength(128)
        call s%widgets%interval_upper%setMaxLength(128)
        call set_status_hint(s%widgets%interval_lower,'Lower bound · Enter applies')
        call set_status_hint(s%widgets%interval_upper,'Upper bound · Enter applies')
        call s%widgets%interval_lower%returnPressed%connect(on_interval_apply,s)
        call s%widgets%interval_upper%returnPressed%connect(on_interval_apply,s)
        call row%addWidget(s%widgets%interval_lower,1)
        dash=QLabel('–')
        call row%addWidget(dash)
        call row%addWidget(s%widgets%interval_upper,1)
        call domain_button(row,s%widgets%interval_apply,s,apply_interval_command,35,'DomainApplyInterval','Apply parameter bounds')
        call layout%addLayout(row)
        s%widgets%interval_validator=QRegularExpressionValidator(parameter_number_pattern,s%window)
        s%widgets%interval_context=0
        s%widgets%interval_revision=-1_int64

        display_row=QWidget(); call display_row%setObjectName('studioPlotControls')
        row=QHBoxLayout(display_row)
        call row%setContentsMargins(panel_inset,0,panel_inset,6)
        call row%setSpacing(4)
        call build_domain_overlays(row,display_row,s)
        call row%addStretch()
        s%widgets%domain_refine_button=QPushButton('Refine…')
        call s%widgets%domain_refine_button%setObjectName('studioRefinePatch')
        call s%widgets%domain_refine_button%setIcon(line_icon(64))
        call s%widgets%domain_refine_button%setIconSize(17,17)
        call set_status_hint(s%widgets%domain_refine_button,'Subdivide spans and elevate degree with one Undo step')
        call connect_button(s%widgets%domain_refine_button,s,refine_command)
        call row%addWidget(s%widgets%domain_refine_button)
        call layout%addWidget(display_row)

        s%domain=QtfWidget()
        call s%domain%setMinimumSize(200,140)
        call s%domain%setAccessibleName('Interactive parameter space')
        call s%domain%setMouseTracking(.true.)
        call s%domain%setFocusPolicy(Qt_StrongFocus)
        call s%domain%setObjectName('domainCanvas')
        call layout%addWidget(s%domain,1)

        s%widgets%domain_knot_tools=QWidget()
        call s%widgets%domain_knot_tools%setObjectName('studioDomainKnotTools')
        row=QHBoxLayout(s%widgets%domain_knot_tools)
        call row%setContentsMargins(panel_inset,3,panel_inset,4)
        call row%setSpacing(4)
        call domain_button(row,s%widgets%domain_knot_buttons(1),s,domain_select_command,18,'DomainSelectKnot','Select or drag knots')
        call domain_button(row,s%widgets%domain_knot_buttons(2),s,domain_insert_command,30,'DomainInsertKnot','Click a ruler gap to insert')
        call domain_button(row,s%widgets%domain_knot_buttons(3),s,domain_remove_command,36,'DomainRemoveKnot','Click an interior knot to remove')
        do i=1,3
            call s%widgets%domain_knot_buttons(i)%setCheckable(.true.)
        end do
        call row%addStretch()
        call domain_button(row,s%widgets%domain_groups_button,s,domain_details_command,67,'DomainGroupsButton','Show exact knot groups')
        call s%widgets%domain_groups_button%setCheckable(.true.)
        call domain_button(row,s%widgets%insert_knot_button,s,insert_knot_command,65,'InsertKnotButton','Insert at an exact value')
        call domain_button(row,s%widgets%remove_knot_button,s,remove_knot_command,66,'RemoveKnotButton','Remove at an exact value')
        call domain_button(row,s%widgets%domain_values_button,s,edit_knot_values_command,68,'DomainKnotValues','Edit exact knot values')
        call layout%addWidget(s%widgets%domain_knot_tools)

        knot_frame=QWidget()
        s%widgets%domain_groups=knot_frame
        groups_layout=QVBoxLayout(knot_frame)
        call groups_layout%setContentsMargins(panel_inset,0,panel_inset,4)
        s%knot_table=QTableWidget(0,3)
        call s%knot_table%setObjectName('knotTable')
        call s%knot_table%setVerticalHeaderVisible(.false.)
        call s%knot_table%setWordWrap(.false.)
        call s%knot_table%setAlternatingRowColors(.true.)
        call s%knot_table%setMouseTracking(.true.)
        knot_rows=s%knot_table%verticalHeader()
        call knot_rows%setSectionResizeMode(QHeaderView_Fixed)
        knot_header=s%knot_table%horizontalHeader()
        call knot_header%setSectionResizeMode(0,QHeaderView_Fixed)
        call knot_header%setSectionResizeMode(1,QHeaderView_Stretch)
        call knot_header%setSectionResizeMode(2,QHeaderView_Fixed)
        call s%knot_table%setProperty('horizontalScrollBarPolicy',Qt_ScrollBarAlwaysOff)
        call s%knot_table%setAccessibleName('Grouped knots of the active parameter axis')
        call s%knot_table%setAccessibleDescription('Indices, knot value and multiplicity for each group.')
        call s%knot_table%setSelectionMode(QAbstractItemView_SingleSelection)
        call s%knot_table%setSelectionBehavior(QAbstractItemView_SelectRows)
        call s%knot_table%setEditTriggers(QAbstractItemView_NoEditTriggers)
        call s%knot_table%cellDoubleClicked%connect(on_knot_double_click,s)
        call s%knot_table%currentCellChanged%connect(on_knot_table_selection,s)
        call groups_layout%addWidget(s%knot_table)
        call layout%addWidget(s%widgets%domain_groups)
        call s%widgets%domain_groups%hide()
        call s%widgets%domain_knot_tools%hide()

    end subroutine build_domain

    subroutine build_inspector(layout,s)
        type(QVBoxLayout),intent(inout)::layout
        type(studio_t),intent(inout),target::s

        call build_inspector_widgets(layout,[s%u,s%v,s%w],s%inspector,s%probe_u,s%probe_v,s%probe_w, &
            s%probe_label,s%diagnostic_label)
        call connect_button(s%inspector%headers(1),s,inspector_section_command_base+1)
        call s%probe_u%valueChanged%connect(on_probe,s)
        call s%probe_v%valueChanged%connect(on_probe,s)
        call s%probe_w%valueChanged%connect(on_probe,s)
        call connect_button(s%inspector%reset_probe,s,reset_probe_command)
        call connect_button(s%inspector%edit_control,s,edit_controls_command)
        call connect_button(s%inspector%headers(2),s,inspector_section_command_base+2)
        call connect_button(s%inspector%repair_diagnostic,s,diagnose_with_ai_command)
    end subroutine build_inspector

    subroutine menu_action(menu,text,shortcut,s,command,icon_kind,menu_item)
        type(QMenu),intent(inout)::menu
        character(len=*),intent(in)::text,shortcut
        type(studio_t),intent(inout),target::s
        integer,intent(in)::command,icon_kind
        type(QAction),intent(out),optional::menu_item
        type(QAction)::action
        call make_command(menu,text,shortcut,s,command,icon_kind,action)
        call menu%addAction(action)
        if(present(menu_item))menu_item=action
    end subroutine menu_action

    subroutine make_command(parent,text,shortcut,s,command,icon_kind,action)
        class(QObject),intent(in)::parent
        character(*),intent(in)::text,shortcut
        type(studio_t),intent(inout),target::s
        integer,intent(in)::command,icon_kind
        type(QAction),intent(out)::action

        associate(visual=>command_visual(command,text,icon_kind))
            action=QAction(visual%caption,parent)
            call set_status_hint(action,studio_command_tooltip(command,visual%caption,shortcut))
            if(visual%glyph>0)call action%setIcon(line_icon(visual%glyph))
        end associate
        if(len_trim(shortcut)>0)then
            call action%setShortcut(shortcut)
            select case(command)
            case(new_project_command,open_project_command,import_iges_command,find_command,canvas_command, &
                 close_project_tab_command,fullscreen_command)
                call s%window%addAction(action)
            case default
                call s%surface%addAction(action)
            end select
        end if
        if(command==edit_controls_command)s%widgets%edit_controls_action=action
        if(command==undo_patch_command)s%widgets%undo_patch_action=action
        if(command==insert_knot_command)s%widgets%insert_knot_action=action
        if(command==remove_knot_command)s%widgets%remove_knot_action=action
        if(command==elevate_degree_command)s%widgets%elevate_degree_action=action
        if(command==reduce_degree_command)s%widgets%reduce_degree_action=action
        if(command==edit_knot_values_command)then
            s%widgets%edit_knot_values_action=action
            call action%setObjectName('studioEditKnotValuesAction')
        end if
        if(command==edit_interval_command)s%widgets%edit_interval_action=action
        select case(command)
        case(transform_command)
            s%widgets%transform_action=action
        case(color_command)
            s%widgets%material_color_action=action
        case(translucent_command)
            s%widgets%translucent_action=action
            call action%setCheckable(.true.)
        case(export_iges_command)
            s%widgets%export_iges_action=action
        case(export_vtk_command)
            s%widgets%export_vtk_action=action
        case(canvas_command)
            s%widgets%canvas_action=action
            call action%setCheckable(.true.)
        case(snap_off_command:snap_closest_command)
            s%widgets%snap_mode_actions(command-snap_off_command+1)=action
            call action%setObjectName('studioSnapMode'//integer_text(command-snap_off_command))
            call action%setCheckable(.true.)
        case(pick_mode_base:pick_mode_base+3,pick_controls_command,pick_multipatch_command)
            s%widgets%pick_actions(pick_action_index(command))=action
            call action%setCheckable(.true.)
        end select
        if(command>=edit_multipatch_command.and.command<=remove_multipatch_command) &
            s%widgets%multipatch_actions(command-edit_multipatch_command+1)=action
        if(command>=select_curve_command.and.command<=select_volume_command)then
            s%widgets%entity_menu_actions(command-select_curve_command+1)=action
        end if
        if(command==0)then
            call action%setEnabled(.false.)
            call set_status_hint(action,studio_tooltip(text,'Planned for a later ForCAD Studio release.'))
        else
            call register_action(s,command)
            call action%triggered%connect(on_action,s%action_bindings(s%action_count))
        end if
    end subroutine make_command

    module procedure sync_controls
        logical::was_updating
        integer::i,preset
        was_updating=s%updating; s%updating=.true.
        preset=0
        if(s%viewport_style==viewport_edges)preset=display_style(s%layers)
        call sync_viewport_style(s)
        call s%probe_u%setValue(s%u); call s%probe_v%setValue(s%v); call s%probe_w%setValue(s%w)
        do i=1,size(s%viewport_layers)
            call s%viewport_layers(i)%setChecked(s%layers(i))
        end do
        do i=1,5
            call s%domain_overlays(i)%setChecked(s%domain_layers(i))
        end do
        call s%basis_buttons(1)%setChecked(.not.s%basis_derivative)
        call s%basis_buttons(2)%setChecked(s%basis_derivative)
        do i=1,size(s%knot_buttons)
            call s%knot_buttons(i)%setChecked(i==s%knot_axis)
        end do
        do i=1,size(s%viewport_presets)
            call s%viewport_presets(i)%setChecked(i==preset)
        end do
        call s%basis_options(1)%setChecked(s%basis_legend)
        call s%basis_options(2)%setChecked(s%basis_probe_guide)
        call s%viewport_grid%setChecked(s%world_grid)
        call s%widgets%ribbon%interactive_quality_action%setChecked(s%interactive_quality)
        call sync_transform_controls(s)
        call s%toolbar_actions(11)%setChecked(s%layers(1))
        call s%toolbar_actions(12)%setChecked(s%layers(2))
        call s%toolbar_actions(15)%setChecked(s%domain_layers(1))
        call s%toolbar_actions(16)%setChecked(s%layers(3))
        call s%toolbar_actions(17)%setChecked(s%domain_layers(4))
        call sync_parameter_panels(s)
        call s%widgets%domain_probe_button%setChecked(.not.s%domain_edit_knots)
        call s%widgets%domain_knots_button%setChecked(s%domain_edit_knots)
        call s%widgets%domain_groups_button%setChecked(s%domain_knot_details)
        if(s%domain_knot_tool==knot_remove_tool)then
            if(.not.can_remove_knot(s))s%domain_knot_tool=knot_select_tool
        end if
        do i=1,3
            call s%widgets%domain_knot_buttons(i)%setChecked(i-1==s%domain_knot_tool)
        end do
        call s%widgets%domain_knot_tools%setVisible(s%domain_edit_knots)
        call s%widgets%domain_groups%setVisible(s%domain_edit_knots.and.s%domain_knot_details)
        call s%widgets%project_inspect_action%setEnabled(s%selected_multipatch==0.and.s%live_scene%valid)
        call s%widgets%canvas_action%setChecked(s%canvas_only)
        call sync_snap_controls(s)
        call sync_selection_filter(s)
        call sync_viewport_orientation(s)
        s%updating=was_updating
        call refresh_workspace_status(s)
    end procedure sync_controls

    integer function display_style(layers) result(index)
        logical,intent(in)::layers(viewport_layer_count)
        integer::i

        index=0
        do i=1,size(display_layer_presets,2)
            if(.not.all(layers.eqv.display_layer_presets(:,i)))cycle
            index=i
            return
        end do
    end function display_style

    subroutine register_action(s,command)
        type(studio_t),intent(inout),target::s
        integer,intent(in)::command
        s%action_count=s%action_count+1
        if(s%action_count>size(s%action_bindings))error stop 'Too many UI actions'
        s%action_bindings(s%action_count)%state=>s; s%action_bindings(s%action_count)%command=command
    end subroutine register_action

    subroutine connect_button(button,s,command)
        type(QPushButton),intent(inout)::button
        type(studio_t),intent(inout),target::s
        integer,intent(in)::command
        call register_action(s,command)
        call button%clicked%connect(on_action,s%action_bindings(s%action_count))
    end subroutine connect_button

    subroutine on_viewport_delete(data)
        class(*),intent(inout),optional::data
        select type(s=>data)
        type is(studio_t)
            if(s%dragging.or.s%updating)return
            if(s%widgets%project_remove_action%isEnabled())call s%widgets%project_remove_action%trigger()
        end select
    end subroutine on_viewport_delete

    subroutine on_action(data)
        class(*),intent(inout),optional::data
        select type(context=>data)
        type is(studio_action_binding_t)
            call invoke_command(context%state,context%command)
        end select
    end subroutine on_action

    module procedure restore_workspace_layout
        integer::i,panel
        logical :: visible

        s%focused_surface=s%selected_multipatch>0
        s%saved_panels=.false.
        s%saved_panels(workspace_panels)=.true.
        do i=1,size(workspace_panels)
            panel=workspace_panels(i)
            visible=panel<5.or..not.s%focused_surface
            call set_panel_visibility(s,panel,visible)
        end do
        call s%columns%setSizes([264,1041,342]); call s%center_rows%setSizes([488,293])
        s%normal_columns=s%columns%sizes(); s%normal_center=s%center_rows%sizes()
        s%assistant_width=default_assistant_width
        call size_assistant_column(s)
        call position_assistant_overlay(s)
        call sync_project_column(s)
    end procedure restore_workspace_layout

    module procedure leave_viewport_focus
        integer::i,panel
        s%focused_surface=.false.
        do i=1,size(workspace_panels)
            panel=workspace_panels(i)
            call set_panel_visibility(s,panel,s%saved_panels(panel))
        end do
        call s%columns%setSizes(s%normal_columns); call s%center_rows%setSizes(s%normal_center)
        call sync_project_column(s)
    end procedure leave_viewport_focus

    module procedure toggle_canvas

        if(.not.s%canvas_only)then
            s%before_canvas_focus=s%focused_surface
            s%before_canvas_project=s%panels(1)%isVisible()
            s%before_canvas_ai=studio_assistant_is_visible(s)
            if(.not.s%focused_surface)call invoke_command(s,toggle_parameter_panels_command)
            s%canvas_only=.true.
            call s%command_ribbon%hide()
            call s%widgets%leave_canvas_button%show()
            call set_panel_visibility(s,1,.false.)
            call hide_studio_assistant(s)
        else
            s%canvas_only=.false.
            call s%command_ribbon%show()
            call s%widgets%leave_canvas_button%hide()
            if(.not.s%before_canvas_focus.and.s%selected_multipatch==0)call invoke_command(s,toggle_parameter_panels_command)
            call set_panel_visibility(s,1,s%before_canvas_project)
            if(s%before_canvas_ai)call invoke_command(s,show_assistant_command)
        end if
        call sync_project_column(s)
        call s%surface%setFocus()
    end procedure toggle_canvas

    module procedure set_panel_visibility
        logical :: in_left_column

        in_left_column=index==1
        if(in_left_column.and.visible.and..not.s%canvas_only)call s%left_column%show()
        if(visible.and.(index==5.or.index==6).and..not.s%canvas_only)s%focused_surface=.false.
        call s%panels(index)%setVisible(visible)
        if(index==6)call s%lower_row%setVisible(visible)
        if(index==4)call s%panel_actions(index)%setChecked(visible)
        if(index==8)s%widgets%inspector_requested=visible
        if(in_left_column)call sync_project_column(s)
    end procedure set_panel_visibility

    module procedure sync_project_column

        if(s%canvas_only)then
            call s%left_column%hide()
        else if(s%workspace%isVisible())then
            call s%left_column%setVisible(s%panels(1)%isVisible())
        end if
    end procedure sync_project_column
end submodule studio_ui_layout
