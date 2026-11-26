submodule(studio_ui) studio_ui_project
    use studio_display_refinement, only: finish_display_refinement
    use studio_document, only: set_document_project_directory
    use studio_fortran_export, only: export_project_fortran
    use studio_status_hints, only: register_status_menu
    use qt, only: QFileDialog_ExistingFiles,QInputDialog_getItem, QHeaderView, QScrollBar, QFontMetrics, Qt_ControlModifier, Qt_Key_A, Qt_ElideRight
    use studio_viewport_scene, only: pick_face, pick_edge, pick_vertex
    use studio_icons, only: icon_measure
    use studio_units, only: unit_labels
    use studio_commands, only: model_ribbon_tab, sweep_command, loft_command, rotate_geometry_command, scale_geometry_command
    use studio_command_presentation, only: command_visual
    use studio_selection, only: selection_set_t, selection_add, change_selection_batch
    use studio_state, only: scene_patch_drawn, scene_mesh_drawn, scene_drawn_patch_count, studio_mesh_source_group
    use studio_project_tree, only: project_row_visibility, scene_visibility_column, reveal_project_row, project_target_row
    use studio_workspace_selection, only: selection_mesh_records
    implicit none(type, external)
contains
    module procedure show_home
        type(QPushButton)::new_button
        if(s%document%document_started)then
            s%resume_assistant=studio_assistant_is_visible(s)
            if(s%resume_assistant)call hide_studio_assistant(s)
        end if
        s%home_active=.true.
        call s%widgets%resume_button%setVisible(s%document%document_started)
        call s%pages%setCurrentWidget(s%home_page)
        call s%widgets%find_command_button%hide()
        if(s%window%isVisible())call s%window%activateWindow()
        if(s%document%document_started)then
            call s%widgets%resume_button%setFocus()
        else
            call s%home_page%findChild('studioHomeAction'//integer_text(new_project_command),new_button)
            call new_button%setFocus()
        end if
        call s%status%clearMessage()
        call s%status%setActivity('')
        call refresh_workspace_status(s)
        call sync_project_title(s)
    end procedure show_home

    module procedure enter_workspace
        logical::was_home,first_document

        was_home=s%home_active
        first_document=.not.s%document%document_started
        s%home_active=.false.
        s%document%document_started=.true.
        call ensure_project_tab(s)
        call s%pages%setCurrentWidget(s%workspace)
        call s%widgets%find_command_button%show()
        call refresh_workspace_status(s)
        call sync_project_title(s)
        if(first_document.and.start_assistant.and..not.s%focused_surface)call invoke_command(s,toggle_parameter_panels_command)
        if(was_home.and.start_assistant)call invoke_command(s,show_assistant_command)
        call s%surface%setFocus()
    end procedure enter_workspace

    module procedure document_patch
        present_in_project=document_has_patch(s%document,kind,variant)
    end procedure document_patch

    module procedure new_studio_project

        call reset_project_document(s%document,.false.)
        if(len_trim(folder)>0)call set_document_project_directory(s%document,folder)
        call prepare_document_view(s)
        call enter_workspace(s,.true.)
        call pick_document_selection(s)
        if(.not.s%focused_surface)call invoke_command(s,toggle_parameter_panels_command)
        call mark_project_saved(s)
        call show_ready(s)
    end procedure new_studio_project

    module procedure reset_document


        call reset_project_document(s%document,templates)
        call prepare_document_view(s)
    end procedure reset_document

    module procedure prepare_document_view
        type(QTabWidget)::tabs

        call cancel_transform_drag(s)
        call cancel_area_drag(s)
        call s%selection%clear()
        call s%isolation%clear()
        if(allocated(s%isolation_mesh_names))deallocate(s%isolation_mesh_names)
        s%selected_multipatch=0
        call clear_mesh_selection(s)
        s%mesh_only=.false.
        s%mesh_pick=.false.
        s%mesh_shaded=.true.
        s%mesh_view=mesh_exterior; s%mesh_layer_axis=3; s%mesh_layer_index=1
        s%mesh_section_axis=0; s%mesh_section_fraction=.5_dp
        s%physical_knots=.false.
        s%selected_control=0
        s%control_dragging=.false.
        s%viewport_feature=0
        s%transform_mode=0
        s%live_scene%valid=.false.
        s%live_scene%dimension=0
        s%live_scene%variant=0
        s%live_probe_axis=0; s%live_plot_axis=0
        call reset_domain_view(s)
        s%zoom=1._dp; s%pan_x=0._dp; s%pan_y=0._dp
        call rebuild_project_tree(s)
        call s%workspace%findChild('studioRibbonTabs',tabs)
        call tabs%setCurrentIndex(model_ribbon_tab)
        call refresh_inspector(s)
        call refresh_knot_table(s)
        call sync_controls(s)
        call repaint_all(s)
    end procedure prepare_document_view

    module procedure pick_document_selection
        integer::kind,variant,pass,visible_members,i

        if(count_multipatches(s%document%multipatches)==1)then
            associate(group=>s%document%multipatches%groups(1))
                if(group%validation%valid)then
                    kind=group%definition%kind
                    visible_members=count([(studio_scene_is_visible(s%document%scene,kind,group%definition%members(i),catalog=s%document%patches), &
                        i=1,size(group%definition%members))])
                    if(visible_members>0.and.visible_members==studio_scene_count(s%document%scene))then
                        call select_multipatch_group(s,1,reveal=.false.)
                        return
                    end if
                end if
            end associate
        end if
        do pass=1,2
            do kind=1,3
                do variant=1,single_patch_count(kind,catalog=s%document%patches)
                    if(.not.document_patch(s,kind,variant))cycle
                    if(pass==1.and..not.studio_scene_is_visible(s%document%scene,kind,variant,catalog=s%document%patches))cycle
                    call select_patch(s,kind,variant,reveal=.false.)
                    return
                end do
            end do
        end do
        call select_patch(s,2,reveal=.false.)
    end procedure pick_document_selection

    module procedure sync_project_title
        character(:),allocatable::caption,tip,title,name

        if(s%widgets%project_caption%isNull())return
        title=studio_product_name
        caption=''; tip=''
        if(.not.s%home_active.and.s%document%document_started)then
            name=document_display_name(s%document)
            tip=name
            if(allocated(s%document%project_path))tip=s%document%project_path
            if(document_has_changes(s%document))then
                name=name//' *'
                tip=tip//new_line('a')//'Unsaved changes'
            end if
            caption='·  '//name
            title=title//' — '//name
        end if
        if(s%widgets%project_caption%text()/=caption)call s%widgets%project_caption%setText(caption)
        call set_status_hint(s%widgets%project_caption,tip)
        if(s%window%windowTitle()/=title)call s%window%setWindowTitle(title)
        call update_project_tabs(s)
        call refresh_project_tab_actions(s)
    end procedure sync_project_title

    module procedure mark_project_saved
        if(s%document%document_started)call ensure_project_tab(s)
        call mark_document_saved(s%document)
        call sync_project_title(s)
    end procedure mark_project_saved

    module procedure confirm_project_change
        integer::choice

        allowed=.true.
        if(.not.document_has_changes(s%document))return
        choice=QMessageBox_question('Save changes?', 'Save changes to '//document_display_name(s%document)//'?', &
            buttons=ior(ior(QMessageBox_Save,QMessageBox_Discard),QMessageBox_Cancel), &
            defaultButton=QMessageBox_Save,parent=s%window)
        allowed=choice==QMessageBox_Discard
        if(choice==QMessageBox_Save)allowed=write_project(s,.false.)
    end procedure confirm_project_change

    module procedure confirm_stored_document
        integer :: choice
        character(:),allocatable :: name

        allowed=.true.
        associate(document=>s%project_tabs%items(index)%document)
            if(.not.document_has_changes(document))return
            name=document_display_name(document)
        end associate
        choice=QMessageBox_question('Save changes?', 'Save changes to '//name//'?', &
            buttons=ior(ior(QMessageBox_Save,QMessageBox_Discard),QMessageBox_Cancel), &
            defaultButton=QMessageBox_Save,parent=s%window)
        allowed=choice==QMessageBox_Discard
        if(choice==QMessageBox_Save)allowed=write_tab_project(s,index,.false.)
    end procedure confirm_stored_document

    module procedure write_project
        saved=write_tab_project(s,s%project_tabs%active,save_as)
    end procedure write_project

    logical function write_tab_project(s,index,save_as) result(saved)
        type(studio_t),intent(inout),target :: s
        integer,intent(in) :: index
        logical,intent(in) :: save_as
        type(studio_document_t),pointer :: document
        character(:),allocatable :: path,error
        integer :: opened

        saved=.false.
        if(index==s%project_tabs%active)then
            document=>s%document
        else
            document=>s%project_tabs%items(index)%document
        end if
        if(.not.document%document_started)return
        path='Untitled.fcad'
        if(allocated(document%sample_title))path=document%sample_title//'.fcad'
        if(allocated(document%project_path))path=document%project_path
        if(save_as.or..not.allocated(document%project_path))then
            path=choose_project_save_path(s%window,document,path,save_as)
            if(len(path)==0)return
        end if
        opened=find_project_tab(s,path)
        if(opened>0.and.opened/=index)then
            call s%status%showMessage('Project already open in another tab · Choose another file',6500,warning=.true.)
            return
        end if
        call save_project_document(document,path,error,remember_path=.true.)
        if(len(error)>0)then
            call s%status%showMessage('Could not save the project',6500,warning=.true.,detail=error)
            return
        end if
        if(index==s%project_tabs%active)call sync_studio_assistant(s)
        call sync_project_title(s)
        call s%status%showMessage('Saved '//document_display_name(document),3500)
        saved=.true.
    end function write_tab_project

    function choose_project_save_path(window,document,suggestion,save_as) result(path)
        type(QMainWindow),intent(in) :: window
        type(studio_document_t),intent(in) :: document
        character(*),intent(in) :: suggestion
        logical,intent(in) :: save_as
        character(:),allocatable :: path,title
        type(QFileDialog) :: dialog
        type(qt_string),allocatable :: files(:)

        title=trim(merge('Save project as','Save project   ',save_as))
        dialog=QFileDialog(title,filter=studio_project_filter,parent=window)
        if(allocated(document%project_directory))call dialog%setDirectory(document%project_directory)
        call dialog%setObjectName('studioSaveProjectDialog')
        call dialog%setAcceptMode(QFileDialog_AcceptSave)
        call dialog%setDefaultSuffix('fcad')
        call dialog%selectFile(suggestion)
        path=''
        if(dialog%exec()==QDialog_Accepted)then
            files=dialog%selectedFiles()
            if(size(files)>0)path=files(1)%s
        end if
        call dialog%delete()
    end function choose_project_save_path

    module procedure choose_project
        type(QFileDialog) :: dialog
        type(qt_string),allocatable :: files(:)
        integer :: i

        dialog=QFileDialog('Open project',filter=studio_project_filter,parent=s%window)
        call dialog%setObjectName('studioOpenProjectDialog')
        call dialog%setFileMode(QFileDialog_ExistingFiles)
        if(allocated(s%document%project_directory))call dialog%setDirectory(s%document%project_directory)
        if(dialog%exec()==QDialog_Accepted)files=dialog%selectedFiles()
        call dialog%delete()
        if(.not.allocated(files))return
        do i=1,size(files)
            call open_project_path(s,files(i)%s)
        end do
    end procedure choose_project

    module procedure open_project_path
        character(:),allocatable :: error

        call open_project_tab(s,path,error)
        if(len(error)>0)call s%status%showMessage('Could not open the project',6500,warning=.true.,detail=error)
    end procedure open_project_path

    module procedure load_project_workspace
        integer :: opened

        if(present(remembered_path))then
            opened=find_project_tab(s,remembered_path)
        else
            opened=find_project_tab(s,path)
        end if
        if(opened>0.and.opened/=s%project_tabs%active)then
            error='This project is already open in another tab.'
            return
        end if
        call open_project_document(s%document,path,error)
        if(len(error)>0)return
        if(present(remembered_path))s%document%project_path=remembered_path
        if(present(project_directory))call set_document_project_directory(s%document,project_directory)
        call prepare_document_view(s)
        call enter_workspace(s,.true.)
        call fit_around_panels(s)
        call sync_project_title(s)
        call s%status%showMessage('Opened '//document_display_name(s%document),3500)
    end procedure load_project_workspace

    module procedure on_studio_close

        select type(s=>data)
        type is(studio_t)
            if(confirm_all_project_tabs(s))then
                call finish_display_refinement(s)
                if(s%remember_workspace.and..not.s%home_active)call save_preferences(s)
            else
                call event%ignore()
            end if
        end select
    end procedure on_studio_close

    module procedure choose_units
        integer :: model_unit,display_unit

        model_unit=s%document%model_unit
        display_unit=s%document%display_unit
        if(.not.edit_project_units(s%window,model_unit,display_unit))return
        call set_workspace_units(s,model_unit,display_unit)
    end procedure choose_units

    module procedure set_workspace_units
        character(*), parameter :: history_error='Could not record the project units in the edit history.'

        if(present(error))error=''
        if(s%document%model_unit==model_unit.and.s%document%display_unit==display_unit)return
        if(.not.record_document_properties(s%document,'Project units',units=.true.))then
            if(present(error))error=history_error
            call s%status%showMessage(history_error,5000,warning=.true.)
            return
        end if
        s%document%model_unit=model_unit
        s%document%display_unit=display_unit
        call reset_transform_values(s)
        call refresh_inspector(s)
        call repaint_all(s)
    end procedure set_workspace_units

    module procedure capture_workspace
        type(QPixmap) :: pixmap

        pixmap=s%window%grab()
        saved=pixmap%save(path)
        call pixmap%delete()
    end procedure capture_workspace

    module procedure active_length_unit
        id=s%document%display_unit
        if(id==unitless)id=s%document%model_unit
    end procedure active_length_unit

    module procedure project_unit_status
        integer::id

        id=active_length_unit(s)
        label='Unitless'
        if(id/=unitless)label=trim(unit_symbol(id))
    end procedure project_unit_status

    module procedure rebuild_project_tree
        type(QTreeWidgetItem)::group,item
        logical::was_updating,restored_current

        was_updating=s%updating
        s%updating=.true.
        call refresh_studio_scene(s%document%scene,catalog=s%document%patches)
        call build_project_rows(s%widgets%project_browser,s,restored_current)
        call sync_project_meshes(s)
        if(.not.restored_current.and.(.not.allocated(s%selection%items).or.s%selection%count()>0))then
            if(s%selected_mesh>0)then
                item=project_target_row(s%widgets%project_browser,selection_target_t(selection_mesh,s%selected_mesh))
                call s%widgets%project_browser%mesh_branch%setExpanded(.true.)
                call s%widgets%project_browser%tree%setCurrentItem(item)
            else if(s%selected_multipatch>0)then
                group=project_target_row(s%widgets%project_browser,selection_target_t(selection_assembly,s%selected_multipatch))
                call s%widgets%project_browser%tree%setCurrentItem(group)
            else
                call select_project_patch_row(s%widgets%project_browser,s,s%entity_kind)
            end if
        end if
        s%updating=was_updating
        call sync_selection_rows(s)
        call sync_patch_visibility_button(s)
    end procedure rebuild_project_tree

    module procedure sync_project_meshes
        integer::id,n

        if(s%widgets%project_browser%tree%isNull())return
        n=studio_mesh_count(s)
        if(s%selected_mesh>n)then
            s%selected_mesh=0
            s%mesh_cell=0
        end if
        call update_project_mesh_rows(s%widgets%project_browser,s)
        if(.not.has_visible_meshes(s))s%mesh_pick=.false.
        id=s%selected_mesh
        if(id<=0)return
        if(studio_mesh_valid(s,id).and.s%document%meshes(id)%visible)return
        s%mesh_cell=0
    end procedure sync_project_meshes

    module procedure select_project_mesh
        type(QTreeWidgetItem)::item
        logical::was_updating

        if(id<1.or.id>studio_mesh_count(s))return
        call forget_workspace_selection(s)
        if(.not.s%focused_surface)call invoke_command(s,toggle_parameter_panels_command)
        if(id/=s%selected_mesh)s%mesh_cell=0
        s%selected_mesh=id
        s%selected_multipatch=0
        s%viewport_feature=0
        s%selected_control=0
        s%control_dragging=.false.
        s%control_snapped=.false.
        s%transform_mode=0
        call cancel_transform_drag(s)
        call sync_project_meshes(s)
        item=project_target_row(s%widgets%project_browser,selection_target_t(selection_mesh,id))
        was_updating=s%updating
        s%updating=.true.
        call s%widgets%project_browser%mesh_branch%setExpanded(.true.)
        call s%widgets%project_browser%tree%setCurrentItem(item)
        s%updating=was_updating
        call sync_patch_visibility_button(s)
        call refresh_inspector(s)
        call refresh_knot_table(s)
        call sync_controls(s)
        call repaint_all(s)
    end procedure select_project_mesh

    module procedure refresh_multipatch_card
        integer::id

        id=s%selected_multipatch
        if(id>0)call refresh_multipatch(s%document%multipatches,id,.false.,catalog=s%document%patches)
        call show_multipatch_inspector_summary(s)
    end procedure refresh_multipatch_card

    module procedure select_multipatch_group
        type(QTreeWidgetItem)::group
        logical::was_updating,show_members
        integer::kind,i
        integer,allocatable::targets(:,:)

        if(id<1.or.id>count_multipatches(s%document%multipatches))return
        call forget_workspace_selection(s)
        show_members=.true.
        if(present(reveal))show_members=reveal
        if(.not.s%focused_surface)call invoke_command(s,toggle_parameter_panels_command)
        call clear_mesh_selection(s)
        if(id/=s%selected_multipatch)then
            s%selected_control=0
            s%canvas_gesture%viewport_shared_hover=0
            if(allocated(s%selected_shared_controls))deallocate(s%selected_shared_controls)
        end if
        s%selected_multipatch=id
        s%viewport_feature=0
        s%control_dragging=.false.
        if(show_members)then
        associate(definition=>s%document%multipatches%groups(id)%definition)
            kind=definition%kind
            do i=1,size(definition%members)
                call studio_scene_set_visible(s%document%scene,kind,definition%members(i),.true.,catalog=s%document%patches)
            end do
            allocate(targets(2,size(definition%members)))
            targets(1,:)=kind
            targets(2,:)=definition%members
            call sync_scene_visibility_batch(s%widgets%project_browser,s,targets)
        end associate
        end if
        group=project_target_row(s%widgets%project_browser,selection_target_t(selection_assembly,id))
        was_updating=s%updating; s%updating=.true.
        call s%widgets%project_browser%tree%setCurrentItem(group)
        if(show_members)call group%setExpanded(.true.)
        s%updating=was_updating
        call sync_patch_visibility_button(s)
        call sync_controls(s)
        call refresh_inspector(s)
        call s%surface%update()
    end procedure select_multipatch_group

    module procedure select_project_patch

        call select_project_patch_row(s%widgets%project_browser,s,kind)
        call sync_patch_visibility_button(s)
    end procedure select_project_patch

    module procedure sync_patch_visibility_button
        integer::i,visible_count,total
        integer,allocatable::mesh_ids(:)
        logical::selected_patch,can_fit,was_updating,isolated
        character(:),allocatable::label,summary
        type(QIcon)::visibility_icon
        type(selection_set_t)::targets
        type(selection_target_t)::focused

        was_updating=s%updating
        s%updating=.true.
        call update_project_summary(s%widgets%project_browser,s,summary)
        call s%status%setSceneSummary(summary)
        s%updating=was_updating
        focused=primary_workspace_target(s)
        selected_patch=focused%category==selection_patch.and.s%selection%count()<=1
        if(allocated(s%selection%items))selected_patch=selected_patch.and.s%selection%count()==1
        call s%widgets%project_duplicate_action%setEnabled(selected_patch)
        call s%widgets%project_rename_action%setEnabled(selected_patch)
        call sync_project_delete_controls(s)
        targets=scene_visibility_targets(s)
        total=targets%count()
        visible_count=count([(scene_target_visible(s,targets%items(i)),i=1,total)])
        can_fit=any([(scene_target_drawn(s,targets%items(i)),i=1,total)])
        if(s%mesh_only.and..not.can_fit)then
            mesh_ids=selection_mesh_records(s)
            can_fit=any([(scene_mesh_drawn(s,mesh_ids(i)),i=1,size(mesh_ids))])
        end if
        call s%widgets%project_fit_action%setEnabled(can_fit)
        call s%widgets%ribbon%view_fit_action%setEnabled(can_fit)
        label='Show selection'
        if(total>0.and.visible_count==total)then
            label='Hide selection'
            visibility_icon=s%widgets%project_browser%visible_icon
        else
            visibility_icon=s%widgets%project_browser%hidden_icon
        end if
        if(s%selection%count()>1)then
            label=label(:4)//' '//integer_text(s%selection%count())//' selected items'
        end if
        call s%widgets%scene_selection_visibility_action%setEnabled(total>0)
        call s%widgets%scene_selection_visibility_action%setText(label)
        call set_status_hint(s%widgets%scene_selection_visibility_action,label)
        call s%widgets%ribbon%view_visibility_action%setEnabled(total>0)
        call s%widgets%ribbon%view_visibility_action%setText(label)
        call s%widgets%ribbon%view_visibility_action%setIconText(label(:4))
        call s%widgets%scene_selection_visibility_action%setIcon(visibility_icon)
        call s%widgets%ribbon%view_visibility_action%setIcon(visibility_icon)
        call set_status_hint(s%widgets%ribbon%view_visibility_action,studio_tooltip(label,shortcut='Ctrl+Shift+H'))
        isolated=s%isolation%count()>0
        if(.not.s%widgets%project_browser%exit_isolation%isNull())then
            call s%widgets%project_browser%exit_isolation%setVisible(isolated)
        end if
        if(.not.s%widgets%viewport_visibility%isNull())then
            call s%widgets%viewport_visibility%setChecked(isolated.or.s%mesh_only)
            call s%widgets%viewport_hide_action%setEnabled(visible_count>0)
            call s%widgets%viewport_isolate_action%setEnabled(visible_count>0)
            call s%widgets%viewport_exit_isolation_action%setProperty('visible',isolated)
            call s%widgets%viewport_exit_isolation_action%setEnabled(isolated)
            if(s%widgets%viewport_visibility%boolProperty('exitsIsolation').neqv.isolated)then
                call s%widgets%viewport_visibility%setProperty('exitsIsolation',isolated)
                if(isolated)then
                    visibility_icon=line_icon(icon_restore)
                    call s%widgets%viewport_visibility%setIcon(visibility_icon)
                    call visibility_icon%delete()
                else
                    call s%widgets%viewport_visibility%setIcon(s%widgets%project_browser%visible_icon)
                end if
            end if
            if(isolated)then
                call set_status_hint(s%widgets%viewport_visibility,'Exit isolation · Restore visibility before isolation')
                call s%widgets%viewport_visibility%setAccessibleName('Exit isolation')
            else
                label='Visibility · Hide or isolate objects, show all, or show meshes only'
                if(s%mesh_only)label='Visibility · Meshes only is active'
                call set_status_hint(s%widgets%viewport_visibility,label)
                call s%widgets%viewport_visibility%setAccessibleName(label)
            end if
        end if
    end procedure sync_patch_visibility_button

    module procedure multipatch_members_shown
        integer::i

        associate(definition=>s%document%multipatches%groups(s%selected_multipatch)%definition)
            visible=all([(studio_scene_is_visible(s%document%scene,definition%kind,definition%members(i),catalog=s%document%patches), &
                i=1,size(definition%members))])
        end associate
    end procedure multipatch_members_shown

    module procedure on_tree_selection
        integer::kind,variant,group,branch_kind,mesh
        logical::was_updating

        select type(s=>data)
        type is(studio_t)
            if(s%updating)return
            mesh=mesh_project_id(s%widgets%project_browser,item)
            if(mesh>0)then
                call select_project_mesh(s,mesh)
                return
            end if
            call project_item_identity(s%widgets%project_browser,s,item,kind,variant,group,branch_kind)
            if(group>0)then
                call select_multipatch_group(s,group,reveal=.false.)
            else if(kind>0)then
                call select_patch(s,kind,variant,reveal=.false.)
                ! Preserve the actual browser row, including a member alias.
                was_updating=s%updating
                s%updating=.true.
                call s%widgets%project_browser%tree%setCurrentItem(item)
                s%updating=was_updating
            else if(branch_kind/=0)then
                call clear_workspace_selection(s)
                return
            end if
            call sync_project_delete_controls(s)
        end select
    end procedure on_tree_selection

    module procedure on_tree_click
        type(selection_target_t) :: target

        select type(s=>data)
        type is(studio_t)
            if(s%updating.or.column==scene_visibility_column)return
            target=tree_selection_target(s,item)
            if(target%category>0)call apply_workspace_selection(s,target,selection_replace)
        end select
    end procedure on_tree_click

    module procedure on_project_key
        integer :: command

        if(event%key()==Qt_Key_Escape)then
            select type(s=>data)
            type is(studio_t)
                call clear_workspace_selection(s)
                call event%accept()
            end select
            return
        end if
        if(event%key()==Qt_Key_A.and.event%modifiers()==Qt_ControlModifier)then
            select type(s=>data)
            type is(studio_t)
                if(s%updating)return
                call select_displayed_scene_rows(s)
                call event%accept()
            end select
            return
        end if
        if(event%modifiers()/=Qt_NoModifier)return
        select case(event%key())
        case(Qt_Key_Delete)
            command=project_remove_command
        case(Qt_Key_F2)
            command=project_rename_command
        case default
            return
        end select
        call event%accept()
        if(event%isAutoRepeat())return
        select type(s=>data)
        type is(studio_t)
            if(s%updating)return
            if(command==project_rename_command)then
                if(.not.s%widgets%project_rename_action%isEnabled())return
            end if
            call invoke_command(s,command)
        end select
    end procedure on_project_key

    subroutine select_displayed_scene_rows(s)
        type(studio_t), intent(inout), target :: s
        type(selection_target_t),allocatable :: targets(:)
        type(selection_target_t) :: focused
        type(QTreeWidgetItem) :: row
        integer :: i,n

        allocate(targets(s%widgets%project_browser%row_count+1))
        n=0
        do i=0,s%widgets%project_browser%tree%topLevelItemCount()-1
            row=s%widgets%project_browser%tree%topLevelItem(i)
            call add_displayed_rows(row)
        end do
        if(n==0)return
        focused=tree_selection_target(s,s%widgets%project_browser%tree%currentItem())
        if(any(targets(:n)%category==focused%category.and.targets(:n)%rank==focused%rank.and.targets(:n)%id==focused%id))then
            n=n+1
            targets(n)=focused
        end if
        call apply_workspace_selection_batch(s,targets(:n),selection_replace)
        call s%status%showMessage(integer_text(s%selection%count())//' selected',3000)
    contains
        recursive subroutine add_displayed_rows(item)
            type(QTreeWidgetItem), intent(in) :: item
            type(QTreeWidgetItem) :: child
            type(selection_target_t) :: target
            integer :: child_index

            if(item%isHidden())return
            target=tree_selection_target(s,item)
            if(target%category>0)then
                n=n+1
                targets(n)=target
            end if
            if(.not.item%isExpanded())return
            do child_index=0,item%childCount()-1
                child=item%child(child_index)
                call add_displayed_rows(child)
            end do
        end subroutine add_displayed_rows
    end subroutine select_displayed_scene_rows

    module procedure on_project_mouse_press
        type(QTreeWidgetItem) :: item
        type(QHeaderView) :: header
        type(selection_target_t) :: target
        integer :: mode,visibility,kind,variant,group,branch

        if(event%button()==Qt_RightButton)then
            call event%accept()
            return
        end if
        if(event%button()/=Qt_LeftButton)return
        select type(s=>data)
        type is(studio_t)
            if(s%updating)return
            mode=pointer_selection_mode(event%modifiers())
            item=s%widgets%project_browser%tree%itemAt(int(event%x()),int(event%y()))
            if(item%isNull())then
                if(mode==selection_replace)call clear_workspace_selection(s)
                return
            end if
            header=s%widgets%project_browser%tree%header()
            if(header%logicalIndexAt(int(event%x()))==scene_visibility_column)then
                call event%accept()
                visibility=project_row_visibility(s%widgets%project_browser,item)
                if(visibility<0)return
                target=tree_selection_target(s,item)
                if(target%category>0)then
                    call set_scene_target_visibility(s,target,visibility/=1)
                    call sync_scene_visibility(s)
                else
                    call project_item_identity(s%widgets%project_browser,s,item,kind,variant,group,branch)
                    call set_project_branch_visibility(s,branch,visibility/=1)
                end if
                return
            end if
            if(mode==selection_replace)return
            target=tree_selection_target(s,item)
            if(target%category==0)return
            call apply_workspace_selection(s,target,mode)
            call event%accept()
        end select
    end procedure on_project_mouse_press

    module procedure on_project_context
        type(QPoint)::at
        type(QRect)::bounds
        type(QTreeWidgetItem)::item
        type(selection_target_t)::target,focused

        select type(s=>data)
        type is(studio_t)
            if(s%updating)return
            at=event%globalPos()
            if(event%reason()==QContextMenuEvent_Keyboard)then
                item=s%widgets%project_browser%tree%currentItem()
                if(.not.item%isNull())then
                    bounds=s%widgets%project_browser%tree%visualItemRect(item)
                    if(.not.bounds%isEmpty())at=s%widgets%project_browser%viewport%mapToGlobal(bounds%center())
                end if
            else
                item=s%widgets%project_browser%tree%itemAt(s%widgets%project_browser%viewport%mapFromGlobal(at))
            end if
            target=tree_selection_target(s,item)
            if(s%selection%count()==1)then
                focused=s%selection%items(1)
            else if(.not.allocated(s%selection%items))then
                focused=primary_workspace_target(s)
            end if
            select case(target%category)
            case(selection_patch)
                if((focused%category==selection_control.or.focused%category==selection_patch).and. &
                    focused%item>0.and.focused%rank==target%rank.and.focused%id==target%id)target=focused
            case(selection_mesh)
                if(focused%category==selection_mesh.and.focused%id==target%id.and.focused%item>0.and. &
                    scene_mesh_drawn(s,target%id).and.s%selected_mesh==target%id)then
                    if(all(s%mesh_cell(:2)==[focused%part,focused%item]))target=focused
                end if
            case default
                continue
            end select
            call event%accept()
            call show_scene_context(s,target,at,item)
        end select
    end procedure on_project_context

    module procedure show_scene_context
        integer, parameter :: select_target=-1,select_owner=-2
        type(QMenu)::menu,submenu
        type(QAction)::choice,heading
        type(QAction),allocatable::items(:)
        type(QFontMetrics)::metrics
        type(selection_set_t)::targets,menu_selection
        integer,allocatable::commands(:)
        integer::i,command,kind,id,group,branch,row_kind,row_id,row_group,profiles,curves
        integer(int64)::target_revision,tree_stamp,mesh_generation,document_revision
        logical::from_tree,multiple,subfeatures,drawn,show,can_restructure,source_available,can_remesh
        character(:),allocatable::title,object_name,action_prefix,feature,cell_label

        from_tree=present(item)
        branch=0
        if(from_tree)call project_item_identity(s%widgets%project_browser,s,item,row_kind,row_id,row_group,branch)
        if(from_tree)then
            menu=QMenu('Scene',s%widgets%project_browser%tree)
            call register_status_menu(menu)
            call menu%setObjectName('studioProjectContextMenu')
            action_prefix='studioProjectContextAction'
        else
            menu=QMenu('Viewport',s%surface)
            call register_status_menu(menu)
            call menu%setObjectName('studioViewportContext')
            action_prefix='studioViewportContext'
        end if
        allocate(items(0),commands(0))
        heading=QAction('',menu)
        call heading%setObjectName('studioContextTarget')
        call heading%setEnabled(.false.)
        call menu%addAction(heading)
        call menu%addSeparator()

        kind=target%rank
        id=target%id
        tree_stamp=s%widgets%project_browser%stamp
        mesh_generation=s%document%mesh_generation
        document_revision=s%document%document_revision
        target_revision=-1_int64
        drawn=scene_target_drawn(s,target)
        subfeatures=.false.
        if(from_tree.and.s%selection%count()>1.and.target%item==0)then
            select case(target%category)
            case(selection_patch)
                subfeatures=all(s%selection%items%category==selection_control.and. &
                    s%selection%items%rank==kind.and.s%selection%items%id==id)
            case(selection_mesh)
                subfeatures=all(s%selection%items%category==selection_mesh.and. &
                    s%selection%items%id==id.and.s%selection%items%item>0)
            case default
                continue
            end select
        end if
        multiple=s%selection%count()>1.and.(s%selection%contains(target).or.subfeatures).and.target%item==0
        show=.false.
        if(multiple)then
            menu_selection=s%selection
            if(subfeatures)then
                if(target%category==selection_patch)then
                    title=integer_text(menu_selection%count())//' control points selected'
                    if(drawn.and.s%layers(4))call append('Fit control points',fit_selection_command,icon_fit_selection)
                    call append('Select owning patch',select_owner,33)
                else
                    title=integer_text(menu_selection%count())//' mesh cells selected'
                    if(drawn)call append('Fit cells',fit_selection_command,icon_fit_selection)
                    call append('Select mesh',select_owner,icon_mesh)
                end if
            else
                targets=scene_visibility_targets(s)
                title=integer_text(menu_selection%count())//' selected'
                if(any([(scene_target_drawn(s,targets%items(i)),i=1,targets%count())]))then
                    call append('Fit selection',fit_selection_command,icon_fit_selection)
                    call append('Isolate selection',project_isolate_command,13)
                end if
                if(targets%count()>0)call append_visibility('selected objects')
                if(len(document_deletion_error(s%document,menu_selection%items))==0)then
                    call menu%addSeparator()
                    call append('Delete selection',project_remove_command,icon_delete)
                end if
            end if
        else
            select case(target%category)
            case(selection_control,selection_patch)
                target_revision=single_patch_revision(kind,id,catalog=s%document%patches)
                object_name=single_patch_name(kind,id,catalog=s%document%patches)
                if(target%category==selection_control)then
                    feature=integer_text(target%item)
                    if(kind>1.and.s%live_scene%valid.and.s%live_scene%variant==id.and.all(s%live_scene%control_counts>0))then
                        feature='('//control_index_text(target%item,s%live_scene%control_counts,kind)//')'
                    end if
                    title='Control point '//feature//' · '//object_name
                    call append('Select control point',select_target,18)
                    call append('Edit control point…',edit_controls_command,22,scoped=.true.)
                    if(s%layers(4).and.scene_patch_drawn(s,kind,id))call append('Fit control point',fit_selection_command,icon_fit_selection)
                    call menu%addSeparator()
                    call append('Select owning patch',select_owner,33)
                    if(.not.from_tree)call append('Reveal in Scene',scene_reveal_command,2)
                else if(target%feature>0.and.target%item>0)then
                    feature='Face'
                    if(target%feature==pick_edge)feature='Edge'
                    if(target%feature==pick_vertex)feature='Vertex'
                    title=feature//' '//integer_text(target%item)//' · '//object_name
                    call append('Select '//feature,select_target,18)
                    if(kind==3.and.target%feature==pick_face)call append('Extract face',extract_face_command,icon_extract_face)
                    if(target%feature/=pick_vertex)call append('Measure '//feature//'…',measure_geometry_command,icon_measure,scoped=.true.)
                    if(drawn)call append('Fit owning patch',fit_selection_command,icon_fit_selection)
                    call menu%addSeparator()
                    call append('Select owning patch',select_owner,33)
                    if(.not.from_tree)call append('Reveal in Scene',scene_reveal_command,2)
                else
                    title=entity_name(kind)//' · '//object_name
                    if(drawn)call append('Fit object',fit_selection_command,icon_fit_selection)
                    if(.not.from_tree)call append('Reveal in Scene',scene_reveal_command,2)
                    call append('Inspect parameters',inspect_patch_command,33)
                    call append('Measure…',measure_geometry_command,icon_measure)
                    call menu%addSeparator()
                    submenu=context_submenu('Transform',9)
                    call append('Move',translate_geometry_command,9,submenu)
                    call append('Rotate',rotate_geometry_command,60,submenu)
                    call append('Scale',scale_geometry_command,61,submenu)
                    submenu=context_submenu('Edit patch',22)
                    call append('Control points and weights…',edit_controls_command,22,submenu)
                    can_restructure=.not.document_patch_in_assembly(s%document,kind,id)
                    if(can_restructure)call append('Split at parameter…',split_patch_command,90,submenu)
                    call append('Extract parameter range…',extract_patch_command,91,submenu)
                    if(kind>=2)call append('Extract section…',extract_section_command,icon_extract_section,submenu)
                    if(kind==3)call append('Extract face…',extract_face_command,icon_extract_face,submenu)
                    if(can_restructure)call append('Reverse parameter direction',reverse_patch_command,92,submenu)
                    if(kind<=2)then
                        submenu=context_submenu('Create from profile',77)
                        curves=context_patch_count(s,1)
                        profiles=curves
                        if(kind==2)profiles=context_patch_count(s,2)
                        do i=1,size(construction_names)
                            if(construction_command+i-1==sweep_command.and.curves<merge(2,1,kind==1))cycle
                            if(construction_command+i-1==loft_command.and.profiles<2)cycle
                            call append(trim(construction_names(i))//'…',construction_command+i-1,76+i,submenu)
                        end do
                    end if
                    call append('Rename patch…',project_rename_command,68)
                    call append('Duplicate patch',project_duplicate_command,11)
                    call append('Generate mesh…',open_mesh_dialog_command,icon_generate_mesh)
                    call menu%addSeparator()
                    call append('Object color…',color_command,37)
                    targets=scene_visibility_targets(s,target)
                    call append_visibility('object')
                    if(.not.show)call append('Isolate object',project_isolate_command,13,scoped=.true.)
                    if(can_delete_document_patch(s%document,kind,id))then
                        call menu%addSeparator()
                        call append('Delete patch',project_remove_command,icon_delete)
                    end if
                end if
            case(selection_mesh)
                if(id<1.or.id>studio_mesh_count(s))then
                    call menu%delete()
                    return
                end if
                object_name=s%document%meshes(id)%name
                if(target%item>0)then
                    cell_label=studio_mesh_cell_label(s,target)
                    if(len(cell_label)==0.or..not.drawn)then
                        call menu%delete()
                        return
                    end if
                    title=cell_label//' · '//object_name
                    call append('Select cell',select_target,18)
                    call append('Fit cell',fit_selection_command,icon_fit_selection)
                    call menu%addSeparator()
                    call append('Select mesh',select_owner,icon_mesh)
                    if(.not.from_tree)call append('Reveal in Scene',scene_reveal_command,2)
                    if(mesh_source_available_in_scene(s,id))call append('Select mesh source',scene_select_mesh_source_command,33)
                else
                    title='Mesh · '//object_name
                    if(drawn)call append('Fit mesh',fit_selection_command,icon_fit_selection)
                    if(studio_mesh_valid(s,id))then
                        call append('Inspector',show_inspector_command,40)
                        if(drawn.and..not.s%mesh_pick)call append('Pick mesh cells',mesh_pick_command,icon_mesh_pick)
                    end if
                    source_available=mesh_source_available_in_scene(s,id)
                    if(source_available)call append('Select mesh source',scene_select_mesh_source_command,33)
                    group=studio_mesh_source_group(s,id)
                    can_remesh=source_available
                    if(group>0)can_remesh=allocated(s%document%multipatches%groups(group)%native)
                    if(studio_mesh_valid(s,id))then
                        call menu%addSeparator()
                        targets=scene_visibility_targets(s,target)
                        call append_visibility('mesh')
                        if(.not.show)call append('Isolate mesh',project_isolate_command,13,scoped=.true.)
                        call append('Export mesh (VTK)…',export_mesh_command,icon_export)
                    end if
                    if(can_remesh)then
                        call menu%addSeparator()
                        call append('Remesh…',open_mesh_dialog_command,icon_generate_mesh,scoped=.true.)
                    end if
                    call menu%addSeparator()
                    call append('Delete mesh',project_remove_command,icon_delete)
                end if
            case(selection_assembly)
                group=id
                if(group<1.or.group>count_multipatches(s%document%multipatches))then
                    call menu%delete()
                    return
                end if
                target_revision=multipatch_undo_stamp(s%document%multipatches,catalog=s%document%patches)
                title='Multipatch · '//s%document%multipatches%groups(group)%definition%name
                if(target%item>0)then
                    title='Control point '//integer_text(target%item)//' · '//title
                    call append('Select control point',select_target,18)
                    call append('Edit control point…',edit_controls_command,22,scoped=.true.)
                else
                    targets=scene_visibility_targets(s,target)
                    if(any([(scene_target_drawn(s,targets%items(i)),i=1,targets%count())]))then
                        call append('Fit members',fit_selection_command,icon_fit_selection)
                    end if
                    call append('Inspector',show_inspector_command,40)
                    if(targets%count()>0)then
                        call append_visibility('members')
                        if(any([(scene_target_visible(s,targets%items(i)),i=1,targets%count())]))then
                            call append('Isolate members',project_isolate_command,13,scoped=.true.)
                        end if
                    end if
                    call menu%addSeparator()
                    call append('Edit multipatch…',edit_multipatch_command,88)
                    call append('Check interfaces',validate_multipatch_command,89)
                    if(allocated(s%document%multipatches%groups(group)%native))then
                        call append('Control points and weights…',edit_controls_command,22)
                        call append('Generate mesh…',open_mesh_dialog_command,icon_generate_mesh)
                    end if
                    call menu%addSeparator()
                    call append('Dissolve multipatch',remove_multipatch_command,36)
                end if
            case default
                if(from_tree)then
                    title='Scene'
                    select case(branch)
                    case(-3)
                        title='Meshes'
                        if(studio_mesh_count(s)>0)then
                            call append('Show all meshes',project_show_branch,31)
                            call append('Hide all meshes',project_hide_branch,32)
                        end if
                    case(-1)
                        title='Multipatches'
                        call append('Create multipatch…',new_multipatch_command,87)
                    case(-4,-2,1:3)
                        title='Geometry folder'
                        if(branch>0)title=entity_name(branch)//' folder'
                        call append('Show this folder',project_show_branch,31)
                        call append('Hide this folder',project_hide_branch,32)
                    case default
                        call append('Show all',project_show_all_command,31)
                    end select
                    if(size(items)>0)call menu%addSeparator()
                    call append('Expand all',project_expand_command,30)
                    call append('Collapse all',project_collapse_command,36)
                else
                    title='Viewport'
                    call append('Create geometry…',primitive_command,30)
                    call menu%addSeparator()
                    call append('Fit scene',fit_scene_command,19)
                    call append('Fit around panels',fit_around_panels_command,icon_fit_panels)
                    call append('Isometric view',isometric_view_command,25)
                end if
            end select
        end if
        if(s%isolation%count()>0)then
            if(size(items)>0)call menu%addSeparator()
            call append('Exit isolation',scene_exit_isolation_command,icon_restore)
        end if
        if(size(items)==0)then
            call menu%delete()
            return
        end if
        call menu%setAccessibleName(title)
        metrics=QFontMetrics(menu%font())
        call heading%setText(menu_title(metrics%elidedText(title,Qt_ElideRight,260)))
        call metrics%delete()
        choice=menu%exec(at)
        command=0
        do i=1,size(items)
            if(.not.c_associated(choice%ptr,items(i)%ptr))cycle
            command=commands(i)
            exit
        end do
        call menu%delete()
        if(command==0)return
        if(document_revision/=s%document%document_revision.or.tree_stamp/=s%widgets%project_browser%stamp.or. &
            mesh_generation/=s%document%mesh_generation)then
            call s%status%showMessage('Scene changed. Open the context menu again.',3500)
            return
        end if
        if(target%category==selection_patch.or.target%category==selection_control)then
            if(.not.multiple.and.single_patch_revision(kind,id,catalog=s%document%patches)/=target_revision)return
        else if(target%category==selection_assembly)then
            if(.not.multiple.and.target_revision/=multipatch_undo_stamp(s%document%multipatches,catalog=s%document%patches))return
        end if
        if(target%category==selection_mesh.and.target%item>0)then
            if(.not.scene_mesh_drawn(s,id))return
        end if
        if(command==scene_exit_isolation_command)then
            call invoke_command(s,command)
            return
        end if
        if(multiple)then
            if(s%selection%count()/=menu_selection%count().or. &
                .not.all([(s%selection%contains(menu_selection%items(i)),i=1,menu_selection%count())]))then
                call s%status%showMessage('Selection changed. Open the context menu again.',3500)
                return
            end if
            select case(command)
            case(toggle_patch_visibility_command)
                call set_scene_visibility(s,targets,show)
                call sync_scene_visibility(s)
            case(select_owner)
                call apply_workspace_selection(s,target,selection_replace)
            case default
                call invoke_command(s,command)
            end select
            return
        end if
        if(command==project_show_branch.or.command==project_hide_branch)then
            call set_project_branch_visibility(s,branch,command==project_show_branch)
            return
        end if
        if(target%category==0)then
            call invoke_command(s,command)
            return
        end if
        if(target%category==selection_assembly.and.target%item>0)then
            call apply_workspace_selection(s,selection_target_t(selection_assembly,id),selection_replace)
            if(.not.s%shared_controls%valid)return
            if(target%item>size(s%shared_controls%weights))return
            s%selected_control=target%item
            s%selected_shared_controls=[target%item]
            call on_control_selection(s)
            call s%surface%update()
            if(command==edit_controls_command)call invoke_command(s,command)
            return
        end if
        select case(command)
        case(fit_selection_command)
            if(fit_workspace_selection(s,[target]))call repaint_all(s)
        case(project_isolate_command)
            call isolate_scene_targets(s,target)
        case(toggle_patch_visibility_command)
            call set_scene_visibility(s,targets,show)
            call sync_scene_visibility(s)
        case(scene_reveal_command)
            call reveal_scene_selection(s,target)
        case(select_target)
            call apply_workspace_selection(s,target,selection_replace)
        case(select_owner)
            if(target%category==selection_mesh)then
                call apply_workspace_selection(s,selection_target_t(selection_mesh,id),selection_replace)
            else
                call apply_workspace_selection(s,selection_target_t(selection_patch,id,kind),selection_replace)
            end if
        case default
            call apply_workspace_selection(s,target,selection_replace)
            if(command==mesh_pick_command.and.s%mesh_pick)return
            call invoke_command(s,command)
        end select

    contains

        function context_submenu(label,icon_number) result(group_menu)
            character(*),intent(in)::label
            integer,intent(in)::icon_number
            type(QMenu)::group_menu
            type(QAction)::group_action

            group_menu=menu%addMenu(label)
            call register_status_menu(group_menu)
            group_action=group_menu%menuAction()
            call group_action%setIcon(line_icon(icon_number))
        end function context_submenu

        subroutine append(label,id,icon_number,host,scoped)
            character(*),intent(in)::label
            integer,intent(in)::id,icon_number
            type(QMenu),intent(in),optional::host
            logical,intent(in),optional::scoped
            type(QMenu)::parent
            type(QAction)::action
            character(:),allocatable::caption
            integer::glyph

            parent=menu
            if(present(host))parent=host
            associate(visual=>command_visual(id,label,icon_number))
                caption=visual%caption
                glyph=visual%glyph
            end associate
            if(present(scoped))then
                if(scoped)caption=label
            end if
            action=QAction(caption,parent)
            call action%setObjectName(action_prefix//integer_text(id))
            call action%setIcon(line_icon(glyph))
            if(id>0)call set_status_hint(action,studio_command_tooltip(id,caption))
            call parent%addAction(action)
            items=[items,action]
            commands=[commands,id]
        end subroutine append

        subroutine append_visibility(noun)
            character(*),intent(in)::noun

            show=.not.all([(scene_target_visible(s,targets%items(i)),i=1,targets%count())])
            if(show)then
                call append('Show '//noun,toggle_patch_visibility_command,31,scoped=.true.)
            else
                call append('Hide '//noun,toggle_patch_visibility_command,32,scoped=.true.)
            end if
        end subroutine append_visibility

    end procedure show_scene_context

    integer function context_patch_count(s,kind) result(counted)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: kind
        integer :: id

        counted=count([(document_patch(s,kind,id),id=1,single_patch_count(kind,catalog=s%document%patches))])
    end function context_patch_count

    function menu_title(title) result(label)
        character(*), intent(in) :: title
        character(:), allocatable :: label
        integer :: i

        label=''
        do i=1,len(title)
            if(title(i:i)=='&')label=label//'&'
            label=label//title(i:i)
        end do
    end function menu_title

    subroutine sync_project_delete_controls(s)
        type(studio_t),intent(in)::s
        type(selection_target_t)::target
        logical::allowed
        character(:),allocatable::label,tip

        if(s%selection%count()>1)then
            label='Delete selection'
            tip=document_deletion_error(s%document,s%selection%items)
            allowed=len(tip)==0
            if(allowed)tip='Del deletes selected objects in one edit. Ctrl+Z restores them.'
            call s%widgets%project_remove_action%setText(label)
            call s%widgets%project_remove_action%setIconText('Delete')
            call sync_delete_glyph(s%widgets%project_remove_action,.false.)
            call set_status_hint(s%widgets%project_remove_action,tip)
            call s%widgets%project_remove_action%setEnabled(allowed)
            return
        end if

        target=primary_workspace_target(s)
        allowed=.true.
        select case(target%category)
        case(selection_mesh)
            label='Delete mesh'
            tip='Delete this mesh and keep its source geometry'
        case(selection_assembly)
            label='Dissolve multipatch'
            tip='Remove this multipatch and keep its members · Ctrl+Z to restore'
        case(selection_patch)
            label='Delete patch'
            allowed=can_delete_document_patch(s%document,target%rank,target%id)
            tip='Delete this patch · Ctrl+Z to restore'
            if(target%id<=2)then
                tip='Built-in templates remain available'
            else if(.not.allowed)then
                tip='Detach this patch from its assembly before deleting it'
            end if
        case default
            label='Delete selection'
            tip='Select a patch, assembly, or mesh to delete it'
            allowed=.false.
        end select
        if(allocated(s%selection%items))allowed=allowed.and.s%selection%count()==1
        if(project_target_is_subfeature(target))then
            allowed=.false.
            tip='Select the whole patch or mesh to delete it'
        end if
        call s%widgets%project_remove_action%setText(label)
        if(target%category==selection_assembly)then
            call s%widgets%project_remove_action%setIconText('Dissolve')
        else
            call s%widgets%project_remove_action%setIconText('Delete')
        end if
        call sync_delete_glyph(s%widgets%project_remove_action,target%category==selection_assembly)
        if(allowed)tip=tip//' · Del'
        call set_status_hint(s%widgets%project_remove_action,tip)
        call s%widgets%project_remove_action%setEnabled(allowed)
    end subroutine sync_project_delete_controls

    subroutine sync_delete_glyph(action,dissolving)
        type(QAction),intent(in) :: action
        logical,intent(in) :: dissolving
        type(QIcon) :: glyph

        if(action%boolProperty('dissolveMultipatch').eqv.dissolving)return
        call action%setProperty('dissolveMultipatch',dissolving)
        if(dissolving)then
            associate(visual=>command_visual(remove_multipatch_command,'Dissolve multipatch',36))
                glyph=line_icon(visual%glyph)
            end associate
        else
            glyph=line_icon(icon_delete)
        end if
        call action%setIcon(glyph)
        call glyph%delete()
    end subroutine sync_delete_glyph

    module procedure run_project_delete
        type(selection_target_t) :: target
        type(selection_target_t), allocatable :: targets(:)
        character(:),allocatable :: error
        integer :: removed

        if(s%selection%count()>1)then
            targets=s%selection%items
            call delete_project_selection(s,targets,removed,error)
            if(len(error)>0)call s%status%showMessage(error,5000,warning=.true.)
            return
        end if
        if(allocated(s%selection%items))then
            if(s%selection%count()==0)return
        end if
        target=primary_workspace_target(s)
        if(project_target_is_subfeature(target))then
            call s%status%showMessage('Select the whole patch or mesh to delete it',3500)
            return
        end if
        select case(target%category)
        case(selection_mesh)
            call delete_project_mesh(s,target%id,error)
        case(selection_assembly)
            call invoke_command(s,remove_multipatch_command)
            return
        case(selection_patch)
            call delete_project_patch(s,target%rank,target%id,error)
        case default
            return
        end select
        if(len(error)>0)call s%status%showMessage(error,5000,warning=.true.)
    end procedure run_project_delete

    module procedure delete_project_selection

        call delete_document_selection(s%document,targets,s%mesh_only,removed,error)
        if(len(error)>0)return
        call refresh_project_merge_view(s)
        call s%status%showMessage('Deleted '//integer_text(removed)//' objects · Ctrl+Z to restore',5000)
    end procedure delete_project_selection

    logical function project_target_is_subfeature(target) result(has_subfeature)
        type(selection_target_t), intent(in) :: target

        has_subfeature=target%category==selection_control.or.target%item>0.or.target%feature>0
    end function project_target_is_subfeature

    module procedure delete_project_patch
        character(:),allocatable :: name
        integer :: i

        call delete_document_patch(s%document,kind,variant,name,error)
        if(len(error)>0)return
        if(.not.has_visible_meshes(s))s%mesh_only=.false.
        call clamp_mesh_layer(s)
        if(s%entity_kind==kind.and.s%live_scene%variant==variant)then
            do i=single_patch_count(kind,catalog=s%document%patches),1,-1
                if(.not.studio_scene_is_visible(s%document%scene,kind,i,catalog=s%document%patches))cycle
                call set_active_patch(kind,i,catalog=s%document%patches)
                exit
            end do
        end if
        s%live_scene%valid=.false.
        call rebuild_project_tree(s)
        if(s%selected_mesh>0)then
            call select_project_mesh(s,int(s%selected_mesh))
        else if(s%selected_multipatch>0)then
            call select_multipatch_group(s,int(s%selected_multipatch),reveal=.false.)
        else
            call select_patch(s,s%entity_kind,reveal=.false.)
        end if
        call s%status%showMessage('Deleted '//name//' · Ctrl+Z to restore',5000)
    end procedure delete_project_patch

    module procedure delete_project_assembly

        call delete_document_assembly(s%document,id,error)
        if(len(error)>0)return
        if(s%selected_multipatch==id)s%selected_multipatch=0
        if(s%selected_multipatch>id)s%selected_multipatch=s%selected_multipatch-1
        if(.not.has_visible_meshes(s))s%mesh_only=.false.
        call clamp_mesh_layer(s)
        call rebuild_project_tree(s)
        if(s%selected_mesh>0)then
            call select_project_mesh(s,int(s%selected_mesh))
        else if(s%selected_multipatch>0)then
            call select_multipatch_group(s,int(s%selected_multipatch),reveal=.false.)
        else
            call select_patch(s,s%entity_kind,reveal=.false.)
        end if
        call s%status%showMessage('Multipatch dissolved · Member patches kept · Ctrl+Z to restore',5000)
    end procedure delete_project_assembly

    module procedure delete_project_mesh
        type(deleted_mesh_t) :: removed
        logical :: was_selected

        was_selected=s%selected_mesh==id
        call delete_document_mesh(s%document,id,s%selected_mesh,s%mesh_only,removed,error)
        if(len(error)>0)return
        call adjust_mesh_selection_after_removal(s,id)
        if(.not.has_visible_meshes(s))s%mesh_only=.false.
        call clamp_mesh_layer(s)
        call rebuild_project_tree(s)
        if(s%selected_mesh>0)then
            call select_project_mesh(s,int(s%selected_mesh))
        else if(.not.was_selected.and.s%selected_multipatch>0)then
            call select_multipatch_group(s,int(s%selected_multipatch),reveal=.false.)
        else if(was_selected.and.removed%source_group>0)then
            call select_multipatch_group(s,removed%source_group,reveal=.false.)
        else if(was_selected.and.document_patch(s,removed%kind,removed%variant))then
            call select_patch(s,removed%kind,removed%variant,reveal=.false.)
        else if(.not.was_selected.and.s%live_scene%valid)then
            call select_patch(s,s%entity_kind,reveal=.false.)
        else
            call pick_document_selection(s)
        end if
        call s%status%showMessage('Deleted '//removed%name//' · Ctrl+Z to restore',5000)
    end procedure delete_project_mesh

    subroutine refresh_restored_mesh(s,selection,message)
        type(studio_t),intent(inout),target :: s
        integer,intent(in) :: selection
        character(*),intent(in) :: message

        if(.not.has_visible_meshes(s))s%mesh_only=.false.
        call clamp_mesh_layer(s)
        call rebuild_project_tree(s)
        if(selection>0)then
            call select_project_mesh(s,selection)
        else if(s%selected_multipatch>0)then
            call select_multipatch_group(s,int(s%selected_multipatch),reveal=.false.)
        else if(s%live_scene%valid)then
            call select_patch(s,s%entity_kind,reveal=.false.)
        else
            call pick_document_selection(s)
        end if
        call s%status%showMessage(message,4500)
    end subroutine refresh_restored_mesh

    module procedure undo_project_edit
        type(document_undo_result_t) :: undone
        type(QScrollBar) :: scrollbar
        integer :: scroll,previous_units(2)

        scrollbar=s%widgets%project_browser%tree%verticalScrollBar()
        scroll=scrollbar%value()
        previous_units=[s%document%model_unit,s%document%display_unit]
        call undo_document_edit(s%document,undone)
        if(.not.undone%changed)then
            if(len(undone%message)>0)call s%status%showMessage(undone%message,5000,warning=.true.)
            return
        end if

        select case(undone%source)
        case(undo_mesh)
            s%mesh_only=undone%mesh_only
            call clear_mesh_selection(s)
            call refresh_restored_mesh(s,undone%mesh_id,undone%message)
        case(undo_material)
            call sync_controls(s)
            call repaint_all(s)
            call s%status%showMessage('Appearance restored',3500)
        case(undo_assembly)
            if(undone%removed_mesh_id>0)call adjust_mesh_selection_after_removal(s,undone%removed_mesh_id)
            if(undone%undid_construction.and.s%mesh_only)then
                if(.not.has_visible_meshes(s))s%mesh_only=.false.
            end if
            s%selected_multipatch=undone%assembly_id
            call rebuild_project_tree(s)
            if(undone%assembly_id>0)then
                call select_multipatch_group(s,undone%assembly_id,reveal=.false.)
            else
                call select_entity(s,s%entity_kind)
            end if
            call s%status%showMessage(undone%message,4500)
        case(undo_patch)
            call refresh_restored_patch(s,undone)
        case(undo_deletion)
            s%mesh_only=undone%mesh_only
            call refresh_project_merge_view(s)
            call s%status%showMessage(undone%message,4500)
        case(undo_project_merge)
            call refresh_project_merge_view(s)
            call s%status%showMessage(undone%message,4500)
        case(undo_properties)
            if(any(previous_units/=[s%document%model_unit,s%document%display_unit]))call reset_transform_values(s)
            call rebuild_project_tree(s)
            call sync_scene_visibility(s)
            call s%status%showMessage(undone%message,4500)
        case default
            return
        end select
        call scrollbar%setValue(scroll)
        call refresh_edit_actions(s)
        call sync_project_title(s)
    end procedure undo_project_edit

    module procedure redo_project_edit
        type(QScrollBar) :: scrollbar
        integer :: assembly,mesh,kind,scroll,previous_units(2)
        character(:), allocatable :: error

        scrollbar=s%widgets%project_browser%tree%verticalScrollBar()
        scroll=scrollbar%value()
        previous_units=[s%document%model_unit,s%document%display_unit]
        if(.not.redo_document_edit(s%document,error))then
            if(len(error)>0)call s%status%showMessage(error,5000,warning=.true.)
            return
        end if
        assembly=s%selected_multipatch
        mesh=s%selected_mesh
        kind=s%entity_kind
        call cancel_transform_drag(s)
        call forget_workspace_selection(s)
        s%live_scene%valid=.false.
        call refresh_studio_scene(s%document%scene,catalog=s%document%patches)
        if(s%mesh_only.and..not.has_visible_meshes(s))s%mesh_only=.false.
        call rebuild_project_tree(s)
        if(assembly>0.and.assembly<=count_multipatches(s%document%multipatches))then
            call select_multipatch_group(s,assembly,reveal=.false.)
        else if(mesh>0.and.mesh<=studio_mesh_count(s))then
            call select_project_mesh(s,mesh)
        else if(document_patch(s,kind,active_patch_id(kind,catalog=s%document%patches)))then
            call select_patch(s,kind,reveal=.false.)
        else
            call pick_document_selection(s)
        end if
        call scrollbar%setValue(scroll)
        if(any(previous_units/=[s%document%model_unit,s%document%display_unit]))call reset_transform_values(s)
        call sync_controls(s)
        call refresh_edit_actions(s)
        call sync_project_title(s)
        call repaint_all(s)
        call s%status%showMessage('Project edit redone',3500)
    end procedure redo_project_edit

    subroutine refresh_restored_patch(s,undone)
        type(studio_t),intent(inout),target :: s
        type(document_undo_result_t),intent(in) :: undone
        integer :: kind,assembly_id

        if(undone%patch_kind<1)return
        s%live_scene%valid=.false.
        kind=undone%patch_kind
        assembly_id=undone%assembly_id
        associate(operation=>undone%operation)
            if(operation=='scene rotation'.or.operation=='scene transform')then
                call refresh_rotated_scene(s)
                call s%status%showMessage('Scene transform undone',4500)
                return
            end if
            if(assembly_id>0.and.assembly_id<=count_multipatches(s%document%multipatches))then
                call select_multipatch_group(s,assembly_id,reveal=.false.)
                if(s%transform_mode==0.and.operation=='multipatch transform')call invoke_command(s,fit_selection_command)
                call s%status%showMessage(s%document%multipatches%groups(assembly_id)%definition%name//' · '//operation//' undone',3500)
                return
            end if
            s%selected_multipatch=0
            select case(operation)
            case('patch creation','primitive creation','construction','geometry import','AI construction','section extraction', &
                 'rename','patch removal','split','parameter extraction')
                call rebuild_project_tree(s)
            end select
            if(undone%patch_axis>0)s%knot_axis=undone%patch_axis
            if(operation=='patch removal'.or.operation=='rename')then
                call select_patch(s,kind,reveal=.false.)
            else
                call select_entity(s,kind)
            end if
            if(operation=='patch creation'.or.operation=='primitive creation'.or.operation=='construction'.or. &
                operation=='geometry import'.or.operation=='AI construction')then
                s%zoom=1._dp; s%pan_x=0._dp; s%pan_y=0._dp
                call frame_workspace_scene(s)
                call s%surface%update()
                call s%status%showMessage(operation//' undone',3500)
            else
                if(s%transform_mode==0.and.(operation=='transform'.or.operation=='multipatch transform')) &
                    call invoke_command(s,fit_selection_command)
                call s%status%showMessage(single_patch_name(kind,active_patch_id(kind,catalog=s%document%patches),catalog=s%document%patches)//' · '//operation//' undone',3500)
            end if
        end associate
    end subroutine refresh_restored_patch

    module procedure run_project_edit
        integer::kind,variant,created,suffix,rgb(3),opacity
        integer(int64)::revision
        type(selection_target_t)::target
        character(:),allocatable::name,edited,error
        logical::accepted

        target=primary_workspace_target(s)
        if(target%category/=selection_patch)return
        kind=target%rank
        variant=target%id
        if(.not.document_patch(s,kind,variant))return
        name=single_patch_name(kind,variant,catalog=s%document%patches)
        select case(command)
        case(project_duplicate_command)
            edited=name//' copy'
            suffix=2
            do while(project_name_exists(s,edited))
                edited=name//' copy '//integer_text(suffix)
                suffix=suffix+1
            end do
            call add_single_patch(kind,variant,edited,created,error,catalog=s%document%patches)
            if(len(error)==0)then
                call get_studio_material(s%document%scene,kind,variant,rgb,opacity)
                call refresh_studio_scene(s%document%scene,catalog=s%document%patches)
                call set_studio_material(s%document%scene,kind,created,rgb,opacity)
            end if
        case(project_rename_command)
            revision=single_patch_revision(kind,variant,catalog=s%document%patches)
            edited=QInputDialog_getText('Rename patch','New patch name',text=name,ok=accepted,parent=s%window)
            if(.not.accepted.or.trim(adjustl(edited))==name)return
            if(single_patch_revision(kind,variant,catalog=s%document%patches)/=revision)return
            call rename_project_patch(s,kind,variant,edited,error)
            if(len(error)>0)call s%status%showMessage(error,5000,warning=.true.)
            return
        end select
        if(len(error)>0)then
            call s%status%showMessage(error,5000,warning=.true.)
            return
        end if
        s%live_scene%valid=.false.
        call rebuild_project_tree(s)
        call select_patch(s,kind,reveal=.false.)
        call s%status%showMessage('Duplicated '//name,3500)
    end procedure run_project_edit

    module procedure rename_project_patch
        logical::renamed

        call rename_document_patch(s%document,kind,variant,name,renamed,error)
        if(present(changed))changed=renamed
        if(.not.renamed)return
        ! Keep the selected object and update linked readouts and tree rows.
        s%live_probe_axis=0
        call rebuild_project_tree(s)
        call refresh_inspector(s)
        call sync_controls(s)
        call repaint_all(s)
        call s%status%showMessage('Renamed patch to '//trim(adjustl(name)),3500)
    end procedure rename_project_patch

    logical function project_name_exists(s,name) result(exists)
        type(studio_t),intent(in)::s
        character(*),intent(in)::name
        integer::kind,variant

        exists=.true.
        do kind=1,3
            do variant=1,single_patch_count(kind,catalog=s%document%patches)
                if(single_patch_name(kind,variant,catalog=s%document%patches)==name)return
            end do
        end do
        exists=.false.
    end function project_name_exists

    module procedure set_patch_visibility
        integer :: targets(2,size(variants))
        character(:), allocatable :: message

        targets(1,:)=kind
        targets(2,:)=variants
        call commit_document_visibility(s%document,visible,changed,message,patches=targets)
        if(present(error))error=message
        if(len(message)>0)call s%status%showMessage(message,5000,warning=.true.)
        if(.not.changed)return
        call sync_scene_visibility_batch(s%widgets%project_browser,s,targets)
        call sync_controls(s)
        call sync_scene_visibility(s)
    end procedure set_patch_visibility

    subroutine set_project_branch_visibility(s,branch,visible)
        type(studio_t), intent(inout) :: s
        integer, intent(in) :: branch
        logical, intent(in) :: visible
        type(selection_set_t) :: targets

        targets=project_branch_visibility_targets(s,branch)
        call set_scene_visibility(s,targets,visible)
        call sync_scene_visibility(s)
    end subroutine set_project_branch_visibility

    function project_branch_visibility_targets(s,branch,include_meshes) result(targets)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: branch
        logical, intent(in), optional :: include_meshes
        type(selection_set_t) :: targets, members
        integer :: i,j,kind,variant
        logical :: meshes

        call targets%clear()
        if(branch==-4.or.branch>0)then
            do kind=1,3
                if(branch>0.and.kind/=branch)cycle
                do variant=1,single_patch_count(kind,catalog=s%document%patches)
                    if(.not.document_patch(s,kind,variant))cycle
                    call targets%change(selection_target_t(selection_patch,variant,kind),selection_add)
                end do
            end do
        else if(branch==-1)then
            do i=1,count_multipatches(s%document%multipatches)
                members=scene_visibility_targets(s,selection_target_t(selection_assembly,i))
                do j=1,members%count()
                    call targets%change(members%items(j),selection_add)
                end do
            end do
        end if
        meshes=branch==-3
        if(present(include_meshes))meshes=meshes.or.include_meshes
        if(meshes)then
            do i=1,studio_mesh_count(s)
                if(.not.studio_mesh_valid(s,i))cycle
                call targets%change(selection_target_t(selection_mesh,i),selection_add)
            end do
        end if
    end function project_branch_visibility_targets

    module procedure update_project_visibility
        type(selection_set_t) :: targets
        logical :: accepted

        if(.not.show_all)then
            call isolate_scene_targets(s)
            return
        end if
        targets=project_branch_visibility_targets(s,-4,include_meshes=.true.)
        call set_scene_visibility(s,targets,.true.,accepted)
        if(.not.accepted)return
        call s%isolation%clear()
        if(allocated(s%isolation_mesh_names))deallocate(s%isolation_mesh_names)
        s%mesh_only=.false.
        call sync_scene_visibility(s)
    end procedure update_project_visibility

    function scene_visibility_targets(s,only) result(targets)
        type(studio_t), intent(in) :: s
        type(selection_target_t), intent(in), optional :: only
        type(selection_set_t) :: targets
        type(selection_target_t),allocatable :: candidates(:)
        type(selection_target_t) :: focused
        integer :: i,n,capacity

        n=0
        if(present(only))then
            allocate(candidates(target_capacity(only)))
            call add_target(only)
        else if(allocated(s%selection%items))then
            capacity=0
            do i=1,s%selection%count()
                capacity=capacity+target_capacity(s%selection%items(i))
            end do
            allocate(candidates(capacity))
            do i=1,s%selection%count()
                call add_target(s%selection%items(i))
            end do
        else
            focused=primary_workspace_target(s)
            allocate(candidates(target_capacity(focused)))
            call add_target(focused)
        end if
        call change_selection_batch(targets,candidates(:n),selection_replace)
    contains
        integer function target_capacity(target) result(length)
            type(selection_target_t), intent(in) :: target

            length=1
            if(target%category/=selection_assembly)return
            if(target%id<1.or.target%id>count_multipatches(s%document%multipatches))return
            length=size(s%document%multipatches%groups(target%id)%definition%members)
        end function target_capacity

        subroutine add_target(target)
            type(selection_target_t), intent(in) :: target
            integer :: i

            select case(target%category)
            case(selection_patch,selection_control)
                if(document_patch(s,target%rank,target%id))then
                    n=n+1
                    candidates(n)=selection_target_t(selection_patch,target%id,target%rank)
                end if
            case(selection_mesh)
                if(studio_mesh_valid(s,target%id))then
                    n=n+1
                    candidates(n)=selection_target_t(selection_mesh,target%id)
                end if
            case(selection_assembly)
                if(target%id<1.or.target%id>count_multipatches(s%document%multipatches))return
                associate(definition=>s%document%multipatches%groups(target%id)%definition)
                    do i=1,size(definition%members)
                        if(.not.document_patch(s,definition%kind,definition%members(i)))cycle
                        n=n+1
                        candidates(n)=selection_target_t(selection_patch,definition%members(i),definition%kind)
                    end do
                end associate
            case default
                continue
            end select
        end subroutine add_target
    end function scene_visibility_targets

    logical function scene_target_visible(s,target) result(visible)
        type(studio_t), intent(in) :: s
        type(selection_target_t), intent(in) :: target

        visible=.false.
        select case(target%category)
        case(selection_patch)
            visible=studio_scene_is_visible(s%document%scene,target%rank,target%id,catalog=s%document%patches)
        case(selection_mesh)
            if(studio_mesh_valid(s,target%id))visible=s%document%meshes(target%id)%visible
        case default
            continue
        end select
    end function scene_target_visible

    logical function scene_target_drawn(s,target) result(drawn)
        type(studio_t), intent(in) :: s
        type(selection_target_t), intent(in) :: target

        drawn=.false.
        select case(target%category)
        case(selection_patch)
            drawn=scene_patch_drawn(s,target%rank,target%id)
        case(selection_mesh)
            drawn=scene_mesh_drawn(s,target%id)
        case default
            continue
        end select
    end function scene_target_drawn

    subroutine set_scene_target_visibility(s,target,visible)
        type(studio_t), intent(inout) :: s
        type(selection_target_t), intent(in) :: target
        logical, intent(in) :: visible
        type(selection_set_t) :: targets

        targets=scene_visibility_targets(s,target)
        call set_scene_visibility(s,targets,visible)
    end subroutine set_scene_target_visibility

    subroutine set_scene_visibility(s,targets,visible,accepted)
        type(studio_t), intent(inout) :: s
        type(selection_set_t), intent(in) :: targets
        logical, intent(in) :: visible
        logical, intent(out), optional :: accepted
        integer, allocatable :: patches(:,:),meshes(:)
        integer :: i,npatches,nmeshes
        logical :: changed
        character(:), allocatable :: error

        allocate(patches(2,targets%count()),meshes(targets%count()))
        npatches=0
        nmeshes=0
        do i=1,targets%count()
            associate(target=>targets%items(i))
                select case(target%category)
                case(selection_patch)
                    npatches=npatches+1
                    patches(:,npatches)=[target%rank,target%id]
                case(selection_mesh)
                    nmeshes=nmeshes+1
                    meshes(nmeshes)=target%id
                case default
                    continue
                end select
            end associate
        end do
        call commit_document_visibility(s%document,visible,changed,error,patches(:,:npatches),meshes(:nmeshes))
        if(present(accepted))accepted=len(error)==0
        if(len(error)>0)call s%status%showMessage(error,5000,warning=.true.)
        if(.not.changed)return
        call sync_scene_visibility_batch(s%widgets%project_browser,s,patches(:,:npatches))
    end subroutine set_scene_visibility

    subroutine sync_scene_visibility(s)
        type(studio_t), intent(inout), target :: s

        if(.not.has_visible_meshes(s))then
            s%mesh_only=.false.
            s%mesh_pick=.false.
        end if
        if(.not.scene_mesh_drawn(s,s%selected_mesh))s%mesh_cell=0
        call refresh_inspector(s)
        call repaint_all(s)
    end subroutine sync_scene_visibility

    module procedure toggle_scene_selection_visibility
        type(selection_set_t) :: targets
        logical :: show
        integer :: i

        targets=scene_visibility_targets(s)
        if(targets%count()==0)return
        show=.not.all([(scene_target_visible(s,targets%items(i)),i=1,targets%count())])
        call set_scene_visibility(s,targets,show)
        call sync_scene_visibility(s)
    end procedure toggle_scene_selection_visibility

    module procedure hide_scene_selection
        type(selection_set_t) :: targets
        integer :: i
        logical :: accepted

        targets=scene_visibility_targets(s)
        if(.not.any([(scene_target_visible(s,targets%items(i)),i=1,targets%count())]))return
        call set_scene_visibility(s,targets,.false.,accepted)
        if(.not.accepted)return
        call sync_scene_visibility(s)
        call s%status%showMessage('Selection hidden · Ctrl+Z to restore',3500)
    end procedure hide_scene_selection

    subroutine isolate_scene_targets(s,only)
        type(studio_t), intent(inout), target :: s
        type(selection_target_t), intent(in), optional :: only
        type(selection_set_t) :: targets
        integer :: id,i,n

        targets=scene_visibility_targets(s,only)
        if(.not.any([(scene_target_visible(s,targets%items(i)),i=1,targets%count())]))then
            call s%status%showMessage('Select a visible scene item to isolate',3500)
            return
        end if
        s%isolation=targets
        if(allocated(s%isolation_mesh_names))deallocate(s%isolation_mesh_names)
        allocate(s%isolation_mesh_names(count(targets%items%category==selection_mesh)))
        n=0
        do i=1,targets%count()
            if(targets%items(i)%category/=selection_mesh)cycle
            id=targets%items(i)%id
            n=n+1
            s%isolation_mesh_names(n)%s=s%document%meshes(id)%name
        end do
        call sync_scene_visibility(s)
        call s%status%showMessage('Scene isolation active',3000)
    end subroutine isolate_scene_targets

    module procedure exit_scene_isolation
        if(s%isolation%count()==0)return
        call s%isolation%clear()
        if(allocated(s%isolation_mesh_names))deallocate(s%isolation_mesh_names)
        call sync_scene_visibility(s)
        call s%status%showMessage('Isolation ended',3000)
    end procedure exit_scene_isolation

    logical function mesh_source_available_in_scene(s,id) result(available)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: id

        available=.false.
        if(id<1.or.id>studio_mesh_count(s))return
        associate(mesh=>s%document%meshes(id))
            if(mesh%group>0)then
                available=studio_mesh_source_group(s,id)>0
            else
                available=document_patch(s,mesh%kind,mesh%variant)
            end if
        end associate
    end function mesh_source_available_in_scene

    module procedure select_scene_mesh_source
        integer :: id,group,kind,variant

        id=s%selected_mesh
        if(.not.mesh_source_available_in_scene(s,id))then
            call s%status%showMessage('Mesh source is unavailable',3500)
            return
        end if
        group=studio_mesh_source_group(s,id)
        kind=s%document%meshes(id)%kind
        variant=s%document%meshes(id)%variant
        if(group>0)then
            call apply_workspace_selection(s,selection_target_t(selection_assembly,group),selection_replace)
        else
            call apply_workspace_selection(s,selection_target_t(selection_patch,variant,kind),selection_replace)
        end if
        call reveal_scene_selection(s)
    end procedure select_scene_mesh_source

    module procedure reveal_scene_selection
        type(selection_target_t) :: focused
        logical :: was_updating,revealed

        if(present(target))then
            focused=target
        else
            if(allocated(s%selection%items).and.s%selection%count()==0)return
            focused=primary_workspace_target(s)
        end if
        if(focused%category==0)return
        if(s%canvas_only)call toggle_canvas(s)
        if(.not.s%panels(1)%isVisible())call invoke_command(s,toggle_panel_command_base+1)
        was_updating=s%updating
        s%updating=.true.
        revealed=reveal_project_row(s%widgets%project_browser,focused)
        s%updating=was_updating
        if(revealed)then
            call sync_selection_rows(s)
            call s%widgets%project_browser%tree%setFocus()
            call s%status%showMessage('Located in Scene',2500)
        end if
    end procedure reveal_scene_selection

    module procedure on_tree_open
        integer::kind,variant,group,branch_kind,mesh
        if(item%isNull().or.column/=0)return
        select type(s=>data)
        type is(studio_t)
            mesh=mesh_project_id(s%widgets%project_browser,item)
            call project_item_identity(s%widgets%project_browser,s,item,kind,variant,group,branch_kind)
            if(mesh>0)then
                call select_project_mesh(s,mesh)
            else if(group>0)then
                call invoke_command(s,edit_multipatch_command)
            else if(kind>0.and.s%focused_surface)then
                call invoke_command(s,inspect_patch_command)
            end if
        end select
    end procedure on_tree_open

    module procedure open_studio_scene
        logical::reset_tool

        call select_entity(s,s%entity_kind,active_patch_id(s%entity_kind,catalog=s%document%patches))
        if(.not.s%focused_surface)call invoke_command(s,toggle_parameter_panels_command)
        reset_tool=.true.
        if(present(keep_tool))reset_tool=.not.keep_tool
        if(reset_tool)call invoke_command(s,select_tool_command)
        if(.not.studio_assistant_is_visible(s))call invoke_command(s,show_assistant_command)
        call s%surface%setFocus()
        call show_ready(s)
    end procedure open_studio_scene

    module procedure show_scene_multipatch
        call rebuild_project_tree(s)
        call select_multipatch_group(s,id)
    end procedure show_scene_multipatch

    module procedure run_geometry_import
        character(:),allocatable::path,error,choice
        character(len(unit_labels)),parameter :: choices(0:ubound(unit_labels,1))= &
            [character(len(unit_labels)) :: 'Keep numeric values',unit_labels(1:)]
        integer::kind,variant,count,skipped,source
        logical::needs_source_unit,accepted

        path=QFileDialog_getOpenFileName('Import IGES',filter='IGES geometry (*.igs *.iges)',parent=s%window)
        if(len(path)==0)return
        call import_document_geometry(s,path,kind,variant,count,skipped,error,needs_source_unit=needs_source_unit)
        if(needs_source_unit)then
            choice=QInputDialog_getItem('IGES source units','Choose how to interpret the coordinates in this file.',choices, &
                current=s%document%model_unit,editable=.false.,ok=accepted,parent=s%window)
            if(.not.accepted)return
            source=findloc(choices,choice,dim=1)-1
            if(source<0)return
            call import_document_geometry(s,path,kind,variant,count,skipped,error,source_unit=source)
        end if
        if(len(error)>0)then
            call record_studio_diagnostic(s,error)
            call refresh_studio_diagnostics(s)
            call s%status%showMessage('Could not import IGES',6500,warning=.true.,detail=error)
        end if
    end procedure run_geometry_import

    module procedure import_document_geometry
        character(:),allocatable::summary

        if(.not.s%document%document_started.and.s%home_active)then
            call reset_document(s,.false.)
            call mark_project_saved(s)
        end if
        call import_document_iges(s%document,path,kind,variant,count,skipped,error,source_unit,needs_source_unit)
        if(len(error)>0)return

        if(s%home_active)call enter_workspace(s,.true.)
        call rebuild_project_tree(s)
        call select_entity(s,kind,variant)
        call invoke_command(s,fit_scene_command)
        summary='Imported '//integer_text(count)//' IGES objects'
        if(count==1)summary='Imported one IGES object'
        if(skipped>0)summary=summary//' · '//integer_text(skipped)//' unsupported entities skipped'
        if(s%document%model_unit==unitless)summary=summary//' · Set units as needed'
        call s%status%showMessage(summary,7000)
    end procedure import_document_geometry

    module procedure choose_export_path
        character(:),allocatable :: chosen_suffix
        integer::confirm
        logical::exists

        path=QFileDialog_getSaveFileName(title,base//suffix,filter,parent=s%window)
        if(len(path)==0)return
        chosen_suffix=studio_exchange_suffix(path)
        if(chosen_suffix==suffix)return
        if(present(alternate))then
            if(chosen_suffix==alternate)return
        end if
        path=path//suffix
        inquire(file=path,exist=exists)
        if(.not.exists)return
        confirm=QMessageBox_question('Replace export file?','Replace '//path//'?', &
            buttons=ior(QMessageBox_Yes,QMessageBox_No),defaultButton=QMessageBox_No,parent=s%window)
        if(confirm/=QMessageBox_Yes)path=''
    end procedure choose_export_path

    module procedure run_geometry_export
        character(:),allocatable::path,error,name
        integer,allocatable::variants(:)
        integer::kind,id

        call selected_geometry(s,kind,variants,name)
        if(size(variants)==0)return
        id=s%selected_multipatch
        if(.not.vtk.and.(id>0.or.kind==3))return
        if(vtk)then
            path=choose_export_path(s,'Export VTK','geometry','VTK visualization (*.vtk)','.vtk')
        else
            path=choose_export_path(s,'Export IGES','geometry','IGES geometry (*.igs *.iges)','.igs','.iges')
        end if
        if(len(path)==0)return
        call export_document_geometry(s%document,kind,variants(1),id,path,error)
        if(len(error)>0)then
            call s%status%showMessage(error,7000,warning=.true.)
        else
            call s%status%showMessage('Exported '//path,6000)
        end if
    end procedure run_geometry_export

    module procedure run_fortran_export
        character(:),allocatable :: base,path,error

        if(s%home_active.or..not.s%document%document_started)return
        base=document_display_name(s%document)
        if(studio_exchange_suffix(base)=='.fcad')base=base(:len(base)-len('.fcad'))
        if(allocated(s%document%project_directory))base=s%document%project_directory//'/'//base
        path=choose_export_path(s,'Export Fortran',base,'Fortran source (*.f90)','.f90')
        if(len(path)==0)return
        call export_project_fortran(s%document,path,error)
        if(len(error)>0)then
            call s%status%showMessage(error,7000,warning=.true.)
        else
            call s%status%showMessage('Exported '//path,6000)
        end if
    end procedure run_fortran_export

    module procedure export_geometry_selection
        integer,allocatable::variants(:)
        character(:),allocatable::name
        integer::kind

        call selected_geometry(s,kind,variants,name)
        error='Select a curve, surface, volume or assembly to export.'
        if(size(variants)==0)return
        call export_document_geometry(s%document,kind,variants(1),s%selected_multipatch,path,error)
    end procedure export_geometry_selection
end submodule studio_ui_project
