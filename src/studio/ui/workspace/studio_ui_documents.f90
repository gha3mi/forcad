submodule(studio_ui) studio_ui_documents
    use qt, only: QShortcut,QTabBar,Qt_ElideMiddle
    use studio_document, only: move_studio_document,set_document_project_directory
    use studio_document_tabs, only: studio_document_tab_t,studio_document_view_t
    use studio_platform_files, only: resolved_file_path,studio_platform,platform_macos
    use studio_theme, only: studio_active_edge,studio_divider,studio_selected_fill
    use studio_project_merge, only: merge_studio_document
    use studio_assistant, only: stop_assistant_for_document_change
    implicit none(type,external)

contains
    module procedure build_document_tabs
        character(4) :: tab_control

        tab_control=merge('Meta','Ctrl',studio_platform()==platform_macos)
        associate(tabs=>s%project_tabs)
            tabs%next_shortcut=QShortcut(tab_control//'+Tab',s%window)
            tabs%previous_shortcut=QShortcut(tab_control//'+Shift+Tab',s%window)
            call tabs%next_shortcut%setEnabled(.false.)
            call tabs%previous_shortcut%setEnabled(.false.)
            call tabs%next_shortcut%setAutoRepeat(.false.)
            call tabs%previous_shortcut%setAutoRepeat(.false.)
            call tabs%next_shortcut%activated%connect(on_next_document_tab,s)
            call tabs%previous_shortcut%activated%connect(on_previous_document_tab,s)
        end associate
        associate(bar=>s%project_tabs%bar)
            bar=QTabBar()
            call bar%setObjectName('studioDocumentTabs')
            call bar%setAccessibleName('Open projects')
            call bar%setStyleSheet( &
                'QTabBar#studioDocumentTabs::tab { background:#1d2933; color:#bbcfde; border:1px solid '//studio_divider// &
                '; padding:4px 8px 4px 10px; margin-right:3px; min-width:132px; max-width:245px; font-size:13px; }'// &
                'QTabBar#studioDocumentTabs::tab:selected { background:'//studio_selected_fill// &
                '; color:#eefaff; border-bottom:2px solid '//studio_active_edge//'; }'// &
                'QTabBar#studioDocumentTabs::tab:!selected:hover { background:#273b49; color:#eefaff; }')
            call bar%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Preferred)
            call bar%setDocumentMode(.true.)
            call bar%setDrawBase(.false.)
            call bar%setExpanding(.false.)
            call bar%setElideMode(Qt_ElideMiddle)
            call bar%setUsesScrollButtons(.true.)
            call bar%setTabsClosable(.true.)
            call bar%setFixedHeight(30)
            call bar%currentChanged%connect(on_document_tab_changed,s)
            call bar%tabBarClicked%connect(on_document_tab_clicked,s)
            call bar%tabCloseRequested%connect(on_document_tab_closed,s)
            call set_status_hint(bar,'Ctrl+Tab switches projects · Ctrl+Shift+Tab goes back · Ctrl+W closes the active tab')
            call row%addWidget(bar,1)
            call bar%hide()
        end associate
    end procedure build_document_tabs

    module procedure ensure_project_tab
        if(s%project_tabs%active>0.or.s%project_tabs%bar%isNull())return
        call append_project_tab(s,document_display_name(s%document))
    end procedure ensure_project_tab

    module procedure start_project_tab
        logical :: continue_assistant

        if(s%document%document_started)call ensure_project_tab(s)
        if(s%project_tabs%active==0)return
        continue_assistant=.false.
        if(present(assistant_folder))continue_assistant=project_tab_in_folder(s,s%project_tabs%active,assistant_folder)
        if(.not.continue_assistant)call stop_assistant_for_document_change(s)
        call park_active_document(s)
        call advance_document_revision(s)
        call append_project_tab(s,'Untitled project')
    end procedure start_project_tab

    subroutine append_project_tab(s,title)
        type(studio_t),intent(inout) :: s
        character(*),intent(in) :: title
        type(studio_document_tab_t),allocatable :: grown(:)
        integer :: i,count

        associate(tabs=>s%project_tabs)
            count=0
            if(allocated(tabs%items))count=size(tabs%items)
            allocate(grown(count+1))
            do i=1,count
                call move_alloc(tabs%items(i)%document,grown(i)%document)
                grown(i)%view=tabs%items(i)%view
                grown(i)%status_log=tabs%items(i)%status_log
            end do
            call move_alloc(grown,tabs%items)
            tabs%changing=.true.
            tabs%active=count+1
            call tabs%bar%addTab(title)
            call tabs%bar%setCurrentIndex(count)
            tabs%changing=.false.
        end associate
    end subroutine append_project_tab

    module procedure update_project_tabs
        character(:),allocatable :: title
        integer :: i
        logical :: can_cycle

        associate(tabs=>s%project_tabs)
            if(tabs%bar%isNull())return
            can_cycle=tabs%active>0.and.tabs%bar%count()>1
            call tabs%next_shortcut%setEnabled(can_cycle)
            call tabs%previous_shortcut%setEnabled(can_cycle)
            call tabs%bar%setVisible(tabs%active>0)
            call s%widgets%project_caption%setVisible(tabs%active==0)
            if(.not.allocated(tabs%items))return
            do i=1,size(tabs%items)
                if(i==tabs%active)then
                    title=tab_title(s%document)
                else
                    if(.not.allocated(tabs%items(i)%document))cycle
                    title=tab_title(tabs%items(i)%document)
                end if
                if(tabs%bar%tabText(i-1)/=title)call tabs%bar%setTabText(i-1,title)
            end do
        end associate
    end procedure update_project_tabs

    function tab_title(document) result(title)
        type(studio_document_t),intent(in) :: document
        character(:),allocatable :: title

        title=document_display_name(document)
        if(document_has_changes(document))title=title//' *'
    end function tab_title

    module procedure refresh_project_tab_actions
        type(selection_target_t) :: focused
        logical :: selected

        associate(tabs=>s%project_tabs)
            if(tabs%copy%isNull())return
            focused=primary_workspace_target(s)
            selected=focused%category>0
            if(allocated(s%selection%items))selected=s%selection%count()>0
            call tabs%copy%setEnabled(.not.s%home_active.and.selected)
            call tabs%paste%setEnabled(.not.s%home_active.and.allocated(tabs%clipboard))
            call tabs%close%setEnabled(tabs%active>0)
        end associate
    end procedure refresh_project_tab_actions

    module procedure find_project_tab
        character(:),allocatable :: wanted
        integer :: i

        tab=0
        wanted=resolved_file_path(path)
        if(.not.allocated(s%project_tabs%items))return
        do i=1,size(s%project_tabs%items)
            if(i==s%project_tabs%active)then
                if(.not.same_project_path(s%document,path,wanted))cycle
            else
                if(.not.allocated(s%project_tabs%items(i)%document))cycle
                if(.not.same_project_path(s%project_tabs%items(i)%document,path,wanted))cycle
            end if
            tab=i
            return
        end do
    end procedure find_project_tab

    logical function same_project_path(document,path,resolved) result(matches)
        type(studio_document_t),intent(in) :: document
        character(*),intent(in) :: path,resolved

        matches=.false.
        if(.not.allocated(document%project_path))return
        matches=document%project_path==path
        if(.not.matches.and.len(resolved)>0)matches=resolved_file_path(document%project_path)==resolved
    end function same_project_path

    module procedure open_project_tab
        type(studio_document_t) :: incoming
        integer :: tab
        integer(int64) :: previous_revision

        error=''
        if(present(project_directory))then
            error='This project tab is outside the active project folder.'
            if(.not.allocated(s%document%project_directory))return
            if(project_directory/=s%document%project_directory.or.len_trim(project_directory)==0)return
            error=''
        end if
        tab=0
        if(.not.present(sample_title))then
            if(present(remembered_path))then
                tab=find_project_tab(s,remembered_path)
            else
                tab=find_project_tab(s,path)
            end if
        end if
        if(tab>0)then
            if(present(project_directory))then
                if(.not.project_tab_in_folder(s,tab,project_directory))then
                    error='The open project tab belongs to a different project folder.'
                    return
                end if
            end if
            call activate_project_tab(s,tab,keep_assistant=present(project_directory))
            return
        end if
        call open_project_document(incoming,path,error)
        if(len(error)>0)return
        if(present(sample_title))then
            incoming%sample_title=sample_title
            deallocate(incoming%project_path)
            call set_document_project_directory(incoming)
        else
            if(present(remembered_path))incoming%project_path=remembered_path
            if(present(project_directory))call set_document_project_directory(incoming,project_directory)
        end if
        call start_project_tab(s,assistant_folder=project_directory)
        previous_revision=s%document%document_revision
        call move_studio_document(incoming,s%document)
        s%document%document_revision=previous_revision
        call advance_document_revision(s)
        call prepare_document_view(s)
        call enter_workspace(s,.true.)
        call fit_around_panels(s)
        call sync_project_title(s)
        call s%status%showMessage('Opened '//document_display_name(s%document),3500)
    end procedure open_project_tab

    logical function project_tab_in_folder(s,index,folder) result(matches)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: index
        character(*),intent(in) :: folder

        matches=.false.
        if(len_trim(folder)==0)return
        if(index==s%project_tabs%active)then
            if(allocated(s%document%project_directory))matches=s%document%project_directory==folder
        else
            associate(tab=>s%project_tabs%items(index))
                if(.not.allocated(tab%document))return
                if(allocated(tab%document%project_directory))matches=tab%document%project_directory==folder
            end associate
        end if
    end function project_tab_in_folder

    subroutine park_active_document(s)
        type(studio_t),intent(inout),target :: s
        type(QTabWidget) :: ribbon

        associate(tab=>s%project_tabs%items(s%project_tabs%active))
            allocate(tab%document)
            associate(view=>tab%view)
                view%selection=s%selection
                if(.not.allocated(view%selection%items))call view%selection%change(primary_workspace_target(s),selection_replace)
                view%isolation=s%isolation
                if(allocated(view%isolated_meshes))deallocate(view%isolated_meshes)
                if(allocated(s%isolation_mesh_names))view%isolated_meshes=s%isolation_mesh_names
                view%camera=[s%yaw,s%pitch,s%zoom,s%pan_x,s%pan_y]
                view%domain=[s%domain_yaw,s%domain_pitch,s%domain_zoom,s%domain_pan_x,s%domain_pan_y]
                view%probe=[s%u,s%v,s%w]
                view%rank=s%entity_kind
                view%axis=s%knot_axis
                view%derivative=s%basis_derivative
                view%feature=s%viewport_feature
                view%mesh_only=s%mesh_only
                view%mesh_pick=s%mesh_pick
                view%mesh_shaded=s%mesh_shaded
                view%mesh_display=[s%mesh_view,s%mesh_layer_axis,s%mesh_layer_index,s%mesh_section_axis]
                view%section_fraction=s%mesh_section_fraction
                view%physical_knots=s%physical_knots
                call s%workspace%findChild('studioRibbonTabs',ribbon)
                view%ribbon=ribbon%currentIndex()
            end associate
            call s%status%takeLog(tab%status_log)
            call move_studio_document(s%document,tab%document)
        end associate
    end subroutine park_active_document

    module procedure activate_project_tab
        integer(int64) :: previous_revision
        logical :: continue_assistant

        associate(tabs=>s%project_tabs)
            if(tabs%changing.or..not.allocated(tabs%items))return
            if(index<1.or.index>size(tabs%items))return
            if(index==tabs%active)then
                if(s%home_active)call invoke_command(s,resume_project_command)
                return
            end if
            if(.not.allocated(tabs%items(index)%document))return
            continue_assistant=.false.
            if(present(keep_assistant))then
                if(keep_assistant.and.allocated(s%document%project_directory))then
                    continue_assistant=project_tab_in_folder(s,index,s%document%project_directory)
                end if
            end if
            if(.not.continue_assistant)call stop_assistant_for_document_change(s)
            if(tabs%active>0)call park_active_document(s)
            previous_revision=s%document%document_revision
            associate(tab=>tabs%items(index))
                call move_studio_document(tab%document,s%document)
                deallocate(tab%document)
                s%document%document_revision=previous_revision
                call advance_document_revision(s)
                tabs%active=index
                tabs%changing=.true.
                call tabs%bar%setCurrentIndex(index-1)
                call s%status%restoreLog(tab%status_log)
                call restore_tab_view(s,tab%view)
                tabs%changing=.false.
            end associate
        end associate
        call sync_studio_assistant(s)
        call s%status%clearMessage()
        call show_ready(s)
    end procedure activate_project_tab

    subroutine advance_document_revision(s)
        type(studio_t),intent(inout) :: s

        s%project_tabs%revision=max(s%project_tabs%revision,s%document%document_revision)+1_int64
        s%document%document_revision=s%project_tabs%revision
    end subroutine advance_document_revision

    subroutine restore_tab_view(s,view)
        type(studio_t),intent(inout),target :: s
        type(studio_document_view_t),intent(in) :: view
        type(QTabWidget) :: ribbon

        s%entity_kind=view%rank
        call prepare_document_view(s)
        s%isolation=view%isolation
        if(allocated(view%isolated_meshes))s%isolation_mesh_names=view%isolated_meshes
        s%yaw=view%camera(1)
        s%pitch=view%camera(2)
        s%zoom=view%camera(3)
        s%pan_x=view%camera(4)
        s%pan_y=view%camera(5)
        s%u=view%probe(1)
        s%v=view%probe(2)
        s%w=view%probe(3)
        s%knot_axis=view%axis
        s%basis_derivative=view%derivative
        s%mesh_only=view%mesh_only
        s%mesh_pick=view%mesh_pick
        s%mesh_shaded=view%mesh_shaded
        s%mesh_view=view%mesh_display(1)
        s%mesh_layer_axis=view%mesh_display(2)
        s%mesh_layer_index=view%mesh_display(3)
        s%mesh_section_axis=view%mesh_display(4)
        s%mesh_section_fraction=view%section_fraction
        s%physical_knots=view%physical_knots
        call enter_workspace(s,s%resume_assistant)
        if(view%selection%count()>0)then
            call apply_workspace_selection_batch(s,view%selection%items,selection_replace)
        else
            call refresh_knot_table(s)
        end if
        s%viewport_feature=view%feature
        s%domain_yaw=view%domain(1)
        s%domain_pitch=view%domain(2)
        s%domain_zoom=view%domain(3)
        s%domain_pan_x=view%domain(4)
        s%domain_pan_y=view%domain(5)
        call s%workspace%findChild('studioRibbonTabs',ribbon)
        call ribbon%setCurrentIndex(view%ribbon)
        call sync_controls(s)
        call repaint_all(s)
    end subroutine restore_tab_view

    subroutine on_next_document_tab(data)
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            call cycle_project_tab(s,1)
        end select
    end subroutine on_next_document_tab

    subroutine on_previous_document_tab(data)
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            call cycle_project_tab(s,-1)
        end select
    end subroutine on_previous_document_tab

    subroutine cycle_project_tab(s,offset)
        type(studio_t),intent(inout),target :: s
        integer,intent(in) :: offset
        integer :: count

        associate(tabs=>s%project_tabs)
            if(tabs%changing.or.tabs%active==0.or..not.allocated(tabs%items))return
            count=size(tabs%items)
            if(count<2)return
            call activate_project_tab(s,modulo(tabs%active-1+offset,count)+1)
        end associate
    end subroutine cycle_project_tab

    subroutine on_document_tab_changed(index,data)
        integer,intent(in) :: index
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            if(s%project_tabs%changing)return
            call activate_project_tab(s,index+1)
        end select
    end subroutine on_document_tab_changed

    subroutine on_document_tab_clicked(index,data)
        integer,intent(in) :: index
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            if(s%home_active.and.index>=0)call activate_project_tab(s,index+1)
        end select
    end subroutine on_document_tab_clicked

    subroutine on_document_tab_closed(index,data)
        integer,intent(in) :: index
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            call close_project_tab(s,index+1)
        end select
    end subroutine on_document_tab_closed

    module procedure close_active_project_tab
        integer :: index

        index=s%project_tabs%active
        call close_project_tab(s,index)
    end procedure close_active_project_tab

    module procedure close_project_tab
        type(studio_document_tab_t),allocatable :: remaining(:)
        logical :: active,allowed
        integer :: i,n

        associate(tabs=>s%project_tabs)
            if(tabs%changing.or..not.allocated(tabs%items))return
            if(index<1.or.index>size(tabs%items))return
            active=index==tabs%active
            tabs%changing=.true.
            if(active)then
                allowed=confirm_project_change(s)
            else
                allowed=confirm_stored_document(s,index)
            end if
            tabs%changing=.false.
            if(.not.allowed)return
            if(active)then
                call stop_assistant_for_document_change(s)
                call s%status%takeLog(tabs%items(index)%status_log)
                if(size(tabs%items)>1)call clear_document_history(s%document)
            else if(allocated(tabs%items(index)%document))then
                call clear_document_history(tabs%items(index)%document)
            end if
            allocate(remaining(size(tabs%items)-1))
            n=0
            do i=1,size(tabs%items)
                if(i==index)cycle
                n=n+1
                call move_alloc(tabs%items(i)%document,remaining(n)%document)
                remaining(n)%view=tabs%items(i)%view
                remaining(n)%status_log=tabs%items(i)%status_log
            end do
            call move_alloc(remaining,tabs%items)
            if(active)then
                tabs%active=0
            else if(index<tabs%active)then
                tabs%active=tabs%active-1
            end if
            tabs%changing=.true.
            call tabs%bar%removeTab(index-1)
            tabs%changing=.false.
            if(active.and.n>0)then
                call activate_project_tab(s,min(index,n))
            else if(active)then
                if(studio_assistant_is_visible(s))call hide_studio_assistant(s)
                s%resume_assistant=.false.
                s%document%document_started=.false.
                call reset_document(s,.false.)
                call advance_document_revision(s)
                call show_home(s)
            end if
        end associate
        call sync_project_title(s)
        call refresh_project_tab_actions(s)
    end procedure close_project_tab

    module procedure confirm_all_project_tabs
        integer :: i

        allowed=confirm_project_change(s)
        if(.not.allowed.or..not.allocated(s%project_tabs%items))return
        do i=1,size(s%project_tabs%items)
            if(i==s%project_tabs%active)cycle
            allowed=confirm_stored_document(s,i)
            if(.not.allowed)return
        end do
    end procedure confirm_all_project_tabs

    module procedure run_project_import
        character(:),allocatable :: path,error

        path=QFileDialog_getOpenFileName('Import FCAD',filter=studio_project_filter,parent=s%window)
        if(len(path)==0)return
        call import_studio_project(s,path,error)
        if(len(error)>0)call s%status%showMessage('Could not import the project',6500,warning=.true.,detail=error)
    end procedure run_project_import

    module procedure import_studio_project
        type(studio_document_t) :: incoming
        type(selection_target_t),allocatable :: added(:)
        character(:),allocatable :: label

        if(present(imported))imported=0
        error='Open a project before importing into it.'
        if(.not.s%document%document_started.or.s%home_active)return
        call open_project_document(incoming,path,error)
        if(len(error)>0)return
        if(.not.merge_studio_document(s%document,incoming,error,added=added))return
        label=document_display_name(incoming)
        if(present(source_name))label=source_name
        if(present(imported))imported=size(added)
        call show_appended_objects(s,added,'Imported '//label)
    end procedure import_studio_project

    module procedure copy_project_selection
        type(studio_document_t),allocatable :: clipboard
        type(selection_target_t),allocatable :: targets(:),added(:)
        type(selection_target_t) :: focused
        character(:),allocatable :: error

        if(present(failure))failure='Select project objects to copy.'
        if(s%home_active)return
        if(allocated(s%selection%items))then
            targets=s%selection%items
        else
            focused=primary_workspace_target(s)
            if(focused%category==0)return
            targets=[focused]
        end if
        if(size(targets)==0)return
        allocate(clipboard)
        call reset_project_document(clipboard,.false.)
        clipboard%model_unit=s%document%model_unit
        clipboard%display_unit=s%document%display_unit
        if(.not.merge_studio_document(clipboard,s%document,error,selection=targets,added=added))then
            if(present(failure))failure=error
            call s%status%showMessage(error,4500,warning=.true.)
            return
        end if
        call clear_document_history(clipboard)
        if(allocated(s%document%project_directory))clipboard%project_directory=s%document%project_directory
        call move_alloc(clipboard,s%project_tabs%clipboard)
        call move_alloc(added,s%project_tabs%clipboard_selection)
        call refresh_project_tab_actions(s)
        call s%status%showMessage('Copied selection · Switch projects and paste',4000)
        if(present(failure))failure=''
    end procedure copy_project_selection

    module procedure paste_project_selection
        type(selection_target_t),allocatable :: added(:)
        character(:),allocatable :: error

        if(present(failure))failure='Copy project objects before pasting.'
        if(s%home_active.or..not.allocated(s%project_tabs%clipboard))return
        if(present(same_folder))then
            if(same_folder.and..not.project_clipboard_in_scope(s))then
                if(present(failure))failure='Copy objects from a project in the assigned folder before pasting.'
                return
            end if
        end if
        if(.not.merge_studio_document(s%document,s%project_tabs%clipboard,error, &
            selection=s%project_tabs%clipboard_selection,added=added))then
            if(present(failure))failure=error
            call s%status%showMessage(error,5000,warning=.true.)
            return
        end if
        call show_appended_objects(s,added,'Pasted objects into this project')
        if(present(failure))failure=''
    end procedure paste_project_selection

    module procedure project_clipboard_in_scope

        current=.false.
        if(s%home_active.or..not.allocated(s%project_tabs%clipboard))return
        if(.not.allocated(s%document%project_directory))return
        associate(clipboard=>s%project_tabs%clipboard)
            if(.not.allocated(clipboard%project_directory))return
            current=clipboard%project_directory==s%document%project_directory
        end associate
    end procedure project_clipboard_in_scope

    subroutine show_appended_objects(s,added,message)
        type(studio_t),intent(inout),target :: s
        type(selection_target_t),intent(in) :: added(:)
        character(*),intent(in) :: message

        call refresh_project_merge_view(s)
        if(size(added)>0)then
            call apply_workspace_selection_batch(s,added,selection_replace)
            if(fit_workspace_selection(s,added))call s%surface%update()
        end if
        call s%status%showMessage(message,5000)
    end subroutine show_appended_objects

    module procedure refresh_project_merge_view
        call cancel_transform_drag(s)
        call forget_workspace_selection(s)
        s%live_scene%valid=.false.
        s%selected_multipatch=0
        call clear_mesh_selection(s)
        if(s%mesh_only.and..not.has_visible_meshes(s))s%mesh_only=.false.
        call rebuild_project_tree(s)
        call pick_document_selection(s)
        call sync_controls(s)
        call refresh_edit_actions(s)
        call sync_project_title(s)
        call repaint_all(s)
    end procedure refresh_project_merge_view
end submodule studio_ui_documents
