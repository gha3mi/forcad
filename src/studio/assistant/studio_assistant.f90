! Local chat with an in-canvas dock and a native detached window.
module studio_assistant
    use, intrinsic :: iso_c_binding, only: c_null_char
    use qt, only: QClipboard, QDir_exists, QGuiApplication_clipboard, QJsonObject, QKeyEvent, QMainWindow, &
        QJsonDocument, QJsonDocument_fromJson, &
        QFileDialog_getSaveFileName, QDir_removeRecursively, QSettings, &
        QRect, QThread, QTimer, QWidget, Qt_RightDockWidgetArea, QtfEventFilter, &
        Qt_Key_Return, Qt_Key_Enter, Qt_NoModifier, Qt_ControlModifier
    use studio_preferences, only: application_settings
    use studio_platform_files, only: last_path_separator,temporary_private_directory,replace_file
    use studio_storage_paths, only: studio_project_directory
    use studio_assistant_panel, only: assistant_panel_t, build_assistant_panel
    use studio_assistant_conversation, only: assistant_conversation_t
    use studio_assistant_transcript, only: assistant_transcript_t, chat_body
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_state, only: studio_t
    use studio_presentation_text, only: count_text=>integer_text, without_nul
    use studio_chat_attachment, only: studio_attachment
    use studio_project_folder_dialog, only: choose_project_folder
    use studio_project_folder, only: canonical_project_folder
    use studio_document, only: set_document_project_directory
    use studio_assistant_transport, only: read_assistant_progress, cleanup_assistant_stream, read_agent_events, &
        read_agent_permission, send_agent_decision, assistant_worker_t, assistant_worker_count, write_text
    use studio_assistant_session, only: assistant_session_t, agent_selection
    use studio_assistant_team_settings, only: edit_assistant_team
    use studio_mcp_bridge, only: start_studio_mcp, stop_studio_mcp, pending_studio_mcp_approval, decide_studio_mcp_approval, &
        rebind_assistant_mcp_project
    use studio_assistant_draft_context, only: capture_assistant_selection
    use studio_live_viewport, only: save_viewport_image
    use studio_diagnostics, only: studio_diagnostic_summary
    use studio_ai_cursor, only: hide_studio_ai_cursor, set_studio_ai_tool_cursor
    use studio_status_hints, only: register_status_window, set_status_hint
    implicit none(type, external)
    private
    public :: show_studio_assistant, sync_studio_assistant, hide_studio_assistant
    public :: studio_assistant_is_docked, studio_assistant_is_visible, finish_studio_assistant
    public :: ask_studio_diagnostic
    public :: begin_studio_assistant_mcp_event, append_studio_assistant_mcp_event
    public :: studio_assistant_text, stop_assistant_for_document_change

    integer, parameter :: stream_interval_ms = 120
    integer, parameter :: cursor_feedback_ms = 4500
    character(*), parameter :: approval_setting='assistant/ask-before-studio-edits'
    character(*), parameter :: team_setting='assistant/team-mode', lead_agent='Lead'
    character(*), parameter :: team_settings_group='assistant/team-worker-'
    character(20), parameter :: team_agents(*)=[character(20)::'Geometry reviewer','Modeling planner',lead_agent]

    abstract interface
        subroutine dock_layout_listener(s)
            import studio_t
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine dock_layout_listener
    end interface

    type :: assistant_t
        type(studio_t), pointer :: workspace => null()
        procedure(dock_layout_listener), pointer, nopass :: layout_changed => null()
        type(assistant_panel_t) :: panel
        type(assistant_transcript_t) :: transcript
        type(assistant_conversation_t) :: conversation
        logical :: updating_agent_options=.false.
        type(assistant_session_t) :: session
        type(QTimer) :: progress_timer, cursor_timer, transcript_tail_timer
        type(studio_attachment),allocatable :: attachments(:)
        integer(int64) :: requested_document=-1
        integer(int64) :: attached_selection_document=-1
        integer(int64) :: elapsed_started=0, elapsed_seconds=-1
        character(:), allocatable :: attachment_folder,attached_selection,attached_selection_label
        logical :: busy = .false.
        logical :: studio_access=.true., ask_before_edits=.true., team_mode=.false.
        type(assistant_worker_t) :: team_workers(assistant_worker_count)
        logical :: active_agents(size(team_agents))=.false., permission_pending=.false.

        integer(int64) :: studio_permission_ticket=0_int64
        integer :: shown_agent_bytes=0, committed_tools=0, failed_tools=0
        character(:),allocatable :: permission_id,last_agent_decision
        logical :: stop_requested = .false.
        logical :: listing_models = .false.
        logical :: model_ready = .false.
        logical :: refresh_after_reply = .false., reply_started = .false.
        logical :: timeline_stream=.false.
        character(len=:), allocatable :: selected_model
        character(len=:), allocatable :: resume_model
    end type assistant_t

    ! Qt slots and the worker retain references to these records until shutdown.
    type(assistant_t), target :: assistant

    interface
        module subroutine set_model_status(state,text,detail)
            implicit none(type, external)
            type(assistant_t), intent(inout) :: state
            character(*), intent(in) :: text,detail
        end subroutine set_model_status

        module subroutine on_model_changed(text,data)
            implicit none(type, external)
            character(len=*), intent(in) :: text
            class(*), intent(inout), optional :: data
        end subroutine on_model_changed

        module subroutine on_refresh_models(data)
            implicit none(type, external)
            class(*), intent(inout), optional :: data
        end subroutine on_refresh_models

        module subroutine finish_model_list(state)
            implicit none(type, external)
            type(assistant_t), intent(inout) :: state
        end subroutine finish_model_list

        module subroutine populate_agent_options(state)
            implicit none(type, external)
            type(assistant_t),intent(inout) :: state
        end subroutine populate_agent_options

        module function current_agent_selection(state) result(choice)
            implicit none(type, external)
            type(assistant_t),intent(in) :: state
            type(agent_selection) :: choice
        end function current_agent_selection

        module subroutine on_agent_model_changed(index,data)
            implicit none(type, external)
            integer,intent(in) :: index
            class(*),intent(inout),optional :: data
        end subroutine on_agent_model_changed

        module subroutine on_agent_effort_changed(index,data)
            implicit none(type, external)
            integer,intent(in) :: index
            class(*),intent(inout),optional :: data
        end subroutine on_agent_effort_changed

        module subroutine check_saved_agent_selection(state)
            implicit none(type, external)
            type(assistant_t),intent(inout) :: state
        end subroutine check_saved_agent_selection

        module subroutine require_replacement_model(state)
            implicit none(type, external)
            type(assistant_t), intent(inout) :: state
        end subroutine require_replacement_model

        module subroutine update_provider_badge(state)
            implicit none(type, external)
            type(assistant_t),intent(inout) :: state
        end subroutine update_provider_badge

        module subroutine on_new_chat(data)
            implicit none(type, external)
            class(*), intent(inout), optional :: data
        end subroutine on_new_chat

        module subroutine reset_chat(state)
            implicit none(type, external)
            type(assistant_t), intent(inout) :: state
        end subroutine reset_chat

        module subroutine update_chat_status(state)
            implicit none(type, external)
            type(assistant_t), intent(inout) :: state
        end subroutine update_chat_status

        module subroutine on_save_chat(data)
            implicit none(type, external)
            class(*), intent(inout), optional :: data
        end subroutine on_save_chat

        module subroutine on_chat_history(data)
            implicit none(type, external)
            class(*), intent(inout), optional :: data
        end subroutine on_chat_history

        module subroutine on_attach(data)
            implicit none(type, external)
            class(*),intent(inout),optional :: data
        end subroutine on_attach

        module subroutine on_attach_selection(data)
            implicit none(type, external)
            class(*), intent(inout), optional :: data
        end subroutine on_attach_selection

        module subroutine on_attach_view(data)
            implicit none(type, external)
            class(*), intent(inout), optional :: data
        end subroutine on_attach_view

        module subroutine on_clear_draft_selection(data)
            implicit none(type, external)
            class(*), intent(inout), optional :: data
        end subroutine on_clear_draft_selection

        module subroutine clear_draft_selection(state)
            implicit none(type, external)
            type(assistant_t), intent(inout) :: state
        end subroutine clear_draft_selection

        module subroutine sync_draft_attachment(state)
            implicit none(type, external)
            type(assistant_t), intent(inout) :: state
        end subroutine sync_draft_attachment

        module subroutine on_detach(data)
            implicit none(type, external)
            class(*),intent(inout),optional :: data
        end subroutine on_detach

        module integer function append_draft_attachments(state,choices,folder) result(added)
            implicit none(type, external)
            type(assistant_t),intent(inout) :: state
            type(studio_attachment),intent(in) :: choices(:)
            character(*),intent(in) :: folder
        end function append_draft_attachments

        module subroutine clear_attachment(state)
            implicit none(type, external)
            type(assistant_t),intent(inout) :: state
        end subroutine clear_attachment
    end interface

contains

    subroutine show_assistant_status(state,text,warning,detail,record)
        type(assistant_t),intent(inout) :: state
        character(*),intent(in) :: text
        logical,intent(in),optional :: warning
        character(*),intent(in),optional :: detail
        logical,intent(in),optional :: record

        if(.not.associated(state%workspace))return
        if(len_trim(text)==0)then
            if(index(state%workspace%status%currentMessage(),'Assistant · ')==1)call state%workspace%status%clearMessage()
            return
        end if
        call state%workspace%status%showMessage('Assistant · '//text,7000,warning,detail=detail,record=record)
    end subroutine show_assistant_status

    subroutine show_assistant_progress(state,text,approval)
        type(assistant_t),intent(inout) :: state
        character(*),intent(in) :: text
        logical,intent(in),optional :: approval

        if(present(approval))then
            state%permission_pending=approval
        else if(state%permission_pending.or.(state%busy.and.state%stop_requested))then
            return
        end if
        if(.not.associated(state%workspace))return
        if(len_trim(text)==0)then
            call state%workspace%status%setAssistantActivity('',approval=approval)
        else
            call state%workspace%status%setAssistantActivity('Assistant · '//text,approval=approval)
        end if
    end subroutine show_assistant_progress

    subroutine show_studio_assistant(s,layout_changed)
        type(studio_t), intent(inout), target :: s
        procedure(dock_layout_listener) :: layout_changed
        logical :: reopening

        if(.not.assistant%panel%dock%isNull())then
            if(.not.assistant_owns_workspace(s))then
                call s%status%showMessage('The Assistant belongs to another workspace.',4500)
                return
            end if
        end if
        assistant%workspace => s
        assistant%layout_changed => layout_changed
        if (assistant%panel%dock%isNull()) call build_assistant(s)
        call sync_studio_assistant(s)
        reopening=.not.assistant%panel%dock%isVisible()
        call assistant%panel%dock%show()
        if (assistant%panel%dock%isFloating()) call assistant%panel%dock%activateWindow()
        call assistant%panel%prompt%setFocus()
        if (reopening) then
            if (assistant%busy) then
                assistant%refresh_after_reply=.true.
            else
                call on_refresh_models(assistant)
            end if
        end if
    end subroutine show_studio_assistant

    subroutine hide_studio_assistant(s)
        type(studio_t), intent(in), optional, target :: s

        if(.not.assistant_owns_workspace(s))return
        if (.not.assistant%panel%dock%isNull()) call assistant%panel%dock%hide()
    end subroutine hide_studio_assistant

    logical function studio_assistant_is_docked(s) result(docked)
        type(studio_t), intent(in), optional, target :: s

        docked=.false.
        if(.not.assistant_owns_workspace(s))return
        if (.not.assistant%panel%dock%isNull()) then
            docked=assistant%panel%dock%isVisible().and..not.assistant%panel%dock%isFloating()
        end if
    end function studio_assistant_is_docked

    logical function studio_assistant_is_visible(s) result(visible)
        type(studio_t), intent(in), optional, target :: s

        visible=.false.
        if(.not.assistant_owns_workspace(s))return
        if(.not.assistant%panel%dock%isNull())visible=assistant%panel%dock%isVisible()
    end function studio_assistant_is_visible

    logical function assistant_owns_workspace(s) result(owned)
        type(studio_t), intent(in), optional, target :: s

        owned=.true.
        if(present(s))owned=associated(assistant%workspace,s)
    end function assistant_owns_workspace

    subroutine stop_assistant_for_document_change(s)
        type(studio_t),intent(inout),target :: s

        if(.not.assistant_owns_workspace(s).or.assistant%panel%dock%isNull())return
        call on_stop(assistant)
        call s%status%setAssistantActivity('',elapsed='',approval=.false.)
        call assistant%cursor_timer%stop()
        call hide_studio_ai_cursor(s)
        call clear_draft_selection(assistant)
        call clear_attachment(assistant)
        call assistant%panel%continue_scene%hide()
        assistant%requested_document=-1_int64
        call update_controls(assistant)
    end subroutine stop_assistant_for_document_change

    subroutine sync_studio_assistant(s)
        type(studio_t), intent(in), target :: s

        if (.not.assistant_owns_workspace(s).or.assistant%panel%dock%isNull()) return
        if(assistant%busy.and..not.assistant%stop_requested)then
            if(assistant%requested_document/=s%document%document_revision)then
                if(rebind_assistant_mcp_project(s))then
                    assistant%requested_document=s%document%document_revision
                else
                    call on_stop(assistant)
                end if
            end if
        end if
        if(allocated(assistant%attachment_folder).and..not.assistant%busy)then
            if(.not.allocated(s%document%project_directory))then
                call clear_attachment(assistant)
                call update_controls(assistant)
            else if(assistant%attachment_folder/=s%document%project_directory)then
                call clear_attachment(assistant)
                call update_controls(assistant)
            end if
        end if
        if(assistant%attached_selection_document/=s%document%document_revision)call clear_draft_selection(assistant)
        if(assistant%requested_document/=s%document%document_revision)call assistant%panel%continue_scene%hide()
    end subroutine sync_studio_assistant

    subroutine build_assistant(s)
        type(studio_t), intent(inout) :: s
        type(QtfEventFilter) :: prompt_events
        type(QMainWindow) :: dock_host
        type(QSettings) :: settings
        type(QWidget) :: content
        character(:),allocatable :: prefix
        integer :: i
        call build_assistant_panel(s%surface,s%assistant_overlay,assistant%panel,dock_host,content)
        call register_status_window(s%window,assistant%panel%dock)
        assistant%selected_model=''
        call assistant%conversation%clear()
        assistant%resume_model=''
        call set_model_status(assistant,'Checking…','')
        call sync_studio_tools_control(assistant,.false.,.false.)
        call assistant%panel%agent_model%currentIndexChanged%connect(on_agent_model_changed,assistant)
        call assistant%panel%agent_effort%currentIndexChanged%connect(on_agent_effort_changed,assistant)
        call assistant%panel%studio_tools%clicked%connect(on_studio_tools_toggle,assistant)
        settings=application_settings()
        assistant%ask_before_edits=settings%boolValue(approval_setting,.true.)
        assistant%team_mode=settings%boolValue(team_setting,.false.)
        do i=1,size(assistant%team_workers)
            prefix=team_settings_group//count_text(i)//'/'
            associate(worker=>assistant%team_workers(i))
                worker%provider=settings%stringValue(prefix//'provider','')
                worker%model=settings%stringValue(prefix//'model','')
                worker%effort=settings%stringValue(prefix//'effort','')
            end associate
        end do
        call assistant%panel%team_mode%clicked%connect(on_team_mode_clicked,assistant)
        call assistant%panel%team_settings%clicked%connect(on_team_settings,assistant)
        call settings%delete()
        call assistant%panel%approval_mode%setCurrentIndex(merge(0,1,assistant%ask_before_edits))
        call assistant%panel%approval_mode%currentIndexChanged%connect(on_approval_mode_changed,assistant)
        call assistant%panel%attach%triggered%connect(on_attach,assistant)
        call assistant%panel%attach_selection%triggered%connect(on_attach_selection,assistant)
        call assistant%panel%attach_view%triggered%connect(on_attach_view,assistant)
        call assistant%panel%detach%clicked%connect(on_detach,assistant)
        call assistant%panel%clear_context%clicked%connect(on_clear_draft_selection,assistant)
        call assistant%panel%deny_request%clicked%connect(on_deny_agent_request,assistant)
        call assistant%panel%allow_once%clicked%connect(on_allow_agent_request,assistant)

        assistant%transcript_tail_timer=QTimer(content)
        call assistant%transcript_tail_timer%setSingleShot(.true.)
        call assistant%transcript%attach(assistant%panel%transcript,assistant%panel%transcript_scrollbar, &
            assistant%transcript_tail_timer,assistant%panel%jump_to_latest)
        call assistant%panel%jump_to_latest%clicked%connect(on_jump_to_latest,assistant)
        call assistant%transcript_tail_timer%timeout%connect(on_transcript_tail,assistant)
        assistant%progress_timer=QTimer(content)

        call assistant%progress_timer%timeout%connect(on_stream_progress,assistant)
        assistant%cursor_timer=QTimer(content)
        call assistant%cursor_timer%setSingleShot(.true.)
        call assistant%cursor_timer%timeout%connect(on_cursor_expired,assistant)
        call assistant%panel%model%currentTextChanged%connect(on_model_changed,assistant)
        call assistant%panel%refresh_models%clicked%connect(on_refresh_models,assistant)
        call assistant%panel%prompt%textChanged%connect(on_input_changed,assistant)
        prompt_events=QtfEventFilter(assistant%panel%prompt)
        call prompt_events%keyPressEvent%connect(on_prompt_key,assistant)
        call assistant%panel%ask%clicked%connect(on_ask,assistant)
        call assistant%panel%stop%clicked%connect(on_stop,assistant)
        call assistant%panel%continue_scene%clicked%connect(on_continue_scene,assistant)
        call assistant%panel%new_chat%clicked%connect(on_new_chat,assistant)
        call assistant%panel%save_button%triggered%connect(on_save_chat,assistant)
        call assistant%panel%history_button%triggered%connect(on_chat_history,assistant)
        call assistant%panel%copy_chat%triggered%connect(on_copy_chat,assistant)
        call assistant%panel%export_chat%triggered%connect(on_export_chat,assistant)
        call update_chat_status(assistant)
        call update_controls(assistant)
        call assistant%panel%dock%setWidget(content)
        call assistant%panel%dock%hide()
        call dock_host%addDockWidget(Qt_RightDockWidgetArea,assistant%panel%dock)
        call s%assistant_overlay%setGeometry(QRect(0,0,0,0))
        call s%assistant_overlay%show()
        call assistant%panel%dock%visibilityChanged%connect(on_dock_changed,assistant)
        call assistant%panel%dock%topLevelChanged%connect(on_dock_changed,assistant)
    end subroutine build_assistant

    subroutine sync_studio_tools_control(state,available,active)
        type(assistant_t),intent(inout) :: state
        logical,intent(in) :: available,active

        if(state%panel%studio_tools%isNull())return
        call state%panel%studio_tools%setChecked(active.and.available)
        call state%panel%studio_tools%setEnabled(available.and..not.state%busy.and..not.state%listing_models)
        if(.not.available)then
            call state%panel%studio_tools%setText('Tools: Setup')
            call set_status_hint(state%panel%studio_tools,'Choose an agent and model to use ForCAD Studio tools.')
        else if(active)then
            call state%panel%studio_tools%setText('Tools: Live')
            call set_status_hint(state%panel%studio_tools,'The agent can inspect and edit the open design. '// &
                'Turn off to chat without live ForCAD Studio access.')
        else
            call state%panel%studio_tools%setText('Tools: Off')
            call set_status_hint(state%panel%studio_tools,'The ForCAD Studio connection is off. The agent still uses its own tools and permissions.')
        end if
    end subroutine sync_studio_tools_control

    recursive subroutine on_dock_changed(changed,data)
        logical, intent(in) :: changed
        class(*), intent(inout), optional :: data
        select type(state=>data)
        type is(assistant_t)
            if (associated(state%layout_changed).and.associated(state%workspace)) then
                call state%layout_changed(state%workspace)
            end if
        end select
    end subroutine on_dock_changed

    subroutine on_input_changed(data)
        class(*), intent(inout), optional :: data
        select type(state=>data)
        type is(assistant_t)
            call update_controls(state)
        end select
    end subroutine on_input_changed

    subroutine on_prompt_key(event,data)
        type(QKeyEvent), intent(in) :: event
        class(*), intent(inout), optional :: data

        if(event%key()/=Qt_Key_Return.and.event%key()/=Qt_Key_Enter)return
        if(event%modifiers()/=Qt_NoModifier.and.event%modifiers()/=Qt_ControlModifier)return
        call event%accept()
        if(.not.event%isAutoRepeat())call on_ask(data)
    end subroutine on_prompt_key

    subroutine on_jump_to_latest(data)
        class(*), intent(inout), optional :: data
        select type(state=>data)
        type is(assistant_t)
            call state%transcript%follow()
        end select
    end subroutine on_jump_to_latest

    subroutine on_copy_chat(data)
        class(*), intent(inout), optional :: data
        type(QClipboard) :: clipboard
        character(:), allocatable :: text

        select type(state=>data)
        type is(assistant_t)
            text=state%transcript%plain_text()
            if(len(text)==0)return
            clipboard=QGuiApplication_clipboard()
            call clipboard%setText(text)
            call show_assistant_status(state,'Chat copied.')
        end select
    end subroutine on_copy_chat

    subroutine on_export_chat(data)
        class(*), intent(inout), optional :: data
        character(:), allocatable :: text,path,folder

        select type(state=>data)
        type is(assistant_t)
            text=state%transcript%plain_text()
            if(len(text)==0)return
            path='Chat.txt'
            folder=studio_project_directory(state%conversation%project)
            if(len(folder)>0)path=folder//'/'//path
            path=QFileDialog_getSaveFileName('Export chat',dir=path,filter='Text files (*.txt)',parent=state%panel%dock)
            if(len(path)==0)return
            if(save_chat_export(path,text))then
                call show_assistant_status(state,'Chat exported.')
            else
                call show_assistant_status(state,'Could not export the chat.',warning=.true.)
            end if
        end select
    end subroutine on_export_chat

    logical function save_chat_export(path,text) result(saved)
        character(*),intent(in) :: path,text
        character(:),allocatable :: folder,staging,temporary
        integer :: separator
        logical :: removed

        saved=.false.
        separator=last_path_separator(path)
        folder='.'
        if(separator>0)folder=path(:separator)
        staging=temporary_private_directory(folder,'.forcad-chat-export.')
        if(len(staging)==0)return
        temporary=staging//'/conversation.txt'
        if(write_text(temporary,text))saved=replace_file(temporary//c_null_char,path//c_null_char)==0
        removed=QDir_removeRecursively(staging)
    end function save_chat_export

    integer function attachment_vision(state) result(vision)
        type(assistant_t),intent(in) :: state
        integer :: model

        vision=-1
        if(.not.state%model_ready.or..not.allocated(state%session%agent_models))return
        model=state%panel%agent_model%currentIndex()
        if(model<0.or.model>ubound(state%session%agent_models,1))return
        vision=state%session%agent_models(model)%vision
    end function attachment_vision

    subroutine update_controls(state)
        type(assistant_t), intent(inout) :: state
        logical :: model_available,active
        active=state%busy
        model_available=state%model_ready.and..not.state%listing_models.and..not.active
        call sync_draft_attachment(state)
        call state%panel%ask%setEnabled(model_available.and.len_trim(state%panel%model%currentText())>0.and. &
                                  len(chat_body(state%panel%prompt%toPlainText()))>0)
        call state%panel%stop%setEnabled(active.and..not.state%stop_requested)
        call state%panel%attach%setEnabled(.not.active)
        call state%panel%attach_view%setEnabled(.not.active)
        call state%panel%attach_selection%setEnabled(.not.active)
        call state%panel%clear_context%setEnabled(.not.active)
        call state%panel%attachment_list%setEnabled(.not.active)
        call state%panel%detach%setEnabled(.not.active.and.state%panel%attachment_list%currentRow()>=0)
        call state%panel%new_chat%setEnabled(.not.active)
        call state%panel%save_button%setEnabled(.not.active.and.state%conversation%dirty)
        call state%panel%history_button%setEnabled(.not.active)
        call state%panel%copy_chat%setEnabled(state%conversation%dirty.or.len(state%conversation%id)>0)
        call state%panel%export_chat%setEnabled(state%conversation%dirty.or.len(state%conversation%id)>0)
        call state%panel%model%setEnabled(model_available)
        call state%panel%agent_model%setEnabled(model_available.and.state%panel%agent_model%count()>1)
        call state%panel%agent_effort%setEnabled(model_available.and.state%panel%agent_effort%count()>1)
        call state%panel%refresh_models%setEnabled(.not.active.and..not.state%listing_models)
        call sync_studio_tools_control(state,state%model_ready,state%studio_access)
        call state%panel%team_mode%setEnabled(model_available.and.state%studio_access)
        call state%panel%team_mode%setChecked(state%team_mode.and.state%studio_access)
        call state%panel%team_settings%setVisible(state%team_mode.and.state%studio_access)
        call state%panel%team_settings%setEnabled(model_available.and.state%studio_access)
        call state%panel%approval_mode%setEnabled(.not.active)

        call state%panel%sync(active,state%model_ready.and.state%panel%agent_effort%count()>1)
    end subroutine update_controls

    subroutine start_worker(list_models)
        logical, intent(in) :: list_models
        type(QThread) :: worker

        worker=assistant%session%prepare_worker(list_models)
        call worker%finished%connect(on_worker_finished,assistant)
        call worker%start()
    end subroutine start_worker

    subroutine on_ask(data)
        class(*), intent(inout), optional :: data
        character(:), allocatable :: activity,question,draft_context,project
        type(agent_selection) :: choice
        integer :: i

        select type(state=>data)
        type is(assistant_t)
            if (state%busy.or.state%listing_models.or..not.state%model_ready) return
            if (len_trim(state%panel%model%currentText())==0) return
            if (.not.associated(state%workspace)) return
            question=without_nul(state%panel%prompt%toPlainText())
            if(len(chat_body(question))==0)return
            if(len(question)>8000)then
                call show_assistant_status(state,'Shorten the message to 8,000 bytes to send it.',warning=.true.)
                return
            end if
            if(.not.ensure_project_folder(state))return
            project=assistant_project_folder(state)
            if(len(project)==0)then
                call show_assistant_status(state,'Choose a valid project folder to use the Assistant.',warning=.true.)
                return
            end if
            if(len(state%conversation%project)>0.and.state%conversation%project/=project)then
                call show_assistant_status(state,'This chat belongs to another project. Start a new chat or return to its folder.',warning=.true.)
                return
            end if
            call state%session%prepare_reply(state%workspace%document%project_directory)

            if (len_trim(state%session%request%error)>0) then
                call show_assistant_status(state,state%session%request%error,warning=.true.)
                return
            end if
            state%session%request%list_only=.false.
            state%session%request%model=state%panel%model%currentText()
            choice=current_agent_selection(state)
            state%session%request%agent_model=choice%model
            state%session%request%agent_effort=choice%effort
            state%session%request%question=trim(question)
            state%shown_agent_bytes=0
            state%timeline_stream=.false.
            state%active_agents=.false.
            state%committed_tools=0
            state%failed_tools=0
            call state%panel%continue_scene%hide()
            state%permission_id=''
            state%last_agent_decision=''
            state%studio_permission_ticket=0_int64
            call state%transcript%clear_reply()
            state%session%request%mcp_session=''
            if(state%studio_access)then
                state%session%request%mcp_session=start_studio_mcp(.true.,require_approval=state%ask_before_edits)
                if(len(state%session%request%mcp_session)==0)then
                    call cleanup_assistant_stream(state%session%request)
                    call show_assistant_status(state,'Could not connect to ForCAD Studio. Retry or switch tools off.',warning=.true.)
                    return
                end if
            end if
            state%session%request%team_mode=state%team_mode.and.len(state%session%request%mcp_session)>0
            if(state%session%request%team_mode)state%session%request%team_workers=state%team_workers
            associate(request=>state%session%request,conversation=>state%conversation)
                request%resume_session=''
                if(request%model==conversation%provider.and.request%working_directory==conversation%provider_workspace)then
                    request%resume_session=conversation%provider_session
                end if
                request%history=''
                if(len(request%resume_session)==0)request%history=state%session%serialize_history()
            end associate
            state%requested_document=state%workspace%document%document_revision

            if(allocated(state%session%request%attachments))deallocate(state%session%request%attachments)
            if(allocated(state%attachments))state%session%request%attachments=state%attachments
            state%session%request%attachment_reference_only=attachment_vision(state)==0

            state%session%request%context=''
            if(allocated(state%attached_selection))state%session%request%context=state%attached_selection
            state%session%request%reply=''
            state%session%request%cancelled=.false.
            state%stop_requested=.false.
            call state%cursor_timer%stop()
            call hide_studio_ai_cursor(state%workspace)
            state%busy=.true.
            call start_elapsed(state)
            associate(request=>state%session%request)
                call state%conversation%begin_turn(request%question,request%model,request%agent_model,request%agent_effort,project)
                state%conversation%provider_session=request%resume_session
                state%conversation%provider_workspace=request%working_directory
            end associate
            call update_chat_status(state)
            draft_context=state%session%request%question
            if(allocated(state%attachments))then
                do i=1,size(state%attachments)
                    draft_context=draft_context//new_line('a')//'Reference · '//state%attachments(i)%name
                end do
            end if
            if(allocated(state%attached_selection_label))draft_context=draft_context//new_line('a')//'Selection · '//state%attached_selection_label
            call state%transcript%append_message('You',draft_context)

            call state%transcript%begin_reply()
            call state%panel%prompt%clear()
            call update_controls(state)
            call state%panel%prompt%setFocus()
            call show_assistant_status(state,'')
            if(state%session%request%team_mode)then
                call show_assistant_progress(state,'Team working…')
            else
                call show_assistant_progress(state,'Waiting for '//state%session%request%model//'…')
            end if
            call sync_studio_assistant(state%workspace)
            state%reply_started=.false.
            activity='Request submitted via '//state%session%request%model
            if(len(choice%model)>0)activity=activity//' · '//choice%model
            if(len(choice%effort)>0)activity=activity//' · effort '//state%panel%agent_effort%currentText()
            if(state%session%request%team_mode)activity=activity//' · Team'
            call log_activity(state,activity)
            if(allocated(state%attachments))then
                do i=1,size(state%attachments)
                    associate(file=>state%attachments(i))
                        if(allocated(file%project_copy))deallocate(file%project_copy)
                        if(allocated(file%copy_directory))deallocate(file%copy_directory)
                    end associate
                end do
            end if
            call start_worker(.false.)

            call state%progress_timer%start(stream_interval_ms)
        end select
    end subroutine on_ask

    subroutine on_studio_tools_toggle(data)
        class(*),intent(inout),optional :: data

        select type(state=>data)
        type is(assistant_t)
            state%studio_access=state%panel%studio_tools%isChecked()
            call update_controls(state)
        end select
    end subroutine on_studio_tools_toggle

    subroutine on_team_mode_clicked(data)
        class(*),intent(inout),optional :: data
        type(QSettings) :: settings

        select type(state=>data)
        type is(assistant_t)
            if(state%busy.or..not.state%studio_access)return
            state%team_mode=state%panel%team_mode%isChecked()
            settings=application_settings()
            call settings%setValue(team_setting,state%team_mode)
            call settings%sync()
            call settings%delete()
            call update_controls(state)
        end select
    end subroutine on_team_mode_clicked

    subroutine on_team_settings(data)
        class(*),intent(inout),optional :: data
        type(QSettings) :: settings
        character(:),allocatable :: lead,prefix
        integer :: i

        select type(state=>data)
        type is(assistant_t)
            if(.not.state%panel%team_settings%isEnabled().or..not.associated(state%workspace))return
            lead='Lead · '//state%panel%model%currentText()//' · '//state%panel%agent_model%currentText()
            if(.not.edit_assistant_team(state%panel%dock,state%session,state%panel%model,lead, &
                team_agents(:assistant_worker_count),state%team_workers))return
            settings=application_settings()
            do i=1,size(state%team_workers)
                prefix=team_settings_group//count_text(i)//'/'
                associate(worker=>state%team_workers(i))
                    call settings%setValue(prefix//'provider',worker%provider)
                    call settings%setValue(prefix//'model',worker%model)
                    call settings%setValue(prefix//'effort',worker%effort)
                end associate
            end do
            call settings%sync()
            call settings%delete()
            call show_assistant_status(state,'Team settings saved.')
        end select
    end subroutine on_team_settings

    subroutine on_approval_mode_changed(index,data)
        integer, intent(in) :: index
        class(*), intent(inout), optional :: data
        type(QSettings) :: settings

        select type(state=>data)
        type is(assistant_t)
            if(state%busy.or.index<0.or.index>1)return
            state%ask_before_edits=index==0
            settings=application_settings()
            call settings%setValue(approval_setting,state%ask_before_edits)
            call settings%sync()
            call settings%delete()
        end select
    end subroutine on_approval_mode_changed

    subroutine collect_agent_events(state,drain)
        type(assistant_t),intent(inout) :: state
        logical,intent(in),optional :: drain
        character(:),allocatable :: progress
        integer :: first,last
        logical :: draining

        draining=.false.
        if(present(drain))draining=drain
        do
            progress=read_agent_events(state%session%request,state%shown_agent_bytes)
            if(len(progress)==0)exit
            first=1
            do while(first<=len(progress))
                last=first-1+index(progress(first:),new_line('a'))
                call append_agent_event(state,progress(first:last-1))
                first=last+1
            end do
            state%shown_agent_bytes=state%shown_agent_bytes+len(progress)
            if(.not.draining)exit
        end do
    end subroutine collect_agent_events

    subroutine append_agent_event(state,line)
        type(assistant_t), intent(inout) :: state
        character(*), intent(in) :: line
        type(QJsonDocument) :: record
        type(QJsonObject) :: payload

        if(len(line)==0)return
        if(line(1:1)=='@')then
            record=QJsonDocument_fromJson(line(2:))
            if(record%isObject())then
                payload=record%object()
                call append_structured_agent_event(state,payload)
                call payload%delete()
            end if
            call record%delete()
        else if(line(1:1)==achar(9))then
            record=QJsonDocument_fromJson('{"text":'//line(2:)//'}')
            if(record%isObject())then
                payload=record%object()
                call state%transcript%append_reasoning(payload%stringValue('text'))
                call payload%delete()
            end if
            call record%delete()
        else
            if(index(line,' · Thinking…')>0)call show_assistant_progress(state,'Thinking…')
            call log_activity(state,line)
        end if
    end subroutine append_agent_event

    subroutine append_structured_agent_event(state,payload)
        type(assistant_t),intent(inout) :: state
        type(QJsonObject),intent(in) :: payload
        character(:),allocatable :: agent,title,status,detail,chunk,kind
        integer :: member

        kind=payload%stringValue('kind')
        agent=''
        member=0
        if(state%session%request%team_mode)then
            member=findloc(team_agents,trim(payload%stringValue('agent')),dim=1)
            if(member>0)agent=trim(team_agents(member))
        end if
        select case(kind)
        case('stream')
            state%timeline_stream=.true.
        case('reply')
            state%timeline_stream=.true.
            chunk=payload%stringValue('text')
            if(len(agent)==0)agent=main_agent(state)
            if(len(agent)==0.or.agent==lead_agent)then
                if(.not.state%reply_started.and.len(chat_body(chunk))>0)call note_reply_started(state)
            end if
            call state%transcript%append_reply_chunk(chunk,agent=agent)
        case('reasoning')
            call state%transcript%append_reasoning(payload%stringValue('text'),agent=agent)
        case('activity')
            call state%transcript%append_activity(payload%stringValue('text'),agent=agent)
        case('team')
            if(member==0)return
            status=trim(payload%stringValue('state'))
            if(len(status)==0)return
            state%active_agents(member)=status=='Running'.or.status=='Queued'
            detail=payload%stringValue('text')
            if(len(detail)>0)status=status//' · '//detail
            call state%transcript%append_activity(status,agent=agent)
        case('plan','tool')
            title=trim(payload%stringValue('title'))
            if(len(title)==0)return
            status=trim(payload%stringValue('status'))
            detail=payload%stringValue('detail')
            if(len(status)>0)title=title//' · '//status
            if(kind=='plan')then
                if(len(detail)>0)title=title//new_line('a')//detail
                call state%transcript%append_plan(title,agent=agent)
            else
                call state%transcript%append_tool(title,detail,agent=agent)
                if(len(agent)==0.or.agent==lead_agent)then
                    if(status=='Running'.or.status=='Queued')then
                        call show_assistant_progress(state,title)
                    else
                        call show_assistant_progress(state,'Working…')
                    end if
                end if
            end if
        case default
            return
        end select
    end subroutine append_structured_agent_event

    function main_agent(state) result(agent)
        type(assistant_t),intent(in) :: state
        character(:),allocatable :: agent

        agent=''
        if(state%session%request%team_mode)agent=lead_agent
    end function main_agent

    subroutine note_reply_started(state)
        type(assistant_t),intent(inout) :: state

        if(state%reply_started)return
        state%reply_started=.true.
        call state%transcript%append_activity('Receiving the reply',agent=main_agent(state))
        call show_assistant_progress(state,'Receiving reply…')
    end subroutine note_reply_started

    subroutine show_agent_permission(state)
        type(assistant_t),intent(inout) :: state
        character(:),allocatable :: id,summary,title
        integer(int64) :: ticket

        summary=pending_studio_mcp_approval(ticket)
        if(ticket/=0_int64)then
            if(state%studio_permission_ticket==ticket)return
            state%studio_permission_ticket=ticket
            state%permission_id=''
            title='ForCAD Studio approval'
        else
            state%studio_permission_ticket=0_int64
            call read_agent_permission(state%session%request,id,summary)
            if(len(id)==0.or.id==state%last_agent_decision)then
                if(state%permission_pending)call show_assistant_progress(state,'Working…',approval=.false.)
                call state%panel%permission_box%hide()
                state%permission_id=''
                return
            end if
            if(id==state%permission_id)return
            state%permission_id=id
            title='Agent permission'
            if(len(summary)==0)summary='The agent is requesting permission to use a tool.'
        end if
        call state%panel%permission_title%setText(title)
        call state%panel%permission_summary%setText(summary)
        call state%panel%permission_box%show()
        call show_assistant_progress(state,'Waiting for approval…',approval=.true.)
        call log_activity(state,title//' requested · '//summary)
    end subroutine show_agent_permission

    subroutine decide_agent_request(data,allowed)
        class(*),intent(inout),optional :: data
        logical,intent(in) :: allowed

        select type(state=>data)
        type is(assistant_t)
            if(.not.state%busy.or.state%stop_requested)return
            if(.not.state%panel%permission_box%isVisible())return
            if(state%studio_permission_ticket/=0_int64)then
                call decide_studio_mcp_approval(allowed,state%studio_permission_ticket)
                state%studio_permission_ticket=0_int64
            else if(send_agent_decision(state%session%request,state%permission_id,allowed))then
                state%last_agent_decision=state%permission_id
            else
                call show_assistant_status(state,'Could not send the permission decision. Try again.',warning=.true.)
                return
            end if
            call state%panel%permission_box%hide()
            call show_assistant_progress(state,'Working…',approval=.false.)
            call log_activity(state,trim(merge('Allowed once','Denied      ',allowed)))
        end select
    end subroutine decide_agent_request

    subroutine on_allow_agent_request(data)
        class(*),intent(inout),optional :: data

        call decide_agent_request(data,.true.)
    end subroutine on_allow_agent_request

    subroutine on_deny_agent_request(data)
        class(*),intent(inout),optional :: data

        call decide_agent_request(data,.false.)
    end subroutine on_deny_agent_request

    subroutine begin_studio_assistant_mcp_event(name,args)
        character(*),intent(in) :: name
        type(QJsonObject),intent(in) :: args

        if(.not.assistant%busy.or.assistant%stop_requested.or..not.associated(assistant%workspace))return
        call assistant%cursor_timer%stop()
        call collect_agent_events(assistant,drain=.true.)
        call assistant%transcript%append_tool(studio_tool_label(name)//' · Running',agent=main_agent(assistant))
        call show_assistant_progress(assistant,'ForCAD Studio · '//studio_tool_label(name))
        call set_studio_ai_tool_cursor(assistant%workspace,name,args)
    end subroutine begin_studio_assistant_mcp_event

    function studio_tool_label(name) result(label)
        character(*),intent(in) :: name
        character(:),allocatable :: label
        integer :: i

        label=name
        do i=1,len(label)
            if(label(i:i)=='_')label(i:i)=' '
        end do
        if(len(label)==0)return
        i=index('abcdefghijklmnopqrstuvwxyz',label(1:1))
        if(i>0)label(1:1)=achar(iachar('A')+i-1)
    end function studio_tool_label

    subroutine append_studio_assistant_mcp_event(name,result,args,from_assistant,committed)
        character(*),intent(in) :: name
        type(QJsonObject),intent(in) :: result,args
        logical,intent(in) :: from_assistant,committed
        character(:),allocatable :: summary,action,detail

        if(.not.from_assistant.or..not.assistant%busy.or..not.associated(assistant%workspace))return
        if(committed)assistant%committed_tools=assistant%committed_tools+1
        action=studio_tool_label(name)
        if(result%boolValue('ok',.false.))then
            summary=action//' · Done'
        else
            assistant%failed_tools=assistant%failed_tools+1
            summary=action//' · Failed'
        end if
        if(committed)summary=summary//' · Scene changed'
        detail=''
        if((name=='preview_refinement'.or.name=='preview_controls').and.result%boolValue('ok',.false.))then
            detail=result%stringValue('summary')
        end if
        if(len(detail)==0)detail=studio_tool_result(result)
        call assistant%transcript%append_tool(summary,detail,agent=main_agent(assistant))
        if(.not.assistant%stop_requested)then
            call set_studio_ai_tool_cursor(assistant%workspace,name,args,result)
            call show_assistant_progress(assistant,'Working…')
        end if
    end subroutine append_studio_assistant_mcp_event

    function studio_tool_result(result) result(detail)
        type(QJsonObject),intent(in) :: result
        character(:),allocatable :: detail
        integer,parameter :: preview_bytes=4800
        integer :: last

        if(result%boolValue('ok',.false.))then
            detail=trim(result%toJson(.false.))
        else
            detail=result%stringValue('error','ForCAD Studio could not complete this tool.')
        end if
        if(len(detail)<=preview_bytes)return
        last=preview_bytes
        do while(iand(iachar(detail(last+1:last+1)),192)==128)
            last=last-1
        end do
        detail=detail(:last)//new_line('a')//'… additional result data omitted'
    end function studio_tool_result

    subroutine offer_scene_continuation(state)
        type(assistant_t), intent(inout) :: state

        if(.not.associated(state%workspace))return
        if(state%requested_document/=state%workspace%document%document_revision)return
        call state%panel%continue_scene%show()
    end subroutine offer_scene_continuation

    subroutine on_continue_scene(data)
        class(*), intent(inout), optional :: data
        character(:), allocatable :: question

        select type(state=>data)
        type is(assistant_t)
            if(state%busy.or..not.associated(state%workspace))return
            if(state%requested_document/=state%workspace%document%document_revision)return
            if(len(chat_body(state%panel%prompt%toPlainText()))==0)then
                question='Inspect the current scene and continue from the work already completed. Reuse existing geometry and meshes.'
                if(allocated(state%session%request%question))then
                    if(len(question)+len(state%session%request%question)+28<=8000)then
                        question=question//new_line('a')//new_line('a')//'Original request: '//state%session%request%question
                    end if
                end if
                call state%panel%prompt%setPlainText(question)
            end if
            call state%panel%prompt%setFocus()
            call show_assistant_status(state,'Review your draft and send when ready.')
        end select
    end subroutine on_continue_scene

    subroutine on_stop(data)
        class(*), intent(inout), optional :: data
        select type(state=>data)
        type is(assistant_t)
            if (.not.state%busy.or.state%stop_requested) return
            state%stop_requested=.true.
            call state%cursor_timer%stop()
            if(associated(state%workspace))call hide_studio_ai_cursor(state%workspace)
            call state%panel%permission_box%hide()
            call stop_studio_mcp(.true.)
            call state%session%interrupt_worker()
            call show_assistant_progress(state,'Stopping…',approval=.false.)
            call update_controls(state)
            if (associated(state%workspace)) call sync_studio_assistant(state%workspace)
        end select
    end subroutine on_stop

    subroutine on_worker_finished(data)
        class(*), intent(inout), optional :: data
        select type(state=>data)
        type is(assistant_t)
            if (state%listing_models) then
                call finish_model_list(state)
            else
                call finish_reply(state)
                if (state%refresh_after_reply) then
                    state%refresh_after_reply=.false.
                    call on_refresh_models(assistant)
                end if
            end if
        end select
    end subroutine on_worker_finished

    subroutine on_stream_progress(data)
        class(*), intent(inout), optional :: data
        character(len=:), allocatable :: progress
        select type(state=>data)
        type is(assistant_t)
            if(state%busy)call update_elapsed(state)
            if (.not.state%busy.or.state%stop_requested) return
            if(.not.state%transcript%defer_update())call collect_agent_events(state)
            call show_agent_permission(state)
            if(state%timeline_stream.or.state%transcript%defer_update())return
            progress=read_assistant_progress(state%session%request)
            if(len(progress)==0)return
            call note_reply_started(state)
            call state%transcript%update_reply(progress)
        end select
    end subroutine on_stream_progress

    subroutine start_elapsed(state)
        type(assistant_t), intent(inout) :: state

        call system_clock(count=state%elapsed_started)
        state%elapsed_seconds=-1
        call update_elapsed(state)
    end subroutine start_elapsed

    subroutine update_elapsed(state)
        type(assistant_t), intent(inout) :: state
        integer(int64) :: now,rate,seconds

        call system_clock(count=now,count_rate=rate)
        seconds=0
        if(rate>0)seconds=max(0_int64,now-state%elapsed_started)/rate
        if(seconds==state%elapsed_seconds)return
        state%elapsed_seconds=seconds
        if(.not.associated(state%workspace).or.state%requested_document==-1_int64)return
        call state%workspace%status%setAssistantActivity(elapsed=elapsed_caption(seconds))
    end subroutine update_elapsed

    function elapsed_caption(seconds) result(text)
        integer(int64),intent(in) :: seconds
        character(:),allocatable :: text
        character(24) :: caption

        write(caption,'(i0,a,i2.2)')seconds/60,':',mod(seconds,60_int64)
        text=trim(caption)
    end function elapsed_caption

    subroutine on_transcript_tail(data)
        class(*), intent(inout), optional :: data

        select type(state=>data)
        type is(assistant_t)
            call state%transcript%finish_follow()
        end select
    end subroutine on_transcript_tail

    subroutine resolve_team_activity(state)
        type(assistant_t),intent(inout) :: state
        character(:),allocatable :: message
        integer :: member

        if(.not.state%session%request%team_mode)return
        message='Ended before reporting a result.'
        if(state%stop_requested.or.state%session%request%cancelled)message='Stopped.'
        do member=1,size(team_agents)
            if(state%active_agents(member))call state%transcript%append_activity(message,agent=trim(team_agents(member)))
        end do
        state%active_agents=.false.
    end subroutine resolve_team_activity

    subroutine finish_reply(state)
        type(assistant_t), intent(inout) :: state
        logical :: completed, stopped, follow, keep_cursor
        character(:),allocatable :: validation_error,answer,turn_note

        stopped=state%stop_requested.or.state%session%request%cancelled
        keep_cursor=len(state%session%request%mcp_session)>0.and..not.stopped
        follow=state%transcript%at_end()
        call stop_studio_mcp(.true.)
        call state%panel%permission_box%hide()
        call collect_agent_events(state,drain=.true.)
        call show_assistant_progress(state,'',approval=.false.)
        call resolve_team_activity(state)
        call update_elapsed(state)
        call state%progress_timer%stop()
        associate(request=>state%session%request)
            if(len(request%provider_session)>0)state%conversation%provider_session=request%provider_session
            completed=.not.stopped.and.len_trim(request%error)==0.and.len(chat_body(request%reply))>0
            if(completed)then
                answer=request%reply
            else
                answer=read_assistant_progress(request)
            end if
            if(.not.state%timeline_stream.or..not.state%reply_started)call state%transcript%update_reply(answer)

            if(completed)then
                call show_assistant_status(state,'Ready for a follow-up.',record=.false.)
            else if(stopped)then
                turn_note='ForCAD Studio: The previous request was stopped before a final answer.'
                call state%transcript%append_message('Studio','Request stopped.')
                if(state%requested_document/=-1_int64)call show_assistant_status(state,'Stopped.')
            else if(len_trim(request%error)>0)then
                turn_note='ForCAD Studio: The provider could not complete the previous request.'
                call show_assistant_status(state,request%error,warning=.true.)
            else
                turn_note='ForCAD Studio: The provider did not finish an answer to the previous request.'
                call show_assistant_status(state,'No answer received. You can retry the request.',warning=.true.)
            end if

            if(.not.completed)then
                if(.not.stopped.and.state%committed_tools==0.and.len_trim(state%panel%prompt%toPlainText())==0)then
                    call state%panel%prompt%setPlainText(request%question)
                end if
                if(state%committed_tools>0)then
                    turn_note=turn_note//' ForCAD Studio confirmed '//count_text(state%committed_tools)// &
                        ' tool calls changed the scene during that turn. Inspect the current scene before resuming those edits.'
                end if
                if(len(chat_body(answer))>0)then
                    turn_note=turn_note//new_line('a')//new_line('a')//'Partial assistant reply:'//new_line('a')//answer
                end if
                answer=turn_note
            end if
            call state%session%remember_turn(request%question,answer)
        end associate
        if(state%committed_tools>0.and.(.not.completed.or.state%failed_tools>0))then
            validation_error=count_text(state%committed_tools)//' tool calls changed the scene. Applied changes are still present.'
            call state%transcript%append_message('Studio',validation_error)
            if(state%requested_document/=-1_int64)call show_assistant_status(state,'Some changes are applied.')
            call offer_scene_continuation(state)
        end if
        call log_activity(state,'Elapsed time · '//elapsed_caption(state%elapsed_seconds))
        if(associated(state%workspace))call state%workspace%status%setAssistantActivity(elapsed='')
        call cleanup_assistant_stream(state%session%request)
        call state%transcript%clear_reply()
        state%studio_permission_ticket=0_int64
        state%busy=.false.
        call update_controls(state)
        if (associated(state%workspace)) then
            if(keep_cursor)then
                call state%cursor_timer%start(cursor_feedback_ms)
            else
                call state%cursor_timer%stop()
                call hide_studio_ai_cursor(state%workspace)
            end if
            call sync_studio_assistant(state%workspace)
            if (follow) call state%transcript%follow()
        end if
    end subroutine finish_reply

    function assistant_project_folder(state) result(folder)
        type(assistant_t),intent(in) :: state
        character(:),allocatable :: folder

        folder=''
        if(.not.associated(state%workspace))return
        if(.not.allocated(state%workspace%document%project_directory))return
        folder=canonical_project_folder(state%workspace%document%project_directory)
    end function assistant_project_folder

    logical function ensure_project_folder(state) result(ready)
        type(assistant_t),intent(inout) :: state
        character(:),allocatable :: folder

        ready=.false.
        if(.not.associated(state%workspace))return
        if(allocated(state%workspace%document%project_directory))then
            if(QDir_exists(state%workspace%document%project_directory).and.len(state%workspace%document%project_directory)>0)then
                ready=.true.
                return
            end if
        end if
        folder=choose_project_folder(state%workspace%window)
        if(len(folder)==0)then
            call show_assistant_status(state,'Choose a project folder to use the assistant.')
            return
        end if
        call set_document_project_directory(state%workspace%document,folder)
        ready=.true.
    end function ensure_project_folder

    subroutine finish_studio_assistant()
        call stop_studio_mcp(.true.)
        if(assistant%panel%dock%isNull())return
        if (.not.assistant%progress_timer%isNull()) call assistant%progress_timer%stop()
        if (.not.assistant%cursor_timer%isNull()) call assistant%cursor_timer%stop()
        if (.not.assistant%transcript_tail_timer%isNull()) call assistant%transcript_tail_timer%stop()
        if (associated(assistant%workspace)) call hide_studio_ai_cursor(assistant%workspace)
        call assistant%session%close_worker()
        call clear_attachment(assistant)
        call clear_draft_selection(assistant)
        nullify(assistant%layout_changed)
        call assistant%transcript%release()
        call assistant%panel%dock%delete()
        if(associated(assistant%workspace))then
            call assistant%workspace%status%setAssistantActivity('',elapsed='',approval=.false.,chat='')
            call assistant%workspace%assistant_overlay%delete()
            if(.not.assistant%workspace%assistant_space%isNull())call assistant%workspace%assistant_space%hide()
        end if
        nullify(assistant%workspace)
        call assistant%session%clear_turns()
        assistant%busy=.false.
        assistant%stop_requested=.false.
        assistant%listing_models=.false.
        assistant%model_ready=.false.
        assistant%refresh_after_reply=.false.
        call assistant%conversation%clear()
        assistant%studio_permission_ticket=0_int64
        assistant%permission_pending=.false.
        assistant%studio_access=.true.
        if (allocated(assistant%session%request%context)) deallocate(assistant%session%request%context)
    end subroutine finish_studio_assistant

    subroutine on_cursor_expired(data)
        class(*), intent(inout), optional :: data

        select type(state=>data)
        type is(assistant_t)
            if(associated(state%workspace))then
                call hide_studio_ai_cursor(state%workspace)
                call sync_studio_assistant(state%workspace)
            end if
        end select
    end subroutine on_cursor_expired

    subroutine ask_studio_diagnostic(s)
        type(studio_t),target,intent(inout) :: s
        character(:),allocatable :: issue,draft

        if(.not.assistant_owns_workspace(s))return
        issue=studio_diagnostic_summary(s)
        if(len(issue)==0.or.assistant%panel%dock%isNull().or.assistant%busy)return
        draft=assistant%panel%prompt%toPlainText()
        if(len_trim(draft)==0)draft='Explain the current ForCAD diagnostic and propose a repair.'
        call assistant%panel%prompt%setPlainText(draft//new_line('a')//new_line('a')//'ForCAD diagnostics'//new_line('a')//issue)
        call show_assistant_status(assistant,'Review the diagnostic in the draft before sending.')
        call assistant%panel%prompt%setFocus()
    end subroutine ask_studio_diagnostic

    subroutine log_activity(state,entry)
        type(assistant_t),intent(inout) :: state
        character(*),intent(in) :: entry

        if(len_trim(entry)==0)return
        call state%transcript%append_activity(chat_body(without_nul(trim(entry))))
    end subroutine log_activity

    function studio_assistant_text() result(text)
        character(:), allocatable :: text

        text=assistant%transcript%plain_text()
    end function studio_assistant_text

end module studio_assistant
