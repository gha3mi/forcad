submodule(studio_assistant) studio_assistant_chat_ui
    use qt, only: QMessageBox_Discard, QMessageBox_Save
    use studio_assistant_panel, only: choose_saved_chat, confirm_chat_replacement
    use studio_chat_attachment, only: choose_studio_attachments,discard_attachment_copy
    use studio_chat_history, only: chat_record_t, chat_summary_t, project_chat_file, load_chat, list_chats
    implicit none(type, external)
contains

    module procedure on_new_chat
        select type(state=>data)
        type is(assistant_t)
            if (state%busy) return
            if (.not.can_replace_chat(state)) return
            call reset_chat(state)
            call show_assistant_status(state,'')
            call state%panel%prompt%setFocus()
        end select
    end procedure on_new_chat

    module procedure reset_chat
        integer :: model_index
        call state%session%clear_turns()
        if (allocated(state%session%request%context)) deallocate(state%session%request%context)
        call clear_attachment(state)
        call clear_draft_selection(state)
        call state%panel%continue_scene%hide()
        call state%cursor_timer%stop()
        state%stop_requested=.false.
        call state%conversation%clear()
        state%resume_model=''
        if (associated(state%workspace)) call state%workspace%status%setAssistantActivity(elapsed='',approval=.false.)
        if (state%model_ready.and.state%panel%model%currentIndex()<0) then
            model_index=state%panel%model%findText(state%selected_model)
            call state%panel%model%setCurrentIndex(max(0,model_index))
        end if
        call state%transcript%restore('')
        call state%panel%prompt%clear()
        call update_chat_status(state)
        call update_controls(state)
        if (associated(state%workspace)) call sync_studio_assistant(state%workspace)
    end procedure reset_chat

    module procedure update_chat_status
        if (.not.associated(state%workspace)) return
        if (state%conversation%dirty) then
            call state%workspace%status%setAssistantActivity(chat='Chat unsaved')
        else if (len(state%conversation%id)>0) then
            call state%workspace%status%setAssistantActivity(chat='Chat saved')
        else
            call state%workspace%status%setAssistantActivity(chat='')
        end if
    end procedure update_chat_status

    module procedure on_save_chat
        logical :: saved
        select type(state=>data)
        type is(assistant_t)
            saved=save_current_chat(state)
        end select
    end procedure on_save_chat

    module procedure on_chat_history
        type(chat_record_t) :: record
        type(chat_summary_t), allocatable :: items(:)
        character(len=:), allocatable :: id,error,project,path,active_id
        logical :: deleted_active,deleted
        select type(state=>data)
        type is(assistant_t)
            if (state%busy) return
            if(.not.ensure_project_folder(state))return
            project=assistant_project_folder(state)
            path=project_chat_file(project)
            call list_chats(items,error,path)
            if (len(error)>0) call show_assistant_status(state,error,warning=.true.)
            if (size(items)==0) then
                if (len(error)==0) call show_assistant_status(state,'No saved chats yet. Use Save chat to keep the current chat.')
                return
            end if
            active_id=''
            if(state%conversation%project==project)active_id=state%conversation%id
            id=choose_saved_chat(state%panel%dock,items,path,active_id,deleted_active,deleted)
            if(deleted_active.and.state%conversation%project==project.and.state%conversation%id==active_id)then
                state%conversation%id=''
                state%conversation%dirty=.true.
                call update_chat_status(state)
                call update_controls(state)
                call show_assistant_status(state,'Saved chat deleted. The open chat is now unsaved.')
            else if (deleted) then
                call show_assistant_status(state,'Saved chat deleted.')
            end if
            if (len(id)==0) return
            if(state%conversation%project==project.and.id==state%conversation%id)then
                call show_assistant_status(state,'This chat is already open.')
                return
            end if
            call load_chat(id,record,error,path)
            if (len(error)>0) then
                call show_assistant_status(state,error,warning=.true.)
                return
            end if
            if (.not.can_replace_chat(state)) return
            if(assistant_project_folder(state)/=project)then
                call show_assistant_status(state,'The project folder changed. Open Saved chats again.',warning=.true.)
                return
            end if
            call restore_chat(state,record,project)
        end select
    end procedure on_chat_history

    subroutine restore_chat(state,record,project)
        type(assistant_t), intent(inout) :: state
        type(chat_record_t), intent(in) :: record
        character(*), intent(in) :: project
        integer :: model_index
        call reset_chat(state)
        call state%conversation%restore(record,project)
        select case(state%conversation%provider)
        case('','Ollama','Codex CLI','Claude Code CLI','OpenCode CLI')
        case default
            state%conversation%agent_model=state%conversation%provider
            state%conversation%provider='Ollama'
        end select
        call state%session%remember_agent_selection(state%conversation%provider,agent_selection(state%conversation%agent_model,state%conversation%agent_effort))
        state%resume_model=state%conversation%provider
        if (state%model_ready.and.len(state%resume_model)>0) then
            model_index=state%panel%model%findText(state%conversation%provider)
            if (model_index>=0) then
                call state%panel%model%setCurrentIndex(model_index)
                state%resume_model=''
                call populate_agent_options(state)
            else
                call require_replacement_model(state)
            end if
        end if
        call state%session%restore_turns(record%history)
        call state%transcript%restore(record%transcript)
        call update_chat_status(state)
        call update_controls(state)
        call show_assistant_status(state,'Saved chat opened. The next message starts a new agent session.')
        if(state%model_ready.and.state%panel%model%currentIndex()>=0)call check_saved_agent_selection(state)
        call state%panel%prompt%setFocus()
    end subroutine restore_chat

    module procedure on_attach
        type(studio_attachment),allocatable :: choices(:)
        character(:),allocatable :: error,folder
        integer :: added,i

        select type(state=>data)
        type is(assistant_t)
            if(state%busy)return
            if(.not.ensure_project_folder(state))return
            folder=state%workspace%document%project_directory
            call choose_studio_attachments(state%workspace%window,choices,error,directory=folder)
            if(state%busy.or.assistant_project_folder(state)/=folder)then
                do i=1,size(choices)
                    call discard_attachment_copy(choices(i))
                end do
                return
            end if

            added=append_draft_attachments(state,choices,folder)
            if(added>0)then
                call update_controls(state)
                call state%panel%attachment_list%setCurrentRow(size(state%attachments)-1)
                call show_assistant_status(state,count_text(added)//trim(merge(' file attached. ',' files attached.',added==1)))
            else if(size(choices)>0)then
                call show_assistant_status(state,'These files are already attached.')
            end if
            if(len(error)>0)call show_assistant_status(state,error,warning=.true.)
        end select
    end procedure on_attach

    module procedure append_draft_attachments
        type(studio_attachment),allocatable :: combined(:)
        logical :: keep(size(choices))
        integer :: i,j,count_before

        keep=.true.
        do i=1,size(choices)
            if(allocated(state%attachments))then
                do j=1,size(state%attachments)
                    if(choices(i)%path/=state%attachments(j)%path)cycle
                    keep(i)=.false.
                    exit
                end do
            end if
            if(keep(i))then
                do j=1,i-1
                    if(.not.keep(j).or.choices(i)%path/=choices(j)%path)cycle
                    keep(i)=.false.
                    exit
                end do
            end if
            if(.not.keep(i))call discard_attachment_copy(choices(i))
        end do
        added=count(keep)
        if(added==0)return
        count_before=0
        if(allocated(state%attachments))count_before=size(state%attachments)
        allocate(combined(count_before+added))
        if(count_before>0)combined(:count_before)=state%attachments
        j=count_before
        do i=1,size(choices)
            if(.not.keep(i))cycle
            j=j+1
            combined(j)=choices(i)
        end do
        call move_alloc(combined,state%attachments)
        state%attachment_folder=folder
    end procedure append_draft_attachments

    module procedure on_detach
        type(studio_attachment),allocatable :: remaining(:)
        integer :: row,count

        select type(state=>data)
        type is(assistant_t)
            if(state%busy.or..not.allocated(state%attachments))return
            row=state%panel%attachment_list%currentRow()+1
            if(row<1.or.row>size(state%attachments))return
            call discard_attachment_copy(state%attachments(row))
            count=size(state%attachments)
            if(count==1)then
                deallocate(state%attachments)
                deallocate(state%attachment_folder)
            else
                allocate(remaining(count-1))
                remaining(:row-1)=state%attachments(:row-1)
                remaining(row:)=state%attachments(row+1:)
                call move_alloc(remaining,state%attachments)
            end if
            call update_controls(state)
        end select
    end procedure on_detach

    module procedure clear_attachment
        integer :: i

        if(allocated(state%attachments))then
            do i=1,size(state%attachments)
                call discard_attachment_copy(state%attachments(i))
            end do
            deallocate(state%attachments)
        end if
        if(allocated(state%attachment_folder))deallocate(state%attachment_folder)
        if(.not.state%panel%attachment_list%isNull())call state%panel%attachment_list%clear()
        if(.not.state%panel%attachment_bar%isNull())call state%panel%attachment_bar%hide()
    end procedure clear_attachment

    logical function save_current_chat(state) result(saved)
        type(assistant_t), intent(inout) :: state
        character(len=:), allocatable :: error
        saved=.false.
        if (state%busy.or..not.state%conversation%dirty) return

        call state%conversation%save(state%transcript%serialized(),state%session%serialize_history(),error)
        if (len(error)>0) then
            call show_assistant_status(state,error,warning=.true.)
            return
        end if
        call update_chat_status(state)
        call update_controls(state)
        call show_assistant_status(state,'Chat saved to its project folder.')
        saved=.true.
    end function save_current_chat
    logical function can_replace_chat(state) result(replace)
        type(assistant_t), intent(inout) :: state
        logical :: unsent
        integer :: choice
        unsent=len_trim(state%panel%prompt%toPlainText())>0
        replace=.not.state%conversation%dirty.and..not.unsent
        if (replace) return

        choice=confirm_chat_replacement(state%panel%dock,state%conversation%dirty,unsent)
        select case(choice)
        case(QMessageBox_Save)
            replace=save_current_chat(state)
        case(QMessageBox_Discard)
            replace=.true.
        case default
            replace=.false.
        end select
    end function can_replace_chat

end submodule studio_assistant_chat_ui
