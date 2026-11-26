submodule(studio_assistant) studio_assistant_models
    implicit none(type, external)

contains

    module procedure on_model_changed
        select type(state=>data)
        type is(assistant_t)
            if (state%model_ready.and.len_trim(text)>0) then
                state%selected_model=text
                state%resume_model=''
                call populate_agent_options(state)
                call set_status_hint(state%panel%model,text)
                call update_provider_badge(state)
            end if
            call update_controls(state)
        end select
    end procedure on_model_changed

    module procedure on_refresh_models
        select type(state=>data)
        type is(assistant_t)
            if (state%busy.or.state%listing_models) return
            state%listing_models=.true.
            call set_model_status(state,'Checking…','')
            call show_assistant_progress(state,'Checking agents…')
            call update_controls(state)
            call state%session%prepare_model_list()
            call start_worker(.true.)
        end select
    end procedure on_refresh_models

    module procedure set_model_status
        state%model_ready=.false.
        call state%panel%model%clear()
        call state%panel%model%addItem(text)
        call state%panel%model%setCurrentIndex(0)
        call set_status_hint(state%panel%model,text)
        call state%panel%agent_model%clear()
        call state%panel%agent_effort%clear()
        if(len_trim(detail)>0)call show_assistant_status(state,detail,warning=.true.)
    end procedure set_model_status

    module procedure finish_model_list
        character(len=:), allocatable :: name
        integer :: start, newline, item
        state%listing_models=.false.
        call show_assistant_progress(state,'')
        if (len_trim(state%session%model_request%error)>0) then
            call set_model_status(state,'Unavailable',state%session%model_request%error)
        else
            call state%panel%model%clear()
            start=1
            do while (start<=len(state%session%model_request%models))
                newline=index(state%session%model_request%models(start:),new_line('a'))
                if (newline==0) newline=len(state%session%model_request%models)-start+2
                name=trim(adjustl(state%session%model_request%models(start:start+newline-2)))
                if (len(name)>0) call state%panel%model%addItem(name)
                start=start+newline
            end do
            if (state%panel%model%count()==0) then
                call set_model_status(state,'No agents','Install or start an agent, then refresh the list.')
            else
                item=state%panel%model%findText(state%resume_model)
                if (item<0.and.len(state%resume_model)>0) then
                    call require_replacement_model(state)
                else
                    if (item<0) item=state%panel%model%findText(state%selected_model)
                    call state%panel%model%setCurrentIndex(max(0,item))
                    state%selected_model=state%panel%model%currentText()
                    state%resume_model=''
                    call set_status_hint(state%panel%model,state%selected_model)
                end if
                state%model_ready=.true.
            end if
        end if
        call state%session%remember_model_catalog()
        call populate_agent_options(state)
        call update_provider_badge(state)
        if(state%model_ready.and.state%selected_model==state%conversation%provider)call check_saved_agent_selection(state)
        call update_controls(state)
        if(associated(state%workspace))call sync_studio_assistant(state%workspace)
    end procedure finish_model_list

    module procedure populate_agent_options
        character(:),allocatable :: provider
        integer :: i,position

        provider=state%panel%model%currentText()
        if(len_trim(provider)==0)then
            state%updating_agent_options=.true.
            call state%panel%agent_model%clear()
            call state%panel%agent_effort%clear()
            state%updating_agent_options=.false.
            return
        end if
        state%updating_agent_options=.true.
        call state%panel%agent_model%clear()
        call state%session%select_agent_models(provider,position)
        do i=0,ubound(state%session%agent_models,1)
            call state%panel%agent_model%addItem(state%session%agent_models(i)%label)
        end do
        call state%panel%agent_model%setCurrentIndex(position)
        call populate_effort_options(state,position)
        state%updating_agent_options=.false.
    end procedure populate_agent_options

    subroutine populate_effort_options(state,model_index)
        type(assistant_t),intent(inout) :: state
        integer,intent(in) :: model_index
        character(:),allocatable :: previous
        type(agent_selection) :: preference
        integer :: effort_index,i

        call state%panel%agent_effort%clear()
        call state%panel%agent_effort%addItem('Default')
        preference=state%session%preferred_agent_selection(state%panel%model%currentText())
        previous=preference%effort
        effort_index=0
        associate(model=>state%session%agent_models(model_index))
            do i=1,size(model%efforts)
                call state%panel%agent_effort%addItem(model%efforts(i)%label)
                if(model%efforts(i)%id==previous)effort_index=i
            end do
        end associate
        call state%panel%agent_effort%setCurrentIndex(effort_index)
        call set_status_hint(state%panel%agent_model,state%panel%agent_model%currentText())
        if(len(state%session%agent_models(model_index)%id)==0)then
            call set_status_hint(state%panel%agent_model,'Uses the provider configuration. Refresh to retry model discovery.')
        end if
        if(state%panel%agent_effort%count()==1)then
            call set_status_hint(state%panel%agent_effort,'This model uses the provider’s default effort.')
        else
            call set_status_hint(state%panel%agent_effort,'Reasoning effort · '//state%panel%agent_effort%currentText())
        end if
    end subroutine populate_effort_options

    module procedure current_agent_selection

        choice=agent_selection('','')
        if(len_trim(state%panel%model%currentText())==0)return
        choice=state%session%agent_choice(state%panel%agent_model%currentIndex(),state%panel%agent_effort%currentIndex())
    end procedure current_agent_selection

    module procedure on_agent_model_changed

        select type(state=>data)
        type is(assistant_t)
            if(state%updating_agent_options.or.index<0)return
            state%updating_agent_options=.true.
            call populate_effort_options(state,index)
            state%updating_agent_options=.false.
            call on_agent_effort_changed(state%panel%agent_effort%currentIndex(),state)
            call update_controls(state)
        end select
    end procedure on_agent_model_changed

    module procedure on_agent_effort_changed

        select type(state=>data)
        type is(assistant_t)
            if(state%updating_agent_options.or.index<0)return
            call state%session%remember_agent_selection(state%panel%model%currentText(),current_agent_selection(state))
            call set_status_hint(state%panel%agent_effort,'Reasoning effort · '//state%panel%agent_effort%currentText())
        end select
    end procedure on_agent_effort_changed

    module procedure check_saved_agent_selection
        type(agent_selection) :: choice

        choice=current_agent_selection(state)
        if(choice%model/=state%conversation%agent_model)then
            call show_assistant_status(state,'Saved model unavailable. Current model: '//state%panel%agent_model%currentText()//'.',warning=.true.)
        else if(choice%effort/=state%conversation%agent_effort)then
            call show_assistant_status(state,'Saved effort unavailable. The model default is selected.',warning=.true.)
        end if
    end procedure check_saved_agent_selection

    module procedure require_replacement_model
        call state%panel%model%setCurrentIndex(-1)
        call set_status_hint(state%panel%model,'Saved agent · '//state%resume_model)
        call show_assistant_status(state,'Saved agent unavailable. Choose an available agent.',warning=.true.)
    end procedure require_replacement_model

    module procedure update_provider_badge

        if(.not.state%model_ready)then
            call set_status_hint(state%panel%provider_label,'Choose an available agent.')
            call state%transcript%set_empty_message('Choose an agent to start a chat.')
        else
            call set_status_hint(state%panel%provider_label,'Uses the selected agent and its configured service.')
            call state%transcript%set_empty_message('Ask a question or describe what you want to build.')
        end if
    end procedure update_provider_badge

end submodule studio_assistant_models
