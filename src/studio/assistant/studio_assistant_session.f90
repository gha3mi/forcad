module studio_assistant_session
    use qt, only: QThread, QThread_create, QJsonDocument, QJsonDocument_fromJson, QJsonObject, QJsonArray
    use studio_assistant_transport, only: assistant_request_t, run_assistant_request, &
        prepare_assistant_stream, cleanup_assistant_stream
    use studio_presentation_text, only: without_nul
    implicit none(type, external)
    private

    public :: assistant_session_t, agent_selection, agent_model_entry

    type :: agent_selection
        character(:), allocatable :: model, effort
    end type agent_selection

    type :: provider_preference
        character(:), allocatable :: provider
        type(agent_selection) :: choice
    end type provider_preference

    type :: agent_option
        character(:), allocatable :: id, label
    end type agent_option

    type :: agent_model_entry
        character(:), allocatable :: id, label
        type(agent_option), allocatable :: efforts(:)
        integer :: vision=-1
    end type agent_model_entry

    type :: assistant_session_t
        type(assistant_request_t) :: request, model_request
        type(agent_model_entry), allocatable :: agent_models(:)
        type(provider_preference), allocatable, private :: preferred_agent(:)
        character(:), allocatable, private :: agent_catalog
        type(QThread), private :: worker
        character(:), allocatable, private :: history
    contains
        procedure :: prepare_model_list
        procedure :: prepare_reply
        procedure :: prepare_worker
        procedure :: interrupt_worker
        procedure :: close_worker
        procedure, private :: release_worker
        procedure :: remember_model_catalog
        procedure :: select_agent_models
        procedure :: read_agent_models
        procedure :: preferred_agent_selection
        procedure :: remember_agent_selection
        procedure :: agent_choice
        procedure :: remember_turn
        procedure :: serialize_history
        procedure :: restore_turns
        procedure :: clear_turns
    end type assistant_session_t

contains

    subroutine prepare_model_list(session)
        class(assistant_session_t), intent(inout) :: session

        session%model_request%list_only=.true.
        session%model_request%models=''
        session%model_request%error=''
    end subroutine prepare_model_list

    subroutine prepare_reply(session,directory)
        class(assistant_session_t), intent(inout) :: session
        character(*), intent(in) :: directory

        session%request%error=''
        session%request%working_directory=directory
        call prepare_assistant_stream(session%request)
    end subroutine prepare_reply

    function prepare_worker(session,list_models) result(worker)
        class(assistant_session_t), intent(inout), target :: session
        logical, intent(in) :: list_models
        type(QThread) :: worker

        call session%release_worker()
        if(list_models)then
            session%worker=QThread_create(run_assistant_request,session%model_request)
        else
            session%worker=QThread_create(run_assistant_request,session%request)
        end if
        worker=session%worker
    end function prepare_worker

    subroutine interrupt_worker(session)
        class(assistant_session_t), intent(inout) :: session

        if(.not.session%worker%isNull())call session%worker%requestInterruption()
    end subroutine interrupt_worker

    subroutine close_worker(session)
        class(assistant_session_t), intent(inout) :: session

        call session%interrupt_worker()
        call session%release_worker()
        call cleanup_assistant_stream(session%request)
    end subroutine close_worker

    subroutine release_worker(session)
        class(assistant_session_t), intent(inout) :: session

        if(session%worker%isNull())return
        call session%worker%wait()
        call session%worker%delete()
    end subroutine release_worker

    subroutine remember_model_catalog(session)
        class(assistant_session_t), intent(inout) :: session

        session%agent_catalog=session%model_request%agent_catalog
    end subroutine remember_model_catalog

    subroutine select_agent_models(session,provider,selected)
        class(assistant_session_t), intent(inout) :: session
        character(*), intent(in) :: provider
        integer, intent(out) :: selected
        type(agent_selection) :: preference
        type(agent_model_entry), allocatable :: entries(:)
        integer :: i

        call session%read_agent_models(provider,entries)
        preference=session%preferred_agent_selection(provider)
        selected=0
        do i=0,ubound(entries,1)
            if(entries(i)%id==preference%model)selected=i
        end do
        call move_alloc(entries,session%agent_models)
    end subroutine select_agent_models

    subroutine read_agent_models(session,provider,entries)
        class(assistant_session_t), intent(in) :: session
        character(*), intent(in) :: provider
        type(agent_model_entry), allocatable, intent(out) :: entries(:)
        type(QJsonDocument) :: document
        type(QJsonObject) :: catalog,entry
        type(QJsonArray) :: models
        character(:), allocatable :: json,id,label
        integer :: i

        json='{}'
        if(allocated(session%agent_catalog))json=session%agent_catalog
        document=QJsonDocument_fromJson(json)
        catalog=document%object()
        models=catalog%arrayValue(provider)
        allocate(entries(0:max(1,models%size())-1))
        if(models%size()==0)then
            entries(0)%id=''
            entries(0)%label='Provider default'
            allocate(entries(0)%efforts(0))
        end if
        do i=0,models%size()-1
            entry=models%objectAt(i)
            id=entry%stringValue('id')
            label=entry%stringValue('label')
            if(len_trim(label)==0)then
                label=id
                if(len_trim(label)==0)label='Provider default'
            end if
            if(len(id)>0.and.entry%boolValue('vision'))label='Vision · '//label
            entries(i)%id=id
            entries(i)%label=label
            entries(i)%efforts=read_effort_options(entry)
            if(entry%contains('vision'))entries(i)%vision=merge(1,0,entry%boolValue('vision'))
            call entry%delete()
        end do
        call models%delete()
        call catalog%delete()
        call document%delete()
    end subroutine read_agent_models

    function read_effort_options(entry) result(efforts)
        type(QJsonObject), intent(in) :: entry
        type(agent_option), allocatable :: efforts(:)
        type(QJsonArray) :: list
        type(QJsonObject) :: effort
        integer :: i

        list=entry%arrayValue('efforts')
        allocate(efforts(list%size()))
        do i=1,size(efforts)
            effort=list%objectAt(i-1)
            efforts(i)%id=effort%stringValue('id')
            efforts(i)%label=effort%stringValue('label')
            if(len(efforts(i)%label)==0)efforts(i)%label=efforts(i)%id
            call effort%delete()
        end do
        call list%delete()
    end function read_effort_options

    function preferred_agent_selection(session,provider) result(choice)
        class(assistant_session_t), intent(in) :: session
        character(*), intent(in) :: provider
        type(agent_selection) :: choice
        integer :: i

        choice=agent_selection('','')
        if(.not.allocated(session%preferred_agent))return
        do i=1,size(session%preferred_agent)
            if(session%preferred_agent(i)%provider/=provider)cycle
            choice=session%preferred_agent(i)%choice
            return
        end do
    end function preferred_agent_selection

    subroutine remember_agent_selection(session,provider,choice)
        class(assistant_session_t), intent(inout) :: session
        character(*), intent(in) :: provider
        type(agent_selection), intent(in) :: choice
        integer :: i

        if(len_trim(provider)==0)return
        if(.not.allocated(session%preferred_agent))allocate(session%preferred_agent(0))
        do i=1,size(session%preferred_agent)
            if(session%preferred_agent(i)%provider/=provider)cycle
            session%preferred_agent(i)%choice=choice
            return
        end do
        session%preferred_agent=[session%preferred_agent,provider_preference(provider,choice)]
    end subroutine remember_agent_selection

    function agent_choice(session,model,effort) result(choice)
        class(assistant_session_t), intent(in) :: session
        integer, intent(in) :: model,effort
        type(agent_selection) :: choice

        choice=agent_selection('','')
        if(.not.allocated(session%agent_models))return
        if(model<0.or.model>ubound(session%agent_models,1))return
        choice%model=session%agent_models(model)%id
        if(effort>0.and.effort<=size(session%agent_models(model)%efforts))then
            choice%effort=session%agent_models(model)%efforts(effort)%id
        end if
    end function agent_choice

    subroutine remember_turn(session,question,answer)
        class(assistant_session_t), intent(inout) :: session
        character(*), intent(in) :: question,answer
        character(:),allocatable :: turn

        turn=without_nul(question)//achar(0)//without_nul(answer)
        if(.not.allocated(session%history))session%history=''
        if(len(session%history)>0)turn=achar(0)//turn
        session%history=session%history//turn
    end subroutine remember_turn

    function serialize_history(session) result(history)
        class(assistant_session_t), intent(in) :: session
        character(:), allocatable :: history

        history=''
        if(allocated(session%history))history=session%history
    end function serialize_history

    subroutine restore_turns(session,history)
        class(assistant_session_t), intent(inout) :: session
        character(*), intent(in) :: history

        session%history=history
    end subroutine restore_turns

    subroutine clear_turns(session)
        class(assistant_session_t), intent(inout) :: session

        if(allocated(session%history))deallocate(session%history)
    end subroutine clear_turns

end module studio_assistant_session
