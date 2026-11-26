module studio_assistant_conversation
    use studio_chat_history, only: chat_record_t, save_chat, project_chat_file
    implicit none(type, external)
    private

    public :: assistant_conversation_t

    type :: assistant_conversation_t
        character(:), allocatable :: id, title, provider, agent_model, agent_effort, provider_session, provider_workspace
        character(:),allocatable :: project
        logical :: dirty=.false.
    contains
        procedure :: clear
        procedure :: begin_turn
        procedure :: save
        procedure :: restore
    end type assistant_conversation_t

contains

    subroutine clear(conversation)
        class(assistant_conversation_t), intent(inout) :: conversation

        conversation%id=''
        conversation%title=''
        conversation%provider=''
        conversation%agent_model=''
        conversation%agent_effort=''
        conversation%provider_session=''
        conversation%provider_workspace=''
        conversation%project=''
        conversation%dirty=.false.
    end subroutine clear

    subroutine begin_turn(conversation,question,provider,agent_model,agent_effort,project)
        class(assistant_conversation_t), intent(inout) :: conversation
        character(*), intent(in) :: question,provider,agent_model,agent_effort,project

        conversation%project=project
        conversation%dirty=.true.
        if(conversation%provider/=provider)then
            conversation%provider_session=''
            conversation%provider_workspace=''
        end if
        conversation%provider=provider
        conversation%agent_model=agent_model
        conversation%agent_effort=agent_effort
        if(len(conversation%title)==0)conversation%title=title_from_prompt(question)
    end subroutine begin_turn

    subroutine save(conversation,transcript,history,error)
        class(assistant_conversation_t), intent(inout) :: conversation
        character(*), intent(in) :: transcript,history
        character(:), allocatable, intent(out) :: error
        type(chat_record_t) :: record
        character(:),allocatable :: path

        error='Could not open the conversation’s project folder.'
        path=project_chat_file(conversation%project)
        if(len(path)==0)return
        record%id=conversation%id
        record%title=conversation%title
        record%transcript=transcript
        record%history=history
        record%model=conversation%provider
        record%agent_model=conversation%agent_model
        record%agent_effort=conversation%agent_effort
        call save_chat(record,error,path)
        if(len(error)>0)return
        conversation%id=record%id
        conversation%dirty=.false.
    end subroutine save

    subroutine restore(conversation,record,project)
        class(assistant_conversation_t), intent(inout) :: conversation
        type(chat_record_t), intent(in) :: record
        character(*),intent(in) :: project

        conversation%project=project
        conversation%id=record%id
        conversation%title=record%title
        conversation%provider=record%model
        conversation%agent_model=record%agent_model
        conversation%agent_effort=record%agent_effort
        conversation%provider_session=''
        conversation%provider_workspace=''
        conversation%dirty=.false.
    end subroutine restore

    function title_from_prompt(question) result(title)
        character(*), intent(in) :: question
        character(:), allocatable :: title
        integer, parameter :: limit=84
        integer :: i,cut,byte

        title=trim(adjustl(question))
        do i=1,len(title)
            if(index(achar(9)//achar(10)//achar(13),title(i:i))>0)title(i:i)=' '
        end do
        if(len_trim(title)==0)title='Conversation'
        if(len(title)<=limit)return
        cut=limit+1
        do while(cut>1)
            byte=iachar(title(cut:cut))
            if(byte<128.or.byte>191)exit
            cut=cut-1
        end do
        title=trim(title(:cut-1))//'…'
    end function title_from_prompt

end module studio_assistant_conversation
