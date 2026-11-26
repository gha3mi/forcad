module studio_assistant_events
    use, intrinsic :: iso_c_binding, only: c_char,c_int,c_size_t,c_null_char
    use studio_assistant_system, only: append_agent_frame
    use qt, only: QJsonObject,QJsonArray,QJsonValue_String,QJsonValue_Array,QJsonValue_Double
    implicit none(type,external)
    private

    public :: assistant_events_t,relay_assistant_event,emit_assistant_activity,emit_assistant_reply
    public :: emit_assistant_thought,append_assistant_event,clip_assistant_text,clean_assistant_text,append_frame
    public :: assistant_failure_hint

    character,parameter :: lf=new_line('a')
    type :: event_note
        character(:),allocatable :: id,title,status,detail
        logical :: hidden=.false.
    end type event_note

    type :: assistant_events_t
        character(:),allocatable :: exchange,agent,failure_hint
        integer :: remaining=15000
        logical :: failed=.false.
        logical,private :: thinking=.false.,claude_partial=.false.,has_reply=.false.
        integer,allocatable,private :: thought_blocks(:)
        type(event_note),allocatable,private :: tools(:),plans(:)
    end type assistant_events_t

contains

    subroutine relay_assistant_event(events,provider,message)
        type(assistant_events_t),intent(inout) :: events
        character(*),intent(in) :: provider
        type(QJsonObject),intent(in) :: message
        type(QJsonObject) :: content

        select case(provider)
        case('Codex CLI')
            call codex_event(events,message)
        case('Claude Code CLI')
            call claude_event(events,message)
        case('OpenCode CLI')
            call opencode_event(events,message)
        case('Ollama')
            content=message%objectValue('message')
            call emit_assistant_thought(events,'Ollama',content%stringValue('thinking'))
            call emit_assistant_reply(events,content%stringValue('content'))
            call content%delete()
        case default
            return
        end select
    end subroutine relay_assistant_event

    subroutine append_assistant_event(events,event)
        type(assistant_events_t),intent(in) :: events
        type(QJsonObject),intent(inout) :: event

        if(allocated(events%agent))then
            if(len(events%agent)>0.and..not.event%contains('agent'))call event%insert('agent',events%agent)
        end if
        if(allocated(events%exchange))call append_frame(events%exchange//'/agent-events','@'//event%toJson(.true.)//lf)
    end subroutine append_assistant_event

    subroutine append_frame(path,frame,written)
        character(*),intent(in) :: path
        character(kind=c_char,len=*),intent(in) :: frame
        logical,optional,intent(out) :: written
        integer(c_int) :: status

        status=0
        if(len(frame)>0)status=append_agent_frame(path//c_null_char,frame,len(frame,c_size_t))
        if(present(written))written=status==0
    end subroutine append_frame

    subroutine emit_assistant_reply(events,text)
        type(assistant_events_t),intent(inout) :: events
        character(*),intent(in) :: text
        character(:),allocatable :: clean,chunk
        type(QJsonObject) :: event
        integer :: cursor

        if(events%remaining<=0)return
        clean=clip_assistant_text(clean_assistant_text(text),events%remaining)
        if(len(clean)==0)return
        events%remaining=events%remaining-len(clean)
        events%has_reply=.true.
        if(allocated(events%exchange))then
            call append_frame(events%exchange//'/progress',clean//achar(0))
            call append_frame(events%exchange//'/reply',clean)
        end if
        cursor=1
        do while(cursor<=len(clean))
            chunk=clip_assistant_text(clean(cursor:),2048)
            event=QJsonObject()
            call event%insert('kind','reply')
            call event%insert('text',chunk)
            call append_assistant_event(events,event)
            call event%delete()
            cursor=cursor+len(chunk)
        end do
    end subroutine emit_assistant_reply

    subroutine emit_assistant_activity(events,text)
        type(assistant_events_t),intent(in) :: events
        character(*),intent(in) :: text
        character(:),allocatable :: clean
        type(QJsonObject) :: event

        clean=activity_label(text)
        if(len_trim(clean)==0)return
        if(allocated(events%agent))then
            if(len(events%agent)>0)then
                event=QJsonObject()
                call event%insert('kind','activity')
                call event%insert('text',clean)
                call append_assistant_event(events,event)
                call event%delete()
                return
            end if
        end if
        if(allocated(events%exchange))call append_frame(events%exchange//'/agent-events',clean//lf)
    end subroutine emit_assistant_activity

    subroutine emit_assistant_thought(events,label,text)
        type(assistant_events_t),intent(inout) :: events
        character(*),intent(in) :: label,text
        character(:),allocatable :: clean,chunk,encoded
        type(QJsonObject) :: event
        type(QJsonArray) :: quoted
        integer :: cursor

        clean=clean_assistant_text(text,.true.)
        if(len(clean)==0)return
        if(.not.events%thinking.and.len_trim(clean)>0)then
            events%thinking=.true.
            call emit_assistant_activity(events,label//' · Thinking…')
        end if
        cursor=1
        do while(cursor<=len(clean))
            chunk=clip_assistant_text(clean(cursor:),2048)
            if(has_agent(events))then
                event=QJsonObject()
                call event%insert('kind','reasoning')
                call event%insert('text',chunk)
                call append_assistant_event(events,event)
                call event%delete()
            else if(allocated(events%exchange))then
                quoted=QJsonArray()
                call quoted%append(chunk)
                encoded=quoted%toJson(.true.)
                call quoted%delete()
                call append_frame(events%exchange//'/agent-events',achar(9)//encoded(2:len(encoded)-1)//lf)
            end if
            cursor=cursor+len(chunk)
        end do
    end subroutine emit_assistant_thought

    logical function has_agent(events)
        type(assistant_events_t),intent(in) :: events

        has_agent=.false.
        if(allocated(events%agent))has_agent=len(events%agent)>0
    end function has_agent

    function clean_assistant_text(text,terminal) result(clean)
        character(*),intent(in) :: text
        logical,intent(in),optional :: terminal
        character(:),allocatable :: clean
        integer :: i,used,code
        logical :: strip_terminal

        strip_terminal=.false.
        if(present(terminal))strip_terminal=terminal
        allocate(character(len(text))::clean)
        i=1
        used=0
        do while(i<=len(text))
            code=iachar(text(i:i))
            if(strip_terminal.and.code==27.and.i<len(text))then
                select case(text(i+1:i+1))
                case('[')
                    i=i+2
                    do while(i<=len(text))
                        code=iachar(text(i:i))
                        i=i+1
                        if(code>=64.and.code<=126)exit
                    end do
                    cycle
                case(']')
                    i=i+2
                    do while(i<=len(text))
                        code=iachar(text(i:i))
                        i=i+1
                        if(code==7)exit
                        if(code==27.and.i<=len(text))then
                            if(text(i:i)==achar(92))then
                                i=i+1
                                exit
                            end if
                        end if
                    end do
                    cycle
                case default
                end select
            end if
            if(code>=32.and.code/=127.or.code==9.or.code==10)then
                used=used+1
                clean(used:used)=text(i:i)
            end if
            i=i+1
        end do
        clean=clean(:used)
    end function clean_assistant_text

    function clip_assistant_text(text,limit) result(clipped)
        character(*),intent(in) :: text
        integer,intent(in) :: limit
        character(:),allocatable :: clipped
        integer :: n,code

        n=max(0,min(limit,len(text)))
        if(n<len(text))then
            do while(n>0)
                code=iachar(text(n+1:n+1))
                if(code<128.or.code>=192)exit
                n=n-1
            end do
        end if
        clipped=text(:n)
    end function clip_assistant_text

    function activity_label(text) result(label)
        character(*),intent(in) :: text
        character(:),allocatable :: label
        integer :: i

        label=clip_assistant_text(clean_assistant_text(text,.true.),320)
        do i=1,len(label)
            if(label(i:i)==lf.or.label(i:i)==achar(9))label(i:i)=' '
        end do
    end function activity_label

    function excerpt(text) result(output)
        character(*),intent(in) :: text
        character(:),allocatable :: output

        output=clean_assistant_text(text,.true.)
        if(len(output)>2400)output=clip_assistant_text(output,2396)//lf//'…'
        if(len_trim(output)==0)output=''
    end function excerpt

    subroutine emit_card(events,kind,id,title,status,detail)
        type(assistant_events_t),intent(in) :: events
        character(*),intent(in) :: kind,id,title,status,detail
        type(QJsonObject) :: event
        character(:),allocatable :: label

        label=activity_label(title)
        if(len_trim(label)==0)label=kind
        event=QJsonObject()
        call event%insert('kind',kind)
        call event%insert('id',clip_assistant_text(id,160))
        call event%insert('title',label)
        call event%insert('status',status)
        call event%insert('detail',excerpt(detail))
        call append_assistant_event(events,event)
        call event%delete()
    end subroutine emit_card

    integer function known_note(notes,id) result(at)
        type(event_note),allocatable,intent(in) :: notes(:)
        character(*),intent(in) :: id
        integer :: i

        at=0
        if(.not.allocated(notes).or.len(id)==0)return
        do i=1,size(notes)
            if(notes(i)%id/=id)cycle
            at=i
            return
        end do
    end function known_note

    subroutine emit_tool(events,id,title,status,detail)
        type(assistant_events_t),intent(inout) :: events
        character(*),intent(in) :: id,title,status,detail
        type(event_note) :: note
        integer :: at

        at=known_note(events%tools,id)
        note%title=activity_label(title)
        note%id=id
        note%status=status
        note%detail=excerpt(detail)
        if(at/=0)then
            if(events%tools(at)%hidden)return
            if(len_trim(note%title)==0)note%title=events%tools(at)%title
            if(note%title==events%tools(at)%title.and.note%status==events%tools(at)%status.and. &
                note%detail==events%tools(at)%detail)return
            events%tools(at)=note
        else if(len(id)>0)then
            if(len_trim(note%title)==0)note%title='Tool call'
            if(allocated(events%tools))then
                if(size(events%tools)<1024)events%tools=[events%tools,note]
            else
                events%tools=[note]
            end if
        end if
        if(len_trim(note%title)==0)note%title='Tool call'
        call emit_card(events,'tool',id,note%title,note%status,note%detail)
    end subroutine emit_tool

    subroutine hide_tool(events,id)
        type(assistant_events_t),intent(inout) :: events
        character(*),intent(in) :: id
        type(event_note) :: note
        integer :: at

        at=known_note(events%tools,id)
        if(at/=0)then
            events%tools(at)%hidden=.true.
        else if(len(id)>0)then
            note%id=id
            note%title=''
            note%detail=''
            note%status=''
            note%hidden=.true.
            if(allocated(events%tools))then
                if(size(events%tools)<1024)events%tools=[events%tools,note]
            else
                events%tools=[note]
            end if
        end if
    end subroutine hide_tool

    subroutine emit_plan(events,provider,id,entries,key,completion_flag)
        type(assistant_events_t),intent(inout) :: events
        character(*),intent(in) :: provider,id,key
        type(QJsonArray),intent(in) :: entries
        logical,intent(in),optional :: completion_flag
        type(event_note) :: note
        type(QJsonObject) :: entry
        character(:),allocatable :: title,status
        integer :: i,at,n,done,queued

        note%id=id
        note%title=provider//' · Plan'
        note%detail=''
        n=0
        done=0
        queued=0
        do i=0,min(entries%size(),128)-1
            entry=entries%objectAt(i)
            title=activity_label(entry%stringValue(key))
            if(present(completion_flag))then
                status='Queued'
                if(entry%boolValue('completed'))status='Done'
            else
                status=event_phase(entry%stringValue('status'),'Queued')
            end if
            call entry%delete()
            if(len_trim(title)==0)cycle
            n=n+1
            if(status=='Done')done=done+1
            if(status=='Queued')queued=queued+1
            call add_line(note%detail,status//' · '//title)
        end do
        at=known_note(events%plans,id)
        if(n==0.and.at==0)return
        note%status='Updated'
        if(n==done)note%status='Done'
        if(n==queued)note%status='Queued'
        if(n==0)note%status='Cleared'
        if(at/=0)then
            if(events%plans(at)%detail==note%detail.and.events%plans(at)%status==note%status)return
            events%plans(at)=note
        else if(allocated(events%plans))then
            if(size(events%plans)<256)events%plans=[events%plans,note]
        else
            events%plans=[note]
        end if
        call emit_card(events,'plan',id,note%title,note%status,note%detail)
    end subroutine emit_plan

    function event_phase(phase,fallback) result(label)
        character(*),intent(in) :: phase,fallback
        character(:),allocatable :: label

        select case(phase)
        case('pending','queued'); label='Queued'
        case('in_progress','running'); label='Running'
        case('completed'); label='Done'
        case('failed'); label='Failed'
        case('cancelled'); label='Cancelled'
        case('declined'); label='Declined'
        case default; label=fallback
        end select
    end function event_phase

    subroutine add_line(text,line)
        character(:),allocatable,intent(inout) :: text
        character(*),intent(in) :: line
        if(len_trim(line)==0)return
        if(len(text)>0)text=text//lf
        text=text//line
    end subroutine add_line

    function number_text(number) result(text)
        integer,intent(in) :: number
        character(:),allocatable :: text
        character(30) :: buffer

        write(buffer,'(i0)')number
        text=trim(buffer)
    end function number_text

    logical function studio_tool(name)
        character(*),intent(in) :: name

        studio_tool=index(name,'forcad_studio_')==1.or.index(name,'mcp__forcad_studio__')==1.or. &
            index(name,': forcad_studio_')>0.or.index(name,': mcp__forcad_studio__')>0
    end function studio_tool

    function visible_block(block) result(text)
        type(QJsonObject),intent(in) :: block
        character(:),allocatable :: text
        type(QJsonObject) :: annotations
        type(QJsonArray) :: audience
        integer :: i
        logical :: visible

        text=''
        if(block%stringValue('type')/='text')return
        annotations=block%objectValue('annotations')
        audience=annotations%arrayValue('audience')
        visible=audience%isEmpty()
        do i=0,audience%size()-1
            if(audience%stringAt(i)=='user')visible=.true.
        end do
        if(visible)text=block%stringValue('text')
        call audience%delete()
        call annotations%delete()
    end function visible_block

    function visible_content(object,key) result(text)
        type(QJsonObject),intent(in) :: object
        character(*),intent(in) :: key
        character(:),allocatable :: text
        type(QJsonArray) :: blocks
        type(QJsonObject) :: block
        integer :: i

        text=object%stringValue(key)
        if(object%valueType(key)/=QJsonValue_Array)return
        blocks=object%arrayValue(key)
        do i=0,blocks%size()-1
            block=blocks%objectAt(i)
            call add_line(text,visible_block(block))
            call block%delete()
            if(len(text)>=2400)exit
        end do
        call blocks%delete()
    end function visible_content

    function assistant_failure_hint(detail) result(hint)
        character(*),intent(in) :: detail
        character(:),allocatable :: hint,text
        integer :: i,byte

        text=clean_assistant_text(detail,.true.)
        do i=1,len(text)
            byte=iachar(text(i:i))
            if(byte>=iachar('A').and.byte<=iachar('Z'))text(i:i)=achar(byte+32)
        end do
        hint=''
        if(index(text,'free usage limit')>0)then
            hint='The selected model reached its provider’s free usage limit.'
        else if(index(text,'rate limit')>0.or.index(text,'too many requests')>0.or.index(text,'quota')>0.or. &
            index(text,'status 429')>0.or.index(text,'status code 429')>0)then
            hint='The provider rate-limited the request or its quota is exhausted.'
        else if(index(text,'authentication')>0.or.index(text,'unauthorized')>0.or.index(text,'invalid api key')>0.or. &
            index(text,'login required')>0.or.index(text,'status 401')>0.or.index(text,'status code 401')>0)then
            hint='The provider rejected authentication. Check the selected agent’s sign-in.'
        else if(index(text,'invalid transport')>0.or.index(text,'mcp handshake')>0.or. &
            index(text,'mcp startup')>0.or.index(text,'mcp client')>0)then
            hint='The ForCAD Studio tool connection could not start. Check the provider’s MCP configuration.'
        else if(index(text,'model not found')>0.or.index(text,'unsupported_model')>0.or.index(text,'model is not supported')>0)then
            hint='The provider does not support the selected model. Refresh the models and choose an available one.'
        else if(index(text,'unexpected argument')>0.or.index(text,'unknown option')>0.or.index(text,'loading config')>0.or. &
            index(text,'invalid configuration')>0)then
            hint='The provider rejected its startup options or configuration.'
        else if(index(text,'econnrefused')>0.or.index(text,'connection refused')>0.or.index(text,'standalone server exited')>0)then
            hint='The provider’s local service could not start or accept connections.'
        else if(index(text,'timed out')>0.or.index(text,'timeout')>0)then
            hint='The provider timed out before completing the request.'
        end if
    end function assistant_failure_hint

    subroutine note_assistant_failure(events,detail)
        type(assistant_events_t),intent(inout) :: events
        character(*),intent(in) :: detail
        character(:),allocatable :: hint

        hint=assistant_failure_hint(detail)
        if(len(hint)>0)events%failure_hint=hint
    end subroutine note_assistant_failure

    subroutine codex_event(events,message)
        type(assistant_events_t),intent(inout) :: events
        type(QJsonObject),intent(in) :: message
        type(QJsonObject) :: item,record,result,problem
        type(QJsonArray) :: entries
        character(:),allocatable :: kind,id,status,detail
        logical :: done
        integer :: i

        kind=message%stringValue('type')
        select case(kind)
        case('turn.failed')
            events%failed=.true.
            problem=message%objectValue('error')
            call note_assistant_failure(events,problem%stringValue('message'))
            call problem%delete()
            return
        case('error')
            call note_assistant_failure(events,message%stringValue('message'))
            return
        case default
        end select
        if(index(kind,'item.')/=1)return
        item=message%objectValue('item')
        done=kind=='item.completed'
        id=item%stringValue('id')
        status=event_phase(item%stringValue('status'),trim(merge('Done   ','Running',done)))
        select case(item%stringValue('type'))
        case('agent_message')
            if(done)then
                events%thinking=.false.
                call emit_assistant_reply(events,item%stringValue('text'))
            end if
        case('reasoning')
            if(done)then
                events%thinking=.false.
                call emit_assistant_thought(events,'Codex',item%stringValue('text'))
            end if
        case('todo_list')
            entries=item%arrayValue('items')
            call emit_plan(events,'Codex','codex-plan:'//id,entries,'text',.true.)
            call entries%delete()
        case('command_execution')
            events%thinking=.false.
            detail=clip_assistant_text(item%stringValue('command'),320)
            if(item%valueType('exit_code')==QJsonValue_Double)call add_line(detail,'Exit code '//number_text(item%intValue('exit_code')))
            call add_line(detail,item%stringValue('aggregated_output'))
            call emit_tool(events,id,'Codex · Command',status,detail)
        case('mcp_tool_call')
            if(item%stringValue('server')/='forcad_studio')then
                events%thinking=.false.
                problem=item%objectValue('error')
                result=item%objectValue('result')
                detail=problem%stringValue('message')
                call add_line(detail,visible_content(result,'content'))
                call emit_tool(events,id,'Codex · '//item%stringValue('tool','Tool call'),status,detail)
                call problem%delete()
                call result%delete()
            end if
        case('web_search')
            events%thinking=.false.
            detail=item%stringValue('query')
            entries=item%arrayValue('results')
            do i=0,entries%size()-1
                record=entries%objectAt(i)
                call add_line(detail,record%stringValue('title'))
                call record%delete()
                if(len(detail)>=2400)exit
            end do
            call emit_tool(events,id,'Codex · Web search',status,detail)
            call entries%delete()
        case('file_change')
            events%thinking=.false.
            detail=''
            entries=item%arrayValue('changes')
            do i=0,entries%size()-1
                record=entries%objectAt(i)
                if(record%valueType('path')==QJsonValue_String) then
                    call add_line(detail,record%stringValue('kind','Update')//' · '//record%stringValue('path'))
                end if
                call record%delete()
                if(len(detail)>=2400)exit
            end do
            call emit_tool(events,id,'Codex · File changes',status,detail)
            call entries%delete()
        case('collab_tool_call')
            events%thinking=.false.
            call emit_tool(events,id,'Codex · '//item%stringValue('tool','Agent collaboration'),status,'')
        case default
        end select
        call item%delete()
    end subroutine codex_event

    subroutine claude_event(events,message)
        type(assistant_events_t),intent(inout) :: events
        type(QJsonObject),intent(in) :: message
        type(QJsonObject) :: stream,block,delta,inner
        type(QJsonArray) :: content
        character(:),allocatable :: id
        integer :: i,at

        select case(message%stringValue('type'))
        case('stream_event')
            stream=message%objectValue('event')
            delta=stream%objectValue('delta')
            if(stream%stringValue('type')=='message_start')then
                call reset_claude_message(events)
            else if(delta%stringValue('type')=='text_delta')then
                events%claude_partial=.true.
                events%thinking=.false.
                call emit_assistant_reply(events,delta%stringValue('text'))
            else if(delta%stringValue('type')=='thinking_delta')then
                i=stream%intValue('index')
                if(len(delta%stringValue('thinking'))>0)then
                    if(allocated(events%thought_blocks))then
                        if(.not.any(events%thought_blocks==i))events%thought_blocks=[events%thought_blocks,i]
                    else
                        events%thought_blocks=[i]
                    end if
                    call emit_assistant_thought(events,'Claude Code',delta%stringValue('thinking'))
                end if
            else if(stream%stringValue('type')=='content_block_start')then
                block=stream%objectValue('content_block')
                if(block%stringValue('type')=='tool_use'.or.block%stringValue('type')=='server_tool_use') then
                    call claude_tool(events,block)
                end if
                call block%delete()
            end if
            call delta%delete()
            call stream%delete()
        case('assistant','user')
            inner=message%objectValue('message')
            content=inner%arrayValue('content')
            do i=0,content%size()-1
                block=content%objectAt(i)
                if(message%stringValue('type')=='assistant')then
                    call claude_block(events,block,i)
                else if(block%stringValue('type')=='tool_result')then
                    call claude_result(events,block)
                end if
                call block%delete()
            end do
            call content%delete()
            call inner%delete()
            if(message%stringValue('type')=='assistant')call reset_claude_message(events)
        case('tool_progress')
            id=message%stringValue('tool_use_id')
            at=known_note(events%tools,id)
            if(at/=0)then
                if(events%tools(at)%status=='Running'.and.message%valueType('elapsed_time_seconds')==QJsonValue_Double) then
                    call emit_tool(events,id,'','Running','Elapsed '//number_text(floor(message%doubleValue('elapsed_time_seconds')))//' s')
                end if
            end if
        case('result')
            if(message%boolValue('is_error'))then
                events%failed=.true.
                call note_assistant_failure(events,message%stringValue('result'))
                call note_assistant_failure(events,visible_content(message,'errors'))
            else if(.not.events%has_reply)then
                call emit_assistant_reply(events,message%stringValue('result'))
            end if
        case default
        end select
    end subroutine claude_event

    subroutine reset_claude_message(events)
        type(assistant_events_t),intent(inout) :: events

        events%claude_partial=.false.
        events%thinking=.false.
        if(allocated(events%thought_blocks))deallocate(events%thought_blocks)
    end subroutine reset_claude_message

    subroutine claude_block(events,block,i)
        type(assistant_events_t),intent(inout) :: events
        type(QJsonObject),intent(in) :: block
        integer,intent(in) :: i
        character(:),allocatable :: kind

        kind=block%stringValue('type')
        select case(kind)
        case('text')
            if(.not.events%claude_partial)then
                events%thinking=.false.
                call emit_assistant_reply(events,block%stringValue('text'))
            end if
        case('thinking')
            if(allocated(events%thought_blocks))then
                if(any(events%thought_blocks==i))return
            end if
            call emit_assistant_thought(events,'Claude Code',block%stringValue('thinking'))
        case('tool_use','server_tool_use')
            call claude_tool(events,block)
        case default
            if(len(kind)>=12)then
                if(kind(len(kind)-11:)=='_tool_result')call claude_result(events,block)
            end if
        end select
    end subroutine claude_block

    subroutine claude_tool(events,block)
        type(assistant_events_t),intent(inout) :: events
        type(QJsonObject),intent(in) :: block
        type(QJsonObject) :: input
        type(QJsonArray) :: todos
        character(:),allocatable :: name,id,status,detail

        name=block%stringValue('name','Tool call')
        if(studio_tool(name))return
        id=block%stringValue('id')
        input=block%objectValue('input')
        select case(name)
        case('TodoWrite')
            if(input%valueType('todos')==QJsonValue_Array)then
                todos=input%arrayValue('todos')
                call emit_plan(events,'Claude Code','claude-todos',todos,'content')
                call todos%delete()
            end if
        case('TaskCreate','TaskUpdate')
            if(.not.input%isEmpty())then
                status=event_phase(input%stringValue('status'),trim(merge('Queued ','Updated',name=='TaskCreate')))
                detail=input%stringValue('subject')
                if(len(detail)==0.and.input%contains('taskId'))then
                    detail=input%stringValue('taskId')
                    if(input%valueType('taskId')==QJsonValue_Double)detail=number_text(input%intValue('taskId'))
                    if(len(detail)>0)detail='Task '//detail
                end if
                call emit_card(events,'plan',id,'Claude Code · Plan',status,detail)
            end if
        case default
            events%thinking=.false.
            detail=input%stringValue('description',input%stringValue('command'))
            call emit_tool(events,id,'Claude Code · '//name,'Running',detail)
        end select
        call input%delete()
    end subroutine claude_tool

    subroutine claude_result(events,block)
        type(assistant_events_t),intent(inout) :: events
        type(QJsonObject),intent(in) :: block
        character(:),allocatable :: id,status

        id=block%stringValue('tool_use_id')
        if(known_note(events%tools,id)==0)return
        status=trim(merge('Failed','Done  ',block%boolValue('is_error')))
        events%thinking=.false.
        call emit_tool(events,id,'',status,visible_content(block,'content'))
    end subroutine claude_result

    subroutine opencode_event(events,message)
        type(assistant_events_t),intent(inout) :: events
        type(QJsonObject),intent(in) :: message
        type(QJsonObject) :: params,update,content,meta,entry
        type(QJsonArray) :: entries

        params=message%objectValue('params')
        update=params%objectValue('update')
        select case(update%stringValue('sessionUpdate'))
        case('agent_message_chunk','agent_thought_chunk')
            content=update%objectValue('content')
            if(content%stringValue('type')=='text')then
                if(update%stringValue('sessionUpdate')=='agent_message_chunk')then
                    events%thinking=.false.
                    call emit_assistant_reply(events,content%stringValue('text'))
                else
                    call emit_assistant_thought(events,'OpenCode',content%stringValue('text'))
                end if
            end if
            call content%delete()
        case('tool_call','tool_call_update')
            events%thinking=.false.
            call opencode_tool(events,update)
        case('plan')
            entries=update%arrayValue('entries')
            call emit_plan(events,'OpenCode','opencode-plan',entries,'content')
            call entries%delete()
        case('session_info_update')
            meta=update%objectValue('_meta')
            entry=meta%objectValue('opencode/retry')
            if(entry%valueType('attempt')==QJsonValue_Double) then
                call emit_assistant_activity(events,'OpenCode · Retry scheduled · Attempt '//number_text(entry%intValue('attempt')))
            end if
            call entry%delete()
            entry=meta%objectValue('opencode/compaction')
            if(entry%valueType('status')==QJsonValue_String) then
                call emit_assistant_activity(events,'OpenCode · Context summary · '//entry%stringValue('status'))
            end if
            call entry%delete()
            call meta%delete()
        case default
        end select
        call update%delete()
        call params%delete()
    end subroutine opencode_event

    subroutine opencode_tool(events,update)
        type(assistant_events_t),intent(inout) :: events
        type(QJsonObject),intent(in) :: update
        type(QJsonArray) :: blocks
        type(QJsonObject) :: block,content
        character(:),allocatable :: name,id,title,status,result
        integer :: i,at

        id=update%stringValue('toolCallId')
        if(len(id)==0)return
        name=activity_label(update%stringValue('name'))
        title=activity_label(update%stringValue('title'))
        if(studio_tool(title).or.studio_tool(name))then
            call hide_tool(events,id)
            return
        end if
        at=known_note(events%tools,id)
        status='Queued'
        if(at/=0)then
            if(events%tools(at)%hidden)return
            status=events%tools(at)%status
        else if(len_trim(title)==0)then
            title='Tool call'
        end if
        if(len_trim(title)>0)title='OpenCode · '//title
        status=event_phase(update%stringValue('status'),status)
        blocks=update%arrayValue('content')
        result=''
        do i=0,blocks%size()-1
            block=blocks%objectAt(i)
            if(block%stringValue('type')=='content')then
                content=block%objectValue('content')
                call add_line(result,visible_block(content))
                call content%delete()
            end if
            call block%delete()
            if(len(result)>=2400)exit
        end do
        call blocks%delete()
        if(len(result)==0.and.at/=0)result=events%tools(at)%detail
        call emit_tool(events,id,title,status,result)
    end subroutine opencode_tool

end module studio_assistant_events
