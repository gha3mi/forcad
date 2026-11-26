module studio_assistant_local
    use, intrinsic :: iso_fortran_env, only: real64
    use qt, only: qt_string,QJsonObject,QJsonArray,QJsonValue,QJsonValue_Object,QJsonValue_Array,QJsonValue_Bool, &
        QJsonValue_String,QJsonValue_Undefined,QJsonValue_Null
    use studio_assistant_events, only: assistant_events_t,relay_assistant_event,emit_assistant_activity
    use studio_assistant_guidance, only: studio_assistant_instructions
    use studio_chat_attachment, only: studio_attachment
    use studio_assistant_local_http, only: local_http_request,local_assistant_models,local_model_info,local_model_efforts, &
        valid_local_model_name,local_model_metadata,has_json_string,parse_json_object
    use studio_assistant_mcp_client, only: assistant_mcp_client_t,connect_assistant_mcp,request_assistant_mcp,close_assistant_mcp
    implicit none(type,external)
    private
    public :: local_assistant_t,prepare_local_assistant,run_local_assistant

    character(*),parameter :: newline=new_line('a'),provider='Ollama'
    integer,parameter :: default_context_size=16384

    type :: local_assistant_t
        logical :: live=.false.
        character(:),allocatable,private :: model,session,think_text
        integer,private :: context_size=default_context_size,think_kind=QJsonValue_Undefined
        logical,private :: think_boolean=.false., image_input=.true.
    end type local_assistant_t

    type :: local_turn_t
        character(:),allocatable :: pending,text,thinking
        type(QJsonArray) :: calls
        type(assistant_events_t),pointer :: events=>null()
        logical :: done=.false.,incomplete=.false.
    end type local_turn_t

contains

    subroutine prepare_local_assistant(local,model,effort,session,error)
        type(local_assistant_t),intent(out) :: local
        character(*),intent(in) :: model,effort,session
        character(:),allocatable,intent(out) :: error
        type(QJsonArray) :: models,capabilities
        type(QJsonObject) :: info

        local%session=session
        local%model=''
        error=''
        if(len(model)>0.and..not.valid_local_model_name(model))then
            error='Choose an installed local model from the provider catalog.'
            return
        end if
        call local_assistant_models(models,error)
        if(len(error)>0)then
            error='Could not load Ollama models from the local service. Check that it is running.'
        else if(len(model)==0)then
            if(.not.models%isEmpty())local%model=models%stringAt(0)
        else if(has_json_string(models,model))then
            local%model=model
        else if(has_json_string(models,model//':latest'))then
            local%model=model//':latest'
        end if
        call models%delete()
        if(len(error)>0)return
        if(len(local%model)==0)then
            error='The selected Ollama model is unavailable. Refresh the model list.'
            return
        end if
        call local_model_info(local%model,5000,info,error)
        if(len(error)>0.or.info%valueType('details')/=QJsonValue_Object.or..not.local_model_metadata(info))then
            error="Could not verify that '"//local%model//"' is a local Ollama model."
            call info%delete()
            return
        end if
        if(len(effort)>0)call select_local_effort(local,info,effort,error)
        capabilities=info%arrayValue('capabilities')
        if(info%valueType('capabilities')==QJsonValue_Array)local%image_input=has_json_string(capabilities,'vision')
        local%live=len(session)>0.and.has_json_string(capabilities,'tools')
        if(local%live)local%context_size=local_context_size(info)
        call capabilities%delete()
        call info%delete()
    end subroutine prepare_local_assistant

    subroutine select_local_effort(local,info,effort,error)
        type(local_assistant_t),intent(inout) :: local
        type(QJsonObject),intent(in) :: info
        character(*),intent(in) :: effort
        character(:),allocatable,intent(out) :: error
        type(QJsonArray) :: efforts
        character(:),allocatable :: text
        integer :: index_effort
        logical :: boolean,wants_boolean

        error=''
        wants_boolean=effort=='boolean:true'.or.effort=='boolean:false'
        boolean=effort=='boolean:true'
        text=effort
        if(index(effort,'string:')==1)text=effort(8:)
        efforts=local_model_efforts(info)
        do index_effort=0,efforts%size()-1
            if(wants_boolean)then
                if(efforts%typeAt(index_effort)/=QJsonValue_Bool)cycle
                if(efforts%boolAt(index_effort).neqv.boolean)cycle
                local%think_kind=QJsonValue_Bool
                local%think_boolean=boolean
            else
                if(efforts%typeAt(index_effort)/=QJsonValue_String)cycle
                if(efforts%stringAt(index_effort)/=text)cycle
                local%think_kind=QJsonValue_String
                local%think_text=text
            end if
            exit
        end do
        call efforts%delete()
        if(local%think_kind==QJsonValue_Undefined) then
          error='The local model does not offer that reasoning effort. Refresh the model list.'
        end if
    end subroutine select_local_effort

    integer function local_context_size(info) result(size_value)
        type(QJsonObject),intent(in) :: info
        type(QJsonObject) :: model_info
        type(qt_string),allocatable :: keys(:)
        character(:),allocatable :: parameters,line,key
        character(*),parameter :: suffix='.context_length',digits='0123456789'
        integer :: position,start,last,separation,configured,status,candidate
        real(real64) :: number,capacity

        configured=0
        parameters=info%stringValue('parameters')//newline
        start=1
        do while(start<=len(parameters))
            last=start+index(parameters(start:),newline)-2
            line=parameters(start:last)
            start=last+2
            if(index(line,'num_ctx')/=1.or.len(line)<=7)cycle
            if(line(8:8)/=' '.and.line(8:8)/=achar(9))cycle
            separation=verify(line(8:),' '//achar(9))
            if(separation==0)cycle
            line=line(7+separation:)
            separation=verify(line,digits)
            if(separation==1)cycle
            if(separation>0)line=line(:separation-1)
            read(line,*,iostat=status)candidate
            if(status==0)configured=max(configured,candidate)
        end do
        if(configured>=default_context_size)then
            size_value=configured
            return
        end if
        model_info=info%objectValue('model_info')
        keys=model_info%keys()
        capacity=0
        do position=1,size(keys)
            key=keys(position)%s
            if(len(key)<len(suffix))cycle
            if(key(len(key)-len(suffix)+1:)/=suffix)cycle
            number=model_info%doubleValue(key)
            if(number>=1)capacity=max(capacity,number)
        end do
        call model_info%delete()
        size_value=default_context_size
        if(capacity>=1)size_value=int(min(real(default_context_size,kind=real64),capacity))
    end function local_context_size

    subroutine run_local_assistant(local,question,history,context,attachments,team_report,team_worker,events,error,executable,role)
        type(local_assistant_t),intent(in) :: local
        character(*),intent(in) :: question,context,team_report
        type(QJsonArray),intent(in) :: history
        type(studio_attachment),intent(in) :: attachments(:)
        logical,intent(in) :: team_worker
        type(assistant_events_t),intent(inout),target :: events
        character(:),allocatable,intent(out) :: error
        character(*),optional,intent(in) :: executable,role
        type(assistant_mcp_client_t) :: client
        type(QJsonObject) :: request,options
        type(QJsonArray) :: tools,messages
        character(:),allocatable :: payload,body,instruction_role
        type(local_turn_t) :: turn
        integer :: round
        logical :: complete

        error=''
        tools=QJsonArray()
        if(local%live)then
            call connect_assistant_mcp(client,local%session,error,executable)
            if(len(error)==0)call load_local_tools(client,tools,error)
            if(len(error)>0)then
                call tools%delete()
                call close_assistant_mcp(client)
                error='The ForCAD Studio connection expired. Send the message again to reconnect.'
                return
            end if
        else if(len(local%session)>0)then
            call emit_assistant_activity(events,provider// &
                ' · This model supports chat only. Choose a model with tool calling for live ForCAD Studio access.')
        end if
        instruction_role=''
        if(len(team_report)>0)instruction_role='Lead'
        if(present(role))instruction_role=role
        messages=initial_local_messages(local,question,history,context,attachments,team_report,team_worker,instruction_role)
        request=QJsonObject()
        call request%insert('model',local%model)
        call request%insert('stream',.true.)
        select case(local%think_kind)
        case(QJsonValue_Bool)
            call request%insert('think',local%think_boolean)
        case(QJsonValue_String)
            call request%insert('think',local%think_text)
        end select
        if(.not.tools%isEmpty())then
            options=QJsonObject()
            call options%insert('num_ctx',local%context_size)
            call request%insert('tools',tools)
            call request%insert('options',options)
            call options%delete()
        end if
        call tools%delete()

        complete=.false.
        do round=1,64
            call request%insert('messages',messages)
            payload=request%toJson(.true.)
            turn%calls=QJsonArray()
            turn%pending=''
            turn%text=''
            turn%thinking=''
            turn%done=.false.
            turn%incomplete=.false.
            turn%events=>events
            call local_http_request('/api/chat',payload,180000,body,error,consume_turn_chunk,turn)
            if(len(error)==0.and.len_trim(turn%pending)>0)then
                if(.not.consume_turn_frame(turn,turn%pending))turn%incomplete=.true.
            end if
            if(turn%incomplete.or.(len(error)==0.and..not.turn%done)) then
              error='Ollama returned an incomplete reply. Check the design before retrying.'
            end if
            if(len(error)>0)then
                call turn%calls%delete()
                exit
            end if
            if(turn%calls%isEmpty())then
                call turn%calls%delete()
                complete=.true.
                exit
            end if
            if(local%live)then
                call run_local_tools(client,turn,messages,error)
            else
                error='The model requested ForCAD Studio tools without an active connection.'
            end if
            call turn%calls%delete()
            if(len(error)>0)exit
        end do
        if(.not.complete.and.len(error)==0) then
          error='The model reached the tool round limit. Continue from the current design in a new message.'
        end if
        call messages%delete()
        call request%delete()
        call close_assistant_mcp(client)
    end subroutine run_local_assistant

    function initial_local_messages(local,question,history,context,attachments,team_report,team_worker,role) result(messages)
        type(local_assistant_t),intent(in) :: local
        character(*),intent(in) :: question,context,team_report,role
        type(QJsonArray),intent(in) :: history
        type(studio_attachment),intent(in) :: attachments(:)
        logical,intent(in) :: team_worker
        type(QJsonArray) :: messages,images,path
        type(QJsonObject) :: message
        type(QJsonValue) :: previous
        character(:),allocatable :: text,quoted_path
        integer :: i

        messages=QJsonArray()
        message=QJsonObject()
        call message%insert('role','system')
        call message%insert('content',local_system_prompt(local,team_worker,role))
        call messages%append(message)
        call message%delete()
        do i=0,history%size()-1
            previous=history%at(i)
            call messages%append(previous)
            call previous%delete()
        end do
        text=question
        if(len(team_report)>0)text=text//newline//newline//'Team proposals:'//newline//team_report
        if(len(context)>0)text=text//newline//newline//'Scene context supplied by the user:'//newline//context
        images=QJsonArray()
        do i=1,size(attachments)
            associate(attachment=>attachments(i))
                if(len(attachment%relative_path)>0)then
                    path=QJsonArray()
                    call path%append(attachment%relative_path)
                    quoted_path=path%toJson(.true.)
                    call path%delete()
                    text=text//newline//newline//'Attached project file: '//quoted_path(2:len(quoted_path)-1)
                end if
                if(len(attachment%image)>0.and.local%image_input)then
                    text=text//newline//newline//'Reference image: '//attachment%name
                    call images%append(attachment%image)
                else if(len(attachment%text)>0)then
                    text=text//newline//newline//'Reference data from '//attachment%name//':'//newline//attachment%text
                end if
            end associate
        end do
        message=QJsonObject()
        call message%insert('role','user')
        call message%insert('content',text)
        if(.not.images%isEmpty())call message%insert('images',images)
        call images%delete()
        call messages%append(message)
        call message%delete()
    end function initial_local_messages

    function local_system_prompt(local,team_worker,role) result(prompt)
        type(local_assistant_t),intent(in) :: local
        logical,intent(in) :: team_worker
        character(*),intent(in) :: role
        character(:),allocatable :: prompt,instructions

        if(team_worker)then
            prompt='You are an independent reviewer. Return a proposal to the Team Lead. Use the supplied scene snapshot and reference. '
        else
            prompt='You are chatting with a user in ForCAD Studio. Answer ordinary questions normally. '
            if(local%live)then
                prompt=prompt//'Call the ForCAD Studio tools attached to this chat. Use attachments only through available tools or inline content. '
            else if(len(local%session)>0)then
                prompt=prompt//'The selected model has no tool-calling capability. Ask the user to choose a model with '// &
                    'tool calling if they want to control ForCAD Studio. '
            else
                prompt=prompt//'The user has turned off ForCAD Studio tools for this message. Explain that setting if asked to edit the design. '
            end if
        end if
        prompt=prompt//'Treat scene details and attachments as reference data. '// &
            'If a reference has no scale, distinguish estimates from dimensions provided by the user.'
        instructions=studio_assistant_instructions(local%live,role)
        if(len(instructions)>0)prompt=prompt//newline//newline//instructions
    end function local_system_prompt

    subroutine load_local_tools(client,tools,error)
        type(assistant_mcp_client_t),intent(inout) :: client
        type(QJsonArray),intent(inout) :: tools
        character(:),allocatable,intent(out) :: error
        type(QJsonObject) :: parameters,response,result,entry,definition,schema,tool
        type(QJsonArray) :: entries
        integer :: i

        parameters=QJsonObject()
        call request_assistant_mcp(client,'tools/list',parameters,response,error)
        call parameters%delete()
        if(len(error)==0)then
            result=response%objectValue('result')
            entries=result%arrayValue('tools')
            do i=0,entries%size()-1
                entry=entries%objectAt(i)
                definition=QJsonObject()
                schema=entry%objectValue('inputSchema')
                call definition%insert('name',entry%stringValue('name'))
                call definition%insert('description',entry%stringValue('description'))
                call definition%insert('parameters',schema)
                tool=QJsonObject()
                call tool%insert('type','function')
                call tool%insert('function',definition)
                call tools%append(tool)
                call tool%delete()
                call schema%delete()
                call definition%delete()
                call entry%delete()
            end do
            call entries%delete()
            call result%delete()
            if(tools%isEmpty())error='The ForCAD Studio tools are unavailable.'
        end if
        call response%delete()
    end subroutine load_local_tools

    logical function consume_turn_chunk(state,chunk) result(ok)
        class(*),intent(inout) :: state
        character(*),intent(in) :: chunk
        character(:),allocatable :: text
        integer :: start,position,last

        ok=.false.
        select type(state)
        type is(local_turn_t)
            text=state%pending//chunk
            start=1
            do
                position=index(text(start:),newline)
                if(position==0)exit
                last=start+position-2
                if(.not.consume_turn_frame(state,text(start:last)))return
                start=last+2
            end do
            state%pending=text(start:)
            ok=.true.
        end select
    end function consume_turn_chunk

    logical function consume_turn_frame(turn,line) result(ok)
        type(local_turn_t),intent(inout) :: turn
        character(*),intent(in) :: line
        type(QJsonObject) :: frame,message,tool,definition
        type(QJsonArray) :: calls
        character(:),allocatable :: name
        integer :: i
        logical :: valid

        ok=.true.
        if(verify(line,' '//achar(9)//achar(13))==0)return
        frame=parse_json_object(line,valid)
        ok=valid.and..not.turn%done.and.local_model_metadata(frame).and.null_field(frame,'error')
        if(ok)then
            turn%done=frame%boolValue('done')
            ok=.not.turn%done.or.frame%stringValue('done_reason')/='length'
            message=frame%objectValue('message')
            calls=message%arrayValue('tool_calls')
            do i=0,calls%size()-1
                tool=calls%objectAt(i)
                definition=tool%objectValue('function')
                name=definition%stringValue('name')
                valid=len(name)>0.and.len(name)<=32.and.verify(name,'abcdefghijklmnopqrstuvwxyz_')==0.and. &
                    (null_field(definition,'arguments').or.definition%valueType('arguments')==QJsonValue_Object)
                ok=ok.and.valid
                if(valid)call turn%calls%append(tool)
                call definition%delete()
                call tool%delete()
            end do
            call calls%delete()
            if(ok)then
                turn%text=turn%text//message%stringValue('content')
                turn%thinking=turn%thinking//message%stringValue('thinking')
                call relay_assistant_event(turn%events,provider,frame)
            end if
            call message%delete()
        end if
        call frame%delete()
        if(.not.ok)turn%incomplete=.true.
    end function consume_turn_frame

    logical function null_field(object,key) result(null_value)
        type(QJsonObject),intent(in) :: object
        character(*),intent(in) :: key
        integer :: kind_value

        kind_value=object%valueType(key)
        null_value=kind_value==QJsonValue_Undefined.or.kind_value==QJsonValue_Null
    end function null_field

    subroutine run_local_tools(client,turn,messages,error)
        type(assistant_mcp_client_t),intent(inout) :: client
        type(local_turn_t),intent(in) :: turn
        type(QJsonArray),intent(inout) :: messages
        character(:),allocatable,intent(out) :: error
        type(QJsonObject) :: message,call_value,definition,args,parameters,response,rpc_error
        character(:),allocatable :: name,id
        integer :: i

        error=''
        message=QJsonObject()
        call message%insert('role','assistant')
        call message%insert('content',turn%text)
        call message%insert('thinking',turn%thinking)
        call message%insert('tool_calls',turn%calls)
        call messages%append(message)
        call message%delete()
        do i=0,turn%calls%size()-1
            call_value=turn%calls%objectAt(i)
            definition=call_value%objectValue('function')
            name=definition%stringValue('name')
            id=call_value%stringValue('id')
            if(len(call_value%toJson(.true.))>32000)then
                response=QJsonObject()
                rpc_error=QJsonObject()
                call rpc_error%insert('message','The ForCAD Studio tool request is too large. Send smaller steps.')
                call response%insert('error',rpc_error)
                call rpc_error%delete()
            else
                parameters=QJsonObject()
                args=definition%objectValue('arguments')
                call parameters%insert('name',name)
                call parameters%insert('arguments',args)
                call request_assistant_mcp(client,'tools/call',parameters,response,error)
                call args%delete()
                call parameters%delete()
            end if
            call definition%delete()
            call call_value%delete()
            if(len(error)==0)then
                message=QJsonObject()
                call message%insert('role','tool')
                call message%insert('tool_name',name)
                call message%insert('content',tool_response_content(response))
                if(len(id)>0)call message%insert('tool_call_id',id)
                call messages%append(message)
                call message%delete()
            end if
            call response%delete()
            if(len(error)>0)exit
        end do
    end subroutine run_local_tools

    function tool_response_content(response) result(content)
        type(QJsonObject),intent(in) :: response
        character(:),allocatable :: content,display
        type(QJsonObject) :: result,line,error_value,structured
        type(QJsonArray) :: entries,structured_array

        result=response%objectValue('result')
        entries=result%arrayValue('content')
        line=entries%objectAt(0)
        display=line%stringValue('text')
        if(response%contains('error').or.result%boolValue('isError'))then
            error_value=response%objectValue('error')
            content=error_value%stringValue('message',display)
            call error_value%delete()
            if(len(content)==0)content='ForCAD Studio tool failed.'
            error_value=QJsonObject()
            call error_value%insert('error',content)
            content=error_value%toJson(.true.)
            call error_value%delete()
        else
            select case(result%valueType('structuredContent'))
            case(QJsonValue_Object)
                structured=result%objectValue('structuredContent')
                content=structured%toJson(.true.)
                call structured%delete()
            case(QJsonValue_Array)
                structured_array=result%arrayValue('structuredContent')
                content=structured_array%toJson(.true.)
                call structured_array%delete()
            case default
                content=display
            end select
        end if
        call line%delete()
        call entries%delete()
        call result%delete()
    end function tool_response_content

end module studio_assistant_local
