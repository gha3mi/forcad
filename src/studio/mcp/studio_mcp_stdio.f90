module studio_mcp_stdio
    use, intrinsic :: iso_c_binding, only: c_char,c_int,c_ptr,c_loc
    use, intrinsic :: iso_fortran_env, only: output_unit,int64
    use qt, only: QJsonArray,QJsonDocument,QJsonDocument_fromJson,QJsonObject,QJsonValue, &
        QJsonValue_Object,QJsonValue_String,QJsonValue_Double,QJsonValue_Bool,QThread_msleep
    use studio_mcp_schema, only: studio_mcp_catalog,studio_mcp_tool_schema,studio_mcp_valid_arguments,studio_rpc_token
    use studio_mcp_mailbox, only: mcp_mailbox,submit_mcp_request,read_mcp_response,close_mcp_mailbox,mcp_mailbox_exists, &
        mcp_request_waiting,mcp_request_claimed,retract_mcp_request,mcp_response_pending,mcp_response_valid,max_mcp_request_bytes
    use studio_presentation_text, only: studio_version
    implicit none(type, external)
    private
    public :: run_studio_mcp_stdio

    character, parameter :: newline=new_line('a')
    integer, parameter :: max_approval_backlog=2*max_mcp_request_bytes
    logical :: output_closed=.false.

    type :: rpc_input
        character(:), allocatable :: pending
        logical :: ended=.false.
    end type rpc_input

    interface
        subroutine open_stdio() bind(C,name='studio_ipc_begin')
            implicit none(type,external)
        end subroutine open_stdio

        subroutine close_stdio() bind(C,name='studio_ipc_end')
            implicit none(type,external)
        end subroutine close_stdio

        integer(c_int) function interrupted() bind(C,name='studio_ipc_interrupted')
            import c_int
            implicit none(type,external)
        end function interrupted

        integer(c_int) function read_input(buffer,capacity,milliseconds) bind(C,name='studio_ipc_read')
            import c_int,c_ptr
            implicit none(type,external)
            type(c_ptr), value :: buffer
            integer(c_int), value :: capacity,milliseconds
        end function read_input
    end interface

contains

    integer function run_studio_mcp_stdio(session) result(status)
        character(*), intent(in) :: session
        type(rpc_input) :: input
        type(QJsonArray) :: catalog
        character(:), allocatable :: frame
        logical :: initialized,oversized

        call open_stdio()
        output_closed=.false.
        input%pending=''
        initialized=.false.
        catalog=studio_mcp_catalog()
        do while(interrupted()==0.and..not.output_closed)
            if(.not.next_frame(input,frame,oversized))exit
            if(oversized)then
                call send_error(-32700,'The MCP message exceeds 4 MiB. Send a smaller request.')
            else
                call handle_frame(frame,session,catalog,input,initialized)
            end if
        end do
        call catalog%delete()
        status=interrupted()
        if(status/=0)status=128+status
        call close_stdio()
    end function run_studio_mcp_stdio

    subroutine collect_input(input,milliseconds)
        type(rpc_input), intent(inout) :: input
        integer(c_int), intent(in) :: milliseconds
        character(kind=c_char,len=65536), target :: bytes
        integer(c_int) :: count

        if(input%ended.or.interrupted()/=0)return
        count=read_input(c_loc(bytes(1:1)),len(bytes,kind=c_int),milliseconds)
        if(count>0)then
            input%pending=input%pending//bytes(:count)
        else if(count<0)then
            input%ended=.true.
        end if
    end subroutine collect_input

    logical function next_frame(input,frame,oversized) result(available)
        type(rpc_input), intent(inout) :: input
        character(:), allocatable, intent(out) :: frame
        logical, intent(out) :: oversized
        integer :: line_end

        available=.false.
        oversized=.false.
        do while(interrupted()==0)
            line_end=index(input%pending,newline)
            if(line_end>0)then
                oversized=oversized.or.line_end-1>max_mcp_request_bytes
                if(.not.oversized)frame=input%pending(:line_end-1)
                input%pending=input%pending(line_end+1:)
                available=.true.
                return
            end if
            if(len(input%pending)>max_mcp_request_bytes)oversized=.true.
            if(oversized)input%pending=''
            if(input%ended)then
                available=oversized.or.len(input%pending)>0
                call move_alloc(input%pending,frame)
                input%pending=''
                return
            end if
            call collect_input(input,250_c_int)
        end do
    end function next_frame

    subroutine parse_frame(payload,object,object_value,valid)
        character(*), intent(in) :: payload
        type(QJsonObject), intent(out) :: object
        logical, intent(out) :: object_value,valid
        type(QJsonDocument) :: document
        type(QJsonArray) :: values

        object_value=.false.
        valid=.false.
        document=QJsonDocument_fromJson('['//payload//']')
        if(document%isArray())then
            values=document%array()
            valid=values%size()==1
            if(valid)then
                object_value=values%typeAt(0)==QJsonValue_Object
                if(object_value)object=values%objectAt(0)
            end if
            call values%delete()
        end if
        call document%delete()
    end subroutine parse_frame

    subroutine handle_frame(frame,session,catalog,input,initialized)
        character(*), intent(in) :: frame,session
        type(QJsonArray), intent(in) :: catalog
        type(rpc_input), intent(inout) :: input
        logical, intent(inout) :: initialized
        type(QJsonObject) :: message,params,response
        type(QJsonValue) :: id
        character(:), allocatable :: method
        logical :: object_value,valid

        call parse_frame(frame,message,object_value,valid)
        if(.not.valid)then
            call send_error(-32700,'Invalid JSON')
            return
        end if
        if(.not.object_value)then
            call send_error(-32600,'Invalid JSON-RPC request')
            return
        end if
        method=message%stringValue('method')
        valid=message%stringValue('jsonrpc')=='2.0'.and.studio_rpc_token(method,.false.)
        if(message%contains('params'))valid=valid.and.message%valueType('params')==QJsonValue_Object
        if(message%contains('id'))valid=valid.and.valid_id(message%valueType('id'))
        if(.not.valid)then
            call send_error(-32600,'Invalid JSON-RPC request')
            call message%delete()
            return
        end if
        if(.not.message%contains('id'))then
            call message%delete()
            return
        end if

        id=message%value('id')
        params=message%objectValue('params')
        select case(method)
        case('initialize')
            call initialize_client(id,params,initialized)
        case('ping')
            response=QJsonObject()
            call send_result(id,response)
            call response%delete()
        case('tools/list','tools/call')
            if(.not.initialized)then
                call send_error(-32002,'Initialize MCP before using ForCAD Studio tools.',id)
            else if(method=='tools/list')then
                response=QJsonObject()
                call response%insert('tools',catalog)
                call send_result(id,response)
                call response%delete()
            else
                call handle_tool(id,params,session,catalog,input)
            end if
        case default
            call send_error(-32601,'Unknown MCP method.',id)
        end select
        call id%delete()
        call params%delete()
        call message%delete()
    end subroutine handle_frame

    logical function valid_id(kind) result(valid)
        integer, intent(in) :: kind

        valid=kind==QJsonValue_Double.or.kind==QJsonValue_String
    end function valid_id

    subroutine initialize_client(id,params,initialized)
        type(QJsonValue), intent(in) :: id
        type(QJsonObject), intent(in) :: params
        logical, intent(inout) :: initialized
        type(QJsonObject) :: response,capabilities,empty,server
        character(:), allocatable :: version

        if(initialized)then
            call send_error(-32600,'MCP is already initialized.',id)
            return
        end if
        version=params%stringValue('protocolVersion')
        if(len(version)==0.or.params%valueType('capabilities')/=QJsonValue_Object.or. &
            params%valueType('clientInfo')/=QJsonValue_Object)then
            call send_error(-32602,'Invalid initialize parameters.',id)
            return
        end if
        select case(version)
        case('2024-11-05','2025-03-26','2025-06-18','2025-11-25')
        case default
            version='2025-11-25'
        end select
        initialized=.true.
        empty=QJsonObject()
        capabilities=QJsonObject()
        call capabilities%insert('tools',empty)
        call empty%delete()
        server=QJsonObject()
        call server%insert('name','forcad-studio')
        call server%insert('version',studio_version)
        response=QJsonObject()
        call response%insert('protocolVersion',version)
        call response%insert('capabilities',capabilities)
        call response%insert('serverInfo',server)
        call response%insert('instructions', &
            'Use inspect for live design IDs, inspect_view for workspace state, and list_commands for available GUI actions.')
        call send_result(id,response)
        call capabilities%delete()
        call server%delete()
        call response%delete()
    end subroutine initialize_client

    subroutine handle_tool(id,params,session,catalog,input)
        type(QJsonValue), intent(in) :: id
        type(QJsonObject), intent(in) :: params
        character(*), intent(in) :: session
        type(QJsonArray), intent(in) :: catalog
        type(rpc_input), intent(inout) :: input
        type(QJsonObject) :: schema,empty,request
        type(QJsonValue) :: arguments
        character(:), allocatable :: name,payload

        name=params%stringValue('name')
        if(.not.studio_rpc_token(name,.true.))name=''
        schema=studio_mcp_tool_schema(catalog,name)
        if(schema%isEmpty())then
            call send_error(-32602,'Unknown ForCAD Studio MCP tool.',id)
            call schema%delete()
            return
        end if
        if(params%contains('arguments'))then
            arguments=params%value('arguments')
        else
            empty=QJsonObject()
            arguments=QJsonValue(empty)
            call empty%delete()
        end if
        if(studio_mcp_valid_arguments(arguments,schema))then
            request=QJsonObject()
            call request%insert('tool',name)
            call request%insert('arguments',arguments)
            payload=request%toJson(.true.)
            call request%delete()
            call relay_tool(id,session,payload,input)
        else
            call send_tool_error(id,'Arguments do not match the input schema for '//name//'.')
        end if
        call schema%delete()
        call arguments%delete()
    end subroutine handle_tool

    subroutine relay_tool(id,session,request,input)
        type(QJsonValue), intent(in) :: id
        character(*), intent(in) :: session,request
        type(rpc_input), intent(inout) :: input
        type(mcp_mailbox) :: mailbox
        character(:), allocatable :: problem,payload
        integer(int64) :: now,clock_rate,deadline
        integer :: response_status,scan_start
        logical :: awaiting_approval,claimed

        problem=submit_mcp_request(mailbox,session,request)
        if(len(problem)>0)then
            call send_tool_error(id,problem)
            return
        end if
        awaiting_approval=.false.
        scan_start=1
        call system_clock(now,clock_rate)
        deadline=now+30*clock_rate
        do while(now<deadline.and.interrupted()==0)
            if(.not.mcp_mailbox_exists(mailbox))then
                call close_mcp_mailbox(mailbox)
                call send_tool_error(id,'The ForCAD Studio session closed before it returned a result. Check the design before retrying.')
                return
            end if
            if(.not.awaiting_approval)then
                awaiting_approval=mcp_request_waiting(mailbox)
                if(awaiting_approval)deadline=now+135*clock_rate
            end if
            if(awaiting_approval)then
                if(cancelled_by_client(id,input,scan_start))then
                    call close_mcp_mailbox(mailbox)
                    return
                end if
                if(len(input%pending)>max_approval_backlog)then
                    call close_mcp_mailbox(mailbox)
                    call send_tool_error(id,'The client queued too much data while waiting for ForCAD Studio approval.')
                    output_closed=.true.
                    return
                end if
            end if
            payload=read_mcp_response(mailbox,response_status)
            if(response_status/=mcp_response_pending)then
                call close_mcp_mailbox(mailbox)
                if(response_status==mcp_response_valid)then
                    call forward_response(id,payload)
                else
                    call send_invalid_response(id)
                end if
                return
            end if
            call QThread_msleep(50)
            call system_clock(now)
        end do
        call retract_mcp_request(mailbox)
        claimed=mcp_request_claimed(mailbox)
        call close_mcp_mailbox(mailbox)
        if(interrupted()/=0)return
        if(claimed)then
            call send_tool_error(id,'ForCAD Studio did not reply in time. The command might have completed. Inspect the design before retrying.')
        else if(awaiting_approval)then
            call send_tool_error(id,'ForCAD Studio did not receive approval in time. The action was not applied.')
        else
            call send_tool_error(id,'ForCAD Studio did not claim the request in time. It may be busy or the MCP session may be closed.')
        end if
    end subroutine relay_tool

    logical function cancelled_by_client(id,input,scan_start) result(cancelled)
        type(QJsonValue), intent(in) :: id
        type(rpc_input), intent(inout) :: input
        integer, intent(inout) :: scan_start
        integer :: line_end,relative_end,length

        cancelled=.false.
        call collect_input(input,0_c_int)
        length=len(input%pending)
        do while(scan_start<=length)
            relative_end=index(input%pending(scan_start:),newline)
            if(relative_end==0.and..not.input%ended)return
            line_end=length+1
            if(relative_end>0)line_end=scan_start+relative_end-1
            if(line_end-scan_start<=max_mcp_request_bytes)then
                cancelled=is_cancellation(input%pending(scan_start:line_end-1),id)
                if(cancelled)then
                    input%pending=input%pending(:scan_start-1)//input%pending(min(length+1,line_end+1):)
                    return
                end if
            end if
            scan_start=line_end+1
        end do
    end function cancelled_by_client

    logical function is_cancellation(payload,id) result(cancelled)
        character(*), intent(in) :: payload
        type(QJsonValue), intent(in) :: id
        type(QJsonObject) :: frame,params
        logical :: object_value,valid

        cancelled=.false.
        call parse_frame(payload,frame,object_value,valid)
        if(.not.object_value)return
        if(frame%stringValue('jsonrpc')=='2.0'.and.frame%stringValue('method')=='notifications/cancelled')then
            params=frame%objectValue('params')
            if(id%isDouble().and.params%valueType('requestId')==QJsonValue_Double)then
                cancelled=id%toDouble()==params%doubleValue('requestId')
            else if(id%isString().and.params%valueType('requestId')==QJsonValue_String)then
                cancelled=id%toString()==params%stringValue('requestId')
            end if
            call params%delete()
        end if
        call frame%delete()
    end function is_cancellation

    subroutine forward_response(id,payload)
        type(QJsonValue), intent(in) :: id
        character(*), intent(in) :: payload
        type(QJsonObject) :: response,data,result
        type(QJsonArray) :: contents
        logical :: object_value,valid

        call parse_frame(payload,response,object_value,valid)
        if(.not.object_value)then
            call send_invalid_response(id)
            return
        end if
        if(response%valueType('ok')/=QJsonValue_Bool)then
            call send_invalid_response(id)
        else if(response%boolValue('ok'))then
            if(.not.response%contains('result'))then
                data=response%copy()
                call data%remove('ok')
            else if(response%valueType('result')==QJsonValue_Object)then
                data=response%objectValue('result')
            else
                call response%delete()
                call send_invalid_response(id)
                return
            end if
            result=QJsonObject()
            contents=text_content(data%toJson(.true.))
            call result%insert('structuredContent',data)
            call result%insert('content',contents)
            call send_result(id,result)
            call data%delete()
            call result%delete()
            call contents%delete()
        else if(response%valueType('error')==QJsonValue_String)then
            call send_tool_error(id,response%stringValue('error'))
        else
            call send_invalid_response(id)
        end if
        call response%delete()
    end subroutine forward_response

    subroutine send_invalid_response(id)
        type(QJsonValue), intent(in) :: id

        call send_tool_error(id,'ForCAD Studio returned an unreadable MCP tool result. Check the design before retrying.')
    end subroutine send_invalid_response

    function text_content(text) result(contents)
        character(*), intent(in) :: text
        type(QJsonArray) :: contents
        type(QJsonObject) :: item

        item=QJsonObject()
        call item%insert('type','text')
        call item%insert('text',text)
        contents=QJsonArray()
        call contents%append(item)
        call item%delete()
    end function text_content

    subroutine send_tool_error(id,message)
        type(QJsonValue), intent(in) :: id
        character(*), intent(in) :: message
        type(QJsonObject) :: response
        type(QJsonArray) :: contents

        contents=text_content(message)
        response=QJsonObject()
        call response%insert('content',contents)
        call response%insert('isError',.true.)
        call send_result(id,response)
        call response%delete()
        call contents%delete()
    end subroutine send_tool_error

    subroutine send_error(code,message,id)
        integer, intent(in) :: code
        character(*), intent(in) :: message
        type(QJsonValue), intent(in), optional :: id
        type(QJsonObject) :: envelope,error

        envelope=QJsonObject()
        call envelope%insert('jsonrpc','2.0')
        if(present(id))then
            call envelope%insert('id',id)
        else
            call envelope%insertNull('id')
        end if
        error=QJsonObject()
        call error%insert('code',code)
        call error%insert('message',message)
        call envelope%insert('error',error)
        call send_frame(envelope)
        call error%delete()
        call envelope%delete()
    end subroutine send_error

    subroutine send_result(id,result)
        type(QJsonValue), intent(in) :: id
        type(QJsonObject), intent(in) :: result
        type(QJsonObject) :: envelope

        envelope=QJsonObject()
        call envelope%insert('jsonrpc','2.0')
        call envelope%insert('id',id)
        call envelope%insert('result',result)
        call send_frame(envelope)
        call envelope%delete()
    end subroutine send_result

    subroutine send_frame(envelope)
        type(QJsonObject), intent(in) :: envelope
        integer :: io

        if(interrupted()/=0.or.output_closed)return
        write(output_unit,'(a)',iostat=io)envelope%toJson(.true.)
        if(io==0)flush(output_unit,iostat=io)
        output_closed=io/=0
    end subroutine send_frame

end module studio_mcp_stdio
