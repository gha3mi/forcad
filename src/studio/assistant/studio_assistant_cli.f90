module studio_assistant_cli
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QProcess,QProcess_NotRunning,QProcess_NormalExit,QProcessEnvironment,QProcessEnvironment_systemEnvironment, &
        QJsonDocument,QJsonDocument_fromJson,QJsonObject,QJsonArray,QJsonValue,QJsonValue_Array,QJsonValue_String, &
        QSaveFile,QThread_isInterruptionRequested,QUrl_fromLocalFile
    use studio_assistant_events, only: assistant_events_t,relay_assistant_event,clip_assistant_text,clean_assistant_text, &
        assistant_failure_hint
    use studio_assistant_process, only: start_provider_process
    use studio_assistant_transport, only: read_text,write_text,valid_provider_session
    use studio_chat_attachment, only: studio_attachment
    use studio_platform_files, only: file_metadata,inspect_path,is_regular_file
    use studio_presentation_text, only: integer_text
    implicit none(type,external)
    private

    public :: run_assistant_cli

    character,parameter :: lf=new_line('a')
    character(*),parameter :: base64_letters='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'
    type :: process_argument
        character(:),allocatable :: text
    end type process_argument

contains

    subroutine run_assistant_cli(provider,model,effort,workspace,session,executable,prompt,attachments, &
            reference_only,worker,events,error,conversation_id)
        character(*),intent(in) :: provider,model,effort,workspace,session,executable,prompt
        type(studio_attachment),intent(in) :: attachments(:)
        logical,intent(in) :: reference_only,worker
        type(assistant_events_t),intent(inout) :: events
        character(:),allocatable,intent(out) :: error
        character(*),optional,intent(in) :: conversation_id
        type(QProcess) :: process
        type(QProcessEnvironment) :: environment
        type(process_argument),allocatable :: args(:)
        character(:),allocatable :: cli,input,saved_session
        character(24) :: exit_code
        logical :: started,completed

        error=''
        if(len_trim(prompt)==0)then
            error='The agent prompt is empty.'
            return
        end if
        if(worker.and.len(session)>0)then
            error='A planning agent cannot use a live ForCAD Studio connection.'
            return
        end if
        select case(provider)
        case('Codex CLI'); cli='codex'
        case('Claude Code CLI'); cli='claude'
        case('OpenCode CLI'); cli='opencode'
        case default
            error='Unknown Assistant provider.'
            return
        end select
        if(.not.allocated(events%exchange))then
            error='The Assistant exchange is unavailable.'
            return
        end if

        saved_session=''
        if(present(conversation_id).and..not.worker)then
            saved_session=conversation_id
            if(len(saved_session)>0.and..not.valid_provider_session(saved_session))then
                error='The saved provider session ID is invalid. Start a new conversation to continue.'
                return
            end if
        end if
        input=prompt
        allocate(args(0))
        select case(cli)
        case('codex')
            call codex_arguments(args,model,effort,executable,session,events%exchange,attachments,reference_only,worker,saved_session,error)
        case('claude')
            call claude_arguments(args,model,effort,executable,session,events%exchange,prompt,attachments,reference_only, &
                worker,input,error,saved_session)
        case('opencode')
            call add_argument(args,'acp')
        case default
            return
        end select
        if(len(error)>0)return

        process=QProcess()
        call process%setWorkingDirectory(workspace)
        if(worker)then
            environment=QProcessEnvironment_systemEnvironment()
            call environment%remove('FORCAD_MCP_SESSION')
            if(cli=='opencode')call restrict_opencode_worker(environment)
            call process%setProcessEnvironment(environment)
            call environment%delete()
        end if
        call start_process(process,cli,args)
        started=process%waitForStarted(5000)
        if(.not.started)then
            error=provider//' is unavailable.'
            call process%delete()
            return
        end if

        if(cli=='opencode')then
            call open_opencode_session(process,model,effort,workspace,session,executable,prompt,attachments, &
                reference_only,worker,events,error,completed,saved_session)
        else
            completed=process%write(input)==len(input)
            call process%closeWriteChannel()
            if(completed)then
                call stream_provider(process,provider,events,completed,worker,saved_session,error)
            else if(len(saved_session)>0)then
                error=resume_failure(provider,'')
            end if
        end if
        call capture_failure_hint(process,events)
        if(.not.completed.or.events%failed.or.provider_stopping(events))then
            if(len(error)==0.and.allocated(events%failure_hint))error=provider//' · '//events%failure_hint
            if(len(error)==0)then
                if(process%state()==QProcess_NotRunning.and.process%exitStatus()==QProcess_NormalExit)then
                    write(exit_code,'(i0)')process%exitCode()
                    error=provider//' ended before completing the turn (exit '//trim(exit_code)//').'
                else
                    error=provider//' disconnected before completing the turn.'
                end if
            end if
        end if
        call close_provider(process)
    end subroutine run_assistant_cli

    subroutine add_argument(args,value)
        type(process_argument),allocatable,intent(inout) :: args(:)
        character(*),intent(in) :: value
        type(process_argument) :: argument

        argument%text=value
        args=[args,argument]
    end subroutine add_argument

    subroutine add_arguments(args,values)
        type(process_argument),allocatable,intent(inout) :: args(:)
        character(*),intent(in) :: values(:)
        integer :: i

        do i=1,size(values)
            call add_argument(args,trim(values(i)))
        end do
    end subroutine add_arguments

    subroutine add_option(args,option,value)
        type(process_argument),allocatable,intent(inout) :: args(:)
        character(*),intent(in) :: option,value

        call add_argument(args,option)
        call add_argument(args,value)
    end subroutine add_option

    subroutine start_process(process,program,args)
        type(QProcess),intent(inout) :: process
        character(*),intent(in) :: program
        type(process_argument),intent(in) :: args(:)
        character(:),allocatable :: parameters(:)
        integer :: i,width

        width=1
        do i=1,size(args)
            width=max(width,len(args(i)%text))
        end do
        allocate(character(width)::parameters(size(args)))
        do i=1,size(args)
            parameters(i)=args(i)%text
        end do
        call start_provider_process(process,program,parameters)
    end subroutine start_process

    subroutine codex_arguments(args,model,effort,executable,session,directory,attachments,reference_only,worker,saved_session,error)
        type(process_argument),allocatable,intent(inout) :: args(:)
        character(*),intent(in) :: model,effort,executable,session,directory,saved_session
        type(studio_attachment),intent(in) :: attachments(:)
        logical,intent(in) :: reference_only,worker
        character(:),allocatable,intent(inout) :: error
        type(QJsonArray) :: bridge

        call add_argument(args,'exec')
        if(len(saved_session)>0)then
            call add_argument(args,'resume')
        else
            call add_option(args,'--color','never')
        end if
        call add_arguments(args,[character(24)::'--skip-git-repo-check','--json'])
        if(worker)call add_arguments(args,[character(48)::'--sandbox','read-only','-c','approval_policy="never"', &
            '-c','mcp_servers={}'])
        if(len(model)>0)call add_option(args,'--model',model)
        if(len(effort)>0)call add_option(args,'-c','model_reasoning_effort='//json_string(effort))
        if(.not.reference_only)call add_codex_images(args,directory,attachments,error)
        if(len(error)>0)return
        if(len(session)>0)then
            call add_option(args,'-c','mcp_servers.forcad_studio.command='//json_string(executable))
            bridge=mcp_arguments(session)
            call add_option(args,'-c','mcp_servers.forcad_studio.args='//bridge%toJson(.true.))
            call bridge%delete()
            call add_arguments(args,[character(48)::'-c','mcp_servers.forcad_studio.tool_timeout_sec=160', &
                '-c','mcp_servers.forcad_studio.enabled=true','-c','mcp_servers.forcad_studio.required=true'])
        end if
        call add_argument(args,'--')
        if(len(saved_session)>0)call add_argument(args,saved_session)
        call add_argument(args,'-')
    end subroutine codex_arguments

    subroutine add_codex_images(args,directory,attachments,error)
        type(process_argument),allocatable,intent(inout) :: args(:)
        character(*),intent(in) :: directory
        type(studio_attachment),intent(in) :: attachments(:)
        character(:),allocatable,intent(inout) :: error
        character(:),allocatable :: media,bytes,path
        integer :: i

        do i=1,size(attachments)
            media=image_media_type(attachments(i)%image)
            if(len(media)==0)cycle
            bytes=decode_base64(attachments(i)%image)
            if(len(bytes)==0)then
                error='ForCAD Studio could not decode the image reference '//attachments(i)%name//'.'
                return
            end if
            path=directory//'/reference-'//integer_text(i)//merge('.png','.jpg',media=='image/png')
            if(.not.write_text(path,bytes))then
                error='ForCAD Studio could not prepare the image reference '//attachments(i)%name//'.'
                return
            end if
            call add_option(args,'--image',path)
        end do
    end subroutine add_codex_images

    subroutine claude_arguments(args,model,effort,executable,session,directory,prompt,attachments,reference_only, &
            worker,input,error,saved_session)
        type(process_argument),allocatable,intent(inout) :: args(:)
        character(*),intent(in) :: model,effort,executable,session,directory,prompt,saved_session
        type(studio_attachment),intent(in) :: attachments(:)
        logical,intent(in) :: reference_only,worker
        character(:),allocatable,intent(inout) :: input,error
        type(QJsonObject) :: config,server,servers,message,frame,source,content_block
        type(QJsonArray) :: bridge,content
        character(:),allocatable :: config_path,media
        integer :: i

        call add_arguments(args,[character(28)::'--print','--output-format','stream-json','--verbose','--include-partial-messages'])
        if(len(saved_session)>0)call add_option(args,'--resume',saved_session)
        if(worker)call add_arguments(args,[character(20)::'--tools','','--disallowedTools','mcp__*','--strict-mcp-config'])
        if(len(model)>0)call add_option(args,'--model',model)
        if(len(effort)>0)call add_option(args,'--effort',effort)
        if(worker.or.len(session)>0)then
            config=QJsonObject()
            servers=QJsonObject()
            if(len(session)>0)then
                server=QJsonObject()
                call server%insert('type','stdio')
                call server%insert('command',executable)
                bridge=mcp_arguments(session)
                call server%insert('args',bridge)
                call bridge%delete()
                call servers%insert('forcad_studio',server)
                call server%delete()
            end if
            call config%insert('mcpServers',servers)
            config_path=directory//'/provider-mcp.json'
            if(.not.write_private_file(config_path,config%toJson(.true.)//lf))error='ForCAD Studio could not prepare the provider tool configuration.'
            call servers%delete()
            call config%delete()
            if(len(error)>0)return
            call add_option(args,'--mcp-config',config_path)
            if(len(session)>0)call add_option(args,'--allowedTools','mcp__forcad_studio__*')
        end if
        if(reference_only)return

        content=text_prompt(prompt)
        do i=1,size(attachments)
            media=image_media_type(attachments(i)%image)
            if(len(media)==0)cycle
            source=QJsonObject()
            call source%insert('type','base64')
            call source%insert('media_type',media)
            call source%insert('data',attachments(i)%image)
            content_block=QJsonObject()
            call content_block%insert('type','image')
            call content_block%insert('source',source)
            call source%delete()
            call content%append(content_block)
            call content_block%delete()
        end do
        if(content%size()==1)then
            call content%delete()
            return
        end if
        call add_option(args,'--input-format','stream-json')
        message=QJsonObject()
        call message%insert('role','user')
        call message%insert('content',content)
        call content%delete()
        frame=QJsonObject()
        call frame%insert('type','user')
        call frame%insertNull('parent_tool_use_id')
        call frame%insert('session_id',saved_session)
        call frame%insert('message',message)
        call message%delete()
        input=frame%toJson(.true.)//lf
        call frame%delete()
    end subroutine claude_arguments

    subroutine restrict_opencode_worker(environment)
        type(QProcessEnvironment),intent(inout) :: environment
        character(20),parameter :: model_fields(*)=[character(20)::'provider','providers','model','small_model','enabled_providers', &
            'disabled_providers']
        character(*),parameter :: restrictions='{"permission":{"*":"deny"},"policies":['// &
            '{"action":"tool.use","resource":"*","effect":"deny"},'// &
            '{"action":"mcp.use","resource":"*","effect":"deny"},'// &
            '{"action":"plugin.use","resource":"*","effect":"deny"},'// &
            '{"action":"integration.use","resource":"mcp:*","effect":"deny"},'// &
            '{"action":"integration.use","resource":"plugin:*","effect":"deny"}]}'
        type(QJsonDocument) :: document
        type(QJsonObject) :: config,original
        type(QJsonValue) :: value
        character(:),allocatable :: key
        integer :: i

        document=QJsonDocument_fromJson(restrictions)
        config=document%object()
        call document%delete()
        if(parse_object(environment%value('OPENCODE_CONFIG_CONTENT'),original))then
            do i=1,size(model_fields)
                key=trim(model_fields(i))
                if(.not.original%contains(key))cycle
                value=original%value(key)
                call config%insert(key,value)
                call value%delete()
            end do
            call original%delete()
        end if
        call environment%insert('OPENCODE_DISABLE_PROJECT_CONFIG','true')
        call environment%insert('OPENCODE_CONFIG_CONTENT',config%toJson(.true.))
        call config%delete()
    end subroutine restrict_opencode_worker

    function mcp_arguments(session) result(args)
        character(*),intent(in) :: session
        type(QJsonArray) :: args

        args=QJsonArray()
        call args%append('--mcp-stdio')
        call args%append('--session')
        call args%append(session)
    end function mcp_arguments

    function json_string(value) result(quoted)
        character(*),intent(in) :: value
        character(:),allocatable :: quoted
        type(QJsonArray) :: array

        array=QJsonArray()
        call array%append(value)
        quoted=array%toJson(.true.)
        call array%delete()
        quoted=quoted(2:len(quoted)-1)
    end function json_string

    function text_prompt(prompt) result(content)
        character(*),intent(in) :: prompt
        type(QJsonArray) :: content
        type(QJsonObject) :: block

        block=QJsonObject()
        call block%insert('type','text')
        call block%insert('text',prompt)
        content=QJsonArray()
        call content%append(block)
        call block%delete()
    end function text_prompt

    subroutine open_opencode_session(process,model,effort,workspace,session,executable,prompt,attachments, &
            reference_only,worker,events,error,complete,saved_session)
        type(QProcess),intent(inout) :: process
        character(*),intent(in) :: model,effort,workspace,session,executable,prompt,saved_session
        type(studio_attachment),intent(in) :: attachments(:)
        logical,intent(in) :: reference_only,worker
        type(assistant_events_t),intent(inout) :: events
        character(:),allocatable,intent(inout) :: error
        logical,intent(out) :: complete
        type(QJsonObject) :: params,empty,client,reply,server,capabilities
        type(QJsonArray) :: servers,args,environment,content
        character(:),allocatable :: id,chosen,requested_effort,effort_id,stop
        integer :: request_id,permission_index
        logical :: supports_load

        complete=.false.
        permission_index=0
        params=QJsonObject()
        empty=QJsonObject()
        client=QJsonObject()
        call client%insert('name','ForCAD Studio')
        call client%insert('version','0.1.0')
        call params%insert('protocolVersion',1)
        call params%insert('clientCapabilities',empty)
        call params%insert('clientInfo',client)
        call empty%delete()
        call client%delete()
        call acp_call(process,1,'initialize',params,25,worker,events,permission_index,reply,error)
        call params%delete()
        capabilities=reply%objectValue('agentCapabilities')
        supports_load=capabilities%boolValue('loadSession')
        call capabilities%delete()
        call reply%delete()
        if(len(error)>0)return
        if(len(saved_session)>0.and..not.supports_load)then
            error='This OpenCode release cannot resume saved sessions over ACP. Update OpenCode or start a new conversation.'
            return
        end if

        params=QJsonObject()
        call params%insert('cwd',workspace)
        servers=QJsonArray()
        if(len(session)>0)then
            server=QJsonObject()
            call server%insert('name','forcad_studio')
            call server%insert('command',executable)
            args=mcp_arguments(session)
            environment=QJsonArray()
            call server%insert('args',args)
            call server%insert('env',environment)
            call args%delete()
            call environment%delete()
            call servers%append(server)
            call server%delete()
        end if
        call params%insert('mcpServers',servers)
        call servers%delete()
        if(len(saved_session)>0)then
            call params%insert('sessionId',saved_session)
            call acp_call(process,2,'session/load',params,45,worker,events,permission_index,reply,error)
            id=reply%stringValue('sessionId')
            if(len(id)==0)id=saved_session
            if(len(error)>0.and..not.provider_stopping(events))then
                error=resume_failure('OpenCode',error)
            end if
        else
            call acp_call(process,2,'session/new',params,45,worker,events,permission_index,reply,error)
            id=reply%stringValue('sessionId')
        end if
        call params%delete()
        if(len(error)==0)call capture_provider_session('OpenCode',events,worker,id,saved_session,error)
        if(len(error)>0.or.provider_stopping(events))then
            call reply%delete()
            return
        end if

        request_id=3
        chosen=model
        if(len(chosen)==0)chosen=environment_value('FORCAD_OPENCODE_MODEL')
        if(len(chosen)>0)then
            call reply%delete()
            call acp_configure(process,request_id,id,'model',chosen,worker,events,permission_index,reply,error)
            request_id=request_id+1
            if(len(error)>0)then
                call reply%delete()
                return
            end if
        end if
        requested_effort=effort
        if(len(effort)==0.and.len(saved_session)>0)requested_effort='default'
        effort_id=''
        if(len(requested_effort)>0)effort_id=opencode_effort(reply,requested_effort)
        if(len(effort_id)==0.and.len(effort)>0)then
            call reply%delete()
            error='The selected OpenCode model does not offer that effort. Refresh the model list and choose an available level.'
            return
        end if
        if(len(effort_id)>0)then
            call reply%delete()
            call acp_configure(process,request_id,id,effort_id,requested_effort,worker,events,permission_index,reply,error)
            request_id=request_id+1
            if(len(error)>0)then
                call reply%delete()
                return
            end if
        end if
        call reply%delete()
        if(provider_stopping(events))return

        content=text_prompt(prompt)
        if(.not.reference_only)call add_opencode_attachments(content,attachments)
        params=QJsonObject()
        call params%insert('sessionId',id)
        call params%insert('prompt',content)
        call content%delete()
        call acp_call(process,request_id,'session/prompt',params,900,worker,events,permission_index,reply,error)
        call params%delete()
        stop=reply%stringValue('stopReason')
        call reply%delete()
        if(provider_stopping(events))call cancel_opencode_prompt(process,id)
        if(len(error)>0)return

        select case(stop)
        case('max_tokens')
            error='OpenCode reached the model output limit before completing the request.'
        case('max_turn_requests','refusal','cancelled')
            error='OpenCode ended the turn without completing the request.'
        case('end_turn')
            call final_opencode_updates(process,events)
            complete=.true.
        case default
            error='OpenCode did not report a completed turn.'
        end select
    end subroutine open_opencode_session

    subroutine add_opencode_attachments(content,attachments)
        type(QJsonArray),intent(inout) :: content
        type(studio_attachment),intent(in) :: attachments(:)
        type(QJsonObject) :: block
        character(:),allocatable :: uri,media
        integer :: i

        do i=1,size(attachments)
            uri=''
            if(len(attachments(i)%path)>0)then
                if(is_inline_reference(attachments(i)%path))uri=file_uri(attachments(i)%path)
            end if
            if(len(uri)>0)then
                block=QJsonObject()
                call block%insert('type','resource_link')
                call block%insert('uri',uri)
                call block%insert('name',attachments(i)%name)
            else
                media=image_media_type(attachments(i)%image)
                if(len(media)==0)cycle
                block=QJsonObject()
                call block%insert('type','image')
                call block%insert('data',attachments(i)%image)
                call block%insert('mimeType',media)
            end if
            call content%append(block)
            call block%delete()
        end do
    end subroutine add_opencode_attachments

    subroutine cancel_opencode_prompt(process,id)
        type(QProcess),intent(inout) :: process
        character(*),intent(in) :: id
        type(QJsonObject) :: frame,params
        logical :: sent

        if(process%state()==QProcess_NotRunning)return
        frame=QJsonObject()
        params=QJsonObject()
        call frame%insert('jsonrpc','2.0')
        call frame%insert('method','session/cancel')
        call params%insert('sessionId',id)
        call frame%insert('params',params)
        sent=send_json(process,frame)
        call params%delete()
        call frame%delete()
        if(sent)sent=process%waitForReadyRead(250)
    end subroutine cancel_opencode_prompt

    subroutine final_opencode_updates(process,events)
        type(QProcess),intent(inout) :: process
        type(assistant_events_t),intent(inout) :: events
        type(QJsonObject) :: frame
        character(:),allocatable :: line
        integer(int64) :: start,now,ticks
        logical :: update

        call system_clock(start,ticks)
        do
            if(provider_stopping(events))exit
            call system_clock(now)
            if(now-start>=2*ticks)exit
            if(.not.process%canReadLine())then
                if(.not.process%waitForReadyRead(250))exit
                if(.not.process%canReadLine())cycle
            end if
            line=process%readLine()
            if(.not.parse_object(line,frame))exit
            update=frame%stringValue('method')=='session/update'
            if(update)call relay_assistant_event(events,'OpenCode CLI',frame)
            call frame%delete()
            if(.not.update)exit
        end do
    end subroutine final_opencode_updates

    subroutine acp_configure(process,id,session,key,value,worker,events,permission_index,reply,error)
        type(QProcess),intent(inout) :: process
        integer,intent(in) :: id
        character(*),intent(in) :: session,key,value
        logical,intent(in) :: worker
        type(assistant_events_t),intent(inout) :: events
        integer,intent(inout) :: permission_index
        type(QJsonObject),intent(out) :: reply
        character(:),allocatable,intent(inout) :: error
        type(QJsonObject) :: params

        params=QJsonObject()
        call params%insert('sessionId',session)
        call params%insert('configId',key)
        call params%insert('value',value)
        call acp_call(process,id,'session/set_config_option',params,15,worker,events,permission_index,reply,error)
        call params%delete()
    end subroutine acp_configure

    subroutine acp_call(process,id,method,params,timeout,worker,events,permission_index,reply,error)
        type(QProcess),intent(inout) :: process
        integer,intent(in) :: id,timeout
        character(*),intent(in) :: method
        type(QJsonObject),intent(in) :: params
        logical,intent(in) :: worker
        type(assistant_events_t),intent(inout) :: events
        integer,intent(inout) :: permission_index
        type(QJsonObject),intent(out) :: reply
        character(:),allocatable,intent(inout) :: error
        type(QJsonObject) :: frame,issue
        character(:),allocatable :: line,route
        integer(int64) :: start,now,ticks
        logical :: ready,handled,interrupted

        reply=QJsonObject()
        frame=QJsonObject()
        call frame%insert('jsonrpc','2.0')
        call frame%insert('id',id)
        call frame%insert('method',method)
        call frame%insert('params',params)
        ready=send_json(process,frame)
        call frame%delete()
        if(.not.ready)then
            error='ForCAD Studio could not send an OpenCode request.'
            return
        end if
        call system_clock(start,ticks)
        do
            interrupted=provider_stopping(events)
            if(interrupted)then
                if(method/='session/new'.or..not.process%canReadLine())exit
            end if
            call capture_failure_hint(process,events)
            if(process%canReadLine())then
                line=process%readLine()
                if(.not.parse_object(line,frame))then
                    error='OpenCode returned an invalid response.'
                    return
                end if
                route=frame%stringValue('method')
                select case(route)
                case('session/update')
                    if(method/='session/load'.and..not.interrupted)call relay_assistant_event(events,'OpenCode CLI',frame)
                case('session/request_permission')
                    if(frame%contains('id').and..not.interrupted)then
                        handled=acp_permission(process,frame,worker,events,permission_index)
                        if(.not.handled)error='ForCAD Studio could not respond to the OpenCode permission request.'
                    end if
                case('')
                    if(frame%intValue('id',-1)==id)then
                        if(frame%contains('result').and..not.frame%contains('error'))then
                            call reply%delete()
                            reply=frame%objectValue('result')
                        else
                            issue=frame%objectValue('error')
                            error=acp_error(issue)
                            call issue%delete()
                        end if
                        call frame%delete()
                        return
                    end if
                case default
                    if(frame%contains('id'))then
                        if(.not.reject_acp_method(process,frame))error='ForCAD Studio could not respond to an OpenCode request.'
                    end if
                end select
                call frame%delete()
                if(len(error)>0)return
                call system_clock(start)
            else
                if(process%state()==QProcess_NotRunning)exit
                call system_clock(now)
                if(now-start>=int(timeout,int64)*ticks)then
                    error='OpenCode did not respond in time.'
                    return
                end if
                ready=process%waitForReadyRead(80)
            end if
        end do
        error='OpenCode closed its connection before completing the request.'
    end subroutine acp_call

    logical function reject_acp_method(process,request) result(ok)
        type(QProcess),intent(inout) :: process
        type(QJsonObject),intent(in) :: request
        type(QJsonObject) :: response,error

        response=acp_response(request)
        error=QJsonObject()
        call error%insert('code',-32601)
        call error%insert('message','This ForCAD Studio host does not implement this client method.')
        call response%insert('error',error)
        call error%delete()
        ok=send_json(process,response)
        call response%delete()
    end function reject_acp_method

    function acp_response(request) result(response)
        type(QJsonObject),intent(in) :: request
        type(QJsonObject) :: response
        type(QJsonValue) :: id

        response=QJsonObject()
        call response%insert('jsonrpc','2.0')
        id=request%value('id')
        call response%insert('id',id)
        call id%delete()
    end function acp_response

    function acp_error(issue) result(error)
        type(QJsonObject),intent(in) :: issue
        type(QJsonObject) :: data
        character(:),allocatable :: error,name
        character(24) :: number
        integer :: code

        data=issue%objectValue('data')
        name=data%stringValue('errorName')
        call data%delete()
        if(name=='provider.free-limit')then
            error='OpenCode free usage limit reached for this model.'
            return
        end if
        error=assistant_failure_hint(issue%stringValue('message'))
        if(len(error)>0)return

        code=issue%intValue('code',0)
        select case(code)
        case(-32000)
            error='OpenCode needs provider authentication. Check your OpenCode sign-in.'
        case(-32602)
            error='OpenCode rejected the request parameters. Check the selected model and effort.'
        case(-32601)
            error='This OpenCode installation does not support the requested session method.'
        case(0)
            error='OpenCode returned an invalid response.'
        case default
            write(number,'(i0)')code
            error='OpenCode returned an error (code '//trim(number)//').'
        end select
    end function acp_error

    logical function acp_permission(process,request,worker,events,sequence) result(ok)
        type(QProcess),intent(inout) :: process
        type(QJsonObject),intent(in) :: request
        logical,intent(in) :: worker
        type(assistant_events_t),intent(in) :: events
        integer,intent(inout) :: sequence
        type(QJsonObject) :: params,entry,response,result,outcome
        type(QJsonArray) :: options
        character(:),allocatable :: decision,option
        integer :: i

        option=''
        if(.not.worker)then
            params=request%objectValue('params')
            sequence=sequence+1
            decision=wait_for_permission(process,events,permission_summary(params),sequence)
            if(len(decision)>0)then
                options=params%arrayValue('options')
                do i=0,options%size()-1
                    entry=options%objectAt(i)
                    if(entry%stringValue('kind')==decision)option=entry%stringValue('optionId')
                    call entry%delete()
                    if(len(option)>0)exit
                end do
                call options%delete()
            end if
            call params%delete()
        end if
        response=acp_response(request)
        result=QJsonObject()
        outcome=QJsonObject()
        if(len(option)>0)then
            call outcome%insert('outcome','selected')
            call outcome%insert('optionId',option)
        else
            call outcome%insert('outcome','cancelled')
        end if
        call result%insert('outcome',outcome)
        call outcome%delete()
        call response%insert('result',result)
        call result%delete()
        ok=send_json(process,response)
        call response%delete()
    end function acp_permission

    function permission_summary(params) result(summary)
        type(QJsonObject),intent(in) :: params
        type(QJsonObject) :: tool,input,location
        type(QJsonArray) :: locations
        character(:),allocatable :: summary,detail
        integer :: i

        tool=params%objectValue('toolCall')
        input=tool%objectValue('rawInput')
        summary=plain_value(tool,'title')
        if(len(summary)==0)summary=plain_value(tool,'kind')
        if(len(summary)==0)summary='OpenCode tool request'
        detail=plain_value(input,'command')
        if(len(detail)==0)detail=plain_value(input,'path')
        if(len(detail)==0)detail=plain_value(input,'filePath')
        if(len(detail)==0)detail=plain_value(input,'file_path')
        if(len(detail)==0)then
            locations=tool%arrayValue('locations')
            location=locations%objectAt(0)
            detail=plain_value(location,'path')
            call location%delete()
            call locations%delete()
        end if
        if(len(detail)>0.and.index(summary,detail)==0)summary=summary//' · '//detail
        summary=clip_assistant_text(clean_assistant_text(summary,.true.),600)
        do i=1,len(summary)
            if(summary(i:i)==lf.or.summary(i:i)==achar(9))summary(i:i)=' '
        end do
        call input%delete()
        call tool%delete()
    end function permission_summary

    function plain_value(object,key) result(text)
        type(QJsonObject),intent(in) :: object
        character(*),intent(in) :: key
        character(:),allocatable :: text
        type(QJsonArray) :: parts
        integer :: i

        text=object%stringValue(key)
        if(object%valueType(key)/=QJsonValue_Array)return
        parts=object%arrayValue(key)
        do i=0,parts%size()-1
            if(parts%typeAt(i)/=QJsonValue_String)then
                text=''
                exit
            end if
            if(i>0)text=text//' '
            text=text//parts%stringAt(i)
        end do
        call parts%delete()
    end function plain_value

    function wait_for_permission(process,events,description,sequence) result(choice)
        type(QProcess),intent(inout) :: process
        type(assistant_events_t),intent(in) :: events
        character(*),intent(in) :: description
        integer,intent(in) :: sequence
        character(:),allocatable :: choice,payload,permission,decision,prefix
        character(24) :: id
        integer :: unit,status,bytes
        logical :: ready

        choice=''
        permission=events%exchange//'/agent-permission'
        decision=events%exchange//'/agent-decision'
        write(id,'(i0)')sequence
        prefix=trim(id)//lf
        if(.not.write_private_file(permission,prefix//description//lf))return
        do while(process%state()/=QProcess_NotRunning)
            if(provider_stopping(events))exit
            open(newunit=unit,file=decision,status='old',access='stream',form='unformatted',iostat=status)
            if(status==0)then
                inquire(unit=unit,size=bytes,iostat=status)
                if(status==0.and.bytes>0.and.bytes<=1024)then
                    allocate(character(bytes)::payload)
                    read(unit,iostat=status)payload
                    if(status==0)then
                        if(payload==prefix//'allow'//lf)choice='allow_once'
                        if(payload==prefix//'deny'//lf)choice='reject_once'
                    end if
                    deallocate(payload)
                end if
                close(unit,status='delete',iostat=status)
                if(len(choice)>0)exit
            end if
            ready=process%waitForFinished(80)
        end do
        open(newunit=unit,file=permission,status='old',iostat=status)
        if(status==0)then
            close(unit,status='delete',iostat=status)
            if(status/=0)choice=''
        end if
    end function wait_for_permission

    function opencode_effort(reply,wanted) result(id)
        type(QJsonObject),intent(in) :: reply
        character(*),intent(in) :: wanted
        character(:),allocatable :: id,candidate
        type(QJsonArray) :: configs,options
        type(QJsonObject) :: config
        integer :: i

        id=''
        configs=reply%arrayValue('configOptions')
        do i=0,configs%size()-1
            config=configs%objectAt(i)
            candidate=config%stringValue('id',config%stringValue('configId'))
            if(config%stringValue('category')=='thought_level'.or.candidate=='effort')then
                options=config%arrayValue('options')
                if(offered(options,wanted,0))id=candidate
                call options%delete()
            end if
            call config%delete()
            if(len(id)>0)exit
        end do
        call configs%delete()
    end function opencode_effort

    recursive logical function offered(options,wanted,depth) result(found)
        type(QJsonArray),intent(in) :: options
        character(*),intent(in) :: wanted
        integer,intent(in) :: depth
        type(QJsonObject) :: option
        type(QJsonArray) :: groups
        integer :: i

        found=.false.
        if(depth>16)return
        do i=0,options%size()-1
            option=options%objectAt(i)
            found=option%stringValue('value')==wanted
            if(.not.found.and.option%valueType('options')==QJsonValue_Array)then
                groups=option%arrayValue('options')
                found=offered(groups,wanted,depth+1)
                call groups%delete()
            end if
            call option%delete()
            if(found)return
        end do
    end function offered

    logical function send_json(process,value) result(ok)
        type(QProcess),intent(inout) :: process
        type(QJsonObject),intent(in) :: value
        character(:),allocatable :: encoded

        encoded=value%toJson(.true.)//lf
        ok=process%write(encoded)==len(encoded)
    end function send_json

    logical function parse_object(line,object) result(ok)
        character(*),intent(in) :: line
        type(QJsonObject),intent(out) :: object
        type(QJsonDocument) :: document

        document=QJsonDocument_fromJson(line)
        ok=document%isObject()
        if(ok)object=document%object()
        call document%delete()
    end function parse_object

    subroutine stream_provider(process,provider,events,completed,worker,saved_session,error)
        type(QProcess),intent(inout) :: process
        character(*),intent(in) :: provider,saved_session
        type(assistant_events_t),intent(inout) :: events
        logical,intent(out) :: completed
        logical,intent(in) :: worker
        character(:),allocatable,intent(inout) :: error
        character(:),allocatable :: line,id
        type(QJsonObject) :: event
        logical :: ready,session_confirmed,interrupted

        completed=.false.
        session_confirmed=.false.
        do
            interrupted=provider_stopping(events)
            if(events%failed)exit
            if(interrupted)then
                if(session_confirmed.or..not.process%canReadLine())return
            end if
            call capture_failure_hint(process,events)
            if(process%canReadLine())then
                line=process%readLine()
            else if(process%state()==QProcess_NotRunning)then
                line=process%readAllStandardOutput()
                if(len(line)==0)exit
            else
                ready=process%waitForReadyRead(80)
                cycle
            end if
            if(parse_object(line,event))then
                if(.not.worker)then
                    if(read_provider_session_event(provider,event,id))then
                        call capture_provider_session(provider,events,worker,id,saved_session,error)
                        session_confirmed=len(error)==0
                    end if
                end if
                if(len(error)==0.and..not.interrupted)then
                    if(len(saved_session)==0.or.session_confirmed)call relay_assistant_event(events,provider,event)
                end if
                call event%delete()
                if(len(error)>0)return
            end if
        end do
        if(len(saved_session)>0.and..not.session_confirmed)then
            call capture_failure_hint(process,events)
            if(allocated(events%failure_hint))then
                error=provider//' · '//events%failure_hint
            else
                error=resume_failure(provider,'')
            end if
            return
        end if
        completed=process%exitStatus()==QProcess_NormalExit.and.process%exitCode()==0
    end subroutine stream_provider

    logical function read_provider_session_event(provider,event,id) result(found)
        character(*),intent(in) :: provider
        type(QJsonObject),intent(in) :: event
        character(:),allocatable,intent(out) :: id

        found=.false.
        select case(provider)
        case('Codex CLI')
            found=event%stringValue('type')=='thread.started'
            if(found)id=event%stringValue('thread_id')
        case('Claude Code CLI')
            found=event%stringValue('type')=='system'.and.event%stringValue('subtype')=='init'
            if(found)id=event%stringValue('session_id')
        case default
            return
        end select
    end function read_provider_session_event

    subroutine capture_provider_session(provider,events,worker,id,expected,error)
        character(*),intent(in) :: provider,id,expected
        type(assistant_events_t),intent(in) :: events
        logical,intent(in) :: worker
        character(:),allocatable,intent(inout) :: error

        if(.not.valid_provider_session(id))then
            error=provider//' did not return a usable session ID.'
        else if(len(expected)>0.and.id/=expected)then
            error=resume_failure(provider,'The provider returned a different session ID.')
        else if(.not.worker)then
            if(.not.write_private_file(events%exchange//'/provider-session',id))then
                error='ForCAD Studio could not save the provider session ID. Check the project working folder.'
            end if
        end if
    end subroutine capture_provider_session

    function resume_failure(provider,detail) result(error)
        character(*),intent(in) :: provider,detail
        character(:),allocatable :: error

        error=provider//' could not resume the saved conversation.'
        if(len(detail)>0)error=error//' '//detail
        error=error//' Verify that the session is available in this working folder, or start a new conversation.'
    end function resume_failure

    logical function provider_stopping(events) result(stopping)
        type(assistant_events_t),intent(in) :: events

        stopping=QThread_isInterruptionRequested()
        if(stopping.or..not.allocated(events%exchange))return
        inquire(file=events%exchange//'/stop-requested',exist=stopping)
    end function provider_stopping

    subroutine capture_failure_hint(process,events)
        type(QProcess),intent(inout) :: process
        type(assistant_events_t),intent(inout) :: events
        character(:),allocatable :: diagnostic,hint

        diagnostic=process%readAllStandardError()
        if(len(diagnostic)==0)return
        hint=assistant_failure_hint(diagnostic(:min(len(diagnostic),8192)))
        if(len(hint)>0.and..not.allocated(events%failure_hint))events%failure_hint=hint
    end subroutine capture_failure_hint

    subroutine close_provider(process)
        type(QProcess),intent(inout) :: process
        logical :: stopped

        call process%closeWriteChannel()
        if(process%state()/=QProcess_NotRunning)then
            if(.not.process%waitForFinished(250))then
                call process%terminate()
                if(.not.process%waitForFinished(250))then
                    call process%kill()
                    stopped=process%waitForFinished(500)
                end if
            end if
        end if
        call process%delete()
    end subroutine close_provider

    function environment_value(name) result(value)
        character(*),intent(in) :: name
        character(:),allocatable :: value
        integer :: length,status

        call get_environment_variable(name,length=length,status=status)
        if(status/=0.or.length<=0)then
            value=''
            return
        end if
        allocate(character(length)::value)
        call get_environment_variable(name,value,status=status)
        if(status/=0)value=''
    end function environment_value

    function image_media_type(image) result(media)
        character(*),intent(in) :: image
        character(:),allocatable :: media

        media=''
        if(len(image)>5592408)return
        if(index(image,'iVBORw0KGgo')==1)media='image/png'
        if(index(image,'/9j/')==1)media='image/jpeg'
    end function image_media_type

    logical function is_inline_reference(path) result(accepted)
        character(*),intent(in) :: path
        character(:),allocatable :: signature
        type(file_metadata) :: metadata

        accepted=.false.
        if(.not.inspect_path(path,metadata))return
        if(.not.is_regular_file(metadata).or.metadata%bytes<12.or.metadata%bytes>20971520)return
        signature=read_text(path,12,0)
        if(len(signature)/=12)return
        accepted=signature(:8)==achar(137)//'PNG'//achar(13)//lf//achar(26)//lf.or. &
            signature(:3)==achar(255)//achar(216)//achar(255).or.signature(:6)=='GIF87a'.or. &
            signature(:6)=='GIF89a'.or.(signature(:4)=='RIFF'.and.signature(9:12)=='WEBP').or.signature(:5)=='%PDF-'
    end function is_inline_reference

    function file_uri(path) result(uri)
        character(*),intent(in) :: path
        character(:),allocatable :: uri,local_url
        character(*),parameter :: hex='0123456789ABCDEF',allowed= &
            'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-._~:/%'
        integer :: i,byte,upper,lower

        local_url=QUrl_fromLocalFile(path)
        uri=''
        do i=1,len(local_url)
            if(index(allowed,local_url(i:i))>0)then
                uri=uri//local_url(i:i)
            else
                byte=iachar(local_url(i:i))
                upper=ishft(byte,-4)+1
                lower=iand(byte,15)+1
                uri=uri//'%'//hex(upper:upper)//hex(lower:lower)
            end if
        end do
    end function file_uri

    function decode_base64(encoded) result(decoded)
        character(*),intent(in) :: encoded
        character(:),allocatable :: decoded
        integer :: source,target,count,values(4),i,padding

        count=len(encoded)
        if(count==0.or.mod(count,4)/=0)then
            decoded=''
            return
        end if
        padding=0
        if(encoded(count:count)=='=')padding=padding+1
        if(encoded(count-1:count-1)=='=')padding=padding+1
        allocate(character(3*(count/4)-padding)::decoded)
        do source=1,count,4
            target=1+3*((source-1)/4)
            do i=1,4
                values(i)=index(base64_letters,encoded(source+i-1:source+i-1))-1
                if(values(i)>=0)cycle
                if(source+i-1>count-padding.and.encoded(source+i-1:source+i-1)=='=')then
                    values(i)=0
                else
                    decoded=''
                    return
                end if
            end do
            decoded(target:target)=achar(4*values(1)+ishft(values(2),-4))
            if(target+1<=len(decoded))decoded(target+1:target+1)=achar(16*iand(values(2),15)+ishft(values(3),-2))
            if(target+2<=len(decoded))decoded(target+2:target+2)=achar(64*iand(values(3),3)+values(4))
        end do
    end function decode_base64

    logical function write_private_file(path,text) result(ok)
        character(*),intent(in) :: path,text
        type(QSaveFile) :: file

        file=QSaveFile(path)
        ok=.false.
        if(file%open())then
            if(file%write(text)==len(text))then
                ok=file%commit()
            else
                call file%cancelWriting()
            end if
        end if
        call file%delete()
    end function write_private_file

end module studio_assistant_cli
