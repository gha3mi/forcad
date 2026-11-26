module studio_assistant_catalog
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: qt_string,QJsonValue,QJsonValue_Object,QProcess,QProcess_NotRunning,QProcess_NormalExit,QThread_isInterruptionRequested, &
        QProcessEnvironment,QProcessEnvironment_systemEnvironment, &
        QJsonObject,QJsonArray,QJsonValue_Array,QJsonValue_Bool,QJsonValue_String
    use studio_assistant_local_http, only: assistant_environment,local_assistant_models,local_model_info, &
        local_model_efforts,local_model_metadata,has_json_string,parse_json_object,assistant_effort_label
    use studio_assistant_process, only: start_provider_process
    implicit none(type,external)
    private
    public :: discover_assistant_catalog

    character(*),parameter :: newline=new_line('a')
    integer,parameter :: metadata_limit=1048576

contains

    subroutine discover_assistant_catalog(workspace,models,catalog,error)
        character(*),intent(in) :: workspace
        character(:),allocatable,intent(out) :: models,catalog,error
        type(QProcess) :: codex,claude,opencode
        type(QJsonObject) :: combined
        type(QJsonArray) :: codex_rows,claude_rows,opencode_rows,local_rows
        integer(int64) :: codex_started,claude_started,opencode_started,rate
        logical :: has_codex,has_claude,has_opencode,codex_ready,claude_ready

        models=''
        catalog='{}'
        error=''
        codex=QProcess()
        claude=QProcess()
        opencode=QProcess()
        call system_clock(codex_started,rate)
        has_codex=start_catalog_process(codex,workspace,'codex',[character(24)::'app-server','--listen','stdio://'])
        codex_ready=.false.
        if(has_codex)codex_ready=write_catalog_frame(codex, &
            '{"id":1,"method":"initialize","params":{"clientInfo":{"name":"forcad_studio","version":"0.1.0"}}}')
        call system_clock(claude_started)
        has_claude=start_catalog_process(claude,workspace,'claude',[character(32)::'--print','--verbose', &
            '--input-format','stream-json','--output-format','stream-json','--no-session-persistence','--strict-mcp-config'])
        claude_ready=.false.
        if(has_claude)claude_ready=write_catalog_frame(claude, &
            '{"type":"control_request","request_id":"forcad-models","request":{"subtype":"initialize"}}')
        call system_clock(opencode_started)
        has_opencode=start_catalog_process(opencode,workspace,'opencode',[character(16)::'--standalone','api','GET','/api/model'])

        if(codex_ready)then
            codex_rows=read_codex_catalog(codex,codex_started,rate)
        else
            codex_rows=QJsonArray()
        end if
        call stop_catalog_process(codex)
        if(claude_ready)then
            claude_rows=read_claude_catalog(claude,remaining_metadata_time(claude_started,rate,7000))
        else
            claude_rows=QJsonArray()
        end if
        call stop_catalog_process(claude)
        if(has_opencode)then
            opencode_rows=read_opencode_catalog(opencode,workspace,remaining_metadata_time(opencode_started,rate,8000))
        else
            opencode_rows=QJsonArray()
        end if
        call stop_catalog_process(opencode)
        local_rows=read_local_catalog()
        call merge_model_vision(opencode_rows,local_rows)

        combined=QJsonObject()
        call combined%insert('Codex CLI',codex_rows)
        call combined%insert('Claude Code CLI',claude_rows)
        call combined%insert('OpenCode CLI',opencode_rows)
        if(has_opencode)models=models//'OpenCode CLI'//newline
        if(has_codex)models=models//'Codex CLI'//newline
        if(has_claude)models=models//'Claude Code CLI'//newline
        if(.not.local_rows%isEmpty())then
            models=models//'Ollama'//newline
            call combined%insert('Ollama',local_rows)
        end if
        catalog=combined%toJson(.true.)
        if(len(models)==0)error='No Assistant providers are available. Start a local service or configure an installed agent.'
        call local_rows%delete()
        call opencode_rows%delete()
        call codex_rows%delete()
        call claude_rows%delete()
        call combined%delete()
    end subroutine discover_assistant_catalog

    integer function remaining_metadata_time(started,rate,limit) result(remaining)
        integer(int64),intent(in) :: started,rate
        integer,intent(in) :: limit
        integer(int64) :: now

        call system_clock(now)
        remaining=max(0,limit-int((now-started)*1000_int64/rate))
    end function remaining_metadata_time

    logical function start_catalog_process(process,workspace,program,arguments) result(started)
        type(QProcess),intent(inout) :: process
        character(*),intent(in) :: workspace,program,arguments(:)
        type(QProcessEnvironment) :: environment

        started=.false.
        if(QThread_isInterruptionRequested())return
        call process%setWorkingDirectory(workspace)
        if(program=='opencode')then
            environment=QProcessEnvironment_systemEnvironment()
            call environment%insert('OPENCODE_DISABLE_PROJECT_CONFIG','true')
            call process%setProcessEnvironment(environment)
            call environment%delete()
        end if
        call start_provider_process(process,program,arguments)
        started=process%waitForStarted(1500)
    end function start_catalog_process

    subroutine stop_catalog_process(process)
        type(QProcess),intent(inout) :: process
        logical :: finished

        if(process%state()/=QProcess_NotRunning)then
            call process%closeWriteChannel()
            call process%terminate()
            finished=process%waitForFinished(200)
        end if
        if(process%state()/=QProcess_NotRunning)then
            call process%kill()
            finished=process%waitForFinished(500)
        end if
        call process%delete()
    end subroutine stop_catalog_process

    logical function write_catalog_frame(process,payload) result(ok)
        type(QProcess),intent(inout) :: process
        character(*),intent(in) :: payload

        ok=process%write(payload//newline)==len(payload)+1
    end function write_catalog_frame

    logical function read_catalog_frame(process,id,frame,timeout_ms) result(ok)
        type(QProcess),intent(inout) :: process
        integer,intent(in) :: id,timeout_ms
        type(QJsonObject),intent(out) :: frame
        type(QJsonObject) :: response
        character(:),allocatable :: line,discard
        integer(int64) :: started,now,rate,total
        logical :: ready,valid

        ok=.false.
        total=0
        frame=QJsonObject()
        call system_clock(started,rate)
        do
            if(QThread_isInterruptionRequested())return
            if(total+process%bytesAvailable()>metadata_limit)return
            if(process%canReadLine())then
                line=process%readLine()
                total=total+len(line)
                call frame%delete()
                frame=parse_json_object(line,valid)
                if(.not.valid)cycle
                if(id>0)then
                    if(frame%intValue('id',-1)/=id)cycle
                else
                    if(frame%stringValue('type')/='control_response')cycle
                    response=frame%objectValue('response')
                    valid=response%stringValue('request_id')=='forcad-models'.and.response%stringValue('subtype')=='success'
                    call response%delete()
                    if(.not.valid)cycle
                end if
                ok=.true.
                return
            end if
            if(process%state()==QProcess_NotRunning)return
            call system_clock(now)
            if(now-started>=int(timeout_ms,int64)*rate/1000_int64)return
            ready=process%waitForReadyRead(25)
            discard=process%readAllStandardError()
        end do
    end function read_catalog_frame

    function read_catalog_output(process,timeout_ms,ok) result(output)
        type(QProcess),intent(inout) :: process
        integer,intent(in) :: timeout_ms
        logical,intent(out) :: ok
        character(:),allocatable :: output,discard
        integer(int64) :: started,now,rate
        logical :: ready

        ok=.false.
        output=''
        call system_clock(started,rate)
        do
            if(QThread_isInterruptionRequested())return
            if(len(output)+process%bytesAvailable()>metadata_limit)return
            output=output//process%readAllStandardOutput()
            discard=process%readAllStandardError()
            if(process%state()==QProcess_NotRunning)exit
            call system_clock(now)
            if(now-started>=int(timeout_ms,int64)*rate/1000_int64)return
            ready=process%waitForReadyRead(25)
        end do
        ok=process%exitStatus()==QProcess_NormalExit.and.process%exitCode()==0
    end function read_catalog_output

    function read_codex_catalog(process,started,rate) result(models)
        type(QProcess),intent(inout) :: process
        integer(int64),intent(in) :: started,rate
        type(QJsonArray) :: models,all_rows,data,modalities,details,effort_list
        type(QJsonObject) :: frame,response,request,params,row,new_row,effort
        character(:),allocatable :: cursor,next_cursor,model,text
        integer :: page,i,j
        logical :: accepted,default_added

        models=QJsonArray()
        accepted=read_catalog_frame(process,1,frame,min(5000,remaining_metadata_time(started,rate,7000)))
        accepted=accepted.and.frame%contains('result')
        call frame%delete()
        if(.not.accepted)return
        if(.not.write_catalog_frame(process,'{"method":"initialized"}'))return

        all_rows=QJsonArray()
        cursor=''
        default_added=.false.
        do page=1,20
            if(remaining_metadata_time(started,rate,7000)==0)exit
            request=QJsonObject()
            params=QJsonObject()
            call request%insert('id',2)
            call request%insert('method','model/list')
            call params%insert('includeHidden',.false.)
            call params%insert('limit',100)
            if(len(cursor)>0)then
                call params%insert('cursor',cursor)
            else
                call params%insertNull('cursor')
            end if
            call request%insert('params',params)
            accepted=write_catalog_frame(process,request%toJson(.true.))
            call params%delete()
            call request%delete()
            if(.not.accepted)exit
            accepted=read_catalog_frame(process,2,frame,min(5000,remaining_metadata_time(started,rate,7000)))
            response=frame%objectValue('result')
            call frame%delete()
            data=response%arrayValue('data')
            next_cursor=response%stringValue('nextCursor')
            call response%delete()
            do i=0,data%size()-1
                row=data%objectAt(i)
                model=row%stringValue('model',row%stringValue('id'))
                if(len(model)>0.and..not.row%boolValue('hidden'))then
                    effort_list=QJsonArray()
                    details=row%arrayValue('supportedReasoningEfforts')
                    do j=0,details%size()-1
                        effort=details%objectAt(j)
                        text=effort%stringValue('reasoningEffort')
                        if(len(text)>0)call add_effort(effort_list,text,assistant_effort_label(text))
                        call effort%delete()
                    end do
                    call details%delete()
                    modalities=row%arrayValue('inputModalities')
                    new_row=model_catalog_row(model,row%stringValue('displayName',model),effort_list)
                    call new_row%insert('vision',has_json_string(modalities,'image'))
                    call modalities%delete()
                    call effort_list%delete()
                    call all_rows%append(new_row)
                    if(.not.default_added.and.row%boolValue('isDefault'))then
                        call new_row%insert('id','')
                        call new_row%insert('label','')
                        call models%append(new_row)
                        default_added=.true.
                    end if
                    call new_row%delete()
                end if
                call row%delete()
            end do
            call data%delete()
            if(.not.accepted.or.len(next_cursor)==0.or.next_cursor==cursor)exit
            cursor=next_cursor
        end do
        do i=0,all_rows%size()-1
            row=all_rows%objectAt(i)
            call models%append(row)
            call row%delete()
        end do
        call all_rows%delete()
    end function read_codex_catalog

    function read_claude_catalog(process,timeout_ms) result(models)
        type(QProcess),intent(inout) :: process
        integer,intent(in) :: timeout_ms
        type(QJsonArray) :: models,data,levels,efforts
        type(QJsonObject) :: frame,inner,response,row,new_row
        character(:),allocatable :: model,effort,label
        integer :: i,j
        logical :: accepted

        models=QJsonArray()
        accepted=read_catalog_frame(process,0,frame,min(5000,timeout_ms))
        if(.not.accepted)then
            call frame%delete()
            return
        end if
        inner=frame%objectValue('response')
        response=inner%objectValue('response')
        data=response%arrayValue('models')
        call response%delete()
        call inner%delete()
        call frame%delete()
        do i=0,data%size()-1
            row=data%objectAt(i)
            if(row%valueType('value')==QJsonValue_String)then
                model=row%stringValue('value')
                label=row%stringValue('displayName',model)
                if(model=='default')model=''
                efforts=QJsonArray()
                if(row%boolValue('supportsEffort',.true.))then
                    levels=row%arrayValue('supportedEffortLevels')
                    do j=0,levels%size()-1
                        effort=levels%stringAt(j)
                        if(len(effort)>0)call add_effort(efforts,effort,assistant_effort_label(effort))
                    end do
                    call levels%delete()
                end if
                new_row=model_catalog_row(model,label,efforts)
                call models%append(new_row)
                call new_row%delete()
                call efforts%delete()
            end if
            call row%delete()
        end do
        call data%delete()
    end function read_claude_catalog

    function read_opencode_catalog(process,workspace,timeout_ms) result(models)
        type(QProcess),intent(inout) :: process
        character(*),intent(in) :: workspace
        integer,intent(in) :: timeout_ms
        type(QJsonArray) :: models,data,rows,modalities,details,efforts
        type(QJsonObject) :: response,item,new_row,capabilities,effort
        character(:),allocatable :: output,model,provider,text,preferred
        integer :: i,j
        logical :: accepted

        models=QJsonArray()
        output=read_catalog_output(process,timeout_ms,accepted)
        response=parse_json_object(output)
        accepted=accepted.and.response%valueType('data')==QJsonValue_Array
        data=response%arrayValue('data')
        call response%delete()
        if(.not.accepted.or.data%isEmpty())then
            call data%delete()
            call stop_catalog_process(process)
            process=QProcess()
            if(start_catalog_process(process,workspace,'opencode',[character(14)::'--standalone','models']))then
                output=read_catalog_output(process,6000,accepted)
                if(accepted)call add_opencode_list(models,output)
            end if
            if(.not.models%isEmpty())call read_opencode_efforts(models,workspace)
            return
        end if
        rows=QJsonArray()
        preferred=assistant_environment('FORCAD_OPENCODE_MODEL')
        do i=0,data%size()-1
            item=data%objectAt(i)
            model=item%stringValue('id')
            provider=item%stringValue('providerID')
            if(len(model)>0.and.len(provider)>0.and.item%boolValue('enabled',.true.))then
                efforts=QJsonArray()
                details=opencode_variants(item)
                do j=0,details%size()-1
                    effort=details%objectAt(j)
                    text=effort%stringValue('id')
                    if(len(text)>0.and.text/='default')call add_effort(efforts,text,effort%stringValue('name',assistant_effort_label(text)))
                    call effort%delete()
                end do
                call details%delete()
                new_row=model_catalog_row(provider//'/'//model,item%stringValue('name',model)//' · '//provider,efforts)
                call efforts%delete()
                capabilities=item%objectValue('capabilities')
                modalities=capabilities%arrayValue('input')
                call new_row%insert('vision',has_json_string(modalities,'image'))
                call modalities%delete()
                call capabilities%delete()
                call rows%append(new_row)
                if(provider//'/'//model==preferred.and.models%isEmpty())then
                    call new_row%insert('id','')
                    call new_row%insert('label','')
                    call models%append(new_row)
                end if
                call new_row%delete()
            end if
            call item%delete()
        end do
        call data%delete()
        do i=0,rows%size()-1
            item=rows%objectAt(i)
            call models%append(item)
            call item%delete()
        end do
        call rows%delete()
    end function read_opencode_catalog

    function opencode_variants(item) result(entries)
        type(QJsonObject),intent(in) :: item
        type(QJsonArray) :: entries
        type(QJsonObject) :: variants
        type(QJsonValue) :: value
        type(qt_string),allocatable :: keys(:)
        integer :: i

        if(item%valueType('variants')/=QJsonValue_Object)then
            entries=item%arrayValue('variants')
            return
        end if
        entries=QJsonArray()
        variants=item%objectValue('variants')
        keys=variants%keys()
        do i=1,size(keys)
            value=variants%value(keys(i)%s)
            call entries%append(value)
            call value%delete()
        end do
        call variants%delete()
    end function opencode_variants

    subroutine add_opencode_list(models,output)
        type(QJsonArray),intent(inout) :: models
        character(*),intent(in) :: output
        character(*),parameter :: allowed='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-'
        type(QJsonArray) :: efforts
        type(QJsonObject) :: row
        character(:),allocatable :: lines,model
        character(240),allocatable :: names(:)
        character(240) :: previous
        integer :: start,last,slash,count,n,position

        lines=output//newline
        n=0
        do position=1,len(lines)
            if(lines(position:position)==newline)n=n+1
        end do
        allocate(names(n))
        efforts=QJsonArray()
        start=1
        count=0
        do while(start<=len(lines))
            last=start+index(lines(start:),newline)-2
            model=lines(start:last)
            start=last+2
            slash=index(model,'/')
            if(slash<=1.or.slash>=len(model).or.len(model)>240)cycle
            if(verify(model(:slash-1),allowed)/=0.or.verify(model(slash+1:),allowed//'/:+@#')/=0)cycle
            position=count+1
            do while(position>1)
                if(names(position-1)<=model)exit
                names(position)=names(position-1)
                position=position-1
            end do
            names(position)=model
            count=count+1
        end do
        previous=''
        do position=1,count
            if(names(position)==previous)cycle
            model=trim(names(position))
            row=model_catalog_row(model,model,efforts)
            call models%append(row)
            call row%delete()
            previous=names(position)
        end do
        call efforts%delete()
    end subroutine add_opencode_list

    subroutine read_opencode_efforts(models,workspace)
        type(QJsonArray),intent(inout) :: models
        character(*),intent(in) :: workspace
        type(QProcess) :: process
        type(QJsonObject) :: params,reply,agent,capabilities,row
        type(QJsonArray) :: servers,rows,efforts
        character(:),allocatable :: session
        integer(int64) :: started,rate
        integer :: i,remaining
        logical :: ok,can_delete

        process=QProcess()
        call system_clock(started,rate)
        probe: block
            if(.not.start_catalog_process(process,workspace,'opencode',['acp']))exit probe
            params=QJsonObject()
            capabilities=QJsonObject()
            call params%insert('protocolVersion',1)
            call params%insert('clientCapabilities',capabilities)
            call capabilities%delete()
            ok=opencode_catalog_call(process,1,'initialize',params,remaining_metadata_time(started,rate,8000),reply)
            call params%delete()
            agent=reply%objectValue('agentCapabilities')
            capabilities=agent%objectValue('sessionCapabilities')
            can_delete=capabilities%valueType('delete')==QJsonValue_Object
            call capabilities%delete()
            call agent%delete()
            call reply%delete()
            if(.not.ok.or..not.can_delete)exit probe

            params=QJsonObject()
            servers=QJsonArray()
            call params%insert('cwd',workspace)
            call params%insert('mcpServers',servers)
            call servers%delete()
            ok=opencode_catalog_call(process,2,'session/new',params,remaining_metadata_time(started,rate,8000),reply)
            call params%delete()
            session=reply%stringValue('sessionId')
            call reply%delete()
            if(.not.ok.or.len(session)==0)exit probe

            rows=QJsonArray()
            params=QJsonObject()
            call params%insert('sessionId',session)
            call params%insert('configId','model')
            do i=0,models%size()-1
                row=models%objectAt(i)
                remaining=remaining_metadata_time(started,rate,8000)
                if(remaining>0.and.process%state()/=QProcess_NotRunning)then
                    call params%insert('value',row%stringValue('id'))
                    ok=opencode_catalog_call(process,i+3,'session/set_config_option',params,remaining,reply)
                    if(ok)then
                        efforts=opencode_config_efforts(reply)
                        call row%insert('efforts',efforts)
                        call efforts%delete()
                    end if
                    call reply%delete()
                end if
                call rows%append(row)
                call row%delete()
            end do
            call params%delete()
            call models%delete()
            models=rows

            params=QJsonObject()
            call params%insert('sessionId',session)
            ok=opencode_catalog_call(process,i+3,'session/delete',params,1200,reply)
            call reply%delete()
            call params%delete()
        end block probe
        call stop_catalog_process(process)
    end subroutine read_opencode_efforts

    logical function opencode_catalog_call(process,id,method,params,timeout_ms,reply) result(ok)
        type(QProcess),intent(inout) :: process
        integer,intent(in) :: id,timeout_ms
        character(*),intent(in) :: method
        type(QJsonObject),intent(in) :: params
        type(QJsonObject),intent(out) :: reply
        type(QJsonObject) :: frame

        ok=.false.
        if(timeout_ms>0)then
            frame=QJsonObject()
            call frame%insert('jsonrpc','2.0')
            call frame%insert('id',id)
            call frame%insert('method',method)
            call frame%insert('params',params)
            ok=write_catalog_frame(process,frame%toJson(.true.))
            call frame%delete()
        end if
        if(ok)then
            ok=read_catalog_frame(process,id,frame,timeout_ms)
            ok=ok.and.frame%valueType('result')==QJsonValue_Object
            reply=frame%objectValue('result')
            call frame%delete()
        else
            reply=QJsonObject()
        end if
    end function opencode_catalog_call

    function opencode_config_efforts(reply) result(efforts)
        type(QJsonObject),intent(in) :: reply
        type(QJsonArray) :: efforts,configs,options
        type(QJsonObject) :: config
        integer :: i

        efforts=QJsonArray()
        configs=reply%arrayValue('configOptions')
        do i=0,configs%size()-1
            config=configs%objectAt(i)
            if(config%stringValue('category')=='thought_level'.or.config%stringValue('id',config%stringValue('configId'))=='effort')then
                options=config%arrayValue('options')
                call add_opencode_effort_options(efforts,options,0)
                call options%delete()
            end if
            call config%delete()
        end do
        call configs%delete()
    end function opencode_config_efforts

    recursive subroutine add_opencode_effort_options(efforts,options,depth)
        type(QJsonArray),intent(inout) :: efforts
        type(QJsonArray),intent(in) :: options
        integer,intent(in) :: depth
        type(QJsonObject) :: option
        type(QJsonArray) :: group
        character(:),allocatable :: id
        integer :: i

        if(depth>16)return
        do i=0,options%size()-1
            option=options%objectAt(i)
            id=option%stringValue('value')
            if(len(id)>0.and.id/='default')call add_effort(efforts,id,option%stringValue('name',assistant_effort_label(id)))
            if(option%valueType('options')==QJsonValue_Array)then
                group=option%arrayValue('options')
                call add_opencode_effort_options(efforts,group,depth+1)
                call group%delete()
            end if
            call option%delete()
        end do
    end subroutine add_opencode_effort_options

    function read_local_catalog() result(rows)
        type(QJsonArray) :: rows,models,unpreferred
        type(QJsonObject) :: row
        character(:),allocatable :: error,preferred,model
        integer :: i
        logical :: valid

        rows=QJsonArray()
        call local_assistant_models(models,error)
        if(len(error)>0)then
            call models%delete()
            return
        end if
        preferred=assistant_environment('FORCAD_OLLAMA_MODEL')
        unpreferred=QJsonArray()
        do i=0,models%size()-1
            if(QThread_isInterruptionRequested())exit
            model=models%stringAt(i)
            row=local_catalog_row(model,valid)
            if(valid)then
                if(model==preferred.or.model==preferred//':latest')then
                    call rows%append(row)
                else
                    call unpreferred%append(row)
                end if
            end if
            call row%delete()
        end do
        call models%delete()
        do i=0,unpreferred%size()-1
            row=unpreferred%objectAt(i)
            call rows%append(row)
            call row%delete()
        end do
        call unpreferred%delete()
    end function read_local_catalog

    function local_catalog_row(model,valid) result(row)
        character(*),intent(in) :: model
        logical,intent(out) :: valid
        type(QJsonObject) :: row,info
        type(QJsonArray) :: capabilities,options,efforts
        character(:),allocatable :: error,effort,id
        integer :: i

        call local_model_info(model,1500,info,error)
        if(len(error)>0)then
            call info%delete()
            info=QJsonObject()
        end if
        valid=local_model_metadata(info)
        options=local_model_efforts(info)
        efforts=QJsonArray()
        do i=0,options%size()-1
            if(options%typeAt(i)==QJsonValue_Bool)then
                if(options%boolAt(i))then
                    call add_effort(efforts,'boolean:true','On')
                else
                    call add_effort(efforts,'boolean:false','Off')
                end if
            else
                effort=options%stringAt(i)
                id=effort
                if(index(id,'boolean:')==1.or.index(id,'string:')==1)id='string:'//id
                call add_effort(efforts,id,assistant_effort_label(effort))
            end if
        end do
        row=model_catalog_row(model,model,efforts)
        call efforts%delete()
        call options%delete()
        if(info%valueType('capabilities')==QJsonValue_Array)then
            capabilities=info%arrayValue('capabilities')
            call row%insert('vision',has_json_string(capabilities,'vision'))
            call row%insert('tools',has_json_string(capabilities,'tools'))
            call capabilities%delete()
        else
            call row%insertNull('vision')
            call row%insertNull('tools')
        end if
        call info%delete()
    end function local_catalog_row

    subroutine merge_model_vision(opencode,local)
        type(QJsonArray),intent(inout) :: opencode
        type(QJsonArray),intent(in) :: local
        type(QJsonArray) :: combined
        type(QJsonObject) :: row,ollama
        character(:),allocatable :: id
        integer :: i,j
        logical :: vision

        if(local%isEmpty().or.opencode%isEmpty())return
        combined=QJsonArray()
        do i=0,opencode%size()-1
            row=opencode%objectAt(i)
            id=row%stringValue('id')
            if(index(id,'ollama/')==1)then
                vision=row%boolValue('vision')
                do j=0,local%size()-1
                    ollama=local%objectAt(j)
                    if(ollama%stringValue('id')==id(8:))vision=vision.or.ollama%boolValue('vision')
                    call ollama%delete()
                end do
                call row%insert('vision',vision)
            end if
            call combined%append(row)
            call row%delete()
        end do
        call opencode%delete()
        opencode=combined
    end subroutine merge_model_vision

    function model_catalog_row(id,label,efforts) result(row)
        character(*),intent(in) :: id,label
        type(QJsonArray),intent(in) :: efforts
        type(QJsonObject) :: row

        row=QJsonObject()
        call row%insert('id',id)
        call row%insert('label',label)
        call row%insert('efforts',efforts)
    end function model_catalog_row

    subroutine add_effort(efforts,id,label)
        type(QJsonArray),intent(inout) :: efforts
        character(*),intent(in) :: id,label
        type(QJsonObject) :: effort

        effort=QJsonObject()
        call effort%insert('id',id)
        call effort%insert('label',label)
        call efforts%append(effort)
        call effort%delete()
    end subroutine add_effort

end module studio_assistant_catalog
