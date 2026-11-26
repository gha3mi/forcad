module studio_assistant_mcp_client
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QProcess,QProcess_NotRunning,QThread_isInterruptionRequested,QCoreApplication_applicationFilePath, &
        QJsonObject,QJsonValue_Object
    use studio_assistant_local_http, only: parse_json_object
    implicit none(type,external)
    private
    public :: assistant_mcp_client_t,connect_assistant_mcp,request_assistant_mcp,close_assistant_mcp

    character(*),parameter :: newline=new_line('a')
    integer,parameter :: connect_timeout=8000,tool_timeout=600000,response_limit=16777216

    type :: assistant_mcp_client_t
        private
        type(QProcess) :: process
        logical :: created=.false.
        integer :: next_id=0
    end type assistant_mcp_client_t

contains

    subroutine connect_assistant_mcp(client,session,error,executable)
        type(assistant_mcp_client_t),intent(inout) :: client
        character(*),intent(in) :: session
        character(:),allocatable,intent(out) :: error
        character(*),optional,intent(in) :: executable
        type(QJsonObject) :: parameters,hello,capabilities,client_info,reply
        character(:),allocatable :: program,frame
        character(:),allocatable :: arguments(:)

        call close_assistant_mcp(client)
        error=''
        if(len(session)==0)then
            error='The ForCAD Studio connection is unavailable.'
            return
        end if
        if(present(executable))then
            program=executable
        else
            program=QCoreApplication_applicationFilePath()
        end if
        if(len(program)==0)then
            error='The ForCAD Studio tool executable could not be located.'
            return
        end if
        allocate(character(max(len(session),len('--mcp-stdio')))::arguments(3))
        arguments(1)='--mcp-stdio'
        arguments(2)='--session'
        arguments(3)=session
        client%process=QProcess()
        client%created=.true.
        call client%process%start(program,arguments)
        if(.not.client%process%waitForStarted(connect_timeout))then
            error='The ForCAD Studio tool process could not be started.'
            call close_assistant_mcp(client)
            return
        end if
        parameters=QJsonObject()
        capabilities=QJsonObject()
        client_info=QJsonObject()
        call client_info%insert('name','ForCAD Studio')
        call client_info%insert('version','1')
        call parameters%insert('clientInfo',client_info)
        call parameters%insert('protocolVersion','2025-11-25')
        call parameters%insert('capabilities',capabilities)
        call client_info%delete()
        call capabilities%delete()
        call request_assistant_mcp(client,'initialize',parameters,reply,error,connect_timeout)
        call parameters%delete()
        if(len(error)==0)then
            hello=reply%objectValue('result')
            capabilities=hello%objectValue('capabilities')
            if(capabilities%valueType('tools')/=QJsonValue_Object)error='The ForCAD Studio tool process did not initialize.'
            call capabilities%delete()
            call hello%delete()
        end if
        call reply%delete()
        if(len(error)==0)then
            frame='{"jsonrpc":"2.0","method":"notifications/initialized"}'//newline
            if(client%process%write(frame)/=len(frame))error='The ForCAD Studio tool process could not be initialized.'
        end if
        if(len(error)>0)call close_assistant_mcp(client)
    end subroutine connect_assistant_mcp

    subroutine request_assistant_mcp(client,method,params,response,error,timeout_ms)
        type(assistant_mcp_client_t),intent(inout) :: client
        character(*),intent(in) :: method
        type(QJsonObject),intent(in) :: params
        type(QJsonObject),intent(out) :: response
        character(:),allocatable,intent(out) :: error
        integer,optional,intent(in) :: timeout_ms
        type(QJsonObject) :: request
        character(:),allocatable :: frame,discard
        integer(int64) :: started,now,rate
        integer :: timeout
        logical :: ready,valid

        response=QJsonObject()
        error=''
        if(.not.client%created)then
            error='The ForCAD Studio connection is unavailable.'
            return
        end if
        timeout=tool_timeout
        if(present(timeout_ms))timeout=timeout_ms
        client%next_id=client%next_id+1
        request=QJsonObject()
        call request%insert('jsonrpc','2.0')
        call request%insert('id',client%next_id)
        call request%insert('method',method)
        call request%insert('params',params)
        frame=request%toJson(.true.)//newline
        call request%delete()
        if(client%process%write(frame)/=len(frame))then
            error='The ForCAD Studio tool process did not accept the request.'
            return
        end if
        call system_clock(started,rate)
        do
            if(QThread_isInterruptionRequested())then
                call cancel_assistant_request(client)
                error='The ForCAD Studio tool request was stopped.'
                return
            end if
            if(client%process%bytesAvailable()>response_limit)then
                error='The ForCAD Studio tool response exceeds the transport size limit.'
                return
            end if
            if(client%process%canReadLine())exit
            if(client%process%state()==QProcess_NotRunning)then
                error='The ForCAD Studio connection closed while a tool was running.'
                return
            end if
            call system_clock(now)
            if(now-started>=int(timeout,int64)*rate/1000_int64)then
                call cancel_assistant_request(client)
                error='The ForCAD Studio tool request timed out.'
                return
            end if
            ready=client%process%waitForReadyRead(25)
            discard=client%process%readAllStandardError()
        end do
        frame=client%process%readLine()
        call response%delete()
        response=parse_json_object(frame,valid)
        if(.not.valid.or.response%stringValue('jsonrpc')/='2.0'.or.response%intValue('id',-1)/=client%next_id.or. &
            (.not.response%contains('result').and..not.response%contains('error'))) then
          error='The ForCAD Studio tool process returned an invalid response.'
        end if
    end subroutine request_assistant_mcp

    subroutine cancel_assistant_request(client)
        type(assistant_mcp_client_t),intent(inout) :: client
        type(QJsonObject) :: notice,params
        character(:),allocatable :: frame
        logical :: sent

        notice=QJsonObject()
        params=QJsonObject()
        call params%insert('requestId',client%next_id)
        call notice%insert('jsonrpc','2.0')
        call notice%insert('method','notifications/cancelled')
        call notice%insert('params',params)
        frame=notice%toJson(.true.)//newline
        call params%delete()
        call notice%delete()
        if(client%process%write(frame)==len(frame))sent=client%process%waitForBytesWritten(100)
        call close_assistant_mcp(client)
    end subroutine cancel_assistant_request

    subroutine close_assistant_mcp(client)
        type(assistant_mcp_client_t),intent(inout) :: client
        logical :: finished

        if(.not.client%created)return
        call client%process%closeWriteChannel()
        finished=client%process%waitForFinished(250)
        if(client%process%state()/=QProcess_NotRunning)then
            call client%process%terminate()
            finished=client%process%waitForFinished(250)
        end if
        if(client%process%state()/=QProcess_NotRunning)then
            call client%process%kill()
            finished=client%process%waitForFinished(500)
        end if
        call client%process%delete()
        client%created=.false.
        client%next_id=0
    end subroutine close_assistant_mcp

end module studio_assistant_mcp_client
