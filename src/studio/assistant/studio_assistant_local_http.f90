module studio_assistant_local_http
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QNetworkAccessManager,QNetworkRequest,QNetworkReply,QNetworkReply_NoError,QThread_isInterruptionRequested, &
        QJsonDocument,QJsonDocument_fromJson,QJsonObject,QJsonArray,QJsonValue,QJsonValue_String, &
        QJsonValue_Array,QJsonValue_Bool,QJsonValue_Undefined,QJsonValue_Null
    implicit none(type,external)
    private
    public :: assistant_environment,local_http_request,local_assistant_models,local_model_info,local_model_efforts
    public :: valid_local_model_name,local_model_metadata,has_json_string,parse_json_object,assistant_effort_label

    character(*),parameter :: loopback='http://127.0.0.1:'
    integer,parameter :: body_limit=1048576

    abstract interface
        logical function consume_http_chunk(state,chunk) result(ok)
            implicit none(type,external)
            class(*),intent(inout) :: state
            character(*),intent(in) :: chunk
        end function consume_http_chunk
    end interface

contains

    function assistant_environment(name) result(value)
        character(*),intent(in) :: name
        character(:),allocatable :: value
        integer :: length,status

        call get_environment_variable(name,length=length,status=status)
        value=''
        if(status/=0.or.length<=0)return
        value=repeat(' ',length)
        call get_environment_variable(name,value,status=status)
        if(status/=0)value=''
    end function assistant_environment

    function local_endpoint() result(endpoint)
        character(:),allocatable :: endpoint,port_text
        integer :: port,status

        endpoint=assistant_environment('FORCAD_OLLAMA_URL')
        if(len(endpoint)==0)then
            endpoint=loopback//'11434'
            return
        end if
        if(index(endpoint,loopback)/=1.or.len(endpoint)<=len(loopback))then
            endpoint=''
            return
        end if
        port_text=endpoint(len(loopback)+1:)
        if(len(port_text)>5.or.verify(port_text,'0123456789')/=0)then
            endpoint=''
            return
        end if
        read(port_text,*,iostat=status)port
        if(status/=0)then
            endpoint=''
        else if(port<1.or.port>65535)then
            endpoint=''
        end if
    end function local_endpoint

    subroutine local_http_request(path,payload,timeout_ms,body,error,consume,state)
        character(*),intent(in) :: path,payload
        integer,intent(in) :: timeout_ms
        character(:),allocatable,intent(out) :: body,error
        procedure(consume_http_chunk),optional :: consume
        class(*),optional,intent(inout) :: state
        type(QNetworkAccessManager) :: manager
        type(QNetworkRequest) :: request
        type(QNetworkReply) :: reply
        character(:),allocatable :: endpoint,chunk
        integer(int64) :: started,now,rate
        integer :: total,http_status
        logical :: finished

        body=''
        error=''
        endpoint=local_endpoint()
        if(len(endpoint)==0)then
            error='The local provider endpoint must use HTTP on 127.0.0.1 and a valid port.'
            return
        end if
        if(QThread_isInterruptionRequested())then
            error='The local provider request was stopped.'
            return
        end if
        manager=QNetworkAccessManager()
        request=QNetworkRequest(endpoint//path)
        call request%setFollowRedirects(.false.)
        call request%setTransferTimeout(timeout_ms)
        if(len(payload)==0)then
            reply=manager%get(request)
        else
            call request%setRawHeader('Content-Type','application/json')
            reply=manager%post(request,payload)
        end if
        total=0
        call system_clock(started,rate)
        do
            finished=reply%waitForFinished(25)
            if(QThread_isInterruptionRequested())then
                error='The local provider request was stopped.'
                exit
            end if
            if(reply%bytesAvailable()>body_limit-total)then
                error='The local provider returned too much data.'
                exit
            end if
            if(reply%bytesAvailable()>0)then
                chunk=reply%readAll()
                total=total+len(chunk)
                if(present(consume).and.present(state))then
                    if(.not.consume(state,chunk))then
                        error='The local provider returned an invalid response.'
                        exit
                    end if
                else
                    body=body//chunk
                end if
                call system_clock(now)
                if(present(consume))started=now
            end if
            if(finished)exit
            call system_clock(now)
            if(now-started>=int(timeout_ms,int64)*rate/1000_int64)then
                error='The local provider took too long to answer.'
                exit
            end if
        end do
        if(len(error)==0)then
            http_status=reply%statusCode()
            if(reply%error()/=QNetworkReply_NoError.or.http_status<200.or.http_status>=300) then
              error='The local provider could not complete the HTTP request.'
            end if
        end if
        if(.not.reply%isFinished())call reply%abort()
        call request%delete()
        call manager%delete()
    end subroutine local_http_request

    function parse_json_object(text,valid) result(object)
        character(*),intent(in) :: text
        logical,optional,intent(out) :: valid
        type(QJsonObject) :: object
        type(QJsonDocument) :: document

        document=QJsonDocument_fromJson(text)
        if(present(valid))valid=document%isObject()
        if(document%isObject())then
            object=document%object()
        else
            object=QJsonObject()
        end if
        call document%delete()
    end function parse_json_object

    logical function has_json_string(values,text) result(found)
        type(QJsonArray),intent(in) :: values
        character(*),intent(in) :: text
        integer :: index_value

        found=.false.
        do index_value=0,values%size()-1
            if(values%typeAt(index_value)/=QJsonValue_String)cycle
            if(values%stringAt(index_value)/=text)cycle
            found=.true.
            exit
        end do
    end function has_json_string

    logical function local_model_metadata(info) result(local)
        type(QJsonObject),intent(in) :: info

        local=empty_remote_field(info,'remote_host').and.empty_remote_field(info,'remote_model')
    end function local_model_metadata

    logical function empty_remote_field(info,name) result(empty)
        type(QJsonObject),intent(in) :: info
        character(*),intent(in) :: name

        select case(info%valueType(name))
        case(QJsonValue_Undefined,QJsonValue_Null)
            empty=.true.
        case(QJsonValue_String)
            empty=len(info%stringValue(name))==0
        case(QJsonValue_Bool)
            empty=.not.info%boolValue(name)
        case default
            empty=.false.
        end select
    end function empty_remote_field

    logical function valid_local_model_name(name) result(valid)
        character(*),intent(in) :: name
        character(*),parameter :: alnum='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
        character(:),allocatable :: folded
        integer :: index_char,code,slash,colon
        logical :: part_start

        valid=.false.
        if(len(name)==0.or.len(name)>160)return
        slash=0
        colon=0
        part_start=.true.
        folded=name
        do index_char=1,len(name)
            code=iachar(name(index_char:index_char))
            if(code>=iachar('A').and.code<=iachar('Z'))folded(index_char:index_char)=achar(code+32)
            if(part_start.and.index(alnum,name(index_char:index_char))==0)return
            part_start=.false.
            select case(name(index_char:index_char))
            case('/')
                if(slash/=0.or.colon/=0)return
                slash=index_char
                part_start=.true.
            case(':')
                if(colon/=0)return
                colon=index_char
                part_start=.true.
            case default
                if(index(alnum//'.-_',name(index_char:index_char))==0)return
            end select
        end do
        if(part_start)return
        if(index(folded//':',':cloud:')>0.or.index(folded//':','-cloud:')>0)return
        if(index(folded,'cloud/')==1.or.index(folded,'ollama.com/')==1)return
        valid=.true.
    end function valid_local_model_name

    subroutine local_assistant_models(models,error)
        type(QJsonArray),intent(out) :: models
        character(:),allocatable,intent(out) :: error
        type(QJsonObject) :: response,entry
        type(QJsonArray) :: entries
        character(:),allocatable :: body,name
        character(160),allocatable :: names(:)
        character(160) :: previous
        integer :: index_entry,count,position
        logical :: valid

        models=QJsonArray()
        call local_http_request('/api/tags','',5000,body,error)
        if(len(error)>0)return
        response=parse_json_object(body,valid)
        valid=valid.and.response%valueType('models')==QJsonValue_Array
        if(.not.valid)then
            call response%delete()
            error='The local provider returned an invalid model list.'
            return
        end if
        entries=response%arrayValue('models')
        allocate(names(entries%size()))
        count=0
        do index_entry=0,entries%size()-1
            entry=entries%objectAt(index_entry)
            name=entry%stringValue('name')
            valid=entry%doubleValue('size')>0.and.local_model_metadata(entry).and.valid_local_model_name(name)
            call entry%delete()
            if(.not.valid)cycle
            position=count+1
            do while(position>1)
                if(names(position-1)<=name)exit
                names(position)=names(position-1)
                position=position-1
            end do
            names(position)=name
            count=count+1
        end do
        previous=''
        do position=1,count
            if(names(position)==previous)cycle
            call models%append(trim(names(position)))
            previous=names(position)
        end do
        call entries%delete()
        call response%delete()
    end subroutine local_assistant_models

    subroutine local_model_info(model,timeout_ms,info,error)
        character(*),intent(in) :: model
        integer,intent(in) :: timeout_ms
        type(QJsonObject),intent(out) :: info
        character(:),allocatable,intent(out) :: error
        type(QJsonObject) :: payload
        character(:),allocatable :: body
        logical :: valid

        payload=QJsonObject()
        call payload%insert('model',model)
        call local_http_request('/api/show',payload%toJson(.true.),timeout_ms,body,error)
        call payload%delete()
        info=parse_json_object(body,valid)
        if(len(error)==0.and..not.valid)error='The local provider returned invalid model metadata.'
    end subroutine local_model_info

    function local_model_efforts(info) result(efforts)
        type(QJsonObject),intent(in) :: info
        type(QJsonArray) :: efforts,capabilities,values
        type(QJsonObject) :: thinking,details,model_info
        type(QJsonValue) :: value
        character(:),allocatable :: family,text
        integer :: position,value_type
        logical :: has_values

        efforts=QJsonArray()
        thinking=info%objectValue('thinking')
        value_type=thinking%valueType('values')
        has_values=value_type/=QJsonValue_Null.and.value_type/=QJsonValue_Undefined
        if(value_type==QJsonValue_Bool)has_values=thinking%boolValue('values')
        if(has_values)then
            if(thinking%valueType('values')==QJsonValue_Array)then
                values=thinking%arrayValue('values')
                do position=0,values%size()-1
                    value_type=values%typeAt(position)
                    select case(value_type)
                    case(QJsonValue_Bool)
                        if(values%size()==1.and..not.values%boolAt(position))cycle
                    case(QJsonValue_String)
                        text=values%stringAt(position)
                        if(len(text)==0.or.len(text)>64)cycle
                    case default
                        cycle
                    end select
                    value=values%at(position)
                    call efforts%append(value)
                    call value%delete()
                end do
                call values%delete()
            end if
        else
            capabilities=info%arrayValue('capabilities')
            if(has_json_string(capabilities,'thinking'))then
                details=info%objectValue('details')
                model_info=info%objectValue('model_info')
                family=details%stringValue('family')
                if(len(family)==0)family=model_info%stringValue('general.architecture')
                if(family=='gptoss')then
                    call efforts%append('low')
                    call efforts%append('medium')
                    call efforts%append('high')
                else
                    call efforts%append(.false.)
                    call efforts%append(.true.)
                end if
                call details%delete()
                call model_info%delete()
            end if
            call capabilities%delete()
        end if
        call thinking%delete()
    end function local_model_efforts

    function assistant_effort_label(value) result(label)
        character(*),intent(in) :: value
        character(:),allocatable :: label
        integer :: first

        label=value
        if(len(label)==0)return
        if(label=='xhigh')then
            label='Extra high'
            return
        end if
        first=iachar(label(1:1))
        if(first>=iachar('a').and.first<=iachar('z'))label(1:1)=achar(first-32)
    end function assistant_effort_label

end module studio_assistant_local_http
