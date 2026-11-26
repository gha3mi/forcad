! Local storage for saved assistant conversations.
module studio_chat_history
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_size_t, c_null_char
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QSettings, QSettings_NoError, QSettings_FormatError, qt_string
    use studio_assistant_system, only: epoch_microseconds
    use studio_platform_files, only: create_private_file,open_regular,close_file,inspect_file,inspect_path,file_metadata,current_directory,restrict_file
    use studio_storage_paths, only: studio_project_directory
    implicit none(type, external)
    private
    public :: chat_record_t, chat_summary_t, project_chat_file, save_chat, load_chat, list_chats, delete_chat

    character(*), parameter :: group='chats-v1', payload_version='v3:', session_payload_version='v2:'
    character(*), parameter :: digits='0123456789', hex_digits=digits//'abcdef'
    integer, parameter :: id_size=32, date_size=19, order_size=20, length_size=10, legacy_length_size=7
    integer, parameter :: title_limit=1024, model_limit=160
    integer, parameter :: agent_model_limit=512, agent_effort_limit=80

    type :: chat_record_t
        character(:), allocatable :: id,title,transcript,history,model,updated
        character(:),allocatable :: agent_model,agent_effort
    end type chat_record_t

    type :: chat_summary_t
        character(:), allocatable :: id,title,updated
        integer(int64), private :: order=0_int64
    end type chat_summary_t

    interface
        integer(c_int) function random_chat_bytes(bytes,count) bind(C,name='studio_random_fill')
            import c_int,c_char,c_size_t
            implicit none(type,external)
            character(kind=c_char),intent(out) :: bytes(*)
            integer(c_size_t),value :: count
        end function random_chat_bytes
    end interface
contains
    subroutine save_chat(record,error,file_name)
        type(chat_record_t), intent(inout) :: record
        character(:), allocatable, intent(out) :: error
        character(*), intent(in) :: file_name
        type(QSettings) :: settings
        character(:), allocatable :: payload
        character(date_size) :: display_time
        character(order_size) :: order_text
        integer(int64) :: order

        error='The chat has missing, invalid, or oversized content.'
        if(.not.allocated(record%agent_model))record%agent_model=''
        if(.not.allocated(record%agent_effort))record%agent_effort=''
        if (.not.valid_content(record)) return
        if (.not.allocated(record%id)) record%id=''
        if (len(record%id)>0.and..not.valid_id(record%id)) return
        error='Could not open the private chat history file for writing.'
        if(.not.secure_file(file_name))return
        settings=QSettings(file_name)
        call settings%beginGroup(group)
        call settings%sync()
        error=store_error(settings)
        if (len(error)==0) then
            if (len(record%id)==0) record%id=new_id(settings)
            if (len(record%id)==0) then
                error='Could not assign a unique identifier to the chat.'
            else if (.not.current_time(display_time,order)) then
                error='Could not read the time for this chat.'
            else
                write(order_text,'(i20.20)') order
                payload=payload_version//order_text//display_time//encode_field(record%title)//encode_field(record%model)// &
                    encode_field(record%history)//encode_field(record%transcript)// &
                    encode_field(record%agent_model)//encode_field(record%agent_effort)
                call settings%setValue(record%id,payload)
                call settings%sync()
                if (settings%status()==QSettings_NoError) then
                    record%updated=display_time
                else
                    error='Could not save the chat history.'
                end if
            end if
        end if
        call settings%delete()
    end subroutine save_chat

    subroutine load_chat(id,record,error,file_name)
        character(*), intent(in) :: id
        type(chat_record_t), intent(out) :: record
        character(:), allocatable, intent(out) :: error
        character(*), intent(in) :: file_name
        type(QSettings) :: settings
        type(chat_record_t) :: decoded
        integer(int64) :: order
        logical :: exists

        error='Invalid chat identifier.'
        if (.not.valid_id(id)) return
        error=check_chat_file(file_name,exists)
        if(len(error)>0)return
        if(.not.exists)then
            error='The saved chat could not be found.'
            return
        end if
        settings=QSettings(file_name)
        call settings%beginGroup(group)
        call settings%sync()
        error=store_error(settings)
        if (len(error)==0) then
            if (.not.settings%contains(id)) then
                error='The saved chat could not be found.'
            else if (unpack_record(settings%stringValue(id),decoded,order)) then
                decoded%id=id
                record=decoded
            else
                error='The saved chat could not be read.'
            end if
            if (settings%status()/=QSettings_NoError) error=store_error(settings)
        end if
        call settings%delete()
    end subroutine load_chat

    subroutine delete_chat(id,error,file_name)
        character(*), intent(in) :: id
        character(:), allocatable, intent(out) :: error
        character(*), intent(in) :: file_name
        type(QSettings) :: settings
        logical :: exists

        error='Invalid chat identifier.'
        if (.not.valid_id(id)) return
        error=check_chat_file(file_name,exists)
        if(len(error)>0.or..not.exists)return
        settings=QSettings(file_name)
        call settings%beginGroup(group)
        call settings%sync()
        error=store_error(settings)
        if (len(error)==0.and.settings%contains(id)) then
            call settings%remove(id)
            call settings%sync()
            if (settings%status()/=QSettings_NoError) error='Could not delete the saved chat.'
        end if
        call settings%delete()
    end subroutine delete_chat

    subroutine list_chats(items,error,file_name)
        type(chat_summary_t), allocatable, intent(out) :: items(:)
        character(:), allocatable, intent(out) :: error
        character(*), intent(in) :: file_name
        type(QSettings) :: settings
        type(qt_string), allocatable :: keys(:)
        type(chat_summary_t), allocatable :: found(:)
        type(chat_record_t) :: record
        integer(int64) :: order
        integer :: i,n,position
        logical :: exists

        error=check_chat_file(file_name,exists)
        if(len(error)>0.or..not.exists)then
            allocate(items(0))
            return
        end if
        settings=QSettings(file_name)
        call settings%beginGroup(group)
        call settings%sync()
        error=store_error(settings)
        n=0
        if (len(error)==0) then
            keys=settings%childKeys()
            allocate(found(size(keys)))
            do i=1,size(keys)
                if (.not.valid_id(keys(i)%s)) cycle
                if (.not.unpack_record(settings%stringValue(keys(i)%s),record,order)) then
                    error='Some saved chats could not be read.'
                    cycle
                end if
                position=n+1
                do while (position>1)
                    if (found(position-1)%order>=order) exit
                    found(position)=found(position-1)
                    position=position-1
                end do
                found(position)=chat_summary_t(keys(i)%s,record%title,record%updated,order)
                n=n+1
            end do
            if (settings%status()/=QSettings_NoError) then
                error=store_error(settings)
                n=0
            end if
        end if
        call settings%delete()
        if (n==0) then
            allocate(items(0))
        else
            items=found(:n)
        end if
    end subroutine list_chats

    function project_chat_file(project) result(path)
        character(*),intent(in) :: project
        character(:),allocatable :: path

        path=studio_project_directory(project)
        if(len(path)>0)path=path//'/chats.ini'
    end function project_chat_file

    function check_chat_file(path,exists) result(error)
        character(*),intent(in) :: path
        logical,intent(out) :: exists
        character(:),allocatable :: error
        type(file_metadata) :: metadata
        integer(c_int) :: fd,inspected,closed

        error='Could not open the private chat history file.'
        exists=.false.
        if(len_trim(path)==0.or.index(path,c_null_char)>0)return
        exists=inspect_path(path,metadata)
        if(exists)then
            fd=open_regular(current_directory,path//c_null_char)
            if(fd<0)return
            inspected=inspect_file(fd,metadata)
            closed=close_file(fd)
            if(inspected/=0.or.closed/=0)return
            if(metadata%owned/=1.or.metadata%links/=1)return
            if(metadata%is_private/=1)then
                if(restrict_file(path//c_null_char)/=0)return
            end if
        end if
        error=''
    end function check_chat_file

    logical function secure_file(path) result(ok)
        character(*), intent(in) :: path
        character(:),allocatable :: error
        integer(c_int) :: fd
        logical :: exists

        ok=.false.
        error=check_chat_file(path,exists)
        if(len(error)>0)return
        if(exists)then
            ok=.true.
            return
        end if
        fd=create_private_file(current_directory,path//c_null_char)
        if(fd<0)return
        ok=close_file(fd)==0
    end function secure_file

    function store_error(settings) result(error)
        type(QSettings), intent(in) :: settings
        character(:), allocatable :: error
        select case (settings%status())
        case (QSettings_NoError)
            error=''
        case (QSettings_FormatError)
            error='The chat history file has an invalid format.'
        case default
            error='Could not open the chat history file.'
        end select
    end function store_error

    function new_id(settings) result(id)
        type(QSettings), intent(in) :: settings
        character(:), allocatable :: id
        character(id_size) :: candidate
        character(kind=c_char,len=id_size/2) :: bytes
        integer :: attempt,i,byte,high,low
        id=''
        do attempt=1,4
            if(random_chat_bytes(bytes,len(bytes,c_size_t))/=0)return
            do i=1,len(bytes)
                byte=iachar(bytes(i:i))
                high=byte/16+1
                low=mod(byte,16)+1
                candidate(2*i-1:2*i)=hex_digits(high:high)//hex_digits(low:low)
            end do
            if (settings%contains(candidate)) cycle
            id=candidate
            exit
        end do
    end function new_id

    logical function current_time(display_time,order) result(ok)
        character(date_size), intent(out) :: display_time
        integer(int64), intent(out) :: order
        integer :: clock(8)
        order=epoch_microseconds()
        call date_and_time(values=clock)
        write(display_time,'(i4.4,2(a,i2.2),a,i2.2,2(a,i2.2))') &
            clock(1),'-',clock(2),'-',clock(3),' ',clock(5),':',clock(6),':',clock(7)
        ok=order>0.and.valid_date(display_time)
    end function current_time

    logical function valid_id(id) result(ok)
        character(*), intent(in) :: id
        ok=len(id)==id_size.and.verify(id,hex_digits)==0
    end function valid_id

    logical function valid_date(value) result(ok)
        character(date_size), intent(in) :: value
        ok=value(5:5)=='-'.and.value(8:8)=='-'.and.value(11:11)==' '.and. &
            value(14:14)==':'.and.value(17:17)==':'.and. &
            verify(value(:4)//value(6:7)//value(9:10)//value(12:13)//value(15:16)//value(18:19),digits)==0.and. &
            value(6:7)>='01'.and.value(6:7)<='12'.and.value(9:10)>='01'.and.value(9:10)<='31'.and. &
            value(12:13)<='23'.and.value(15:16)<='59'.and.value(18:19)<='60'
    end function valid_date

    logical function valid_content(record) result(ok)
        type(chat_record_t), intent(in) :: record
        ok=.false.
        if (.not.allocated(record%title).or..not.allocated(record%transcript).or. &
            .not.allocated(record%history).or..not.allocated(record%model)) return
        if (len_trim(record%title)==0.or.len(record%title)>title_limit.or. &
            scan(record%title,achar(0)//achar(10)//achar(13))/=0) return
        if (len_trim(record%transcript)==0.or.index(record%transcript,achar(0))/=0) return
        if (len(record%model)>model_limit.or.scan(record%model,achar(0)//achar(10)//achar(13))/=0) return
        if(len(record%agent_model)>agent_model_limit.or.scan(record%agent_model,achar(0)//achar(10)//achar(13))/=0)return
        if(len(record%agent_effort)>agent_effort_limit.or.scan(record%agent_effort,achar(0)//achar(10)//achar(13))/=0)return
        ok=valid_history(record%history)
    end function valid_content

    logical function valid_history(history) result(ok)
        character(*), intent(in) :: history
        integer :: i,start,parts
        ok=.false.
        if (len(history)==0) then
            ok=.true.
            return
        end if
        start=1
        parts=1
        do i=1,len(history)
            if (history(i:i)/=achar(0)) cycle
            if (verify(history(start:i-1),' '//achar(9)//achar(10)//achar(13))==0) return
            start=i+1
            parts=parts+1
        end do
        ok=mod(parts,2)==0.and.verify(history(start:),' '//achar(9)//achar(10)//achar(13))/=0
    end function valid_history

    ! The Qt string binding passes explicit UTF-8 byte lengths both ways, so
    ! history's embedded NUL delimiters survive along with Unicode and newlines.
    ! A complete chat lives in one settings value with a size prefix per field.
    function encode_field(value) result(field)
        character(*), intent(in) :: value
        character(:), allocatable :: field
        character(length_size) :: length_text
        write(length_text,'(i10.10)') len(value)
        field=length_text//value
    end function encode_field

    logical function unpack_record(payload,record,order) result(ok)
        character(*), intent(in) :: payload
        type(chat_record_t), intent(out) :: record
        integer(int64), intent(out) :: order
        integer(int64) :: position,last
        integer :: width,io
        logical :: modern_format,session_format

        ok=.false.
        session_format=index(payload,session_payload_version)==1
        modern_format=session_format.or.index(payload,payload_version)==1
        width=legacy_length_size
        position=1
        if(modern_format)then
            width=length_size
            position=position+len(payload_version)
        end if
        if(len(payload,kind=int64)-position+1<order_size+date_size+4*width)return
        last=position+order_size-1
        if(verify(payload(position:last),digits)/=0)return
        read(payload(position:last),'(i20)',iostat=io)order
        if(io/=0)return
        if(order<=0)return
        position=last+date_size+1
        record%updated=payload(last+1:position-1)
        if(.not.valid_date(record%updated))return
        if(.not.take_field(payload,position,width,record%title,title_limit))return
        if(.not.take_field(payload,position,width,record%model,model_limit))return
        if(.not.take_field(payload,position,width,record%history))return
        if(.not.take_field(payload,position,width,record%transcript))return
        record%agent_model=''
        record%agent_effort=''
        if(modern_format.or.position<=len(payload,kind=int64))then
            if(.not.take_field(payload,position,width,record%agent_model,agent_model_limit))return
            if(.not.take_field(payload,position,width,record%agent_effort,agent_effort_limit))return
        end if
        if(session_format)then
            if(.not.take_field(payload,position,width))return
            if(.not.take_field(payload,position,width))return
        end if
        if(position/=len(payload,kind=int64)+1)return
        ok=valid_content(record)
    end function unpack_record

    logical function take_field(payload,position,width,value,limit) result(ok)
        character(*), intent(in) :: payload
        integer(int64), intent(inout) :: position
        integer, intent(in) :: width
        character(:), allocatable, optional, intent(out) :: value
        integer, optional, intent(in) :: limit
        integer(int64) :: last,length
        integer :: io

        ok=.false.
        last=position+width-1
        if(last>len(payload,kind=int64))return
        if(verify(payload(position:last),digits)/=0)return
        read(payload(position:last),'(i10)',iostat=io)length
        if(io/=0)return
        if(length>len(payload,kind=int64)-last)return
        if(present(limit))then
            if(length>limit)return
        end if
        if(present(value))value=payload(last+1:last+length)
        position=last+length+1
        ok=.true.
    end function take_field
end module studio_chat_history
