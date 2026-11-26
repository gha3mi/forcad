module studio_mcp_mailbox
    use, intrinsic :: iso_c_binding, only: c_char,c_int,c_int64_t,c_null_char,c_loc,c_size_t
    use studio_platform_files, only: file_metadata,current_directory,open_directory,open_regular,close_file, &
        inspect_file,inspect_entry,is_regular_file,same_file,create_private_file,rename_child,remove_child, &
        read_file,write_file,make_private_child
    use studio_project_folder, only: canonical_project_folder
    implicit none(type, external)
    private
    public :: mcp_mailbox, submit_mcp_request, read_mcp_response, close_mcp_mailbox
    public :: mcp_mailbox_exists, mcp_request_waiting, mcp_request_claimed, retract_mcp_request
    public :: mcp_response_pending, mcp_response_valid, mcp_response_invalid
    public :: max_mcp_request_bytes

    integer, parameter :: max_mcp_request_bytes=4*1024*1024
    integer, parameter :: max_response=262144
    integer, parameter :: mcp_response_pending=0, mcp_response_valid=1, mcp_response_invalid=2
    character(*), parameter :: request_file='request.json', temporary_file='.request.tmp', running_file='running.json', &
        approval_file='awaiting-approval', response_file='response.json'

    type :: mcp_mailbox
        integer(c_int) :: root_fd=-1, directory_fd=-1
        type(file_metadata) :: identity=file_metadata(0,0,0,0,0_c_int64_t,0_c_int64_t,0_c_int64_t)
        character(:), allocatable :: name
    end type mcp_mailbox

contains

    function submit_mcp_request(mailbox,session,payload) result(problem)
        type(mcp_mailbox), intent(out) :: mailbox
        character(*), intent(in) :: session,payload
        character(:), allocatable :: problem,canonical
        type(file_metadata) :: metadata
        integer(c_int) :: status,file
        logical :: written

        if(len(payload)>max_mcp_request_bytes)then
            problem='The ForCAD Studio tool request exceeds 4 MiB. Break the edit into smaller requests.'
            return
        end if
        problem='ForCAD Studio MCP is unavailable. Send a new Assistant message with tools enabled, or reconnect an external client to a live session.'
        if(len(session)<=1)return
        if(index(session,c_null_char)>0)return
        canonical=canonical_project_folder(session)
        if(len(canonical)/=len(session).or.canonical/=session)return

        mailbox%root_fd=open_directory(current_directory,session//c_null_char)
        if(mailbox%root_fd<0)return
        status=inspect_file(mailbox%root_fd,metadata)
        if(status/=0.or.metadata%is_private/=1)then
            call close_mcp_mailbox(mailbox)
            return
        end if

        problem='Could not open a ForCAD Studio mailbox. Check whether the session is still active.'
        mailbox%name=make_private_child(mailbox%root_fd,'mcp-')
        if(len(mailbox%name)==0)then
            call close_mcp_mailbox(mailbox)
            return
        end if
        mailbox%directory_fd=open_directory(mailbox%root_fd,mailbox%name//c_null_char)
        status=inspect_file(mailbox%directory_fd,mailbox%identity)
        if(status/=0.or.mailbox%identity%is_private/=1)then
            status=remove_child(mailbox%root_fd,mailbox%name//c_null_char,1_c_int)
            call close_mcp_mailbox(mailbox)
            return
        end if

        problem='Could not submit the request to ForCAD Studio.'
        file=create_private_file(mailbox%directory_fd,temporary_file//c_null_char)
        if(file<0)then
            call close_mcp_mailbox(mailbox)
            return
        end if
        written=write_payload(file,payload)
        status=close_file(file)
        if(written.and.status==0)then
            status=rename_child(mailbox%directory_fd,temporary_file//c_null_char,mailbox%directory_fd,request_file//c_null_char,0_c_int)
            if(status==0)then
                problem=''
                return
            end if
        end if
        call close_mcp_mailbox(mailbox)
    end function submit_mcp_request

    logical function write_payload(descriptor,payload) result(written)
        integer(c_int), intent(in) :: descriptor
        character(kind=c_char,len=*), intent(in), target :: payload
        integer(c_size_t) :: offset,length
        integer(c_int64_t) :: count

        offset=0
        length=len(payload,kind=c_size_t)
        do while(offset<length)
            count=write_file(descriptor,c_loc(payload(offset+1:offset+1)),length-offset)
            if(count<=0)exit
            offset=offset+count
        end do
        written=offset==length
    end function write_payload

    logical function mcp_mailbox_exists(mailbox) result(exists)
        type(mcp_mailbox), intent(in) :: mailbox
        type(file_metadata) :: metadata

        exists=.false.
        if(mailbox%directory_fd<0)return
        if(inspect_entry(mailbox%root_fd,mailbox%name//c_null_char,metadata)/=0)return
        exists=metadata%is_private==1.and.same_file(metadata,mailbox%identity)
    end function mcp_mailbox_exists

    logical function mcp_request_waiting(mailbox) result(waiting)
        type(mcp_mailbox), intent(in) :: mailbox
        type(file_metadata) :: metadata

        waiting=.false.
        if(mailbox%directory_fd<0)return
        if(inspect_entry(mailbox%directory_fd,approval_file//c_null_char,metadata)/=0)return
        waiting=is_regular_file(metadata).and.metadata%is_private==1.and.metadata%links==1
    end function mcp_request_waiting

    logical function file_present(directory,name) result(present)
        integer(c_int), intent(in) :: directory
        character(*), intent(in) :: name
        type(file_metadata) :: metadata

        present=inspect_entry(directory,name//c_null_char,metadata)==0
    end function file_present

    logical function mcp_request_claimed(mailbox) result(claimed)
        type(mcp_mailbox), intent(in) :: mailbox

        claimed=file_present(mailbox%directory_fd,running_file).or.file_present(mailbox%directory_fd,response_file)
    end function mcp_request_claimed

    function read_mcp_response(mailbox,response_status) result(response)
        type(mcp_mailbox), intent(in) :: mailbox
        integer, intent(out) :: response_status
        character(kind=c_char,len=:), allocatable, target :: response
        type(file_metadata) :: metadata
        integer(c_int) :: file,status
        integer(c_size_t) :: length,offset
        integer(c_int64_t) :: count

        response_status=mcp_response_pending
        if(.not.file_present(mailbox%directory_fd,response_file))return
        response_status=mcp_response_invalid
        file=open_regular(mailbox%directory_fd,response_file//c_null_char)
        if(file<0)return
        status=inspect_file(file,metadata)
        if(status==0)then
            if(is_regular_file(metadata).and.metadata%is_private==1.and.metadata%links==1.and. &
                metadata%bytes>0.and.metadata%bytes<=max_response)then
                length=int(metadata%bytes,c_size_t)
                allocate(character(kind=c_char,len=length)::response)
                offset=0
                do while(offset<length)
                    count=read_file(file,c_loc(response(offset+1:offset+1)),length-offset)
                    if(count<=0)exit
                    offset=offset+count
                end do
                if(offset==length)response_status=mcp_response_valid
            end if
        end if
        status=close_file(file)
    end function read_mcp_response

    subroutine retract_mcp_request(mailbox)
        type(mcp_mailbox), intent(in) :: mailbox
        integer(c_int) :: status

        if(mailbox%directory_fd>=0)status=remove_child(mailbox%directory_fd,request_file//c_null_char,0_c_int)
    end subroutine retract_mcp_request

    subroutine close_mcp_mailbox(mailbox)
        type(mcp_mailbox), intent(inout) :: mailbox
        integer(c_int) :: status

        if(mailbox%directory_fd>=0)then
            call retract_mcp_request(mailbox)
            status=remove_child(mailbox%directory_fd,temporary_file//c_null_char,0_c_int)
            if(file_present(mailbox%directory_fd,response_file).or..not.file_present(mailbox%directory_fd,running_file))then
                status=remove_child(mailbox%directory_fd,running_file//c_null_char,0_c_int)
                status=remove_child(mailbox%directory_fd,response_file//c_null_char,0_c_int)
                status=remove_child(mailbox%directory_fd,approval_file//c_null_char,0_c_int)
                if(mcp_mailbox_exists(mailbox))status=remove_child(mailbox%root_fd,mailbox%name//c_null_char,1_c_int)
            end if
            status=close_file(mailbox%directory_fd)
        end if
        if(mailbox%root_fd>=0)status=close_file(mailbox%root_fd)
        mailbox=mcp_mailbox()
    end subroutine close_mcp_mailbox

end module studio_mcp_mailbox
