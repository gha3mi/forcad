! Transport for the optional local workspace assistant. Call the worker on a QThread.
module studio_assistant_transport
    use, intrinsic :: iso_c_binding, only: c_int,c_int64_t,c_size_t,c_loc,c_ptr,c_null_char,c_associated
    use studio_assistant_system, only: start_agent_worker,poll_agent_worker,close_agent_worker,seek_agent_file
    use qt, only: QThread_isInterruptionRequested, QThread_msleep, QSaveFile, &
        QDir_removeRecursively, QCoreApplication_applicationFilePath, QDir_entryList, QDir_Dirs, qt_string
    use studio_chat_attachment, only: studio_attachment
    use studio_presentation_text, only: integer_text
    use studio_project_folder, only: canonical_project_file
    use studio_storage_paths, only: studio_work_directory
    use studio_platform_files, only: temporary_private_directory,create_private_file,open_regular,inspect_file, &
        read_file,write_file,close_file,remove_child,current_directory,file_metadata,create_private_directory, &
        inspect_path,is_regular_file,private_directory

    implicit none(type, external)
    private
    public :: assistant_request_t, run_assistant_request
    public :: prepare_assistant_stream, read_assistant_progress, cleanup_assistant_stream, read_agent_events
    public :: read_agent_permission, send_agent_decision
    public :: register_assistant_worker
    public :: assistant_worker_t, assistant_worker_count, read_text, write_text
    public :: valid_provider_session, write_assistant_attachments, read_assistant_attachments
    character(:),allocatable :: worker_executable
    integer, parameter :: assistant_worker_count=2
    integer, parameter :: attachment_name_limit=255, attachment_path_limit=4096
    integer, parameter :: attachment_image_limit=5592408, attachment_text_limit=24576
    character(20), parameter :: worker_files(assistant_worker_count)=[character(20)::'team-worker-1','team-worker-2']
    type :: assistant_worker_t
        character(:), allocatable :: provider, model, effort
    end type assistant_worker_t

    type :: assistant_request_t
        character(:), allocatable :: model, question, context
        character(:), allocatable :: history, mcp_session, agent_model, agent_effort
        character(:), allocatable :: resume_session, provider_session
        character(:), allocatable :: working_directory
        type(studio_attachment),allocatable :: attachments(:)
        character(:), allocatable :: reply, error
        character(:), allocatable :: models, agent_catalog
        character(:), allocatable, private :: stream_directory
        logical :: list_only=.false.
        logical :: attachment_reference_only=.false., team_mode=.false.
        logical :: cancelled=.false.
        type(assistant_worker_t) :: team_workers(assistant_worker_count)
    end type assistant_request_t

contains
    logical function valid_provider_session(value) result(valid)
        character(*), intent(in) :: value
        character(*), parameter :: allowed='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_-'

        valid=.false.
        if(len(value)==0.or.len(value)>256)return
        valid=value(1:1)/='-'.and.verify(value,allowed)==0
    end function valid_provider_session

    subroutine register_assistant_worker()
        worker_executable=QCoreApplication_applicationFilePath()
    end subroutine register_assistant_worker

    subroutine run_assistant_request(data)
        class(*), optional, intent(inout) :: data
        if (.not.present(data)) return
        select type(request => data)
        type is(assistant_request_t)
            call exchange(request)
        end select
    end subroutine run_assistant_request

    ! The GUI prepares the request before starting the worker and releases the
    ! directory only after the worker has stopped. Its path never changes while
    ! the worker is active.
    subroutine prepare_assistant_stream(request)
        type(assistant_request_t), intent(inout) :: request
        call cleanup_assistant_stream(request)
        request%error=''
        request%stream_directory=make_exchange(request)
        if (len(request%stream_directory)==0) then
            deallocate(request%stream_directory)
            request%error='Could not prepare a private workspace for the assistant.'
        end if
    end subroutine prepare_assistant_stream

    function read_assistant_progress(request) result(text)
        type(assistant_request_t), intent(in) :: request
        character(:), allocatable :: text, frames
        integer :: i,n,last
        text=''
        if (.not.allocated(request%stream_directory)) return
        frames=read_text(request%stream_directory//'/progress',32768)
        last=index(frames,achar(0),back=.true.)
        if (last==0) return
        ! Frames contain complete UTF-8 deltas. Ignore a partial last write.
        text=repeat(' ',last)
        n=0
        do i=1,last
            if (frames(i:i)==achar(0)) cycle
            n=n+1
            text(n:n)=frames(i:i)
        end do
        text=text(:n)
    end function read_assistant_progress

    function read_agent_events(request,after) result(text)
        type(assistant_request_t),intent(in) :: request
        integer,intent(in) :: after
        character(:),allocatable :: text
        integer :: complete

        text=''
        if(.not.allocated(request%stream_directory))return
        text=read_text(request%stream_directory//'/agent-events',16000,after)
        complete=index(text,new_line('a'),back=.true.)
        text=text(:complete)
    end function read_agent_events

    subroutine read_agent_permission(request,id,summary)
        type(assistant_request_t),intent(in) :: request
        character(:),allocatable,intent(out) :: id,summary
        character(:),allocatable :: text
        integer :: split

        id=''
        summary=''
        if(.not.allocated(request%stream_directory))return
        text=read_text(request%stream_directory//'/agent-permission',1024)
        split=index(text,new_line('a'))
        if(split<2.or.len(text)<=split)return
        if(text(len(text):)/=new_line('a'))return
        if(verify(text(:split-1),'0123456789')/=0)return
        id=text(:split-1)
        summary=text(split+1:len(text)-1)
    end subroutine read_agent_permission

    logical function send_agent_decision(request,id,allowed) result(ok)
        type(assistant_request_t),intent(in) :: request
        character(*),intent(in) :: id
        logical,intent(in) :: allowed
        type(QSaveFile) :: file
        character(:),allocatable :: payload

        ok=.false.
        if(.not.allocated(request%stream_directory))return
        payload=id//new_line('a')//trim(merge('allow','deny ',allowed))//new_line('a')
        file=QSaveFile(request%stream_directory//'/agent-decision')
        if(file%open())then
            if(file%write(payload)==len(payload))then
                ok=file%commit()
            else
                call file%cancelWriting()
            end if
        end if
        call file%delete()
    end function send_agent_decision

    subroutine cleanup_assistant_stream(request)
        type(assistant_request_t), intent(inout) :: request
        if (.not.allocated(request%stream_directory)) return
        if (len(request%stream_directory)>0) call remove_exchange(request%stream_directory)
        deallocate(request%stream_directory)
    end subroutine cleanup_assistant_stream

    function make_exchange(request) result(directory)
        type(assistant_request_t), intent(in) :: request
        character(:),allocatable :: directory,parent

        if(request%list_only)then
            parent=studio_work_directory()
        else
            parent=studio_work_directory(request%working_directory)
        end if
        directory=temporary_private_directory(parent,'assistant.')
    end function make_exchange

    subroutine exchange(request)
        type(assistant_request_t), intent(inout) :: request
        character(:), allocatable :: directory,session_id
        logical :: successful, own_directory, prepared

        request%reply=''
        request%error=''
        request%models=''
        request%agent_catalog='{}'
        request%provider_session=''
        request%cancelled=QThread_isInterruptionRequested()
        if (request%cancelled) return
        if (.not.allocated(worker_executable)) then
            if (.not.request%list_only) request%error='The Assistant worker is unavailable.'
            return
        end if
        own_directory=.not.allocated(request%stream_directory)
        if (own_directory) then
            directory=make_exchange(request)
        else
            directory=request%stream_directory
        end if
        if (len(directory)==0) then
            request%error='Could not prepare a private workspace for the assistant.'
            return
        end if

        if (request%list_only) then
            prepared=write_text(directory//'/list-only','')
            if (.not.prepared) request%error='Could not prepare the provider catalog request.'
        else
            prepared=write_chat_input(request,directory)
        end if
        if (prepared) then
            call run_worker(worker_executable,directory,successful,request%cancelled)
            if (.not.request%list_only) then
                session_id=read_text(directory//'/provider-session',256)
                if(valid_provider_session(session_id))request%provider_session=session_id
            end if
            ! Stop discards any reply that arrived as cancellation began.
            if (.not.request%cancelled) then
                if (successful) then
                    if (request%list_only) then
                        request%models=read_text(directory//'/models',1048576)
                        request%agent_catalog=read_text(directory//'/agent-catalog.json',4194304)
                    else
                        request%reply=read_text(directory//'/reply',16000)
                        if (len_trim(request%reply)==0) request%error='The assistant did not return an answer.'
                    end if
                else
                    request%error=read_text(directory//'/error',2000)
                    if (len_trim(request%error)==0) request%error='The Assistant worker could not run.'
                end if
            end if
        end if
        if (own_directory) call remove_exchange(directory)
    end subroutine exchange

    logical function write_chat_input(request,directory) result(ok)
        type(assistant_request_t), intent(inout) :: request
        character(*), intent(in) :: directory
        character(:), allocatable :: history,model,session_id
        ok=.false.
        if (.not.allocated(request%question).or..not.allocated(request%context).or..not.allocated(request%model)) then
            request%error='The assistant request is incomplete.'
            return
        end if
        if (len_trim(request%question)==0.or.len(request%question)>8000) then
            request%error='Enter a question of up to 8,000 bytes.'
            return
        end if
        if (len(request%context)>8000) then
            request%error='The workspace context is too large.'
            return
        end if
        session_id=''
        if (allocated(request%resume_session)) session_id=request%resume_session
        history=''
        if(len(session_id)>0)then
            if(.not.valid_provider_session(session_id))then
                request%error='ForCAD Studio could not read the saved provider session.'
                return
            end if
        else if(allocated(request%history))then
            history=request%history
        end if
        model=trim(adjustl(request%model))
        if (len(model)==0.or.len(model)>160) then
            request%error='Choose an available Assistant provider.'
            return
        end if
        ok=write_text(directory//'/model',model).and. &
            write_text(directory//'/question',request%question).and. &
            write_text(directory//'/context',request%context).and. &
            write_text(directory//'/history',history)
        if(ok.and.allocated(request%working_directory)) &
            ok=write_text(directory//'/project-directory',request%working_directory)
        if(ok.and.allocated(request%mcp_session))then
            if(len(request%mcp_session)>0)ok=write_text(directory//'/mcp-session',request%mcp_session)
        end if
        if(ok.and.len(session_id)>0)ok=write_text(directory//'/resume-session',session_id)
        if(ok.and.request%team_mode)ok=write_team_settings(directory,request%team_workers)
        if(ok.and.allocated(request%agent_model))then
            if(len(request%agent_model)>0)ok=write_text(directory//'/agent-model',request%agent_model)
        end if
        if(ok.and.allocated(request%agent_effort))then
            if(len(request%agent_effort)>0)ok=write_text(directory//'/agent-effort',request%agent_effort)
        end if
        if(ok.and.allocated(request%attachments))then
            call write_assistant_attachments(directory,request%attachments,request%working_directory,request%error)
            ok=len(request%error)==0
        end if
        if(ok.and.request%attachment_reference_only)ok=write_text(directory//'/attachment-mode','reference')
        if (.not.ok.and.len(request%error)==0) request%error='Could not prepare the local assistant request.'
    end function write_chat_input

    logical function write_team_settings(directory,workers) result(ok)
        character(*),intent(in) :: directory
        type(assistant_worker_t),intent(in) :: workers(assistant_worker_count)
        character,parameter :: line_break=new_line('a')
        integer :: i

        ok=write_text(directory//'/team-enabled','1')
        do i=1,size(workers)
            if(.not.ok)return
            if(.not.allocated(workers(i)%provider))cycle
            if(len(workers(i)%provider)==0)cycle
            ok=write_text(directory//'/'//trim(worker_files(i)),workers(i)%provider//line_break// &
                workers(i)%model//line_break//workers(i)%effort//line_break)
        end do
    end function write_team_settings

    subroutine write_assistant_attachments(directory,attachments,project,error)
        character(*),intent(in) :: directory
        type(studio_attachment),intent(in) :: attachments(:)
        character(*),intent(in),optional :: project
        character(:),allocatable,intent(out) :: error
        character(:),allocatable :: root,item,path
        integer :: i

        error=''
        if(size(attachments)==0)return
        error='Could not prepare the attached files.'
        root=directory//'/attachments'
        if(.not.create_private_directory(root))return
        do i=1,size(attachments)
            associate(attachment=>attachments(i))
                if(.not.valid_attachment(attachment))return
                item=root//'/'//integer_text(i)
                if(.not.create_private_directory(item))return
                if(.not.write_text(item//'/name',attachment%name))return
                if(allocated(attachment%path))then
                    if(len(attachment%path)>0)then
                        if(.not.present(project))return
                        path=canonical_project_file(project,attachment%path)
                        if(.not.regular_attachment_path(path))then
                            error='An attached file is no longer available inside the project folder.'
                            return
                        end if
                        if(.not.write_text(item//'/path',path))return
                    end if
                end if
                if(.not.write_attachment_field(item//'/image',attachment%image))return
                if(.not.write_attachment_field(item//'/text',attachment%text))return
            end associate
        end do
        if(write_text(directory//'/attachment-count',integer_text(size(attachments))))error=''
    end subroutine write_assistant_attachments

    subroutine read_assistant_attachments(directory,project,attachments,error)
        character(*),intent(in) :: directory,project
        type(studio_attachment),allocatable,intent(out) :: attachments(:)
        character(:),allocatable,intent(out) :: error
        type(qt_string),allocatable :: entries(:)
        character(:),allocatable :: root,count_text,item,path
        integer :: count,i,io

        error='Could not read the attached files.'
        root=directory//'/attachments'
        count_text=read_text(directory//'/attachment-count',24)
        if(len(count_text)==0)then
            if(attachment_entry_exists(root).or.attachment_entry_exists(directory//'/attachment-count'))return
            allocate(attachments(0))
            error=''
            return
        end if
        if(verify(count_text,'0123456789')/=0.or..not.private_directory(root))return
        read(count_text,*,iostat=io)count
        if(io/=0)return
        if(count<1)return
        entries=QDir_entryList(root,QDir_Dirs)
        if(count>size(entries))return
        allocate(attachments(count),stat=io)
        if(io/=0)return
        do i=1,count
            item=root//'/'//integer_text(i)
            if(.not.private_directory(item))return
            associate(attachment=>attachments(i))
                attachment%name=read_text(item//'/name',attachment_name_limit)
                attachment%path=read_text(item//'/path',attachment_path_limit)
                attachment%image=read_text(item//'/image',attachment_image_limit)
                attachment%text=read_text(item//'/text',attachment_text_limit)
                attachment%relative_path=''
                if(.not.valid_attachment(attachment))return
                if(len(attachment%path)>0)then
                    path=canonical_project_file(project,attachment%path)
                    if(path/=attachment%path.or..not.regular_attachment_path(path))then
                        error='An attached file is no longer readable inside the project folder.'
                        return
                    end if
                    attachment%relative_path=path(len(project)+2:)
                end if
            end associate
        end do
        error=''
    end subroutine read_assistant_attachments

    logical function valid_attachment(attachment) result(valid)
        type(studio_attachment),intent(in) :: attachment
        integer :: path_bytes,image_bytes,text_bytes

        valid=.false.
        if(.not.allocated(attachment%name))return
        if(len_trim(attachment%name)==0.or.len(attachment%name)>attachment_name_limit)return
        path_bytes=0
        image_bytes=0
        text_bytes=0
        if(allocated(attachment%path))path_bytes=len(attachment%path)
        if(allocated(attachment%image))image_bytes=len(attachment%image)
        if(allocated(attachment%text))text_bytes=len(attachment%text)
        if(path_bytes>attachment_path_limit.or.image_bytes>attachment_image_limit.or.text_bytes>attachment_text_limit)return
        if(image_bytes>0.and.text_bytes>0)return
        valid=path_bytes+image_bytes+text_bytes>0
    end function valid_attachment

    logical function write_attachment_field(path,value) result(ok)
        character(*),intent(in) :: path
        character(:),allocatable,intent(in) :: value

        ok=.true.
        if(allocated(value))then
            if(len(value)>0)ok=write_text(path,value)
        end if
    end function write_attachment_field

    logical function attachment_entry_exists(path) result(exists)
        character(*),intent(in) :: path
        type(file_metadata) :: metadata

        exists=inspect_path(path,metadata)
    end function attachment_entry_exists

    logical function regular_attachment_path(path) result(valid)
        character(*),intent(in) :: path
        type(file_metadata) :: metadata

        valid=.false.
        if(len(path)==0.or.len(path)>attachment_path_limit)return
        if(inspect_path(path,metadata))valid=is_regular_file(metadata)
    end function regular_attachment_path

    subroutine run_worker(program,directory,successful,cancelled)
        character(*),intent(in) :: program,directory
        logical,intent(out) :: successful,cancelled
        type(c_ptr) :: worker
        integer(c_int) :: status,completed
        integer :: retry

        successful=.false.
        cancelled=QThread_isInterruptionRequested()
        if(cancelled)return
        worker=start_agent_worker(program//c_null_char,directory//c_null_char)
        if(.not.c_associated(worker))return
        do
            if(QThread_isInterruptionRequested())exit
            status=poll_agent_worker(worker,completed)
            if(status/=0)then
                successful=status>0.and.completed==1
                exit
            end if
            call QThread_msleep(25)
        end do
        cancelled=QThread_isInterruptionRequested()
        if(cancelled)then
            successful=.false.
            if(write_text(directory//'/stop-requested',''))then
                do retry=1,40
                    status=poll_agent_worker(worker,completed)
                    if(status/=0)exit
                    call QThread_msleep(25)
                end do
            end if
        end if
        call close_agent_worker(worker)
    end subroutine run_worker

    logical function write_text(path,value) result(ok)
        character(*),intent(in) :: path
        character(*),target,intent(in) :: value
        integer(c_int) :: output,status
        integer(c_int64_t) :: written
        integer :: offset

        ok=.false.
        output=create_private_file(current_directory,path//c_null_char)
        if(output<0)return
        offset=0
        do while(offset<len(value))
            written=write_file(output,c_loc(value(offset+1:offset+1)),int(len(value)-offset,c_size_t))
            if(written<=0)exit
            offset=offset+int(written)
        end do
        status=close_file(output)
        ok=status==0.and.offset==len(value)
        if(.not.ok)status=remove_child(current_directory,path//c_null_char,0_c_int)
    end function write_text

    function read_text(path,limit,after,readable) result(value)
        character(*),intent(in) :: path
        integer,intent(in),optional :: limit,after
        logical,intent(out),optional :: readable
        character(:),allocatable :: value
        character(:),allocatable,target :: buffer
        type(file_metadata) :: metadata
        integer(c_int) :: input,status
        integer(c_int64_t) :: offset,length,received,maximum
        integer :: cursor,allocation
        logical :: complete

        value=''
        if(present(readable))readable=.false.
        maximum=huge(0)-1
        if(present(limit))maximum=limit
        if(maximum<=0)return
        offset=0
        if(present(after))offset=max(0,after)
        input=open_regular(current_directory,path//c_null_char)
        if(input<0)return
        status=inspect_file(input,metadata)
        complete=.false.
        if(status==0)then
            length=metadata%bytes-offset
            if(present(after))length=min(length,maximum)
            complete=length==0
            if(length>0.and.length<=maximum)then
                if(offset>0)status=seek_agent_file(input,offset)
                if(status==0)then
                    allocate(character(int(length))::buffer,stat=allocation)
                    if(allocation==0)then
                        cursor=0
                        do while(cursor<len(buffer))
                            received=read_file(input,c_loc(buffer(cursor+1:cursor+1)),int(len(buffer)-cursor,c_size_t))
                            if(received<=0)exit
                            cursor=cursor+int(received)
                        end do
                        complete=cursor==len(buffer)
                    end if
                end if
            end if
        end if
        status=close_file(input)
        if(status/=0.or..not.complete)return
        if(allocated(buffer))call move_alloc(buffer,value)
        if(present(readable))readable=.true.
    end function read_text

    subroutine remove_exchange(directory)
        character(*), intent(in) :: directory
        logical :: removed

        removed=QDir_removeRecursively(directory)
    end subroutine remove_exchange
end module studio_assistant_transport
