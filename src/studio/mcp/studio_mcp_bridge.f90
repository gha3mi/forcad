! Dispatch local MCP mailbox requests on ForCAD Studio's Qt event thread.
module studio_mcp_bridge
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: iso_c_binding, only: c_char,c_int,c_int64_t,c_size_t,c_loc,c_null_char
    use qt, only: QApplication_activeModalWidget, QDir_Dirs, QDir_entryList, QDir_NoDotAndDotDot, &
        QJsonDocument, QJsonDocument_fromJson, QJsonObject, QJsonValue_Object, QJsonValue_String, qt_string, &
        QDir_removeRecursively, QTimer, QWidget
    use studio_ai_cursor, only: set_studio_ai_tool_cursor, hide_studio_ai_cursor
    use studio_state, only: studio_t
    use studio_storage_paths, only: studio_work_directory
    use studio_mcp_refinement, only: refinement_approval_summary
    use studio_mcp_controls, only: control_approval_summary
    use studio_mcp_mailbox, only: max_mcp_request_bytes
    use studio_edit_history, only: document_undo_stamp
    use studio_knot_interaction, only: knot_edit_in_progress
    use studio_live_viewport, only: transform_dragging
    use studio_selection, only: selection_set_t
    use studio_exchange, only: studio_exchange_suffix
    use studio_platform_files, only: file_metadata,current_directory,open_directory,open_regular,close_file, &
        inspect_file,inspect_entry,inspect_path,is_regular_file,same_file,create_private_file,rename_child,remove_child, &
        read_file,write_file,temporary_private_directory
    implicit none(type, external)
    private
    public :: bind_studio_mcp, start_studio_mcp, stop_studio_mcp, finish_studio_mcp
    public :: pending_studio_mcp_approval, decide_studio_mcp_approval, rebind_assistant_mcp_project

    integer,parameter :: calls_per_tick=2
    integer,parameter :: approval_waiting=0, approval_allowed=1, approval_denied=2, approval_seconds=120
    character(*),parameter :: request_file='request.json', running_file='running.json', &
        response_file='response.json', temporary_file='.response.tmp', approval_file='awaiting-approval'
    abstract interface
        subroutine tool_executor(s,name,args,result)
            import studio_t,QJsonObject
            implicit none(type, external)
            type(studio_t),intent(inout),target :: s
            character(*),intent(in) :: name
            type(QJsonObject),intent(in) :: args
            type(QJsonObject),intent(out) :: result
        end subroutine tool_executor

        subroutine tool_start_observer(name,args)
            import QJsonObject
            implicit none(type, external)
            character(*),intent(in) :: name
            type(QJsonObject),intent(in) :: args
        end subroutine tool_start_observer

        subroutine tool_observer(name,result,args,from_assistant,committed)
            import QJsonObject
            implicit none(type, external)
            character(*),intent(in) :: name
            type(QJsonObject),intent(in) :: result
            type(QJsonObject),intent(in) :: args
            logical,intent(in) :: from_assistant,committed
        end subroutine tool_observer
    end interface

    type(studio_t),pointer :: workspace=>null()
    procedure(tool_executor),pointer :: execute=>null()
    procedure(tool_observer),pointer :: observe=>null()
    procedure(tool_start_observer),pointer :: observe_start=>null()
    type(QTimer) :: timer
    type :: mailbox_session
        character(:),allocatable :: path, project_directory
        integer(int64) :: document_revision=-1
        integer(c_int) :: root_fd=-1
    end type mailbox_session

    type :: studio_approval
        integer(c_int) :: directory=-1
        integer :: decision=approval_waiting
        integer(int64) :: started=0, stamp=0, ticket=0, document_revision=-1
        integer :: units(2)=0, focus(6)=0
        logical :: selection_dependent=.false.
        type(selection_set_t) :: selection
        character(:),allocatable :: summary,payload,project_path
    end type studio_approval

    type(mailbox_session) :: external_session,assistant_session
    type(file_metadata) :: external_project_identity
    type(studio_approval) :: pending
    integer(int64) :: next_ticket=0
    logical :: dispatching=.false.,assistant_requires_approval=.false.,assistant_switching_project=.false.

contains

    subroutine bind_studio_mcp(s,executor,observer,started)
        type(studio_t),intent(inout),target :: s
        procedure(tool_executor) :: executor
        procedure(tool_observer) :: observer
        procedure(tool_start_observer),optional :: started

        workspace=>s
        execute=>executor
        observe=>observer
        nullify(observe_start)
        if(present(started))observe_start=>started
        timer=QTimer(s%window)
        call timer%timeout%connect(poll_mailbox)
    end subroutine bind_studio_mcp

    function start_studio_mcp(for_assistant,require_approval) result(path)
        logical,intent(in),optional :: for_assistant,require_approval
        character(:),allocatable :: path

        path=''
        if(.not.associated(workspace))return
        if(present(for_assistant))then
            if(for_assistant)then
                if(workspace%home_active.or..not.allocated(workspace%document%project_directory))return
                if(.not.assistant_session_current())call stop_studio_mcp(.true.)
                path=open_session(assistant_session,workspace%document%project_directory)
                if(len(path)>0)then
                    assistant_session%document_revision=workspace%document%document_revision
                    assistant_session%project_directory=workspace%document%project_directory
                    assistant_requires_approval=.false.
                    if(present(require_approval))assistant_requires_approval=require_approval
                end if
                return
            end if
        end if
        if(allocated(external_session%path))then
            path=external_session%path
            return
        end if
        path=open_session(external_session,workspace%document%project_directory)
        if(len(path)>0.and.allocated(workspace%document%project_directory))then
            if(inspect_path(workspace%document%project_directory,external_project_identity))then
                external_session%project_directory=workspace%document%project_directory
            end if
        end if
    end function start_studio_mcp

    function open_session(session,project) result(path)
        type(mailbox_session),intent(inout) :: session
        character(*),intent(in),optional :: project
        character(:),allocatable :: path
        type(file_metadata) :: metadata
        integer(c_int) :: status

        if(allocated(session%path))then
            path=session%path
            return
        end if
        path=temporary_private_directory(studio_work_directory(project),'mcp.')
        if(len(path)==0)return
        session%path=path
        session%root_fd=open_directory(current_directory,path//c_null_char)
        status=inspect_file(session%root_fd,metadata)
        if(status/=0)then
            call close_session(session)
            path=''
            return
        end if
        if(metadata%is_private/=1)then
            call close_session(session)
            path=''
            return
        end if
        if(.not.timer%isActive())call timer%start(60)
    end function open_session

    subroutine stop_studio_mcp(for_assistant)
        logical,intent(in),optional :: for_assistant
        logical :: embedded

        embedded=.false.
        if(present(for_assistant))embedded=for_assistant
        if(embedded)then
            assistant_requires_approval=.false.
            call clear_pending_approval()
            call close_session(assistant_session)
        else
            call close_session(external_session)
            if(associated(workspace).and.assistant_session%root_fd<0)call hide_studio_ai_cursor(workspace)
        end if
        if(external_session%root_fd<0.and.assistant_session%root_fd<0)then
            if(.not.timer%isNull())call timer%stop()
        end if
    end subroutine stop_studio_mcp

    subroutine close_session(session)
        type(mailbox_session),intent(inout) :: session
        integer(c_int) :: status
        logical :: removed

        if(session%root_fd>=0)then
            if(session_path_current(session))removed=QDir_removeRecursively(session%path)
            status=close_file(session%root_fd)
        end if
        session%root_fd=-1
        if(allocated(session%path))deallocate(session%path)
        if(allocated(session%project_directory))deallocate(session%project_directory)
        session%document_revision=-1
    end subroutine close_session

    logical function session_path_current(session) result(current)
        type(mailbox_session),intent(in) :: session
        type(file_metadata) :: original,entry

        current=.false.
        if(.not.allocated(session%path))return
        if(inspect_file(session%root_fd,original)/=0)return
        if(.not.inspect_path(session%path,entry))return
        current=same_file(original,entry).and.original%is_private==1
    end function session_path_current

    logical function assistant_session_current() result(current)
        current=.false.
        if(workspace%home_active.or.assistant_session%root_fd<0)return
        if(workspace%document%document_revision/=assistant_session%document_revision)return
        if(.not.allocated(workspace%document%project_directory).or..not.allocated(assistant_session%project_directory))return
        current=workspace%document%project_directory==assistant_session%project_directory
    end function assistant_session_current

    logical function external_project_current() result(current)
        type(file_metadata) :: identity

        current=.false.
        if(.not.allocated(external_session%project_directory).or..not.allocated(workspace%document%project_directory))return
        if(workspace%document%project_directory/=external_session%project_directory)return
        if(.not.inspect_path(workspace%document%project_directory,identity))return
        current=same_file(identity,external_project_identity)
    end function external_project_current

    subroutine finish_studio_mcp()
        call stop_studio_mcp(.true.)
        call stop_studio_mcp()
        if(.not.timer%isNull())call timer%delete()
        nullify(workspace,execute,observe,observe_start)
    end subroutine finish_studio_mcp

    logical function rebind_assistant_mcp_project(s) result(rebound)
        type(studio_t),intent(in),target :: s

        rebound=.false.
        if(.not.assistant_switching_project.or..not.associated(workspace,s).or.assistant_session%root_fd<0)return
        if(.not.allocated(s%document%project_directory).or..not.allocated(assistant_session%project_directory))return
        if(s%document%project_directory/=assistant_session%project_directory)return
        assistant_session%document_revision=s%document%document_revision
        rebound=.true.
    end function rebind_assistant_mcp_project

    function current_project_path() result(path)
        character(:),allocatable :: path

        path=''
        if(allocated(workspace%document%project_path))path=workspace%document%project_path
    end function current_project_path

    subroutine run_studio_mcp_tool(name,args,result,from_assistant)
        character(*),intent(in) :: name
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(out) :: result
        logical,intent(in) :: from_assistant
        logical :: external_cursor

        if(.not.associated(execute).or..not.associated(workspace))then
            result=failed('ForCAD Studio is still starting.')
        else if(dispatching)then
            result=failed('ForCAD Studio is running another command.')
        else
            dispatching=.true.
            assistant_switching_project=from_assistant.and.(name=='open_project'.or.name=='new_project'.or.name=='switch_project_tab')
            external_cursor=.not.from_assistant.and.assistant_session%root_fd<0
            if(from_assistant.and.associated(observe_start))call observe_start(name,args)
            if(external_cursor)call set_studio_ai_tool_cursor(workspace,name,args)
            call execute(workspace,name,args,result)
            if(assistant_switching_project.and.result%boolValue('ok'))then
                call result%insert('assistant_session_current',rebind_assistant_mcp_project(workspace))
            end if
            if(external_cursor)call set_studio_ai_tool_cursor(workspace,name,args,result)
            assistant_switching_project=.false.
            dispatching=.false.
        end if
    end subroutine run_studio_mcp_tool

    subroutine poll_mailbox(data)
        class(*),intent(inout),optional :: data
        type(QWidget) :: modal

        if(dispatching.or..not.associated(workspace))return
        if(assistant_session%root_fd>=0)then
            if(.not.assistant_session_current())call stop_studio_mcp(.true.)
        end if
        if(workspace%updating.or.workspace%dragging.or.workspace%control_dragging.or. &
            knot_edit_in_progress(workspace).or.transform_dragging(workspace))return
        modal=QApplication_activeModalWidget()
        if(.not.modal%isNull())return
        call process_pending_approval()
        if(pending%directory<0)call poll_session(assistant_session,.true.)
        call poll_session(external_session,.false.)
    end subroutine poll_mailbox

    subroutine poll_session(session,from_assistant)
        type(mailbox_session),intent(inout) :: session
        logical,intent(in) :: from_assistant
        type(qt_string),allocatable :: entries(:)
        integer(c_int) :: fd,status
        type(file_metadata) :: metadata
        integer :: i,count
        logical :: require_approval
        character(:),allocatable :: name,payload,rejection

        if(.not.session_path_current(session))return
        require_approval=from_assistant.and.assistant_requires_approval
        entries=QDir_entryList(session%path,ior(QDir_Dirs,QDir_NoDotAndDotDot))
        count=0
        do i=1,size(entries)
            name=entries(i)%s
            if(.not.mailbox_name(name))cycle
            fd=open_directory(session%root_fd,name//c_null_char)
            if(fd<0)cycle
            status=inspect_file(fd,metadata)
            if(status==0)then
                if(metadata%is_private/=1)then
                    status=close_file(fd)
                    cycle
                end if
                if(require_approval)then
                    payload=read_request(fd,request_file)
                    if(queue_approval(fd,payload,rejection))return
                end if
                status=rename_child(fd,request_file//c_null_char,fd,running_file//c_null_char,0_c_int)
                if(status==0)then
                    if(require_approval)then
                        call deliver_request(fd,from_assistant,payload,rejection)
                    else
                        call deliver_request(fd,from_assistant)
                    end if
                    count=count+1
                end if
            end if
            status=close_file(fd)
            if(count>=calls_per_tick.or.session%root_fd<0)exit
        end do
    end subroutine poll_session

    function pending_studio_mcp_approval(ticket) result(summary)
        integer(int64),intent(out),optional :: ticket
        character(:),allocatable :: summary

        summary=''
        if(present(ticket))ticket=0
        if(pending%directory<0.or.pending%decision/=approval_waiting)return
        summary=pending%summary
        if(present(ticket))ticket=pending%ticket
    end function pending_studio_mcp_approval

    subroutine decide_studio_mcp_approval(allowed,ticket)
        logical,intent(in) :: allowed
        integer(int64),intent(in),optional :: ticket

        if(pending%directory<0.or.pending%decision/=approval_waiting)return
        if(present(ticket))then
            if(ticket/=pending%ticket)return
        end if
        pending%decision=merge(approval_allowed,approval_denied,allowed)
    end subroutine decide_studio_mcp_approval

    logical function queue_approval(directory,payload,rejection) result(queued)
        integer(c_int),intent(in) :: directory
        character(*),intent(in) :: payload
        character(:),allocatable,intent(out) :: rejection
        type(QJsonObject) :: args
        character(:),allocatable :: name,details,summary
        integer(c_int) :: marker,status
        integer :: i,last

        queued=.false.
        rejection=''
        call decode_request(payload,name,args)
        if(len(name)==0)return
        if(.not.requires_approval(name,args))then
            call args%delete()
            return
        end if

        summary=''
        if(name=='edit_patch')then
            if(args%stringValue('op')=='control')then
                call control_approval_summary(workspace,args,summary,rejection)
            else
                call refinement_approval_summary(workspace,args,summary,rejection)
            end if
        end if
        if(len(rejection)>0)then
            call args%delete()
            return
        end if
        call system_clock(pending%started)
        next_ticket=next_ticket+1
        pending%ticket=next_ticket
        pending%directory=directory
        pending%payload=payload
        pending%summary=summary
        if(len(pending%summary)==0)then
            pending%summary='Allow ForCAD Studio to run '//name//'?'
            if(.not.args%isEmpty())then
                details=args%toJson(.true.)
                last=min(480,len(details))
                do i=last+1,2,-1
                    if(i>len(details))exit
                    if(iand(iachar(details(i:i)),192)/=128)exit
                    last=i-2
                end do
                pending%summary=pending%summary//new_line('a')//details(:last)
                if(last<len(details))pending%summary=pending%summary//'…'
            end if
        end if
        pending%stamp=document_undo_stamp(workspace%document)
        pending%document_revision=workspace%document%document_revision
        pending%project_path=current_project_path()
        pending%units=[workspace%document%model_unit,workspace%document%display_unit]
        pending%selection_dependent=uses_selection(name,args)
        if(pending%selection_dependent)then
            pending%focus=selection_focus()
            pending%selection=workspace%selection
        end if
        call args%delete()
        marker=create_private_file(directory,approval_file//c_null_char)
        if(marker>=0)status=close_file(marker)
        queued=.true.
    end function queue_approval

    logical function requires_approval(name,args) result(required)
        character(*),intent(in) :: name
        type(QJsonObject),intent(in) :: args

        required=.true.
        select case(name)
        case('inspect','inspect_patch','inspect_assembly','inspect_view','preview_refinement','preview_controls','preview_mesh','measure', &
             'list_commands','list_project_tabs','studio_skill', &
             'select','set_selection','select_region','set_probe','set_camera','set_parameter_view','set_workspace','set_mesh_view','pick_mesh_cell')
            required=.false.
        case('execute')
            select case(args%stringValue('command'))
            case('workspace.curve','workspace.surface','workspace.volume','workspace.reset', &
                 'camera.isometric','camera.fit','camera.fit-selection','camera.fit-panels','camera.front','camera.back', &
                 'camera.left','camera.right','camera.top','camera.bottom','camera.zoom-in','camera.zoom-out', &
                 'layout.toggle-focus','layout.canvas','layout.restore','layout.scene','layout.domain','layout.basis','layout.inspector', &
                 'patch.inspect','scene.reveal','scene.expand','scene.collapse','scene.exit-isolation','selection.clear', &
                 'display.shaded-faces','display.shaded-edges','display.fast-shaded','display.wireframe', &
                 'display.clean','display.structure','display.shaded','display.wire-layer','display.control-net', &
                 'display.control-points','display.parameter-lines','display.ground-grid','display.isolate','display.physical-knots','display.interactive', &
                 'view.select','view.box','view.lasso','view.orbit','view.move-gizmo','view.rotate-gizmo','view.scale-gizmo', &
                 'view.end-transform','view.snap','view.snap-off','view.snap-vertex','view.snap-closest', &
                 'pick.automatic','pick.faces','pick.edges','pick.vertices','pick.controls','pick.multipatches', &
                 'mesh.shaded','mesh.only','mesh.pick','mesh.select-source','multipatch.validate','inspector.show','probe.reset', &
                 'domain.probe','domain.reset-view','domain.zoom-in','domain.zoom-out','domain.grid','domain.support', &
                 'knots.show','knots.groups','knots.select-tool','knots.insert-tool','knots.remove-tool','knots.u','knots.v','knots.w', &
                 'basis.show','basis.functions','basis.derivative','basis.legend','basis.probe','basis.u','basis.v','basis.w', &
                 'basis.du','basis.dv','basis.dw')
                required=.false.
            case default
                return
            end select
        case default
            return
        end select
    end function requires_approval

    logical function uses_selection(name,args) result(dependent)
        character(*),intent(in) :: name
        type(QJsonObject),intent(in) :: args
        character(:),allocatable :: suffix

        dependent=.false.
        select case(name)
        case('transform_selection','set_material')
            dependent=.true.
        case('delete')
            dependent=args%boolValue('selection')
        case('execute')
            dependent=args%stringValue('command')/='display.show-all'
        case('edit_patch','mesh')
            dependent=.not.args%contains('target')
        case('export_file')
            suffix=studio_exchange_suffix(args%stringValue('path'))
            dependent=suffix=='.iges'.or.suffix=='.igs'.or.suffix=='.vtk'
        case default
            return
        end select
    end function uses_selection

    function selection_focus() result(focus)
        integer :: focus(6)

        focus=[workspace%entity_kind,merge(workspace%live_scene%variant,0,workspace%live_scene%valid), &
            workspace%selected_multipatch,workspace%selected_mesh,workspace%selected_control,workspace%knot_axis]
    end function selection_focus

    logical function selection_current() result(current)
        integer :: i

        current=.false.
        if(any(selection_focus()/=pending%focus))return
        if(workspace%selection%count()/=pending%selection%count())return
        do i=1,pending%selection%count()
            associate(a=>workspace%selection%items(i),b=>pending%selection%items(i))
                if(any([a%category,a%id,a%rank,a%part,a%item,a%feature]/=[b%category,b%id,b%rank,b%part,b%item,b%feature]))return
            end associate
        end do
        current=.true.
    end function selection_current

    subroutine process_pending_approval()
        type(file_metadata) :: metadata
        integer(c_int) :: directory,status
        integer(int64) :: now,rate
        character(:),allocatable :: rejection
        logical :: expired

        directory=pending%directory
        if(directory<0)return
        status=inspect_entry(directory,request_file//c_null_char,metadata)
        if(status/=0)then
            call clear_pending_approval()
            return
        end if
        call system_clock(now,rate)
        expired=rate>0.and.now-pending%started>=approval_seconds*rate
        if(.not.expired.and.pending%decision==approval_waiting)return

        rejection=''
        if(expired)then
            rejection='ForCAD Studio approval expired. Ask again to make this change.'
        else if(pending%decision==approval_denied)then
            rejection='The user denied this ForCAD Studio action.'
        else if(workspace%document%document_revision/=pending%document_revision.or. &
                current_project_path()/=pending%project_path.or.document_undo_stamp(workspace%document)/=pending%stamp.or. &
                any([workspace%document%model_unit,workspace%document%display_unit]/=pending%units))then
            rejection='The design changed while waiting for approval. Inspect the design and ask again.'
        else if(pending%selection_dependent)then
            if(.not.selection_current())rejection='The selection changed while waiting for approval. Inspect the selection and ask again.'
        end if
        status=rename_child(directory,request_file//c_null_char,directory,running_file//c_null_char,0_c_int)
        if(status==0)call deliver_request(directory,.true.,pending%payload,rejection)
        call clear_pending_approval()
    end subroutine process_pending_approval

    subroutine clear_pending_approval()
        integer(c_int) :: status

        if(pending%directory>=0)then
            status=remove_child(pending%directory,approval_file//c_null_char,0_c_int)
            status=close_file(pending%directory)
        end if
        pending=studio_approval()
    end subroutine clear_pending_approval

    logical function mailbox_name(name) result(valid)
        character(*),intent(in) :: name
        character(*),parameter :: allowed='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'

        valid=.false.
        if(len(name)/=20)return
        if(name(:4)/='mcp-')return
        valid=verify(name(5:),allowed)==0
    end function mailbox_name

    subroutine deliver_request(directory,from_assistant,reviewed_payload,rejection)
        integer(c_int),intent(in) :: directory
        logical,intent(in) :: from_assistant
        character(*),intent(in),optional :: reviewed_payload,rejection
        type(QJsonObject) :: args,result
        character(:),allocatable :: payload,name,problem
        integer(int64) :: before
        logical :: committed

        payload=read_request(directory,running_file)
        committed=.false.
        call decode_request(payload,name,args)
        problem=''
        if(present(rejection))problem=rejection
        if(.not.from_assistant)then
            if(.not.external_project_current())problem='This MCP session does not have access to the active project.'// &
                ' Start ForCAD Studio with --project and --mcp for the project you want to use.'
        end if
        if(present(reviewed_payload))then
            if(payload/=reviewed_payload)problem='The ForCAD Studio request changed while waiting.'
        end if
        if(len(problem)>0)then
            result=failed(problem)
        else if(len(name)==0)then
            result=failed('Invalid ForCAD Studio tool request.')
        else
            before=document_undo_stamp(workspace%document)
            call run_studio_mcp_tool(name,args,result,from_assistant)
            committed=document_undo_stamp(workspace%document)/=before
        end if
        if(len(name)>0)then
            if(associated(observe))call observe(name,result,args,from_assistant,committed)
            call args%delete()
        end if
        call publish_response(directory,result%toJson(.true.))
        call result%delete()
    end subroutine deliver_request

    subroutine decode_request(payload,name,args)
        character(*),intent(in) :: payload
        character(:),allocatable,intent(out) :: name
        type(QJsonObject),intent(out) :: args
        type(QJsonDocument) :: document
        type(QJsonObject) :: request

        name=''
        document=QJsonDocument_fromJson(payload)
        if(document%isObject())then
            request=document%object()
            if(request%size()==2.and.request%valueType('tool')==QJsonValue_String.and. &
                request%valueType('arguments')==QJsonValue_Object)then
                name=request%stringValue('tool')
                if(len(name)>0)args=request%objectValue('arguments')
            end if
            call request%delete()
        end if
        call document%delete()
    end subroutine decode_request

    function read_request(directory,file) result(payload)
        integer(c_int),intent(in) :: directory
        character(*),intent(in) :: file
        character(kind=c_char,len=:),allocatable,target :: payload
        type(file_metadata) :: metadata
        integer(c_int) :: fd,status
        integer(c_int64_t) :: count
        integer :: offset

        payload=''
        fd=open_regular(directory,file//c_null_char)
        if(fd<0)return
        status=inspect_file(fd,metadata)
        if(status==0)then
            if(metadata%is_private==1.and. &
                is_regular_file(metadata).and. &
                metadata%links==1.and.metadata%bytes>0.and.metadata%bytes<=max_mcp_request_bytes)then
                payload=repeat(' ',int(metadata%bytes))
                offset=0
                do while(offset<len(payload))
                    count=read_file(fd,c_loc(payload(offset+1:offset+1)),int(len(payload)-offset,c_size_t))
                    if(count<=0)exit
                    offset=offset+int(count)
                end do
                if(offset/=len(payload))payload=''
            end if
        end if
        status=close_file(fd)
    end function read_request

    subroutine publish_response(directory,text)
        integer(c_int),intent(in) :: directory
        character(kind=c_char,len=*),intent(in),target :: text
        integer(c_int) :: fd,status
        integer(c_int64_t) :: count

        fd=create_private_file(directory,temporary_file//c_null_char)
        if(fd<0)return
        count=write_file(fd,c_loc(text(1:1)),len(text,kind=c_size_t))
        status=close_file(fd)
        if(status==0.and.count==len(text,kind=c_size_t))then
            status=rename_child(directory,temporary_file//c_null_char,directory,response_file//c_null_char,0_c_int)
            if(status==0)return
        end if
        status=remove_child(directory,temporary_file//c_null_char,0_c_int)
    end subroutine publish_response

    function failed(message) result(result)
        character(*),intent(in) :: message
        type(QJsonObject) :: result

        result=QJsonObject()
        call result%insert('ok',.false.)
        call result%insert('error',message)
    end function failed
end module studio_mcp_bridge
