module studio_process_mode
    use, intrinsic :: iso_c_binding, only: c_int,c_null_char
    use, intrinsic :: iso_fortran_env, only: error_unit
    use qt, only: QApplication
    use studio_assistant_worker, only: run_studio_assistant_worker
    use studio_assistant_system, only: set_process_environment,unset_process_environment
    use studio_mcp_stdio, only: run_studio_mcp_stdio
    use studio_exchange, only: dispatch_studio_iges_worker
    implicit none(type,external)
    private
    public :: dispatch_studio_process_mode,studio_argument

contains

    function studio_argument(index) result(value)
        integer,intent(in) :: index
        character(:),allocatable :: value
        integer :: length,status

        call get_command_argument(index,length=length,status=status)
        if(status/=0)then
            value=''
            return
        end if
        allocate(character(length)::value)
        call get_command_argument(index,value,status=status)
        if(status/=0)value=''
    end function studio_argument

    subroutine dispatch_studio_process_mode()
        character(:),allocatable :: mode,platform,role,session
        type(QApplication) :: app
        integer :: length,status,environment_status
        integer(c_int) :: ignored

        call dispatch_studio_iges_worker()
        mode=studio_argument(1)
        if(mode/='--mcp-stdio'.and.mode/='--assistant-worker')return

        call get_environment_variable('QT_QPA_PLATFORM',length=length,status=environment_status)
        if(environment_status==0)then
            allocate(character(length)::platform)
            call get_environment_variable('QT_QPA_PLATFORM',platform)
        end if
        ignored=set_process_environment('QT_QPA_PLATFORM'//c_null_char,'offscreen'//c_null_char)
        app=QApplication()
        if(allocated(platform))then
            ignored=set_process_environment('QT_QPA_PLATFORM'//c_null_char,platform//c_null_char)
        else
            ignored=unset_process_environment('QT_QPA_PLATFORM'//c_null_char)
        end if

        status=2
        select case(mode)
        case('--mcp-stdio')
            if(command_argument_count()==3.and.studio_argument(2)=='--session')then
                session=studio_argument(3)
                status=run_studio_mcp_stdio(session)
            else if(command_argument_count()==1)then
                call get_environment_variable('FORCAD_MCP_SESSION',length=length,status=environment_status)
                session=''
                if(environment_status==0)then
                    deallocate(session)
                    allocate(character(length)::session)
                    call get_environment_variable('FORCAD_MCP_SESSION',session)
                end if
                status=run_studio_mcp_stdio(session)
            else
                write(error_unit,'(a)') 'Usage: forcad-studio --mcp-stdio [--session DIRECTORY]'
            end if
        case('--assistant-worker')
            if(command_argument_count()==2.or.command_argument_count()==3)then
                role=studio_argument(3)
                if(len(role)>0)ignored=unset_process_environment('FORCAD_MCP_SESSION'//c_null_char)
                status=run_studio_assistant_worker(studio_argument(2),role)
            end if
        end select
        call app%delete()
        stop status,quiet=.true.
    end subroutine dispatch_studio_process_mode

end module studio_process_mode
