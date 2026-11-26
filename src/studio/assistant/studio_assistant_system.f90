module studio_assistant_system
    use, intrinsic :: iso_c_binding, only: c_char,c_int,c_int64_t,c_ptr,c_size_t
    use studio_platform_files, only: assistant_platform=>studio_platform,platform_macos,platform_windows
    implicit none(type,external)
    private
    public :: start_agent_worker,poll_agent_worker,close_agent_worker
    public :: append_agent_frame,seek_agent_file,epoch_microseconds,assistant_platform,platform_macos,platform_windows
    public :: set_process_environment,unset_process_environment

    interface
        type(c_ptr) function start_agent_worker(program,directory) bind(C,name='studio_agent_start')
            import c_ptr,c_char
            implicit none(type,external)
            character(kind=c_char),intent(in) :: program(*),directory(*)
        end function start_agent_worker

        integer(c_int) function poll_agent_worker(worker,successful) bind(C,name='studio_agent_poll')
            import c_ptr,c_int
            implicit none(type,external)
            type(c_ptr),value :: worker
            integer(c_int),intent(out) :: successful
        end function poll_agent_worker

        subroutine close_agent_worker(worker) bind(C,name='studio_agent_close')
            import c_ptr
            implicit none(type,external)
            type(c_ptr),value :: worker
        end subroutine close_agent_worker

        integer(c_int) function append_agent_frame(path,data,length) bind(C,name='studio_append_event')
            import c_char,c_int,c_size_t
            implicit none(type,external)
            character(kind=c_char),intent(in) :: path(*),data(*)
            integer(c_size_t),value :: length
        end function append_agent_frame

        integer(c_int) function seek_agent_file(descriptor,offset) bind(C,name='studio_agent_seek')
            import c_int,c_int64_t
            implicit none(type,external)
            integer(c_int),value :: descriptor
            integer(c_int64_t),value :: offset
        end function seek_agent_file

        integer(c_int64_t) function epoch_microseconds() bind(C,name='studio_epoch_microseconds')
            import c_int64_t
            implicit none(type,external)
        end function epoch_microseconds

        integer(c_int) function set_process_environment(name,value) bind(C,name='studio_set_environment')
            import c_char,c_int
            implicit none(type,external)
            character(kind=c_char),intent(in) :: name(*),value(*)
        end function set_process_environment

        integer(c_int) function unset_process_environment(name) bind(C,name='studio_unset_environment')
            import c_char,c_int
            implicit none(type,external)
            character(kind=c_char),intent(in) :: name(*)
        end function unset_process_environment
    end interface

end module studio_assistant_system
