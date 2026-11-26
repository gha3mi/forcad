module studio_assistant_process
    use qt, only: QProcess,QProcessEnvironment,QProcessEnvironment_systemEnvironment
    use studio_assistant_system, only: assistant_platform,platform_windows
    use studio_platform_files, only: file_metadata,inspect_path,is_regular_file,resolved_file_path
    implicit none(type,external)
    private
    public :: start_provider_process

contains

    subroutine start_provider_process(process,program,arguments)
        type(QProcess),intent(inout) :: process
        character(*),intent(in) :: program,arguments(:)
        type(QProcessEnvironment) :: environment
        character(:),allocatable :: paths,executable,script,node_script,package_arguments(:)

        if(assistant_platform()/=platform_windows)then
            call process%start(program,arguments)
            return
        end if

        environment=QProcessEnvironment_systemEnvironment()
        paths=environment%value('PATH')
        call environment%delete()
        executable=windows_file_on_path(paths,program//'.exe')
        if(len(executable)==0)then
            select case(program)
            case('codex')
                node_script='@openai/codex/bin/codex.js'
            case('claude')
                executable=windows_file_on_path(paths,'node_modules/@anthropic-ai/claude-code/bin/claude.exe')
                node_script='@anthropic-ai/claude-code/cli.js'
            case('opencode')
                executable=windows_file_on_path(paths,'node_modules/opencode-ai/bin/opencode.exe')
                node_script='opencode-ai/bin/opencode'
            case default
                return
            end select
            if(len(executable)==0)then
                script=windows_file_on_path(paths,'node_modules/'//node_script)
                if(len(script)==0)return
                executable=windows_file_on_path(paths,'node.exe')
                if(len(executable)==0)return
                allocate(character(max(len(arguments),len(script)))::package_arguments(size(arguments)+1))
                package_arguments(1)=script
                package_arguments(2:)=arguments
                call process%start(executable,package_arguments)
                return
            end if
        end if

        call process%start(executable,arguments)
    end subroutine start_provider_process

    function windows_file_on_path(paths,leaf) result(found)
        character(*),intent(in) :: paths,leaf
        character(:),allocatable :: found,folder,candidate
        type(file_metadata) :: metadata
        integer :: start,last,separator,i

        found=''
        start=1
        do while(start<=len(paths))
            separator=index(paths(start:),';')
            last=len(paths)
            if(separator>0)last=start+separator-2
            folder=trim(adjustl(paths(start:last)))
            if(len(folder)>2)then
                if(folder(1:1)=='"'.and.folder(len(folder):)=='"')folder=folder(2:len(folder)-1)
            end if
            if(len(folder)>2)then
                do i=1,len(folder)
                    if(folder(i:i)==achar(92))folder(i:i)='/'
                end do
                if(folder(2:3)==':/'.or.folder(1:2)=='//')then
                    candidate=resolved_file_path(folder//'/'//leaf)
                    if(len(candidate)>0)then
                        if(inspect_path(candidate,metadata))then
                            if(is_regular_file(metadata))then
                                found=candidate
                                return
                            end if
                        end if
                    end if
                end if
            end if
            start=last+2
        end do
    end function windows_file_on_path

end module studio_assistant_process
