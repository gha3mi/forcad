module studio_storage_paths
    use, intrinsic :: iso_c_binding, only: c_int,c_null_char
    use qt, only: QCoreApplication_applicationFilePath,QDir_mkpath
    use studio_platform_files, only: resolved_file_path,last_path_separator,current_directory, &
        open_directory,close_file,inspect_file,file_metadata,create_private_child,private_directory
    implicit none(type,external)
    private
    public :: studio_application_directory,studio_project_directory,studio_work_directory

contains

    function studio_application_directory(area) result(folder)
        character(*),intent(in),optional :: area
        character(:),allocatable :: folder,root,executable
        integer :: length,status,separator

        folder=''
        call get_environment_variable('FORCAD_STUDIO_CONFIG_DIR',length=length,status=status)
        if(status==0.and.length>0)then
            allocate(character(length)::root)
            call get_environment_variable('FORCAD_STUDIO_CONFIG_DIR',value=root,status=status)
            if(status/=0.or.len_trim(root)==0)return
            if(index(root,c_null_char)>0)return
            root=portable_override(trim(root))
        else
            call get_command_argument(0,length=length,status=status)
            executable=''
            if(status==0.and.length>0)then
                deallocate(executable)
                allocate(character(length)::executable)
                call get_command_argument(0,value=executable,status=status)
                if(status/=0)executable=''
            end if
            if(last_path_separator(executable)==0)executable=QCoreApplication_applicationFilePath()
            executable=resolved_file_path(executable)
            separator=last_path_separator(executable)
            if(separator==0)return
            root=storage_child(executable(:separator),'.forcad-studio')
        end if
        if(present(area))then
            folder=storage_child(root,area)
        else
            folder=root
        end if
    end function studio_application_directory

    function portable_override(path) result(folder)
        character(*),intent(in) :: path
        character(:),allocatable :: folder,parent
        integer :: last,separator

        folder=resolved_file_path(path)
        if(len(folder)>0)then
            if(.not.private_directory(folder))folder=''
            return
        end if
        last=len(path)
        do while(last>1)
            if(last_path_separator(path(:last))/=last)exit
            last=last-1
        end do
        separator=last_path_separator(path(:last))
        parent='.'
        if(separator>0)parent=path(:separator)
        if(.not.QDir_mkpath(parent))return
        folder=storage_child(parent,path(separator+1:last))
    end function portable_override

    function studio_project_directory(project,area) result(folder)
        character(*),intent(in) :: project
        character(*),intent(in),optional :: area
        character(:),allocatable :: folder

        character(:),allocatable :: root

        folder=''
        root=resolved_file_path(project)
        if(len(root)==0.or.root=='/')return
        if(len(root)==3)then
            if(root(2:)==':/')return
        end if
        folder=storage_child(root,'.forcad')
        if(present(area))folder=storage_child(folder,area)
    end function studio_project_directory

    function studio_work_directory(project) result(folder)
        character(*),intent(in),optional :: project
        character(:),allocatable :: folder

        if(present(project))then
            if(len_trim(project)>0)then
                folder=studio_project_directory(project,'work')
                return
            end if
        end if
        folder=studio_application_directory('work')
    end function studio_work_directory

    function storage_child(parent,name) result(folder)
        character(*),intent(in) :: parent,name
        character(:),allocatable :: folder,root
        type(file_metadata) :: metadata
        integer(c_int) :: parent_fd,child_fd,status,closed
        logical :: created

        folder=''
        if(len(parent)==0.or.len(name)==0.or.name=='.'.or.name=='..')return
        if(scan(name,'/:'//achar(92)//c_null_char)>0)return
        root=resolved_file_path(parent)
        if(len(root)==0)return
        parent_fd=open_directory(current_directory,root//c_null_char)
        if(parent_fd<0)return
        created=create_private_child(parent_fd,name)
        child_fd=open_directory(parent_fd,name//c_null_char)
        closed=close_file(parent_fd)
        if(child_fd<0)return
        status=inspect_file(child_fd,metadata)
        closed=close_file(child_fd)
        if(status/=0)return
        if(metadata%is_private/=1.or.metadata%owned/=1)return
        if(root(len(root):)/='/')root=root//'/'
        folder=root//name
    end function storage_child
end module studio_storage_paths
