! The assistant only reads references from an explicitly assigned project folder.
module studio_project_folder
    use, intrinsic :: iso_c_binding, only: c_null_char,c_int
    use, intrinsic :: iso_fortran_env, only: error_unit
    use qt_core_classes, only: QDir_exists
    use studio_platform_files, only: resolved_file_path,remove_child,current_directory
    implicit none(type, external)
    private
    public :: canonical_project_folder, project_file_folder, canonical_project_file, save_project_plan

contains

    function canonical_project_folder(path) result(folder)
        character(*), intent(in) :: path
        character(:), allocatable :: folder

        folder=resolved_file_path(path)
        if(len(folder)==3)then
            if(folder(2:)==':/')folder=''
        end if
        if(len(folder)==0.or.folder=='/'.or..not.QDir_exists(folder))then
            folder=''
        end if
    end function canonical_project_folder

    function project_file_folder(path) result(folder)
        character(*),intent(in) :: path
        character(:),allocatable :: folder,file
        integer :: split

        folder=''
        file=resolved_file_path(path)
        split=index(file,'/',back=.true.)
        if(split<=0)return
        folder=file(:split-1)
        if(split==3)then
            if(file(2:3)==':/')folder=file(:3)
        end if
    end function project_file_folder

    function canonical_project_file(folder,path) result(file)
        character(*),intent(in) :: folder,path
        character(:),allocatable :: file

        file=resolved_file_path(path)
        if(len_trim(folder)==0.or.index(file,trim(folder)//'/')/=1)file=''
    end function canonical_project_file

    function save_project_plan(folder,json) result(filename)
        character(*),intent(in) :: folder,json
        character(:),allocatable :: filename
        character(16) :: serial
        integer :: i,io,unit,closed,looked_up
        logical :: exists

        filename=''
        if(len_trim(folder)==0.or.len_trim(json)==0)return
        if(.not.QDir_exists(folder))return
        do i=1,9999
            write(serial,'(i0.4)')i
            filename=folder//'/forcad-plan-'//trim(serial)//'.json'
            open(newunit=unit,file=filename,status='new',action='write',form='formatted',iostat=io)
            if(io==0)exit
            inquire(file=filename,exist=exists,iostat=looked_up)
            if(looked_up/=0)exit
            if(.not.exists)exit
        end do
        if(io/=0)then
            filename=''
            return
        end if
        write(unit,'(a)',iostat=io)json
        if(io==0)flush(unit,iostat=io)
        close(unit,iostat=closed)
        if(io/=0.or.closed/=0)then
            if(remove_child(current_directory,filename//c_null_char,0_c_int)/=0)then
                write(error_unit,'(a)')'Could not remove incomplete plan: '//filename
            end if
            filename=''
        end if
    end function save_project_plan
end module studio_project_folder
