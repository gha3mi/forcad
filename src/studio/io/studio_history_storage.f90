module studio_history_storage
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QDir_removeRecursively
    use studio_platform_files, only: temporary_private_directory,file_metadata,inspect_path,same_file
    use studio_storage_paths, only: studio_application_directory,studio_project_directory
    implicit none(type, external)
    private
    public :: history_blob,create_history_blob,open_history_blob,release_history_blob,history_blob_bytes,finish_history_storage

    type :: history_blob
        private
        integer(int64) :: id=0
        integer :: folder=0
    contains
        procedure :: recorded => history_blob_recorded
    end type history_blob

    type :: history_directory
        character(:), allocatable :: root,path
        type(file_metadata) :: identity
    end type history_directory

    type(history_directory), allocatable :: directories(:)
    integer(int64) :: next_id=0

contains

    pure logical function history_blob_recorded(blob) result(recorded)
        class(history_blob), intent(in) :: blob

        recorded=blob%id/=0
    end function history_blob_recorded

    subroutine create_history_blob(blob,unit,status,project)
        type(history_blob), intent(out) :: blob
        integer, intent(out) :: unit,status
        character(*), intent(in), optional :: project

        if(present(project))then
            if(len_trim(project)>0)then
                call create_blob_in(studio_project_directory(project,'undo'),blob,unit,status)
                if(status==0)return
            end if
        end if
        call create_blob_in(studio_application_directory('undo'),blob,unit,status)
    end subroutine create_history_blob

    subroutine create_blob_in(root,blob,unit,status)
        character(*), intent(in) :: root
        type(history_blob), intent(out) :: blob
        integer, intent(out) :: unit,status
        integer :: folder

        status=1
        folder=history_directory_index(root)
        if(folder==0)return
        next_id=next_id+1
        blob=history_blob(id=next_id,folder=folder)
        open(newunit=unit,file=blob_path(blob),status='new',access='stream',form='unformatted',action='write',iostat=status)
        if(status/=0)blob=history_blob()
    end subroutine create_blob_in

    integer function history_directory_index(root) result(folder)
        character(*), intent(in) :: root
        type(history_directory) :: fresh
        integer :: i

        folder=0
        if(len(root)==0)return
        if(.not.allocated(directories))allocate(directories(0))
        do i=1,size(directories)
            if(directories(i)%root/=root)cycle
            if(.not.directory_is_current(directories(i)))cycle
            folder=i
            return
        end do
        fresh%root=root
        fresh%path=temporary_private_directory(root,'session-')
        if(len(fresh%path)==0)return
        if(.not.inspect_path(fresh%path,fresh%identity))return
        directories=[directories,fresh]
        folder=size(directories)
    end function history_directory_index

    logical function directory_is_current(directory) result(current)
        type(history_directory), intent(in) :: directory
        type(file_metadata) :: identity

        current=inspect_path(directory%path,identity)
        if(current)current=same_file(directory%identity,identity).and.identity%owned==1.and.identity%is_private==1
    end function directory_is_current

    subroutine open_history_blob(blob,unit,status)
        type(history_blob), intent(in) :: blob
        integer, intent(out) :: unit,status
        integer :: io
        character(:), allocatable :: path

        status=1
        path=blob_path(blob)
        if(len(path)==0)return
        open(newunit=unit,file=path,status='old',access='stream',form='unformatted',action='read',iostat=io)
        status=merge(0,1,io==0)
    end subroutine open_history_blob

    subroutine release_history_blob(blob)
        type(history_blob), intent(inout) :: blob
        integer :: unit,status

        call open_history_blob(blob,unit,status)
        if(status==0)then
            close(unit,status='delete',iostat=status)
            if(status/=0)return
        end if
        blob=history_blob()
    end subroutine release_history_blob

    integer(int64) function history_blob_bytes(blob) result(bytes)
        type(history_blob), intent(in) :: blob
        integer :: status
        character(:), allocatable :: path

        bytes=0
        path=blob_path(blob)
        if(len(path)==0)return
        inquire(file=path,size=bytes,iostat=status)
        if(status/=0.or.bytes<0)bytes=0
    end function history_blob_bytes

    function blob_path(blob) result(path)
        type(history_blob), intent(in) :: blob
        character(:), allocatable :: path
        character(24) :: name

        path=''
        if(blob%id==0.or..not.allocated(directories))return
        if(blob%folder<1.or.blob%folder>size(directories))return
        if(.not.directory_is_current(directories(blob%folder)))return
        write(name,'(i0)')blob%id
        path=directories(blob%folder)%path//'/'//trim(name)
    end function blob_path

    subroutine finish_history_storage()
        integer :: i
        logical :: removed

        if(.not.allocated(directories))return
        do i=1,size(directories)
            if(.not.directory_is_current(directories(i)))cycle
            removed=QDir_removeRecursively(directories(i)%path)
        end do
        deallocate(directories)
    end subroutine finish_history_storage

end module studio_history_storage
