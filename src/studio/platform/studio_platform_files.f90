module studio_platform_files
    use, intrinsic :: iso_c_binding, only: c_char,c_int,c_int64_t,c_size_t,c_ptr,c_null_char,c_associated,c_f_pointer
    implicit none(type, external)
    private

    integer(c_int), parameter, public :: current_directory=-100_c_int
    integer(c_int), parameter, public :: platform_macos=1,platform_windows=2
    public :: studio_platform
    integer(c_int), parameter :: regular_file=1_c_int

    type,bind(C),public :: file_metadata
        integer(c_int) :: kind,owned,is_private,links
        integer(c_int64_t) :: volume,inode,bytes
    end type file_metadata

    public :: open_directory,open_regular,create_private_file,create_private_directory,private_directory,restrict_file
    public :: close_file,inspect_path,inspect_file,inspect_entry,is_regular_file,same_file
    public :: read_file,write_file,rename_child,remove_child,replace_file,resolved_file_path
    public :: temporary_private_directory,make_private_child,create_private_child,copy_open_file,publish_file,last_path_separator

    interface
        integer(c_int) function studio_platform() bind(C,name='studio_assistant_platform')
            import c_int
            implicit none(type,external)
        end function studio_platform

        integer(c_int) function open_directory(parent,path) bind(C,name='studio_open_directory')
            import c_char,c_int
            implicit none(type, external)
            integer(c_int),value :: parent
            character(kind=c_char),intent(in) :: path(*)
        end function open_directory

        integer(c_int) function open_regular(parent,path) bind(C,name='studio_open_regular')
            import c_char,c_int
            implicit none(type, external)
            integer(c_int),value :: parent
            character(kind=c_char),intent(in) :: path(*)
        end function open_regular

        integer(c_int) function create_private_file(parent,path) bind(C,name='studio_create_file')
            import c_char,c_int
            implicit none(type, external)
            integer(c_int),value :: parent
            character(kind=c_char),intent(in) :: path(*)
        end function create_private_file

        integer(c_int) function mkdir_private(parent,path) bind(C,name='studio_mkdir_private')
            import c_char,c_int
            implicit none(type, external)
            integer(c_int),value :: parent
            character(kind=c_char),intent(in) :: path(*)
        end function mkdir_private

        integer(c_int) function restrict_file(path) bind(C,name='studio_restrict_file')
            import c_int,c_char
            implicit none(type,external)
            character(kind=c_char),intent(in) :: path(*)
        end function restrict_file

        integer(c_int) function is_private_directory(path) bind(C,name='studio_is_private_directory')
            import c_char,c_int
            implicit none(type, external)
            character(kind=c_char),intent(in) :: path(*)
        end function is_private_directory

        integer(c_int) function close_file(descriptor) bind(C,name='studio_close_fd')
            import c_int
            implicit none(type, external)
            integer(c_int),value :: descriptor
        end function close_file

        integer(c_int) function inspect_file(descriptor,metadata) bind(C,name='studio_inspect_fd')
            import c_int,file_metadata
            implicit none(type, external)
            integer(c_int),value :: descriptor
            type(file_metadata),intent(out) :: metadata
        end function inspect_file

        integer(c_int) function inspect_entry(parent,path,metadata) bind(C,name='studio_inspect_entry')
            import c_char,c_int,file_metadata
            implicit none(type, external)
            integer(c_int),value :: parent
            character(kind=c_char),intent(in) :: path(*)
            type(file_metadata),intent(out) :: metadata
        end function inspect_entry

        integer(c_int64_t) function read_file(descriptor,buffer,bytes) bind(C,name='studio_read_fd')
            import c_int,c_int64_t,c_size_t,c_ptr
            implicit none(type, external)
            integer(c_int),value :: descriptor
            type(c_ptr),value :: buffer
            integer(c_size_t),value :: bytes
        end function read_file

        integer(c_int64_t) function write_file(descriptor,buffer,bytes) bind(C,name='studio_write_fd')
            import c_int,c_int64_t,c_size_t,c_ptr
            implicit none(type, external)
            integer(c_int),value :: descriptor
            type(c_ptr),value :: buffer
            integer(c_size_t),value :: bytes
        end function write_file

        integer(c_int) function rename_child(parent,old_name,target,new_name,replace) bind(C,name='studio_rename_child')
            import c_int,c_char
            implicit none(type, external)
            integer(c_int),value :: parent,target,replace
            character(kind=c_char),intent(in) :: old_name(*),new_name(*)
        end function rename_child

        integer(c_int) function remove_child(parent,path,directory) bind(C,name='studio_remove_child')
            import c_int,c_char
            implicit none(type, external)
            integer(c_int),value :: parent,directory
            character(kind=c_char),intent(in) :: path(*)
        end function remove_child

        integer(c_int) function replace_file(old_path,new_path) bind(C,name='studio_replace_path')
            import c_int,c_char
            implicit none(type, external)
            character(kind=c_char),intent(in) :: old_path(*),new_path(*)
        end function replace_file

        type(c_ptr) function native_realpath(path) bind(C,name='studio_realpath')
            import c_char,c_ptr
            implicit none(type, external)
            character(kind=c_char),intent(in) :: path(*)
        end function native_realpath

        integer(c_size_t) function string_length(address) bind(C,name='strlen')
            import c_ptr,c_size_t
            implicit none(type, external)
            type(c_ptr),value :: address
        end function string_length

        subroutine release(address) bind(C,name='free')
            import c_ptr
            implicit none(type, external)
            type(c_ptr),value :: address
        end subroutine release

        pure integer(c_size_t) function find_separator(path,length) bind(C,name='studio_path_separator')
            import c_char,c_size_t
            implicit none(type,external)
            character(kind=c_char),intent(in) :: path(*)
            integer(c_size_t),value :: length
        end function find_separator

        integer(c_int) function fill_suffix(name,length) bind(C,name='studio_fill_suffix')
            import c_char,c_int,c_size_t
            implicit none(type, external)
            character(kind=c_char),intent(out) :: name(*)
            integer(c_size_t),value :: length
        end function fill_suffix

        integer(c_int) function copy_open_file(descriptor,path) bind(C,name='studio_copy_from_fd')
            import c_char,c_int
            implicit none(type, external)
            integer(c_int),value :: descriptor
            character(kind=c_char),intent(in) :: path(*)
        end function copy_open_file

        integer(c_int) function publish_file(path,parent,leaf,replace) bind(C,name='studio_publish_file')
            import c_char,c_int
            implicit none(type, external)
            character(kind=c_char),intent(in) :: path(*),leaf(*)
            integer(c_int),value :: parent,replace
        end function publish_file
    end interface

contains

    pure integer function last_path_separator(path) result(position)
        character(*),intent(in) :: path

        position=int(find_separator(path,len(path,kind=c_size_t)))
    end function last_path_separator

    logical function inspect_path(path,metadata) result(found)
        character(*),intent(in) :: path
        type(file_metadata),intent(out) :: metadata

        found=inspect_entry(current_directory,path//c_null_char,metadata)==0
    end function inspect_path

    pure logical function is_regular_file(metadata) result(regular)
        type(file_metadata),intent(in) :: metadata

        regular=metadata%kind==regular_file
    end function is_regular_file

    pure logical function same_file(first,second) result(matches)
        type(file_metadata),intent(in) :: first,second

        matches=first%kind==second%kind.and.first%volume==second%volume.and.first%inode==second%inode
    end function same_file

    function resolved_file_path(path) result(resolved)
        character(*),intent(in) :: path
        character(:),allocatable :: resolved
        character(kind=c_char),pointer :: bytes(:)
        type(c_ptr) :: address
        integer :: i,n

        resolved=''
        if(len_trim(path)==0.or.index(path,c_null_char)>0)return
        address=native_realpath(trim(path)//c_null_char)
        if(.not.c_associated(address))return
        n=int(string_length(address))
        call c_f_pointer(address,bytes,[n])
        resolved=repeat(' ',n)
        do i=1,n
            resolved(i:i)=bytes(i)
        end do
        call release(address)
    end function resolved_file_path

    logical function private_directory(path) result(valid)
        character(*),intent(in) :: path

        valid=is_private_directory(path//c_null_char)==1_c_int
    end function private_directory

    logical function create_private_directory(path) result(created)
        character(*),intent(in) :: path

        created=mkdir_private(current_directory,path//c_null_char)==0
    end function create_private_directory

    logical function create_private_child(parent,name) result(created)
        integer(c_int),intent(in) :: parent
        character(*),intent(in) :: name

        created=mkdir_private(parent,name//c_null_char)==0
    end function create_private_child

    function make_private_child(parent,prefix) result(name)
        integer(c_int),intent(in) :: parent
        character(*),intent(in) :: prefix
        character(:),allocatable :: name
        character(kind=c_char,len=16) :: suffix
        type(file_metadata) :: existing
        integer :: attempt

        name=''
        do attempt=1,8
            if(fill_suffix(suffix,len(suffix,kind=c_size_t))/=0)return
            if(mkdir_private(parent,prefix//suffix//c_null_char)==0)then
                name=prefix//suffix
                return
            end if
            if(inspect_entry(parent,prefix//suffix//c_null_char,existing)/=0)return
        end do
    end function make_private_child

    function temporary_private_directory(parent,prefix) result(path)
        character(*),intent(in) :: parent,prefix
        character(:),allocatable :: path,root,leaf
        integer(c_int) :: descriptor,status

        path=''
        root=resolved_file_path(parent)
        if(len(root)==0)return
        descriptor=open_directory(current_directory,root//c_null_char)
        if(descriptor<0)return
        leaf=make_private_child(descriptor,prefix)
        status=close_file(descriptor)
        if(len(leaf)==0)return
        if(root(len(root):)=='/')then
            path=root//leaf
        else
            path=root//'/'//leaf
        end if
    end function temporary_private_directory
end module studio_platform_files
