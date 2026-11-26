module studio_chat_attachment
    use, intrinsic :: iso_c_binding, only: c_int,c_null_char,c_size_t,c_loc
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QWidget, qt_string, QFileDialog_getOpenFileNames, QFileDialog_DontResolveSymlinks, QDir_removeRecursively
    use studio_project_folder, only: canonical_project_folder
    use studio_platform_files, only: file_metadata,current_directory,resolved_file_path,inspect_file,open_regular,close_file, &
        copy_open_file,read_file,temporary_private_directory,restrict_file,last_path_separator, &
        inspect_path,same_file,open_directory,inspect_entry,remove_child
    use studio_storage_paths, only: studio_project_directory
    implicit none(type, external)
    private
    public :: studio_attachment, choose_studio_attachments, read_studio_attachment, prepare_studio_attachment, discard_attachment_copy

    integer, parameter :: max_image_bytes=4*1024*1024, max_text_bytes=24*1024
    integer, parameter :: max_image_pixels=16*1024*1024
    character(*), parameter :: png_signature=char(137)//'PNG'//achar(13)//achar(10)//achar(26)//achar(10)
    character(*), parameter :: jpeg_signature=char(255)//char(216)//char(255)

    type :: studio_attachment
        character(:), allocatable :: name, path, image, text, relative_path, project_copy
        type(file_metadata),allocatable :: copy_directory
    end type studio_attachment
contains
    subroutine choose_studio_attachments(parent,attachments,error,directory)
        class(QWidget),intent(in) :: parent
        type(studio_attachment),allocatable,intent(out) :: attachments(:)
        character(:),allocatable,intent(out) :: error
        character(*),intent(in) :: directory
        type(qt_string),allocatable :: paths(:)
        type(studio_attachment),allocatable :: staged(:)
        character(:),allocatable :: problem
        integer :: i,count

        error=''
        paths=QFileDialog_getOpenFileNames('Attach files',dir=directory,filter='All files (*)',parent=parent, &
            options=QFileDialog_DontResolveSymlinks)
        allocate(staged(size(paths)))
        count=0
        do i=1,size(paths)
            call prepare_studio_attachment(paths(i)%s,directory,staged(count+1),problem)
            if(allocated(staged(count+1)%name))then
                count=count+1
            else if(len(error)==0)then
                error=problem
            end if
        end do
        if(count==size(staged))then
            call move_alloc(staged,attachments)
        else
            attachments=staged(:count)
        end if
    end subroutine choose_studio_attachments

    subroutine prepare_studio_attachment(path,directory,attachment,error)
        character(*), intent(in) :: path,directory
        type(studio_attachment), intent(out) :: attachment
        character(:), allocatable, intent(out) :: error
        character(:), allocatable :: name,project,staged,copy_folder,resolved,references,managed
        type(file_metadata) :: metadata,copy_directory
        integer(c_int) :: input,status
        logical :: reusable

        project=canonical_project_folder(directory)
        if(len(project)==0)then
            error='The project folder is unavailable.'
            return
        end if
        error='Choose a readable file.'
        input=open_reference(path,name,resolved,metadata,selected=.true.)
        if(input<0)return

        references=studio_project_directory(project,'references')
        managed=references
        if(index(resolved,project//'/.forcad/screenshots/')==1)managed=studio_project_directory(project,'screenshots')
        reusable=len(managed)>0.and.index(resolved,managed//'/')==1.and.metadata%owned==1.and.metadata%links==1
        if(reusable.and.metadata%is_private/=1)reusable=restrict_file(resolved//c_null_char)==0
        if(reusable)then
            call read_reference(input,name,resolved,metadata%bytes,attachment,error)
            return
        end if

        error='Could not copy the selected file into the project folder.'
        copy_folder=temporary_private_directory(references,'reference-')
        if(len(copy_folder)==0)then
            status=close_file(input)
            return
        end if
        staged=copy_folder//'/'//name
        status=-1
        if(inspect_path(copy_folder,copy_directory))status=copy_open_file(input,staged//c_null_char)
        if(close_file(input)/=0)status=-1
        if(status==0)then
            call read_studio_attachment(staged,attachment,error)
            if(allocated(attachment%path))then
                attachment%project_copy=copy_folder
                attachment%copy_directory=copy_directory
                return
            end if
        end if
        if(.not.QDir_removeRecursively(copy_folder))error='Could not remove the temporary attachment folder.'
    end subroutine prepare_studio_attachment

    subroutine discard_attachment_copy(attachment)
        type(studio_attachment),intent(in) :: attachment
        character(:),allocatable :: path,leaf
        type(file_metadata) :: metadata
        integer(c_int) :: folder,parent,status
        integer :: separator

        if(.not.allocated(attachment%project_copy).or..not.allocated(attachment%copy_directory))return
        if(.not.allocated(attachment%name))return
        path=attachment%project_copy
        separator=last_path_separator(path)
        if(separator<1.or.separator==len(path))return
        leaf=path(separator+1:)
        parent=open_directory(current_directory,path(:separator-1)//c_null_char)
        if(parent<0)return
        folder=open_directory(parent,leaf//c_null_char)
        if(folder>=0)then
            status=inspect_file(folder,metadata)
            if(status==0.and.same_file(metadata,attachment%copy_directory))then
                status=remove_child(folder,attachment%name//c_null_char,0_c_int)
            end if
            status=close_file(folder)
            status=inspect_entry(parent,leaf//c_null_char,metadata)
            if(status==0.and.same_file(metadata,attachment%copy_directory))then
                status=remove_child(parent,leaf//c_null_char,1_c_int)
            end if
        end if
        status=close_file(parent)
    end subroutine discard_attachment_copy

    integer(c_int) function open_reference(path,name,resolved,metadata,selected) result(input)
        character(*),intent(in) :: path
        character(:),allocatable,intent(out) :: name,resolved
        type(file_metadata),intent(out) :: metadata
        logical,intent(in) :: selected
        type(file_metadata) :: resolved_metadata
        character(:),allocatable :: source
        integer :: position
        integer(c_int) :: status

        input=-1
        if(len_trim(path)==0.or.index(path,achar(0))/=0)return
        position=last_path_separator(path)
        name=path(position+1:)
        if(len_trim(name)==0.or.len(name)>255)return
        resolved=resolved_file_path(path)
        if(len(resolved)==0)return
        source=path
        if(selected)source=resolved
        input=open_regular(current_directory,source//c_null_char)
        if(input<0)return
        status=inspect_file(input,metadata)
        if(status==0)then
            if(.not.inspect_path(resolved,resolved_metadata))then
                status=-1
            else if(.not.same_file(metadata,resolved_metadata).or.metadata%bytes<0)then
                status=-1
            end if
        end if
        if(status/=0)then
            status=close_file(input)
            input=-1
            return
        end if

        do position=1,len(name)
            if(iachar(name(position:position))<32.or.iachar(name(position:position))==127)name(position:position)='_'
        end do
    end function open_reference

    subroutine read_studio_attachment(path,attachment,error)
        character(*), intent(in) :: path
        type(studio_attachment), intent(out) :: attachment
        character(:), allocatable, intent(out) :: error
        character(:),allocatable :: name,resolved
        integer(c_int) :: input
        type(file_metadata) :: metadata

        error='Choose a readable file.'
        input=open_reference(path,name,resolved,metadata,selected=.false.)
        if(input<0)return
        call read_reference(input,name,path,metadata%bytes,attachment,error)
    end subroutine read_studio_attachment

    subroutine read_reference(input,name,path,bytes,attachment,error)
        integer(c_int),intent(in) :: input
        character(*),intent(in) :: name,path
        integer(int64),intent(in) :: bytes
        type(studio_attachment),intent(out) :: attachment
        character(:),allocatable,intent(out) :: error
        character(:),allocatable,target :: content
        integer(int64) :: read_bytes

        error='Could not read the selected attachment.'
        if(bytes>0.and.bytes<=max_image_bytes)then
            allocate(character(len=int(bytes))::content)
            read_bytes=read_file(input,c_loc(content(1:1)),int(bytes,c_size_t))
        else
            read_bytes=bytes
        end if
        if(close_file(input)/=0.or.read_bytes/=bytes)return
        attachment%name=name
        attachment%path=path
        error=''
        if(.not.allocated(content))return

        if(index(content,png_signature)==1.or.index(content,jpeg_signature)==1)then
            if(inline_image_size(content))attachment%image=encode_base64(content)
        else if(bytes<=max_text_bytes.and.valid_reference_text(content))then
            if(index(content,char(239)//char(187)//char(191))==1)content=content(4:)
            if(verify(content,' '//achar(9)//achar(10)//achar(13))/=0)attachment%text=content
        end if
    end subroutine read_reference

    logical function inline_image_size(content) result(valid)
        character(*),intent(in) :: content
        integer :: offset,marker,segment
        integer(int64) :: columns,rows

        valid=.false.
        if(index(content,png_signature)==1)then
            if(len(content)<45)return
            if(network_integer(content(9:12))/=13.or.content(13:16)/='IHDR')return
            if(content(len(content)-11:len(content)-4)/=repeat(achar(0),4)//'IEND')return
            columns=network_integer(content(17:20))
            rows=network_integer(content(21:24))
        else
            if(len(content)<15.or.index(content,char(255)//char(217),back=.true.)==0)return
            columns=0
            rows=0
            offset=3
            do
                if(offset+3>len(content))return
                if(content(offset:offset)/=char(255))return
                do
                    offset=offset+1
                    if(offset+2>len(content))return
                    if(content(offset:offset)/=char(255))exit
                end do
                marker=iachar(content(offset:offset))
                if(marker<192.or.(marker>=208.and.marker<=217))return
                segment=int(network_integer(content(offset+1:offset+2)))
                if(segment<2.or.offset+segment>len(content))return
                if(marker==218)exit
                select case(marker)
                case(192:195,197:199,201:203,205:207)
                    if(segment<8.or.columns>0)return
                    rows=network_integer(content(offset+4:offset+5))
                    columns=network_integer(content(offset+6:offset+7))
                case default
                end select
                offset=offset+segment+1
            end do
        end if
        if(min(columns,rows)<1.or.max(columns,rows)>max_image_pixels)return
        valid=columns*rows<=max_image_pixels
    end function inline_image_size

    integer(int64) function network_integer(bytes) result(value)
        character(*),intent(in) :: bytes
        integer :: i

        value=0
        do i=1,len(bytes)
            value=256*value+iachar(bytes(i:i))
        end do
    end function network_integer

    logical function valid_reference_text(input) result(valid)
        character(*), intent(in) :: input
        integer :: position,first,next,extra,code,least,i

        valid=.false.
        position=1
        do while(position<=len(input))
            first=iachar(input(position:position))
            select case(first)
            case(9,10,13,32:126)
                extra=0
                code=first
                least=0
            case(194:223)
                extra=1
                code=iand(first,31)
                least=128
            case(224:239)
                extra=2
                code=iand(first,15)
                least=2048
            case(240:244)
                extra=3
                code=iand(first,7)
                least=65536
            case default
                return
            end select
            if(position+extra>len(input))return
            do i=1,extra
                next=iachar(input(position+i:position+i))
                if(next<128.or.next>191)return
                code=64*code+iand(next,63)
            end do
            if(code<least.or.code>1114111.or.(code>=55296.and.code<=57343))return
            position=position+extra+1
        end do
        valid=.true.
    end function valid_reference_text

    function encode_base64(bytes) result(encoded)
        character(*), intent(in) :: bytes
        character(:), allocatable :: encoded
        character(*), parameter :: letters='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'
        integer :: source,target,b1,b2,b3,indices(4),i,remaining

        allocate(character(len=4*((len(bytes)+2)/3))::encoded)
        do source=1,len(bytes),3
            target=1+4*((source-1)/3)
            remaining=min(len(bytes)-source+1,3)
            b1=iachar(bytes(source:source))
            b2=0
            b3=0
            if(remaining>=2)b2=iachar(bytes(source+1:source+1))
            if(remaining==3)b3=iachar(bytes(source+2:source+2))
            indices=[ishft(b1,-2),ishft(iand(b1,3),4)+ishft(b2,-4),ishft(iand(b2,15),2)+ishft(b3,-6),iand(b3,63)]+1
            do i=1,remaining+1
                encoded(target+i-1:target+i-1)=letters(indices(i):indices(i))
            end do
            encoded(target+remaining+1:target+3)=repeat('=',3-remaining)
        end do
    end function encode_base64
end module studio_chat_attachment
