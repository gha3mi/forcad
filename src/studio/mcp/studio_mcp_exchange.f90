! Project-scoped file operations for the MCP tools. ForCAD Studio's usual ForCAD
! exchange routines read/write the geometry and track the import in Undo.
module studio_mcp_exchange
    use, intrinsic :: iso_c_binding, only: c_int,c_null_char
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QDir_exists,QDir_removeRecursively,QJsonObject,QJsonArray,QJsonValue_String
    use studio_state, only: studio_t,studio_mesh_valid
    use studio_storage_paths, only: studio_project_directory
    use studio_platform_files, only: file_metadata,current_directory,open_directory,open_regular,close_file, &
        is_regular_file,inspect_entry,inspect_file,same_file,copy_open_file,publish_file,temporary_private_directory
    use studio_exchange, only: studio_exchange_suffix,export_studio_mesh
    use studio_project_file, only: save_project_document
    use studio_fortran_export, only: export_project_fortran
    use studio_document, only: studio_document_t,document_has_changes,document_display_name
    use studio_live_viewport, only: save_viewport_image
    use studio_ui, only: import_document_geometry,export_geometry_selection, &
        load_project_workspace,mark_project_saved,new_studio_project,import_studio_project,activate_project_tab,close_project_tab,open_project_tab,start_project_tab
    use studio_units, only: unitless, unit_mm, unit_ft, unit_symbol
    use studio_mcp_arguments, only: only_fields,bounded_integer,optional_boolean
    implicit none(type, external)
    private
    public :: import_mcp_file,export_mcp_file
    public :: open_mcp_project,save_mcp_project,new_mcp_project
    public :: list_mcp_project_tabs,switch_mcp_project_tab,close_mcp_project_tab

    integer,parameter :: maximum_path=1024
    character(7),parameter :: geometry_names(*)=[character(7)::'curve','surface','volume']

contains

    subroutine new_mcp_project(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        character(:),allocatable :: folder
        logical :: new_tab

        error='New project accepts only the optional Boolean new_tab.'
        if(.not.only_fields(args,'|new_tab|').or..not.optional_boolean(args,'new_tab'))return
        new_tab=args%boolValue('new_tab')
        error='Assign an existing project folder before starting a new project.'
        if(.not.allocated(s%document%project_directory))return
        folder=s%document%project_directory
        if(len_trim(folder)==0.or..not.QDir_exists(folder))return
        error='Save the current project before starting a new one.'
        if(.not.new_tab.and.document_has_changes(s%document))return
        if(new_tab)call start_project_tab(s,assistant_folder=folder)
        call new_studio_project(s,folder)
        call result%insert('new_project',.true.)
        if(new_tab)call result%insert('tab',s%project_tabs%active)
        error=''
    end subroutine new_mcp_project

    subroutine open_mcp_project(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        character(:),allocatable :: path,leaf,folder,staging
        logical :: removed,new_tab

        path=file_argument(args,error,allow_new_tab=.true.)
        if(len(error)>0)return
        error='new_tab must be a Boolean.'
        if(.not.optional_boolean(args,'new_tab'))return
        new_tab=args%boolValue('new_tab')
        error='Choose a .fcad project file.'
        if(studio_exchange_suffix(path)/='.fcad')return
        error='Save the current project before opening another project.'
        if(.not.new_tab.and.document_has_changes(s%document))return

        call stage_project_input(s,path,'source.fcad',staging,leaf,error)
        if(len(error)>0)return
        folder=s%document%project_directory
        if(new_tab)then
            call open_project_tab(s,staging//'/source.fcad',error, &
                remembered_path=project_prefix(folder)//path,project_directory=folder)
        else
            call load_project_workspace(s,staging//'/source.fcad',error, &
                remembered_path=project_prefix(folder)//path,project_directory=folder)
        end if
        removed=QDir_removeRecursively(staging)
        if(len(error)>0)return
        call result%insert('path',path)
        if(new_tab)call result%insert('tab',s%project_tabs%active)
    end subroutine open_mcp_project

    subroutine save_mcp_project(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        character(:),allocatable :: path,destination,current,leaf,source
        integer(c_int) :: parent,ignored,status
        type(file_metadata) :: metadata

        error='Open or create a ForCAD Studio project before saving.'
        if(.not.s%document%document_started)return
        current=''
        if(allocated(s%document%project_path))current=s%document%project_path
        if(args%isEmpty())then
            path=current_project_relative_path(s,current,error)
        else
            path=file_argument(args,error)
        end if
        if(len(error)>0)return
        error='Choose a .fcad project file.'
        if(studio_exchange_suffix(path)/='.fcad')return

        call open_project_parent(s,path,parent,leaf,error)
        if(len(error)>0)return
        destination=project_prefix(s%document%project_directory)//path
        status=inspect_entry(parent,leaf//c_null_char,metadata)
        if(status/=0)then
            call publish_staged_export(s,parent,leaf,'.fcad',source,error)
        else if(destination/=current)then
            error='The save destination already exists. Choose a new filename.'
        else if(.not.is_regular_file(metadata))then
            error='The current project path must point to a regular file inside the project folder.'
        else
            call publish_staged_export(s,parent,leaf,'.fcad',source,error,replace=.true.)
        end if
        ignored=close_file(parent)
        if(len(error)>0)return

        s%document%project_path=destination
        call mark_project_saved(s)
        call s%status%showMessage('Saved '//document_display_name(s%document),3500)
        call result%insert('path',path)
    end subroutine save_mcp_project

    subroutine list_mcp_project_tabs(s,args,result,error)
        type(studio_t),intent(in),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(QJsonArray) :: tabs
        type(QJsonObject) :: entry
        type(studio_document_t),pointer :: document
        character(:),allocatable :: path,path_error
        integer :: i

        error='List project tabs takes an empty object.'
        if(.not.args%isEmpty())return
        tabs=QJsonArray()
        if(allocated(s%project_tabs%items))then
            do i=1,size(s%project_tabs%items)
                document=>scoped_project_tab(s,i)
                if(.not.associated(document))cycle
                entry=QJsonObject()
                call entry%insert('tab',i)
                call entry%insert('revision',project_tab_revision(document%document_revision))
                call entry%insert('title',document_display_name(document))
                call entry%insert('active',i==s%project_tabs%active)
                call entry%insert('dirty',document_has_changes(document))
                if(allocated(document%project_path))then
                    path=current_project_relative_path(s,document%project_path,path_error)
                    if(len(path_error)==0.and.index('/'//path//'/', '/../')==0.and.index('/'//path//'/', '/./')==0)then
                        call entry%insert('path',path)
                    end if
                end if
                call tabs%append(entry)
                call entry%delete()
            end do
        end if
        call result%insert('tabs',tabs)
        call tabs%delete()
        error=''
    end subroutine list_mcp_project_tabs

    subroutine switch_mcp_project_tab(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(studio_document_t),pointer :: document
        integer :: tab

        call require_project_tab(s,args,tab,document,error)
        if(len(error)>0)return
        nullify(document)
        call activate_project_tab(s,tab,keep_assistant=.true.)
        error='Could not activate the project tab.'
        if(s%project_tabs%active/=tab)return
        call result%insert('tab',tab)
        call result%insert('revision',project_tab_revision(s%document%document_revision))
        call result%insert('title',document_display_name(s%document))
        error=''
    end subroutine switch_mcp_project_tab

    subroutine close_mcp_project_tab(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(studio_document_t),pointer :: document
        integer :: tab,count

        call require_project_tab(s,args,tab,document,error)
        if(len(error)>0)return
        error='Switch to another project tab before closing this one.'
        if(tab==s%project_tabs%active)return
        error='Save this project before closing its tab.'
        if(document_has_changes(document))return
        nullify(document)
        count=size(s%project_tabs%items)
        call close_project_tab(s,tab)
        error='Could not close the project tab.'
        if(size(s%project_tabs%items)/=count-1)return
        call result%insert('closed_tab',tab)
        error=''
    end subroutine close_mcp_project_tab

    subroutine require_project_tab(s,args,tab,document,error)
        type(studio_t),intent(in),target :: s
        type(QJsonObject),intent(in) :: args
        integer,intent(out) :: tab
        type(studio_document_t),pointer,intent(out) :: document
        character(:),allocatable,intent(out) :: error

        nullify(document)
        error='Supply a tab and revision from list_project_tabs.'
        if(.not.only_fields(args,'|tab|revision|').or.args%valueType('revision')/=QJsonValue_String)return
        tab=bounded_integer(args,'tab',1,huge(1))
        if(tab<1)return
        document=>scoped_project_tab(s,tab)
        error='That project tab is outside the active project folder or unavailable.'
        if(.not.associated(document))return
        error='The project tabs changed. List them again before proceeding.'
        if(args%stringValue('revision')/=project_tab_revision(document%document_revision))return
        error=''
    end subroutine require_project_tab

    function scoped_project_tab(s,index) result(document)
        type(studio_t),intent(in),target :: s
        integer,intent(in) :: index
        type(studio_document_t),pointer :: document

        nullify(document)
        if(.not.allocated(s%document%project_directory).or..not.allocated(s%project_tabs%items))return
        if(len_trim(s%document%project_directory)==0)return
        if(index<1.or.index>size(s%project_tabs%items))return
        if(index==s%project_tabs%active)then
            document=>s%document
            return
        end if
        if(.not.allocated(s%project_tabs%items(index)%document))return
        associate(other=>s%project_tabs%items(index)%document)
            if(.not.allocated(other%project_directory))return
            if(other%project_directory==s%document%project_directory)document=>other
        end associate
    end function scoped_project_tab

    function project_tab_revision(revision) result(text)
        integer(int64),intent(in) :: revision
        character(:),allocatable :: text
        character(24) :: buffer

        write(buffer,'(i0)')revision
        text=trim(buffer)
    end function project_tab_revision

    function current_project_relative_path(s,current,error) result(path)
        type(studio_t),intent(in) :: s
        character(*),intent(in) :: current
        character(:),allocatable,intent(out) :: error
        character(:),allocatable :: path,prefix

        path=''
        error='Supply a new project-relative .fcad path to save this project.'
        if(len(current)==0.or..not.allocated(s%document%project_directory))return
        if(len_trim(s%document%project_directory)==0)return
        prefix=project_prefix(s%document%project_directory)
        error='The current file is outside the project folder. Supply a new project-relative .fcad path.'
        if(index(current,prefix)/=1.or.len(current)<=len(prefix))return
        path=current(len(prefix)+1:)
        if(.not.valid_file_path(path))return
        error=''
    end function current_project_relative_path

    function project_prefix(folder) result(prefix)
        character(*),intent(in) :: folder
        character(:),allocatable :: prefix
        integer :: last

        last=len_trim(folder)
        do while(last>0)
            if(folder(last:last)/='/')exit
            last=last-1
        end do
        prefix=folder(:last)//'/'
    end function project_prefix

    subroutine import_mcp_file(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        character(:),allocatable :: path,leaf,suffix,staging
        integer :: kind,variant,count,skipped,source_unit,i
        logical :: needs_source_unit,removed

        path=file_argument(args,error,allow_source_unit=.true.)
        if(len(error)>0)return
        suffix=studio_exchange_suffix(path)
        if(suffix=='.fcad')then
            if(args%contains('source_unit'))then
                error='Source units apply only to IGES imports. FCAD projects store their units.'
                return
            end if
            call stage_project_input(s,path,'source.fcad',staging,leaf,error,managed_area='references')
            if(len(error)>0)return
            call import_studio_project(s,staging//'/source.fcad',error,source_name=leaf,imported=count)
            removed=QDir_removeRecursively(staging)
            if(len(error)>0)return
            call result%insert('path',path)
            call result%insert('format','fcad')
            call result%insert('imported',count)
            return
        end if
        if(suffix/='.igs'.and.suffix/='.iges')then
            error='Use an .fcad, .igs or .iges file for geometry import.'
            return
        end if

        source_unit=-1
        if(args%contains('source_unit'))then
            error='Choose mm, cm, m, in, ft, or keep for the IGES source unit.'
            if(args%valueType('source_unit')/=QJsonValue_String)return
            leaf=args%stringValue('source_unit')
            if(leaf=='keep')source_unit=unitless
            do i=unit_mm,unit_ft
                if(leaf==unit_symbol(i))source_unit=i
            end do
            if(source_unit<0)return
            error=''
        end if
        call stage_project_input(s,path,'source.igs',staging,leaf,error,managed_area='references')
        if(len(error)>0)return

        needs_source_unit=.false.
        if(source_unit>=0)then
            call import_document_geometry(s,staging//'/source.igs',kind,variant,count,skipped,error,source_unit=source_unit)
        else
            call import_document_geometry(s,staging//'/source.igs',kind,variant,count,skipped,error,needs_source_unit=needs_source_unit)
        end if
        removed=QDir_removeRecursively(staging)
        if(needs_source_unit)call result%insert('needs_source_unit',.true.)
        if(len(error)>0)return

        call result%insert('path',path)
        call result%insert('format','iges')
        call result%insert('imported',count)
        call result%insert('skipped',skipped)
        call result%insert('selected_kind',trim(geometry_names(kind)))
        call result%insert('selected_id',variant)
    end subroutine import_mcp_file

    subroutine export_mcp_file(s,args,result,error)
        type(studio_t),intent(inout) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        character(*),parameter :: screenshot_prefix='.forcad/screenshots/'
        character(:),allocatable :: path,leaf,suffix,source,folder
        integer(c_int) :: parent,ignored
        type(file_metadata) :: metadata
        integer(c_int) :: status

        path=file_argument(args,error)
        if(len(error)>0)return
        suffix=studio_exchange_suffix(path)
        select case(suffix)
        case('.igs','.iges','.vtk','.fcad','.f90','.png')
        case default
            error='Use .igs, .iges, .vtk, .fcad, .f90 or .png for the export destination.'
            return
        end select
        if(suffix=='.png'.and.index(path,screenshot_prefix)==1)then
            error='Assign a project folder before taking a screenshot.'
            if(.not.allocated(s%document%project_directory))return
            folder=studio_project_directory(s%document%project_directory,'screenshots')
            error='Could not open the project screenshot folder.'
            if(len(folder)==0)return
        end if
        if(suffix=='.png')then
            call open_project_parent(s,path,parent,leaf,error,managed_area='screenshots')
        else
            call open_project_parent(s,path,parent,leaf,error)
        end if
        if(len(error)>0)return

        status=inspect_entry(parent,leaf//c_null_char,metadata)
        if(status==0)then
            ignored=close_file(parent)
            error='The export destination already exists. Choose a new filename.'
            return
        end if

        call publish_staged_export(s,parent,leaf,suffix,source,error)
        ignored=close_file(parent)
        if(len(error)>0)return

        if(suffix=='.igs')suffix='.iges'
        call result%insert('path',path)
        call result%insert('format',suffix(2:))
        call result%insert('source',source)
        call s%status%showMessage('Exported '//path,6000)
    end subroutine export_mcp_file

    subroutine stage_project_input(s,path,stage_name,staging,leaf,error,managed_area)
        type(studio_t),intent(in) :: s
        character(*),intent(in) :: path,stage_name
        character(*),intent(in),optional :: managed_area
        character(:),allocatable,intent(out) :: staging,leaf,error
        type(file_metadata) :: metadata
        integer(c_int) :: parent,handle,ignored,status

        call open_project_parent(s,path,parent,leaf,error,managed_area)
        if(len(error)>0)return
        handle=open_regular(parent,leaf//c_null_char)
        ignored=close_file(parent)
        error='Could not open a regular file inside the project folder.'
        if(handle<0)return
        status=inspect_file(handle,metadata)
        if(status==0)then
            if(metadata%links/=1)status=-1
        end if
        if(status/=0)then
            ignored=close_file(handle)
            error='ForCAD Studio can only import or open files with one filesystem link.'
            return
        end if
        staging=stage_open_file(handle,stage_name,s%document%project_directory)
        ignored=close_file(handle)
        if(len(staging)==0)then
            error='Could not stage the input file for ForCAD Studio.'
            return
        end if
        error=''
    end subroutine stage_project_input

    function stage_open_file(handle,leaf,project) result(staging)
        integer(c_int),intent(in) :: handle
        character(*),intent(in) :: leaf,project
        character(:),allocatable :: staging
        logical :: removed

        staging=temporary_private_directory(studio_project_directory(project,'work'),'read.')
        if(len(staging)==0)return
        if(copy_open_file(handle,staging//'/'//leaf//c_null_char)==0)return
        removed=QDir_removeRecursively(staging)
        staging=''
    end function stage_open_file

    subroutine publish_staged_export(s,parent,leaf,suffix,source,error,replace)
        type(studio_t),intent(inout) :: s
        integer(c_int),intent(in) :: parent
        character(*),intent(in) :: leaf,suffix
        character(:),allocatable,intent(out) :: source,error
        logical,intent(in),optional :: replace
        character(:),allocatable :: output,staging
        integer(c_int) :: replacement
        logical :: removed

        staging=temporary_private_directory(studio_project_directory(s%document%project_directory,'work'),'export.')
        error='Could not prepare a private export workspace.'
        if(len(staging)==0)return
        output=staging//'/export'//suffix
        call write_staged_export(s,suffix,output,source,error)
        if(len(error)==0)then
            replacement=0_c_int
            if(present(replace))replacement=merge(1_c_int,0_c_int,replace)
            if(publish_file(output//c_null_char,parent,leaf//c_null_char,replacement)/=0)then
                error='Could not publish the export. Check the filename and project folder permissions.'
            end if
        end if
        removed=QDir_removeRecursively(staging)
    end subroutine publish_staged_export

    subroutine write_staged_export(s,suffix,output,source,error)
        type(studio_t),intent(inout) :: s
        character(*),intent(in) :: suffix,output
        character(:),allocatable,intent(out) :: source,error
        integer :: mesh_id

        if(suffix=='.fcad')then
            source='project'
            call save_project_document(s%document,output,error)
            return
        end if
        if(suffix=='.f90')then
            source='project'
            call export_project_fortran(s%document,output,error)
            return
        end if
        if(suffix=='.png')then
            source='viewport'
            error=''
            if(.not.save_viewport_image(s,output))error='Could not save the 3D viewport image.'
            return
        end if
        mesh_id=s%selected_mesh
        if(mesh_id>0)then
            source='mesh'
            if(suffix/='.vtk')then
                error='Use .vtk to export a generated mesh.'
            else if(.not.studio_mesh_valid(s,mesh_id))then
                error='The selected mesh is out of date. Generate a new mesh before exporting.'
            else
                call export_studio_mesh(s%document%meshes(mesh_id),output,error)
            end if
            return
        end if
        source='selection'
        call export_geometry_selection(s,output,error)
    end subroutine write_staged_export

    function file_argument(args,error,allow_source_unit,allow_new_tab) result(path)
        type(QJsonObject),intent(in) :: args
        character(:),allocatable,intent(out) :: error
        character(:),allocatable :: path
        logical,optional,intent(in) :: allow_source_unit,allow_new_tab
        character(:),allocatable :: fields

        path=''
        error='Supply a relative path inside the active project folder.'
        fields='|path|'
        if(present(allow_source_unit))then
            if(allow_source_unit)fields='|path|source_unit|'
        end if
        if(present(allow_new_tab))then
            if(allow_new_tab)fields=fields//'new_tab|'
        end if
        if(.not.only_fields(args,fields).or.args%valueType('path')/=QJsonValue_String)return
        path=args%stringValue('path')
        if(valid_file_path(path))error=''
    end function file_argument

    logical function valid_file_path(path) result(valid)
        character(*),intent(in) :: path
        integer :: i,code

        valid=.false.
        if(len(path)<1.or.len(path)>maximum_path.or.len_trim(path)/=len(path))return
        if(path(1:1)=='/'.or.path(len(path):)=='/'.or.scan(path,achar(92)//':')>0)return
        do i=1,len(path)
            code=iachar(path(i:i))
            if(code<32.or.code==127)return
        end do
        valid=.true.
    end function valid_file_path

    subroutine open_project_parent(s,path,parent,leaf,error,managed_area)
        type(studio_t),intent(in) :: s
        character(*),intent(in) :: path
        character(*),intent(in),optional :: managed_area
        integer(c_int),intent(out) :: parent
        character(:),allocatable,intent(out) :: leaf,error
        character(:),allocatable :: component
        integer(c_int) :: child,ignored
        integer :: first,next

        leaf=''
        parent=-1
        error='Assign a project folder before importing or exporting files.'
        if(.not.allocated(s%document%project_directory))return
        if(len_trim(s%document%project_directory)==0)return
        parent=open_directory(current_directory,trim(s%document%project_directory)//c_null_char)
        if(parent<0)return

        error='Choose an existing directory inside the project folder and use a path without links or traversal.'
        first=1
        do
            next=index(path(first:),'/')
            if(next==0)then
                leaf=path(first:)
                if(len(leaf)>0.and.leaf/='.'.and.leaf/='..')then
                    error=''
                    return
                end if
                exit
            end if
            component=path(first:first+next-2)
            if(len(component)==0.or.component=='.'.or.component=='..')exit
            child=open_directory(parent,component//c_null_char)
            if(child<0)exit
            if(private_exchange_directory(parent,child))then
                if(.not.allowed_managed_area(path,first,managed_area))then
                    ignored=close_file(child)
                    error='Choose a location outside ForCAD Studio’s private storage.'
                    exit
                end if
            end if
            ignored=close_file(parent)
            parent=child
            first=first+next
        end do
        ignored=close_file(parent)
        parent=-1
    end subroutine open_project_parent

    logical function private_exchange_directory(parent,child) result(private)
        integer(c_int),intent(in) :: parent,child
        type(file_metadata) :: directory,entry

        private=.true.
        if(inspect_file(child,directory)/=0)return
        if(inspect_entry(parent,'.forcad'//c_null_char,entry)==0)then
            if(same_file(directory,entry))return
        end if
        if(inspect_entry(parent,'.forcad-studio'//c_null_char,entry)==0)then
            if(same_file(directory,entry))return
        end if
        private=.false.
    end function private_exchange_directory

    logical function allowed_managed_area(path,first,area) result(allowed)
        character(*),intent(in) :: path
        integer,intent(in) :: first
        character(*),intent(in),optional :: area

        allowed=.false.
        if(first/=1.or..not.present(area))return
        allowed=index(path,'.forcad/'//area//'/')==1
    end function allowed_managed_area

end module studio_mcp_exchange
