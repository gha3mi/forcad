module studio_project_file
    use, intrinsic :: iso_c_binding, only: c_null_char,c_double,c_int
    use, intrinsic :: iso_fortran_env, only: int64,iostat_end
    use studio_forcad, only: rk
    use studio_patch_stream, only: write_patch,read_patch,indexable_shape
    use studio_single_patch, only: studio_patch_catalog,allocate_studio_patch_catalog,stage_studio_patch, &
        install_studio_patch_catalog,studio_patch_is_staged,studio_patch_resolution,single_patch_count, &
        single_patch_available,active_patch_id,single_patch_name,copy_single_patch_geometry,single_patch_revision
    use studio_scene, only: studio_scene_t,frame_studio_scene,clear_material_history,move_studio_scene
    use studio_multipatch, only: multipatch_draft
    use studio_multipatch_stream, only: write_group,read_group
    use studio_multipatch_collection, only: multipatch_collection,count_multipatches,restore_multipatch_definitions,clear_multipatch_history
    use studio_units, only: valid_unit_id,unitless
    use studio_document, only: studio_document_t,reset_document_session,mark_document_saved,set_document_project_directory
    use studio_mesh, only: studio_mesh_t
    use studio_project_mesh, only: write_project_meshes,read_project_meshes,rebind_project_meshes
    use studio_project_stream, only: write_text,read_text,remaining_file_units,integer_units
    use studio_project_folder, only: project_file_folder
    use qt, only: QDir_removeRecursively
    use studio_platform_files, only: temporary_private_directory,replace_file,last_path_separator, &
        create_private_file,current_directory,close_file
    implicit none(type, external)
    private
    public :: save_studio_project,open_studio_project,reset_studio_project
    public :: save_project_document,open_project_document,reset_project_document
    public :: write_patch,read_patch

    character(*), parameter :: signature='ForCAD Studio project'
    integer, parameter :: version=2,endian_marker=int(z'01020304')

contains
    subroutine reset_project_document(document,include_templates)
        type(studio_document_t), intent(inout) :: document
        logical, intent(in) :: include_templates

        call reset_studio_project(document%scene,document%multipatches,catalog=document%patches)
        if(allocated(document%project_path))deallocate(document%project_path)
        call set_document_project_directory(document)
        document%show_templates=include_templates
        document%model_unit=unitless
        document%display_unit=unitless
        call reset_document_session(document)
    end subroutine reset_project_document

    subroutine open_project_document(document,path,error)
        type(studio_document_t), intent(inout) :: document
        character(*), intent(in) :: path
        character(:), allocatable, intent(out) :: error
        type(studio_mesh_t), allocatable :: meshes(:)

        call open_studio_project(path,document%scene,document%multipatches,document%model_unit,document%display_unit, &
            error,document%show_templates,catalog=document%patches,meshes=meshes)
        if(len(error)/=0)return
        call reset_document_session(document)
        call move_alloc(meshes,document%meshes)
        call remember_project_path(document,path)
    end subroutine open_project_document

    subroutine save_project_document(document,path,error,remember_path)
        type(studio_document_t), intent(inout) :: document
        character(*), intent(in) :: path
        character(:), allocatable, intent(out) :: error
        logical, intent(in), optional :: remember_path

        call save_studio_project(path,document%scene,document%multipatches,document%model_unit,document%display_unit, &
            error,document%show_templates,catalog=document%patches,meshes=document%meshes)
        if(len(error)/=0.or..not.present(remember_path))return
        if(.not.remember_path)return
        if(allocated(document%project_path).and.allocated(document%project_directory))then
            if(document%project_path==path)then
                call mark_document_saved(document)
                return
            end if
        end if
        call remember_project_path(document,path)
    end subroutine save_project_document

    subroutine remember_project_path(document,path)
        type(studio_document_t), intent(inout) :: document
        character(*), intent(in) :: path

        document%project_path=path
        call set_document_project_directory(document,project_file_folder(path))
        call mark_document_saved(document)
    end subroutine remember_project_path

    subroutine reset_studio_project(scene,collection,catalog)
        type(studio_scene_t), intent(inout) :: scene
        type(multipatch_collection), intent(inout) :: collection
        type(studio_patch_catalog), intent(inout), optional :: catalog
        type(studio_patch_catalog) :: empty
        type(studio_scene_t) :: empty_scene
        type(multipatch_collection) :: empty_collection
        integer :: kind

        call clear_material_history(scene)
        call clear_multipatch_history(collection)
        scene=empty_scene
        collection=empty_collection
        call install_studio_patch_catalog(empty,catalog=catalog)
        do kind=1,3
            ! Keep the internal templates hidden until the user chooses one.
            allocate(scene%catalog(kind)%entries(single_patch_count(kind,catalog=catalog)))
        end do
    end subroutine reset_studio_project

    subroutine save_studio_project(path,scene,collection,units_id,display_unit,error,include_templates,catalog,meshes)
        character(*), intent(in) :: path
        type(studio_scene_t), intent(in) :: scene
        type(multipatch_collection), intent(in) :: collection
        integer, intent(in) :: units_id,display_unit
        character(:), allocatable, intent(out) :: error
        logical, intent(in), optional :: include_templates
        type(studio_patch_catalog), intent(in), optional :: catalog
        type(studio_mesh_t), allocatable, intent(in), optional :: meshes(:)
        character(:),allocatable :: staging,temporary,parent
        integer :: unit,status,closed,templates,separator
        integer(c_int) :: file
        logical :: removed

        error='Could not create a temporary file next to the project.'
        if(len_trim(path)==0.or.index(path,c_null_char)/=0)return
        templates=0
        if(present(include_templates))templates=merge(1,0,include_templates)
        separator=last_path_separator(path)
        parent='.'
        if(separator>0)parent=path(:separator)
        staging=temporary_private_directory(parent,'.fcad-save.')
        if(len(staging)==0)return
        temporary=staging//'/project.fcad'
        file=create_private_file(current_directory,temporary//c_null_char)
        status=1
        if(file>=0)status=close_file(file)
        if(status==0)then
            open(newunit=unit,file=temporary,status='old',action='write', &
                access='stream',form='unformatted',iostat=status)
        end if
        if(status==0)then
            error='Could not write this project.'
            call write_document(unit,scene,collection,units_id,display_unit,templates,status,catalog)
            if(status==0)call write_project_meshes(unit,meshes,catalog,status)
            close(unit,iostat=closed)
            if(status==0.and.closed==0)then
                error='Could not replace the destination project file.'
                if(replace_file(temporary//c_null_char,path//c_null_char)==0)error=''
            end if
        end if
        removed=QDir_removeRecursively(staging)
    end subroutine save_studio_project

    subroutine write_document(unit,scene,collection,units_id,display_unit,templates,status,catalog)
        integer, intent(in) :: unit,units_id,display_unit,templates
        type(studio_scene_t), intent(in) :: scene
        type(multipatch_collection), intent(in) :: collection
        integer, intent(out) :: status
        type(studio_patch_catalog), intent(in), optional :: catalog
        integer :: counts(3),selection(3),kind,variant,index,present_flag,material(5)
        class(*), allocatable :: geometry
        character(:), allocatable :: name

        status=1
        if(.not.valid_unit_id(units_id).or..not.valid_unit_id(display_unit))return
        counts=[(single_patch_count(kind,catalog=catalog),kind=1,3)]
        selection=[(active_patch_id(kind,catalog=catalog),kind=1,3)]
        write(unit,iostat=status)signature,version,endian_marker,storage_size(0),storage_size(0._rk), &
            storage_size(0._c_double),units_id,display_unit,templates,counts,selection,count_multipatches(collection)
        if(status/=0)return
        do kind=1,3
            do variant=1,counts(kind)
                present_flag=merge(1,0,single_patch_available(kind,variant,catalog=catalog))
                write(unit,iostat=status)present_flag
                if(status/=0)return
                if(present_flag==0)cycle
                call copy_single_patch_geometry(kind,variant,geometry,catalog=catalog)
                status=1
                if(.not.allocated(geometry))return
                name=single_patch_name(kind,variant,catalog=catalog)
                material=[0,-1,-1,-1,255]
                if(allocated(scene%catalog(kind)%entries))then
                    if(variant<=size(scene%catalog(kind)%entries))then
                        associate(entry=>scene%catalog(kind)%entries(variant))
                            material=[merge(1,0,entry%visible),entry%rgb,entry%opacity]
                        end associate
                    end if
                end if
                call write_text(unit,name,status)
                if(status/=0)return
                write(unit,iostat=status)material,studio_patch_resolution(kind,variant,catalog=catalog)
                if(status/=0)return
                call write_patch(unit,kind,geometry,status)
                if(status/=0)return
            end do
        end do
        do index=1,count_multipatches(collection)
            call write_group(unit,collection%groups(index)%definition,status)
            if(status/=0)return
        end do
    end subroutine write_document

    subroutine open_studio_project(path,scene,collection,units_id,display_unit,error,include_templates,catalog,meshes)
        character(*), intent(in) :: path
        type(studio_scene_t), intent(inout) :: scene
        type(multipatch_collection), intent(inout) :: collection
        integer, intent(inout) :: units_id,display_unit
        character(:), allocatable, intent(out) :: error
        logical, intent(inout), optional :: include_templates
        type(studio_patch_catalog), intent(inout), optional :: catalog
        type(studio_mesh_t), allocatable, intent(inout), optional :: meshes(:)
        type(studio_mesh_t), allocatable :: stored_meshes(:)
        type(studio_patch_catalog) :: staged
        type(studio_scene_t) :: replacement
        type(multipatch_draft), allocatable :: definitions(:)
        integer :: unit,status,closed,stored_units,stored_display,templates,kind,variant
        integer(int64) :: length
        character :: extra

        error='Could not read the selected project file.'
        if(len_trim(path)==0.or.index(path,c_null_char)/=0)return
        length=0
        inquire(file=path,size=length,iostat=status)
        if(status/=0)return
        if(length<=0)return
        open(newunit=unit,file=path,status='old',action='read',access='stream',form='unformatted',iostat=status)
        if(status/=0)return
        call read_document(unit,staged,replacement,definitions,stored_units,stored_display,templates,error)
        if(len(error)==0)then
            call read_project_meshes(unit,staged,stored_meshes,status)
            if(status/=0)error='The project mesh data is incomplete or invalid.'
        end if
        if(len(error)==0)then
            read(unit,iostat=status)extra
            if(status==0)then
                error='Unexpected data at the end of the project file.'
            else if(status/=iostat_end)then
                error='Could not finish reading the project file.'
            end if
        end if
        close(unit,iostat=closed)
        if(len(error)>0)return
        if(closed/=0)then
            error='Could not close the project file.'
            return
        end if
        call install_studio_patch_catalog(staged,catalog=catalog)
        if(allocated(stored_meshes))call rebind_project_meshes(stored_meshes,catalog)
        if(present(meshes))call move_alloc(stored_meshes,meshes)
        call clear_multipatch_history(collection)
        call restore_multipatch_definitions(collection,definitions,catalog=catalog)
        do kind=1,size(replacement%catalog)
            do variant=1,size(replacement%catalog(kind)%entries)
                associate(entry=>replacement%catalog(kind)%entries(variant))
                    if(entry%patch%valid)entry%revision=single_patch_revision(kind,variant,catalog)
                end associate
            end do
        end do
        call frame_studio_scene(replacement)
        call clear_material_history(scene)
        call move_studio_scene(replacement,scene)
        units_id=stored_units
        display_unit=stored_display
        if(present(include_templates))include_templates=templates==1
    end subroutine open_studio_project

    subroutine read_document(unit,staged,scene,definitions,units_id,display_unit,templates,error)
        integer, intent(in) :: unit
        type(studio_patch_catalog), intent(out) :: staged
        type(studio_scene_t), intent(out) :: scene
        type(multipatch_draft), allocatable, intent(out) :: definitions(:)
        integer, intent(out) :: units_id,display_unit,templates
        character(:), allocatable, intent(out) :: error
        character(len(signature)) :: header
        integer :: file_version,endian,integer_bits,kernel_bits,ui_bits,counts(3),selection(3),groups,kind,variant
        integer :: status,present_flag,material(5),resolution(3),index,member
        character(:), allocatable :: name,patch_error
        class(*), allocatable :: geometry

        error='This file is not a ForCAD Studio project.'
        read(unit,iostat=status)header,file_version,endian,integer_bits,kernel_bits,ui_bits
        if(status/=0)return
        if(header/=signature)return
        error='This ForCAD Studio project version or numeric encoding is unsupported.'
        if(file_version/=version.or.endian/=endian_marker)return
        if(integer_bits/=storage_size(0).or.kernel_bits/=storage_size(0._rk).or.ui_bits/=storage_size(0._c_double))return
        error='The project file is incomplete or contains invalid document data.'
        read(unit,iostat=status)units_id,display_unit,templates,counts,selection,groups
        if(status/=0)return
        if(.not.valid_unit_id(units_id).or..not.valid_unit_id(display_unit))return
        if(any(counts<2).or.any(selection<1).or.any(selection>counts))return
        if(templates<0.or.templates>1)return
        if(groups<0)return
        if(sum(int(counts,int64))+5_int64*groups>remaining_file_units(unit)/integer_units)return
        call allocate_studio_patch_catalog(staged,counts,selection,status)
        if(status/=0)return
        do kind=1,3
            allocate(scene%catalog(kind)%entries(counts(kind)),stat=status)
            if(status/=0)return
            do variant=1,counts(kind)
                read(unit,iostat=status)present_flag
                if(status/=0)return
                if(present_flag==0.and.variant>2)cycle
                if(present_flag/=1)return
                call read_text(unit,name,status)
                if(status/=0)return
                read(unit,iostat=status)material,resolution
                if(status/=0)return
                if(material(1)<0.or.material(1)>1.or.material(5)<0.or.material(5)>255)return
                if(any(material(2:4)>255).or.(any(material(2:4)<0).and..not.all(material(2:4)==-1)))return
                if(any(resolution<0))return
                if(.not.indexable_shape(max(17,resolution(:kind))))return
                call read_patch(unit,kind,geometry,status)
                if(status/=0)return
                associate(entry=>scene%catalog(kind)%entries(variant))
                    entry%visible=material(1)==1
                    entry%rgb=material(2:4)
                    entry%opacity=material(5)
                    if(entry%visible)then
                        call stage_studio_patch(staged,kind,variant,name,geometry,resolution,patch_error,entry%patch)
                    else
                        call stage_studio_patch(staged,kind,variant,name,geometry,resolution,patch_error)
                    end if
                end associate
                if(len(patch_error)>0)then
                    error=patch_error
                    return
                end if
            end do
            if(.not.studio_patch_is_staged(staged,kind,selection(kind)))return
        end do
        allocate(definitions(groups),stat=status)
        if(status/=0)return
        do index=1,groups
            call read_group(unit,definitions(index),status)
            if(status/=0)return
            associate(group=>definitions(index))
                do member=1,size(group%members)
                    variant=group%members(member)
                    if(.not.studio_patch_is_staged(staged,group%kind,variant))return
                    if(any(group%members(:member-1)==variant))return
                end do
            end associate
        end do
        error=''
    end subroutine read_document

end module studio_project_file
