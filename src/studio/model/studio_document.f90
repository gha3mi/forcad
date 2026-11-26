module studio_document
    use, intrinsic :: iso_fortran_env, only: int64,iostat_end
    use studio_single_patch, only: studio_patch_catalog,next_patch_revision,single_patch_mesh_revision,usable_undo_stamp, &
        single_patch_undo_stamp,undo_single_patch,redo_single_patch,studio_patch_replay,release_patch_replay,clear_patch_history, &
        remove_single_patches,move_studio_patch_catalog
    use studio_scene, only: studio_scene_t,project_patch_exists,project_geometry_targets,studio_scene_is_visible,get_studio_material, &
        studio_material_replay,clear_material_history,snapshot_scene_properties,restore_scene_properties,move_studio_scene
    use studio_multipatch_collection, only: multipatch_collection,multipatch_undo_stamp,count_multipatches,multipatch_contains, &
        multipatch_replay,release_multipatch_replay,clear_multipatch_history,append_multipatch_definitions,truncate_multipatches, &
        erase_multipatch_definitions,insert_multipatch_definitions,move_multipatch_collection
    use studio_multipatch, only: multipatch_draft
    use studio_multipatch_stream, only: save_group_history,load_group_history
    use studio_property_stream, only: studio_property_values,save_property_history,load_property_history
    use studio_mesh, only: studio_mesh_t, move_studio_mesh
    use studio_project_mesh, only: write_session_mesh,read_session_mesh
    use studio_history_storage, only: history_blob,create_history_blob,open_history_blob,release_history_blob
    use studio_units, only: unitless
    use studio_platform_files, only: last_path_separator
    implicit none(type, external)
    private

    type, public :: studio_mesh_edit
        integer(int64) :: stamp=0_int64
        integer :: selection=0
        logical :: mesh_only=.false.
        integer, allocatable :: prior_ids(:),published_ids(:)
        type(history_blob) :: payload
    end type studio_mesh_edit

    type, public :: studio_mesh_restore
        type(studio_mesh_t), allocatable, private :: meshes(:)
        integer, allocatable, private :: source_ids(:),target_ids(:)
    end type studio_mesh_restore

    type :: studio_merge_edit
        integer(int64) :: stamp=0_int64, patch_stamp=0_int64
        integer :: assemblies=0,meshes=0
        character(:),allocatable :: operation
    end type studio_merge_edit

    type :: studio_properties_edit
        integer(int64) :: stamp=0
        type(history_blob) :: payload
        character(:), allocatable :: operation
    end type studio_properties_edit

    type :: studio_merge_replay
        type(studio_merge_edit) :: edit
        type(studio_patch_replay) :: patches
        type(studio_mesh_edit) :: meshes
        type(history_blob) :: groups
    end type studio_merge_replay

    type :: studio_deletion_edit
        integer(int64) :: stamp=0_int64,patch_stamp=0_int64
        integer, allocatable :: groups(:)
        type(history_blob) :: definitions
        type(studio_mesh_edit) :: meshes
        logical :: mesh_only=.false.
    end type studio_deletion_edit

    type :: studio_deletion_replay
        type(studio_deletion_edit) :: edit
        type(studio_patch_replay) :: patches
    end type studio_deletion_replay

    type :: saved_patch_state
        integer :: kind,id,rgb(3),opacity
        integer(int64) :: revision
        logical :: visible
        character(:), allocatable :: name
    end type saved_patch_state

    type :: saved_project_state
        type(saved_patch_state), allocatable :: patches(:)
        integer(int64), allocatable :: mesh_revisions(:)
        logical, allocatable :: mesh_visibility(:)
        integer(int64) :: multipatch_stamp
        integer :: units(2)
    end type saved_project_state

    type, public :: studio_redo_entry_t
        type(studio_patch_replay) :: patches
        type(studio_mesh_edit) :: meshes
        type(multipatch_replay) :: assembly
        type(studio_material_replay) :: material
        type(studio_merge_replay) :: merge
        type(studio_deletion_replay) :: deletion
        type(studio_properties_edit) :: properties
        type(history_blob) :: appearance
        type(studio_redo_entry_t), allocatable :: next
        integer(int64) :: undo_stamp=0
        integer :: source=0
    end type studio_redo_entry_t

    type, public :: studio_document_t
        type(studio_patch_catalog) :: patches
        type(studio_scene_t) :: scene
        type(multipatch_collection) :: multipatches
        type(studio_mesh_t), allocatable :: meshes(:)
        type(studio_mesh_edit), allocatable, private :: mesh_history(:)
        integer, private :: mesh_undo_count=0
        type(studio_merge_edit), allocatable, private :: merge_history(:)
        type(studio_deletion_edit), allocatable, private :: deletion_history(:)
        type(studio_properties_edit), allocatable, private :: properties_history(:)
        integer(int64) :: mesh_generation=0
        type(studio_redo_entry_t), allocatable, private :: redo
        character(:), allocatable :: project_path,project_directory,sample_title
        type(saved_project_state), allocatable, private :: saved_project
        integer(int64) :: document_revision=0
        integer :: model_unit=unitless,display_unit=unitless
        logical :: document_started=.false.,show_templates=.true.
    end type studio_document_t

    public :: reset_document_session, mark_document_saved, document_has_changes, document_display_name, move_studio_document
    public :: set_document_project_directory

    public :: document_has_patch, document_patch_in_assembly, document_mesh_count, document_mesh_source_group
    public :: append_document_mesh, record_document_mesh_edit, remove_document_mesh, document_mesh_valid
    public :: document_mesh_undo_stamp,undo_document_mesh_edit,redo_document_mesh_edit,remove_construction_mesh
    public :: capture_document_mesh_edit,stage_document_mesh_edit,apply_document_mesh_edit
    public :: swap_document_mesh_edit,release_document_mesh_edit,clear_document_mesh_history
    public :: remember_document_redo,take_document_redo,release_document_redo,clear_document_history
    public :: discard_stale_redo,can_restore_document_redo
    public :: record_document_merge,document_merge_undo_stamp,undo_document_merge,redo_document_merge
    public :: record_document_properties,document_properties_undo_stamp,undo_document_properties,redo_document_properties
    public :: stage_document_mesh_removal
    public :: commit_document_deletion,document_deletion_undo_stamp,undo_document_deletion,redo_document_deletion

contains

    logical function document_has_patch(document,kind,variant) result(available)
        type(studio_document_t),intent(in)::document
        integer,intent(in)::kind,variant

        available=project_patch_exists(document%show_templates,kind,variant,patches=document%patches)
    end function document_has_patch

    logical function document_patch_in_assembly(document,kind,variant) result(member)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: kind,variant
        integer :: id

        member=.true.
        do id=1,count_multipatches(document%multipatches)
            if(multipatch_contains(document%multipatches,id,kind,variant))return
        end do
        member=.false.
    end function document_patch_in_assembly

    subroutine set_document_project_directory(document,folder)
        type(studio_document_t), intent(inout) :: document
        character(*), intent(in), optional :: folder
        character(:), allocatable :: assigned

        if(present(folder))then
            if(len_trim(folder)>0)assigned=folder
        end if
        call move_alloc(assigned,document%project_directory)
        call document%patches%set_history_project(document%project_directory)
    end subroutine set_document_project_directory

    subroutine reset_document_session(document)
        type(studio_document_t), intent(inout) :: document

        document%document_revision=document%document_revision+1
        call clear_document_history(document)
        if(allocated(document%meshes))deallocate(document%meshes)
        if(allocated(document%sample_title))deallocate(document%sample_title)
        if(allocated(document%saved_project))deallocate(document%saved_project)
    end subroutine reset_document_session

    subroutine move_studio_document(source,destination)
        type(studio_document_t),intent(inout) :: source
        type(studio_document_t),intent(out) :: destination

        call move_studio_patch_catalog(source%patches,destination%patches)
        call move_studio_scene(source%scene,destination%scene)
        call move_multipatch_collection(source%multipatches,destination%multipatches)
        call move_alloc(source%meshes,destination%meshes)
        call move_alloc(source%mesh_history,destination%mesh_history)
        call move_alloc(source%merge_history,destination%merge_history)
        call move_alloc(source%deletion_history,destination%deletion_history)
        call move_alloc(source%properties_history,destination%properties_history)
        destination%mesh_undo_count=source%mesh_undo_count
        destination%mesh_generation=source%mesh_generation
        source%mesh_undo_count=0
        source%mesh_generation=0
        call move_alloc(source%redo,destination%redo)
        call move_alloc(source%sample_title,destination%sample_title)
        call move_alloc(source%project_path,destination%project_path)
        call move_alloc(source%project_directory,destination%project_directory)
        call move_alloc(source%saved_project,destination%saved_project)
        destination%document_revision=source%document_revision
        destination%model_unit=source%model_unit
        destination%display_unit=source%display_unit
        destination%document_started=source%document_started
        destination%show_templates=source%show_templates
        source%model_unit=unitless
        source%display_unit=unitless
        source%document_started=.false.
        source%show_templates=.true.
    end subroutine move_studio_document

    subroutine remember_document_redo(document,entry,undo_stamp)
        type(studio_document_t), intent(inout) :: document
        type(studio_redo_entry_t), allocatable, intent(inout) :: entry
        integer(int64), intent(in) :: undo_stamp

        call move_alloc(document%redo,entry%next)
        entry%undo_stamp=undo_stamp
        call move_alloc(entry,document%redo)
    end subroutine remember_document_redo

    subroutine take_document_redo(document,entry)
        type(studio_document_t), intent(inout) :: document
        type(studio_redo_entry_t), allocatable, intent(out) :: entry

        call move_alloc(document%redo,entry)
        if(allocated(entry))call move_alloc(entry%next,document%redo)
    end subroutine take_document_redo

    logical function can_restore_document_redo(document,undo_stamp) result(available)
        type(studio_document_t), intent(in) :: document
        integer(int64), intent(in) :: undo_stamp

        available=.false.
        if(allocated(document%redo))available=undo_stamp==document%redo%undo_stamp
    end function can_restore_document_redo

    subroutine discard_stale_redo(document,undo_stamp)
        type(studio_document_t), intent(inout) :: document
        integer(int64), intent(in) :: undo_stamp

        if(.not.can_restore_document_redo(document,undo_stamp))call clear_document_redo(document)
    end subroutine discard_stale_redo

    subroutine clear_document_redo(document)
        type(studio_document_t), intent(inout) :: document
        type(studio_redo_entry_t), allocatable :: entry

        do while(allocated(document%redo))
            call take_document_redo(document,entry)
            call release_document_redo(entry)
        end do
    end subroutine clear_document_redo

    subroutine release_document_redo(entry)
        type(studio_redo_entry_t), allocatable, intent(inout) :: entry

        if(.not.allocated(entry))return
        call release_patch_replay(entry%patches)
        call release_document_mesh_edit(entry%meshes)
        call release_multipatch_replay(entry%assembly)
        call release_history_blob(entry%material%payload)
        call release_patch_replay(entry%merge%patches)
        call release_document_mesh_edit(entry%merge%meshes)
        call release_history_blob(entry%merge%groups)
        call release_patch_replay(entry%deletion%patches)
        call release_deletion_edit(entry%deletion%edit)
        call release_history_blob(entry%properties%payload)
        call release_history_blob(entry%appearance)
        deallocate(entry)
    end subroutine release_document_redo

    subroutine clear_document_history(document)
        type(studio_document_t), intent(inout) :: document
        integer :: i

        call clear_document_redo(document)
        call clear_patch_history(document%patches)
        call clear_material_history(document%scene)
        call clear_multipatch_history(document%multipatches)
        call clear_document_mesh_history(document)
        if(allocated(document%merge_history))deallocate(document%merge_history)
        if(allocated(document%deletion_history))then
            do i=1,size(document%deletion_history)
                call release_deletion_edit(document%deletion_history(i))
            end do
            deallocate(document%deletion_history)
        end if
        if(.not.allocated(document%properties_history))return
        do i=1,size(document%properties_history)
            call release_history_blob(document%properties_history(i)%payload)
        end do
        deallocate(document%properties_history)
    end subroutine clear_document_history

    logical function record_document_properties(document,operation,patches,meshes,units) result(saved)
        type(studio_document_t), intent(inout) :: document
        character(*), intent(in) :: operation
        integer, intent(in), optional :: patches(:,:),meshes(:)
        logical, intent(in), optional :: units
        type(studio_properties_edit) :: edit
        type(studio_property_values) :: values

        call collect_properties(document,values,patches,meshes,units)
        saved=save_property_history(values,edit%payload,document%project_directory)
        if(.not.saved)return
        edit%stamp=next_patch_revision()
        edit%operation=operation
        call append_property_edit(document,edit)
        call clear_document_redo(document)
    end function record_document_properties

    subroutine collect_properties(document,values,patches,meshes,units)
        type(studio_document_t), intent(in) :: document
        type(studio_property_values), intent(out) :: values
        integer, intent(in), optional :: patches(:,:),meshes(:)
        logical, intent(in), optional :: units
        integer :: i

        if(present(units))then
            if(units)values%units=[document%model_unit,document%display_unit]
        end if
        if(present(patches))values%patches=snapshot_scene_properties(document%scene,patches)
        if(present(meshes))then
            allocate(values%meshes(2,size(meshes)))
            do i=1,size(meshes)
                values%meshes(:,i)=[meshes(i),merge(1,0,document%meshes(meshes(i))%visible)]
            end do
        end if
    end subroutine collect_properties

    subroutine append_property_edit(document,edit)
        type(studio_document_t), intent(inout) :: document
        type(studio_properties_edit), intent(in) :: edit

        if(.not.allocated(document%properties_history))allocate(document%properties_history(0))
        document%properties_history=[document%properties_history,edit]
    end subroutine append_property_edit

    integer(int64) function document_properties_undo_stamp(document) result(stamp)
        type(studio_document_t), intent(in) :: document
        integer :: n

        stamp=0
        if(.not.allocated(document%properties_history))return
        n=size(document%properties_history)
        if(n>0)stamp=usable_undo_stamp(document%properties_history(n)%stamp,catalog=document%patches)
    end function document_properties_undo_stamp

    subroutine undo_document_properties(document,inverse,changed,message)
        type(studio_document_t), intent(inout) :: document
        type(studio_properties_edit), intent(out) :: inverse
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: message
        type(studio_properties_edit) :: edit
        integer :: n

        n=size(document%properties_history)
        edit=document%properties_history(n)
        call swap_properties(document,edit,inverse,changed,message)
        if(.not.changed)return
        call release_history_blob(edit%payload)
        document%properties_history=document%properties_history(:n-1)
        message=edit%operation//' undone'
    end subroutine undo_document_properties

    subroutine redo_document_properties(document,edit,changed,message)
        type(studio_document_t), intent(inout) :: document
        type(studio_properties_edit), intent(inout) :: edit
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: message
        type(studio_properties_edit) :: inverse

        call swap_properties(document,edit,inverse,changed,message)
        if(.not.changed)return
        call append_property_edit(document,inverse)
        call release_history_blob(edit%payload)
        message=edit%operation//' redone'
    end subroutine redo_document_properties

    subroutine swap_properties(document,edit,inverse,changed,message)
        type(studio_document_t), intent(inout) :: document
        type(studio_properties_edit), intent(in) :: edit
        type(studio_properties_edit), intent(out) :: inverse
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: message
        type(studio_property_values) :: values,current
        integer :: i

        changed=.false.
        message='Could not read the project property history.'
        if(.not.load_property_history(edit%payload,values))return
        if(any(values%meshes(1,:)<1).or.any(values%meshes(1,:)>document_mesh_count(document)))return
        call collect_properties(document,current,values%patches(:2,:),values%meshes(1,:),all(values%units>=0))
        message='Could not write the project property history.'
        if(.not.save_property_history(current,inverse%payload,document%project_directory))return
        inverse%stamp=edit%stamp
        inverse%operation=edit%operation
        if(all(values%units>=0))then
            document%model_unit=values%units(1)
            document%display_unit=values%units(2)
        end if
        call restore_scene_properties(document%scene,values%patches,document%patches)
        do i=1,size(values%meshes,2)
            document%meshes(values%meshes(1,i))%visible=values%meshes(2,i)==1
        end do
        changed=.true.
        message=''
    end subroutine swap_properties

    subroutine mark_document_saved(document)
        type(studio_document_t), intent(inout) :: document
        integer, allocatable :: targets(:,:)
        integer :: i,meshes

        if(allocated(document%saved_project))deallocate(document%saved_project)
        allocate(document%saved_project)
        targets=project_geometry_targets(document%show_templates,document%patches)
        meshes=document_mesh_count(document)
        associate(saved=>document%saved_project)
            saved%units=[document%model_unit,document%display_unit]
            saved%multipatch_stamp=multipatch_undo_stamp(document%multipatches,document%patches)
            allocate(saved%patches(size(targets,2)),saved%mesh_revisions(meshes),saved%mesh_visibility(meshes))
            do i=1,size(targets,2)
                call read_patch_save_state(document,targets(1,i),targets(2,i),saved%patches(i))
            end do
            do i=1,meshes
                saved%mesh_revisions(i)=document%meshes(i)%content_revision
                saved%mesh_visibility(i)=document%meshes(i)%visible
            end do
        end associate
    end subroutine mark_document_saved

    logical function document_has_changes(document) result(changed)
        type(studio_document_t), intent(in) :: document
        type(saved_patch_state) :: current
        integer, allocatable :: targets(:,:)
        integer :: i,meshes

        changed=.false.
        if(.not.allocated(document%saved_project))return
        changed=.true.
        associate(saved=>document%saved_project)
            if(any(saved%units/=[document%model_unit,document%display_unit]))return
            if(saved%multipatch_stamp/=multipatch_undo_stamp(document%multipatches,document%patches))return
            meshes=document_mesh_count(document)
            if(meshes/=size(saved%mesh_revisions))return
            do i=1,meshes
                if(saved%mesh_revisions(i)/=document%meshes(i)%content_revision)return
                if(saved%mesh_visibility(i).neqv.document%meshes(i)%visible)return
            end do

            targets=project_geometry_targets(document%show_templates,document%patches)
            if(size(targets,2)/=size(saved%patches))return
            do i=1,size(targets,2)
                call read_patch_save_state(document,targets(1,i),targets(2,i),current)
                associate(previous=>saved%patches(i))
                    if(current%kind/=previous%kind.or.current%id/=previous%id.or.current%revision/=previous%revision)return
                    if(current%visible.neqv.previous%visible)return
                    if(current%opacity/=previous%opacity.or.any(current%rgb/=previous%rgb))return
                    if(len(current%name)/=len(previous%name).or.current%name/=previous%name)return
                end associate
            end do
        end associate
        changed=.false.
    end function document_has_changes

    subroutine read_patch_save_state(document,kind,id,state)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: kind,id
        type(saved_patch_state), intent(out) :: state

        state%kind=kind
        state%id=id
        state%revision=document%patches%patch_revision(kind,id)
        state%name=document%patches%patch_name(kind,id)
        state%visible=studio_scene_is_visible(document%scene,kind,id,document%patches)
        call get_studio_material(document%scene,kind,id,state%rgb,state%opacity)
    end subroutine read_patch_save_state

    function document_display_name(document) result(name)
        type(studio_document_t),intent(in) :: document
        character(:),allocatable :: name

        name='Untitled project'
        if(allocated(document%project_path))then
            name=document%project_path(last_path_separator(document%project_path)+1:)
        else if(allocated(document%sample_title))then
            name=document%sample_title
        end if
    end function document_display_name

    subroutine record_document_merge(document,assemblies,meshes,operation)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: assemblies,meshes
        character(*), intent(in) :: operation
        type(studio_merge_edit) :: edit

        edit%stamp=next_patch_revision()
        edit%patch_stamp=single_patch_undo_stamp(document%patches)
        edit%assemblies=assemblies
        edit%meshes=meshes
        edit%operation=operation
        call append_merge_edit(document,edit)
        call clear_document_redo(document)
    end subroutine record_document_merge

    subroutine append_merge_edit(document,edit)
        type(studio_document_t), intent(inout) :: document
        type(studio_merge_edit), intent(in) :: edit

        if(.not.allocated(document%merge_history))allocate(document%merge_history(0))
        document%merge_history=[document%merge_history,edit]
    end subroutine append_merge_edit

    integer(int64) function document_merge_undo_stamp(document) result(stamp)
        type(studio_document_t), intent(in) :: document
        integer :: n

        stamp=0
        if(.not.allocated(document%merge_history))return
        n=size(document%merge_history)
        if(n==0)return
        associate(edit=>document%merge_history(n))
            if(usable_undo_stamp(edit%patch_stamp,catalog=document%patches)==0)return
            stamp=usable_undo_stamp(edit%stamp,catalog=document%patches)
        end associate
    end function document_merge_undo_stamp

    subroutine stage_document_mesh_removal(document,ids,inverse,staged,saved)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: ids(:)
        type(studio_mesh_edit), intent(out) :: inverse
        type(studio_mesh_restore), intent(out) :: staged
        logical, intent(out) :: saved
        type(studio_mesh_edit) :: removal

        saved=.true.
        if(size(ids)==0)return
        call capture_document_mesh_edit(document,ids,[integer::],inverse,saved)
        if(.not.saved)return
        call capture_document_mesh_edit(document,[integer::],ids,removal,saved)
        if(saved)call stage_document_mesh_edit(document,removal,staged,saved)
        if(.not.saved)call release_document_mesh_edit(inverse)
    end subroutine stage_document_mesh_removal

    subroutine undo_document_merge(document,message,replay)
        type(studio_document_t), intent(inout) :: document
        character(:), allocatable, intent(out) :: message
        type(studio_merge_replay), intent(out) :: replay
        type(studio_mesh_restore) :: staged
        type(studio_merge_edit) :: edit
        integer :: n,rank,group_count,mesh_count,i
        logical :: saved

        message=''
        if(document_merge_undo_stamp(document)==0_int64)return
        n=size(document%merge_history)
        edit=document%merge_history(n)
        if(single_patch_undo_stamp(document%patches)/=edit%patch_stamp)return
        group_count=count_multipatches(document%multipatches)
        mesh_count=document_mesh_count(document)
        if(group_count<edit%assemblies.or.mesh_count<edit%meshes)return
        if(group_count>edit%assemblies)then
            if(.not.save_group_history(document%multipatches%groups(edit%assemblies+1:)%definition,replay%groups,document%project_directory))then
                message='Could not record the imported assemblies for Redo.'
                return
            end if
        end if
        call stage_document_mesh_removal(document,[(i,i=edit%meshes+1,mesh_count)],replay%meshes,staged,saved)
        if(.not.saved)then
            message='Could not record the imported meshes for Redo.'
            return
        end if
        call undo_single_patch(rank,message,catalog=document%patches,replay=replay%patches)
        if(len(message)>0.or.rank==0)return
        call truncate_multipatches(document%multipatches,edit%assemblies)
        if(mesh_count>edit%meshes)call apply_document_mesh_edit(document,staged)
        document%merge_history=document%merge_history(:n-1)
        replay%edit=edit
        message=edit%operation//' undone'
    end subroutine undo_document_merge

    subroutine redo_document_merge(document,replay,changed,message)
        type(studio_document_t), intent(inout) :: document
        type(studio_merge_replay), intent(inout) :: replay
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: message
        type(multipatch_draft), allocatable :: definitions(:)
        type(studio_mesh_restore) :: staged
        logical :: saved

        changed=.false.
        message='Could not read the project import history.'
        if(replay%groups%recorded())then
            if(.not.load_group_history(replay%groups,definitions))return
        end if
        if(allocated(replay%meshes%prior_ids))then
            call stage_document_mesh_edit(document,replay%meshes,staged,saved)
            if(.not.saved)return
        end if
        call redo_single_patch(replay%patches,message,document%patches)
        if(len(message)>0)return
        if(allocated(definitions))call append_multipatch_definitions(document%multipatches,definitions,document%patches)
        if(allocated(replay%meshes%prior_ids))call apply_document_mesh_edit(document,staged)
        call append_merge_edit(document,replay%edit)
        message=replay%edit%operation//' redone'
        changed=.true.
    end subroutine redo_document_merge

    subroutine commit_document_deletion(document,patches,groups,meshes,mesh_only,error)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: patches(:,:),groups(:),meshes(:)
        logical, intent(in) :: mesh_only
        character(:), allocatable, intent(out) :: error
        type(studio_deletion_edit) :: edit
        type(studio_mesh_restore) :: removal
        type(multipatch_draft), allocatable :: definitions(:)
        integer :: i
        logical :: saved

        edit%groups=groups
        edit%mesh_only=mesh_only
        if(size(groups)>0)then
            allocate(definitions(size(groups)))
            do i=1,size(groups)
                definitions(i)=document%multipatches%groups(groups(i))%definition
            end do
            if(.not.save_group_history(definitions,edit%definitions,document%project_directory))then
                error='Could not record the selected multipatches for Undo.'
                return
            end if
        end if
        call stage_document_mesh_removal(document,meshes,edit%meshes,removal,saved)
        if(.not.saved)then
            call release_deletion_edit(edit)
            error='Could not record the affected meshes for Undo.'
            return
        end if
        call remove_single_patches(patches,error,document%patches)
        if(len(error)>0)then
            call release_deletion_edit(edit)
            return
        end if
        if(size(patches,2)>0)edit%patch_stamp=single_patch_undo_stamp(document%patches)
        call erase_multipatch_definitions(document%multipatches,groups)
        if(size(meshes)>0)call apply_document_mesh_edit(document,removal)
        edit%stamp=next_patch_revision()
        call append_deletion_edit(document,edit)
        call clear_document_redo(document)
    end subroutine commit_document_deletion

    subroutine append_deletion_edit(document,edit)
        type(studio_document_t), intent(inout) :: document
        type(studio_deletion_edit), intent(inout) :: edit

        if(.not.allocated(document%deletion_history))allocate(document%deletion_history(0))
        document%deletion_history=[document%deletion_history,edit]
        edit=studio_deletion_edit()
    end subroutine append_deletion_edit

    subroutine release_deletion_edit(edit)
        type(studio_deletion_edit), intent(inout) :: edit

        call release_history_blob(edit%definitions)
        call release_document_mesh_edit(edit%meshes)
        edit=studio_deletion_edit()
    end subroutine release_deletion_edit

    integer(int64) function document_deletion_undo_stamp(document) result(stamp)
        type(studio_document_t), intent(in) :: document
        integer :: n

        stamp=0_int64
        if(.not.allocated(document%deletion_history))return
        n=size(document%deletion_history)
        if(n==0)return
        associate(edit=>document%deletion_history(n))
            if(edit%patch_stamp>0.and.usable_undo_stamp(edit%patch_stamp,catalog=document%patches)==0)return
            stamp=usable_undo_stamp(edit%stamp,catalog=document%patches)
        end associate
    end function document_deletion_undo_stamp

    subroutine undo_document_deletion(document,mesh_only,message,replay)
        type(studio_document_t), intent(inout) :: document
        logical, intent(out) :: mesh_only
        character(:), allocatable, intent(out) :: message
        type(studio_deletion_replay), intent(out) :: replay
        type(studio_deletion_edit) :: edit
        type(studio_mesh_restore) :: staged
        type(multipatch_draft), allocatable :: definitions(:)
        integer :: n,rank
        logical :: saved

        mesh_only=.false.
        message='Could not read the deletion history.'
        if(document_deletion_undo_stamp(document)==0_int64)return
        n=size(document%deletion_history)
        edit=document%deletion_history(n)
        if(edit%definitions%recorded())then
            if(.not.load_group_history(edit%definitions,definitions))return
            if(size(definitions)/=size(edit%groups))return
        end if
        if(allocated(edit%meshes%prior_ids))then
            call stage_document_mesh_edit(document,edit%meshes,staged,saved)
            if(.not.saved)return
        end if
        if(edit%patch_stamp>0)then
            if(single_patch_undo_stamp(document%patches)/=edit%patch_stamp)return
            call undo_single_patch(rank,message,catalog=document%patches,replay=replay%patches)
            if(len(message)>0.or.rank==0)return
        end if
        if(allocated(definitions))call insert_multipatch_definitions(document%multipatches,edit%groups,definitions,document%patches)
        if(allocated(edit%meshes%prior_ids))call apply_document_mesh_edit(document,staged)
        document%deletion_history=document%deletion_history(:n-1)
        replay%edit=edit
        mesh_only=edit%mesh_only
        message='Deleted objects restored'
    end subroutine undo_document_deletion

    subroutine redo_document_deletion(document,replay,changed,message)
        type(studio_document_t), intent(inout) :: document
        type(studio_deletion_replay), intent(inout) :: replay
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: message
        type(studio_mesh_restore) :: staged
        type(studio_mesh_edit) :: removal
        logical :: saved

        changed=.false.
        message='Could not read the deletion history.'
        associate(edit=>replay%edit)
            if(allocated(edit%meshes%prior_ids))then
                removal%prior_ids=[integer::]
                removal%published_ids=edit%meshes%prior_ids
                call stage_document_mesh_edit(document,removal,staged,saved)
                if(.not.saved)return
            end if
            if(edit%patch_stamp>0)then
                call redo_single_patch(replay%patches,message,document%patches)
                if(len(message)>0)return
            end if
            call erase_multipatch_definitions(document%multipatches,edit%groups)
            if(allocated(edit%meshes%prior_ids))call apply_document_mesh_edit(document,staged)
            call append_deletion_edit(document,edit)
        end associate
        changed=.true.
        message='Deletion redone'
    end subroutine redo_document_deletion

    integer function document_mesh_count(document) result(n)
        type(studio_document_t), intent(in) :: document
        n=0
        if(allocated(document%meshes))n=size(document%meshes)
    end function document_mesh_count

    subroutine append_document_mesh(document,mesh)
        type(studio_document_t), intent(inout) :: document
        type(studio_mesh_t), intent(inout) :: mesh
        type(studio_mesh_t), allocatable :: grown(:)
        integer :: i,n

        n=document_mesh_count(document)
        mesh%content_revision=next_patch_revision()
        allocate(grown(n+1))
        do i=1,n
            call move_studio_mesh(document%meshes(i),grown(i))
        end do
        call move_studio_mesh(mesh,grown(n+1))
        call move_alloc(grown,document%meshes)
    end subroutine append_document_mesh

    subroutine remove_document_mesh(document,id)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: id
        type(studio_mesh_t), allocatable :: retained(:)
        integer :: i,j

        if(id<1.or.id>document_mesh_count(document))return
        allocate(retained(size(document%meshes)-1))
        j=0
        do i=1,size(document%meshes)
            if(i==id)cycle
            j=j+1
            call move_studio_mesh(document%meshes(i),retained(j))
        end do
        call move_alloc(retained,document%meshes)
    end subroutine remove_document_mesh

    subroutine record_document_mesh_edit(document,prior_ids,published_id,selection,mesh_only,saved)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: prior_ids(:),published_id,selection
        logical, intent(in) :: mesh_only
        logical, intent(out), optional :: saved
        type(studio_mesh_edit) :: edit
        logical :: recorded

        call capture_document_mesh_edit(document,prior_ids,pack([published_id],[published_id>0]),edit,recorded)
        if(present(saved))saved=recorded
        if(.not.recorded)return
        edit%stamp=next_patch_revision()
        edit%selection=selection
        edit%mesh_only=mesh_only
        call push_mesh_edit(document,edit)
        call clear_document_redo(document)
    end subroutine record_document_mesh_edit

    subroutine push_mesh_edit(document,edit)
        type(studio_document_t), intent(inout) :: document
        type(studio_mesh_edit), intent(inout) :: edit
        type(studio_mesh_edit), allocatable :: grown(:)
        integer :: count,capacity

        count=document%mesh_undo_count
        capacity=0
        if(allocated(document%mesh_history))capacity=size(document%mesh_history)
        if(count==capacity)then
            allocate(grown(max(8,2*capacity)))
            if(count>0)grown(:count)=document%mesh_history(:count)
            call move_alloc(grown,document%mesh_history)
        end if
        document%mesh_undo_count=count+1
        document%mesh_history(count+1)=edit
        edit=studio_mesh_edit()
    end subroutine push_mesh_edit

    subroutine clear_document_mesh_history(document)
        type(studio_document_t), intent(inout) :: document
        integer :: i

        do i=1,document%mesh_undo_count
            call release_document_mesh_edit(document%mesh_history(i))
        end do
        if(allocated(document%mesh_history))deallocate(document%mesh_history)
        document%mesh_undo_count=0
    end subroutine clear_document_mesh_history

    subroutine release_document_mesh_edit(edit)
        type(studio_mesh_edit), intent(inout) :: edit

        call release_history_blob(edit%payload)
        edit=studio_mesh_edit()
    end subroutine release_document_mesh_edit

    subroutine capture_document_mesh_edit(document,prior_ids,published_ids,edit,saved)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: prior_ids(:),published_ids(:)
        type(studio_mesh_edit), intent(out) :: edit
        logical, intent(out) :: saved
        integer :: unit,status,closed,i

        saved=.false.
        if(.not.valid_mesh_slots(prior_ids,document_mesh_count(document)))return
        edit%prior_ids=prior_ids
        edit%published_ids=published_ids
        if(size(prior_ids)==0)then
            saved=.true.
            return
        end if
        call create_history_blob(edit%payload,unit,status,document%project_directory)
        if(status/=0)return
        do i=1,size(prior_ids)
            call write_session_mesh(unit,document%meshes(prior_ids(i)),status)
            if(status/=0)exit
        end do
        close(unit,iostat=closed)
        saved=status==0.and.closed==0
        if(.not.saved)call release_document_mesh_edit(edit)
    end subroutine capture_document_mesh_edit

    subroutine stage_document_mesh_edit(document,edit,staged,saved)
        type(studio_document_t), intent(in) :: document
        type(studio_mesh_edit), intent(in) :: edit
        type(studio_mesh_restore), intent(out) :: staged
        logical, intent(out) :: saved
        integer :: unit,status,closed,i,mesh_count,retained,total,trailer_status
        character :: extra

        saved=.false.
        if(.not.allocated(edit%prior_ids).or..not.allocated(edit%published_ids))return
        mesh_count=document_mesh_count(document)
        if(.not.valid_mesh_slots(edit%published_ids,mesh_count))return
        retained=mesh_count-size(edit%published_ids)
        if(size(edit%prior_ids)>huge(total)-retained)return
        total=retained+size(edit%prior_ids)
        if(.not.valid_mesh_slots(edit%prior_ids,total))return
        allocate(staged%meshes(total),stat=status)
        if(status/=0)return
        allocate(staged%source_ids(retained),stat=status)
        if(status/=0)return
        allocate(staged%target_ids(retained),stat=status)
        if(status/=0)return
        call available_mesh_slots(mesh_count,edit%published_ids,staged%source_ids)
        call available_mesh_slots(total,edit%prior_ids,staged%target_ids)
        if(size(edit%prior_ids)==0)then
            saved=.true.
            return
        end if
        call open_history_blob(edit%payload,unit,status)
        if(status/=0)return
        do i=1,size(edit%prior_ids)
            call read_session_mesh(unit,staged%meshes(edit%prior_ids(i)),status)
            if(status/=0)exit
        end do
        if(status==0)then
            read(unit,iostat=trailer_status)extra
            status=merge(0,1,trailer_status==iostat_end)
        end if
        close(unit,iostat=closed)
        saved=status==0.and.closed==0
    end subroutine stage_document_mesh_edit

    subroutine apply_document_mesh_edit(document,staged)
        type(studio_document_t), intent(inout) :: document
        type(studio_mesh_restore), intent(inout) :: staged
        integer :: i

        do i=1,size(staged%source_ids)
            call move_studio_mesh(document%meshes(staged%source_ids(i)),staged%meshes(staged%target_ids(i)))
        end do
        call move_alloc(staged%meshes,document%meshes)
        document%mesh_generation=document%mesh_generation+1
    end subroutine apply_document_mesh_edit

    logical function valid_mesh_slots(slots,last) result(valid)
        integer, intent(in) :: slots(:),last
        integer :: i,previous

        valid=.false.
        previous=0
        do i=1,size(slots)
            if(slots(i)<=previous.or.slots(i)>last)return
            previous=slots(i)
        end do
        valid=.true.
    end function valid_mesh_slots

    subroutine available_mesh_slots(count,excluded,remaining)
        integer, intent(in) :: count,excluded(:)
        integer, intent(out) :: remaining(:)
        integer :: id,i,next

        next=1
        i=0
        do id=1,count
            if(next<=size(excluded))then
                if(id==excluded(next))then
                    next=next+1
                    cycle
                end if
            end if
            i=i+1
            remaining(i)=id
        end do
    end subroutine available_mesh_slots

    subroutine swap_document_mesh_edit(document,edit,inverse,changed,message)
        type(studio_document_t), intent(inout) :: document
        type(studio_mesh_edit), intent(in) :: edit
        type(studio_mesh_edit), intent(out) :: inverse
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: message
        type(studio_mesh_restore) :: staged

        message='Could not read the mesh history.'
        call stage_document_mesh_edit(document,edit,staged,changed)
        if(.not.changed)return
        message='Could not write the mesh history.'
        call capture_document_mesh_edit(document,edit%published_ids,edit%prior_ids,inverse,changed)
        if(.not.changed)return
        inverse%stamp=edit%stamp
        inverse%selection=edit%selection
        inverse%mesh_only=edit%mesh_only
        call apply_document_mesh_edit(document,staged)
        message=''
    end subroutine swap_document_mesh_edit

    integer function document_mesh_source_group(document,id) result(group)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: id
        integer :: i

        group=0
        if(id<1.or.id>document_mesh_count(document))return
        if(.not.allocated(document%meshes(id)%source_group))return
        do i=1,count_multipatches(document%multipatches)
            if(document%multipatches%groups(i)%definition%name/=document%meshes(id)%source_group)cycle
            group=i
            return
        end do
    end function document_mesh_source_group

    integer(int64) function document_mesh_undo_stamp(document) result(stamp)
        type(studio_document_t), intent(in) :: document

        stamp=0_int64
        if(document%mesh_undo_count==0)return
        stamp=usable_undo_stamp(document%mesh_history(document%mesh_undo_count)%stamp,catalog=document%patches)
    end function document_mesh_undo_stamp

    subroutine undo_document_mesh_edit(document,changed,selection,mesh_only,message,inverse)
        type(studio_document_t), intent(inout) :: document
        logical, intent(out) :: changed
        integer, intent(out) :: selection
        logical, intent(out) :: mesh_only
        character(:), allocatable, intent(out) :: message
        type(studio_mesh_edit), intent(out), optional :: inverse
        type(studio_mesh_edit) :: edit
        type(studio_mesh_restore) :: staged
        integer :: index

        changed=.false.
        selection=0
        mesh_only=.false.
        message=''
        if(document_mesh_undo_stamp(document)==0_int64)return
        index=document%mesh_undo_count
        edit=document%mesh_history(index)
        if(present(inverse))then
            call swap_document_mesh_edit(document,edit,inverse,changed,message)
            if(.not.changed)return
        else
            call stage_document_mesh_edit(document,edit,staged,changed)
            if(.not.changed)then
                message='Could not read the mesh history.'
                return
            end if
            call apply_document_mesh_edit(document,staged)
        end if

        selection=edit%selection
        mesh_only=edit%mesh_only
        if(size(edit%published_ids)==0)then
            message='Mesh deletion undone'
            if(selection==0)selection=edit%prior_ids(1)
        else if(size(edit%prior_ids)==0)then
            message='Mesh generation undone'
        else
            message='Remeshing undone'
        end if
        call release_document_mesh_edit(document%mesh_history(index))
        document%mesh_undo_count=index-1
    end subroutine undo_document_mesh_edit

    subroutine redo_document_mesh_edit(document,edit,changed,message)
        type(studio_document_t), intent(inout) :: document
        type(studio_mesh_edit), intent(inout) :: edit
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: message
        type(studio_mesh_edit) :: inverse

        call swap_document_mesh_edit(document,edit,inverse,changed,message)
        if(.not.changed)return
        call push_mesh_edit(document,inverse)
        call release_document_mesh_edit(edit)
    end subroutine redo_document_mesh_edit

    subroutine remove_construction_mesh(document,stamp,id)
        type(studio_document_t), intent(inout) :: document
        integer(int64), intent(in) :: stamp
        integer, intent(out) :: id

        id=0

        if(stamp<=0.or..not.allocated(document%meshes))return
        id=findloc(document%meshes%construction_stamp,stamp,dim=1)
        call remove_document_mesh(document,id)
    end subroutine remove_construction_mesh

    logical function document_mesh_valid(document,id) result(valid)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: id
        integer :: i,group

        valid=.false.
        if(id<1.or.id>document_mesh_count(document))return
        associate(mesh=>document%meshes(id))
            if(.not.allocated(mesh%parts))return
            if(mesh%group>0)then
                group=document_mesh_source_group(document,id)
                if(group==0)return
                associate(source=>document%multipatches%groups(group))
                    if(.not.source%validation%valid.or.source%definition%kind/=mesh%kind)return
                    if(size(source%definition%members)/=size(mesh%parts))return
                    if(any(source%definition%members/=mesh%parts%variant))return
                end associate
            end if
            do i=1,size(mesh%parts)
                if(single_patch_mesh_revision(mesh%kind,mesh%parts(i)%variant,catalog=document%patches)/=mesh%parts(i)%revision)return
            end do
        end associate
        valid=.true.
    end function document_mesh_valid

end module studio_document
