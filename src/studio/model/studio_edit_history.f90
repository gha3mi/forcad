module studio_edit_history
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_multipatch_collection, only: multipatch_undo_stamp,multipatch_construction_stamp,undo_multipatch,redo_multipatch
    use studio_scene, only: studio_material_undo_stamp,undo_studio_material,redo_studio_material, &
        snapshot_scene_properties,restore_scene_properties
    use studio_single_patch, only: single_patch_undo_stamp,undo_single_patch,redo_single_patch,get_patch_undo_targets
    use studio_document, only: studio_document_t,studio_redo_entry_t,studio_mesh_restore, &
        document_mesh_undo_stamp,undo_document_mesh_edit,redo_document_mesh_edit, &
        stage_document_mesh_edit,apply_document_mesh_edit,stage_document_mesh_removal, &
        document_merge_undo_stamp,undo_document_merge,redo_document_merge, &
        document_deletion_undo_stamp,undo_document_deletion,redo_document_deletion, &
        document_properties_undo_stamp,undo_document_properties,redo_document_properties, &
        remember_document_redo,take_document_redo,release_document_redo,discard_stale_redo,can_restore_document_redo
    use studio_property_stream, only: studio_property_values,save_property_history,load_property_history
    implicit none(type, external)
    private

    integer, parameter, public :: undo_mesh=1,undo_material=2,undo_patch=3,undo_assembly=4,undo_project_merge=5,undo_properties=6,undo_deletion=7

    type, public :: document_undo_result_t
        integer :: source=0,patch_kind=0,patch_axis=0,assembly_id=0,mesh_id=0
        integer :: removed_mesh_id=0
        logical :: changed=.false.,undid_construction=.false.,mesh_only=.false.
        character(:), allocatable :: operation,message
    end type document_undo_result_t

    public :: document_undo_stamp,undo_document_edit,document_can_redo,redo_document_edit,prune_document_redo

contains

    function document_history_stamps(document) result(stamps)
        type(studio_document_t), intent(in) :: document
        integer(int64) :: stamps(undo_deletion)

        stamps(undo_mesh)=document_mesh_undo_stamp(document)
        stamps(undo_material)=studio_material_undo_stamp(document%scene,document%patches)
        stamps(undo_assembly)=multipatch_undo_stamp(document%multipatches,document%patches)
        stamps(undo_patch)=single_patch_undo_stamp(document%patches)
        stamps(undo_project_merge)=document_merge_undo_stamp(document)
        stamps(undo_properties)=document_properties_undo_stamp(document)
        stamps(undo_deletion)=document_deletion_undo_stamp(document)
    end function document_history_stamps

    integer(int64) function document_undo_stamp(document) result(stamp)
        type(studio_document_t), intent(in) :: document

        stamp=maxval(document_history_stamps(document))
    end function document_undo_stamp

    logical function document_can_redo(document) result(available)
        type(studio_document_t), intent(in) :: document

        available=can_restore_document_redo(document,document_undo_stamp(document))
    end function document_can_redo

    subroutine prune_document_redo(document)
        type(studio_document_t), intent(inout) :: document

        call discard_stale_redo(document,document_undo_stamp(document))
    end subroutine prune_document_redo

    logical function redo_document_edit(document,error) result(changed)
        type(studio_document_t), intent(inout) :: document
        character(:), allocatable, optional, intent(out) :: error
        type(studio_redo_entry_t), allocatable :: replay
        integer(int64) :: stamp
        character(:), allocatable :: message

        changed=.false.
        if(present(error))error=''
        stamp=document_undo_stamp(document)
        call discard_stale_redo(document,stamp)
        call take_document_redo(document,replay)
        if(.not.allocated(replay))return

        call replay_document_edit(document,replay,changed,message)
        if(changed)then
            call release_document_redo(replay)
        else
            call remember_document_redo(document,replay,stamp)
            if(present(error))error=message
        end if
    end function redo_document_edit

    subroutine replay_document_edit(document,replay,changed,message)
        type(studio_document_t), intent(inout) :: document
        type(studio_redo_entry_t), intent(inout) :: replay
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: message
        type(studio_property_values) :: appearance
        type(studio_mesh_restore) :: meshes
        logical :: saved

        changed=.false.
        message='Could not read appearance history.'
        if(replay%appearance%recorded())then
            if(.not.load_property_history(replay%appearance,appearance))return
        end if
        message=''
        select case(replay%source)
        case(undo_mesh)
            call redo_document_mesh_edit(document,replay%meshes,changed,message)
        case(undo_material)
            changed=redo_studio_material(document%scene,replay%material,document%patches)
            if(.not.changed)message='Could not restore appearance history.'
        case(undo_patch)
            call redo_single_patch(replay%patches,message,document%patches)
            changed=len(message)==0
        case(undo_assembly)
            if(allocated(replay%meshes%prior_ids))then
                call stage_document_mesh_edit(document,replay%meshes,meshes,saved)
                if(.not.saved)then
                    message='Could not read the assembly mesh history.'
                    return
                end if
            end if
            call redo_multipatch(document%multipatches,replay%assembly,changed,message,document%patches)
            if(changed.and.allocated(replay%meshes%prior_ids))call apply_document_mesh_edit(document,meshes)
        case(undo_project_merge)
            call redo_document_merge(document,replay%merge,changed,message)
        case(undo_deletion)
            call redo_document_deletion(document,replay%deletion,changed,message)
        case(undo_properties)
            call redo_document_properties(document,replay%properties,changed,message)
        case default
            message='The Redo entry has an unknown action.'
        end select
        if(changed.and.allocated(appearance%patches))call restore_scene_properties(document%scene,appearance%patches,document%patches)
    end subroutine replay_document_edit

    subroutine undo_document_edit(document,undone)
        type(studio_document_t), intent(inout) :: document
        type(document_undo_result_t), intent(out) :: undone
        integer(int64) :: stamps(undo_deletion),stamp,construction_stamp
        type(studio_redo_entry_t), allocatable :: replay
        type(studio_mesh_restore) :: meshes
        logical :: saved

        undone%operation=''
        undone%message=''
        stamps=document_history_stamps(document)
        stamp=maxval(stamps)
        call discard_stale_redo(document,stamp)
        if(stamp==0_int64)return

        allocate(replay)
        undone%source=maxloc(stamps,dim=1)
        replay%source=undone%source
        construction_stamp=0
        if(undone%source==undo_assembly)construction_stamp=multipatch_construction_stamp(document%multipatches,document%patches)
        if(undone%source==undo_patch.or.undone%source==undo_project_merge.or.construction_stamp>0)then
            if(.not.record_redo_appearance(document,replay))then
                undone%message='Could not record the scene appearance for Redo.'
                return
            end if
        end if

        select case(undone%source)
        case(undo_mesh)
            call undo_document_mesh_edit(document,undone%changed,undone%mesh_id,undone%mesh_only,undone%message,replay%meshes)
        case(undo_material)
            call undo_studio_material(document%scene,document%patches,replay%material)
            if(replay%material%stamp==0)undone%message='Could not undo the appearance edit.'
        case(undo_assembly)
            if(construction_stamp>0.and.allocated(document%meshes))then
                undone%removed_mesh_id=findloc(document%meshes%construction_stamp,construction_stamp,dim=1)
            end if
            if(undone%removed_mesh_id>0)then
                call stage_document_mesh_removal(document,[undone%removed_mesh_id],replay%meshes,meshes,saved)
                if(.not.saved)then
                    undone%message='Could not record the assembly mesh for Redo.'
                    call release_document_redo(replay)
                    return
                end if
            end if
            call undo_multipatch(document%multipatches,undone%assembly_id,undone%message, &
                creation_stamp=construction_stamp,catalog=document%patches,replay=replay%assembly)
            undone%undid_construction=construction_stamp>0
            if(undone%undid_construction.and.undone%removed_mesh_id>0)call apply_document_mesh_edit(document,meshes)
        case(undo_patch)
            call undo_single_patch(undone%patch_kind,undone%message,undone%operation,undone%patch_axis, &
                undone%assembly_id,catalog=document%patches,replay=replay%patches)
        case(undo_project_merge)
            call undo_document_merge(document,undone%message,replay%merge)
        case(undo_deletion)
            call undo_document_deletion(document,undone%mesh_only,undone%message,replay%deletion)
        case(undo_properties)
            call undo_document_properties(document,replay%properties,undone%changed,undone%message)
        case default
            undone%message='The Undo entry has an unknown action.'
        end select
        stamps=document_history_stamps(document)
        undone%changed=stamps(undone%source)/=stamp
        if(undone%changed)then
            call remember_document_redo(document,replay,maxval(stamps))
        else
            call release_document_redo(replay)
        end if
    end subroutine undo_document_edit

    logical function record_redo_appearance(document,replay) result(saved)
        type(studio_document_t), intent(in) :: document
        type(studio_redo_entry_t), intent(inout) :: replay
        type(studio_property_values) :: appearance
        integer, allocatable :: targets(:,:)

        call get_patch_undo_targets(document%patches,targets)
        saved=.true.
        if(size(targets,2)==0)return
        appearance%patches=snapshot_scene_properties(document%scene,targets)
        saved=save_property_history(appearance,replay%appearance,document%project_directory)
    end function record_redo_appearance

end module studio_edit_history
