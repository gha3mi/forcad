! Scene links, native validation snapshots, and Undo for multipatch definitions.
module studio_multipatch_collection
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_multipatch, only: multipatch_draft, multipatch_join, multipatch_result, inspect_multipatch
    use studio_single_patch, only: single_patch_revision, single_patch_mesh_revision, next_patch_revision, usable_undo_stamp, &
        studio_native_patch, publish_native_patches, single_patch_undo_stamp, undo_single_patch, studio_patch_catalog, &
        studio_patch_replay,redo_single_patch,release_patch_replay,commit_native_patch_edits,patch_history_project
    use studio_presentation_text, only: integer_text
    use studio_mesh, only: studio_mesh_t, build_planned_volume_mesh
    use studio_history_storage, only: history_blob,release_history_blob
    use studio_multipatch_stream, only: save_group_history,load_group_history
    implicit none(type, external)
    private
    integer, parameter :: created=1, revised=2, removed=3

    type :: multipatch_entry
        type(multipatch_draft) :: definition
        type(multipatch_result) :: validation
        class(*), allocatable :: native
        integer(int64), allocatable :: revisions(:)
    end type multipatch_entry

    type :: multipatch_edit
        integer :: action=0,index=0
        integer(int64) :: stamp=0_int64, batch_stamp=0_int64
        type(history_blob) :: previous
    end type multipatch_edit

    type, public :: multipatch_replay
        type(multipatch_edit), private :: edit
        type(history_blob), private :: displaced
        type(studio_patch_replay), private :: patches
    end type multipatch_replay

    type, public :: multipatch_collection
        type(multipatch_entry), allocatable :: groups(:)
        type(multipatch_edit), allocatable, private :: history(:)
    end type multipatch_collection

    public :: count_multipatches, save_multipatch, remove_multipatch, refresh_multipatch
    public :: multipatch_undo_stamp, undo_multipatch, multipatch_contains
    public :: restore_multipatch_definitions, publish_multipatch_volume, move_multipatch_collection
    public :: redo_multipatch,release_multipatch_replay,clear_multipatch_history,multipatch_construction_stamp
    public :: append_multipatch_definitions,truncate_multipatches
    public :: erase_multipatch_definitions,insert_multipatch_definitions
    public :: commit_multipatch_patch_edits, validate_multipatch_patch_edits
contains
    subroutine commit_multipatch_patch_edits(collection,id,kind,targets,candidates,operation,catalog,changed,error,preserve_mesh)
        type(multipatch_collection), intent(inout) :: collection
        integer, intent(in) :: id,kind,targets(:)
        type(studio_native_patch), intent(inout) :: candidates(:)
        character(*), intent(in) :: operation
        type(studio_patch_catalog), intent(inout), target :: catalog
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: error
        logical, intent(in), optional :: preserve_mesh
        integer :: group,i

        changed=.false.
        call validate_multipatch_patch_edits(collection,id,kind,targets,candidates,catalog,error)
        if(len(error)>0)return
        call commit_native_patch_edits(kind,targets,candidates,id,changed,error,catalog,operation,preserve_mesh)
        if(.not.changed)return
        do group=1,count_multipatches(collection)
            associate(definition=>collection%groups(group)%definition)
                if(definition%kind/=kind)cycle
                if(.not.any([(any(definition%members==targets(i)),i=1,size(targets))]))cycle
            end associate
            call refresh_multipatch(collection,group,.false.,catalog=catalog)
        end do
    end subroutine commit_multipatch_patch_edits

    subroutine validate_multipatch_patch_edits(collection,id,kind,targets,candidates,catalog,error)
        type(multipatch_collection), intent(in) :: collection
        integer, intent(in) :: id,kind,targets(:)
        type(studio_native_patch), intent(in) :: candidates(:)
        type(studio_patch_catalog), intent(in) :: catalog
        character(:), allocatable, intent(out) :: error
        type(multipatch_result) :: checked
        integer :: group,join

        do group=1,count_multipatches(collection)
            associate(definition=>collection%groups(group)%definition)
                if(definition%kind/=kind)cycle
                if(group/=id)then
                    if(.not.allocated(definition%joins))cycle
                    if(.not.any([(touches_join(definition%joins(join),targets),join=1,size(definition%joins))]))cycle
                end if
                call inspect_multipatch(definition,checked,patches=candidates,catalog=catalog,patch_ids=targets)
                if(.not.checked%valid)then
                    error='The edit would violate a join in '//definition%name//'.'
                    do join=1,size(checked%joins)
                        if(checked%joins(join)%valid)cycle
                        error=error//' Join '//integer_text(join)//' · '//checked%joins(join)%message
                        exit
                    end do
                    return
                end if
            end associate
        end do

        error=''
    end subroutine validate_multipatch_patch_edits

    logical function touches_join(join,targets) result(touches)
        type(multipatch_join), intent(in) :: join
        integer, intent(in) :: targets(:)

        touches=join%continuity>=0.and.any(targets==join%patch_a.or.targets==join%patch_b)
    end function touches_join

    subroutine move_multipatch_collection(source,destination)
        type(multipatch_collection), intent(inout) :: source
        type(multipatch_collection), intent(out) :: destination

        call move_alloc(source%groups,destination%groups)
        call move_alloc(source%history,destination%history)
    end subroutine move_multipatch_collection

    subroutine restore_multipatch_definitions(collection,definitions,catalog)
        type(multipatch_collection), intent(out) :: collection
        type(multipatch_draft), intent(in) :: definitions(:)
        type(studio_patch_catalog), intent(in), optional :: catalog
        integer :: i

        allocate(collection%groups(size(definitions)))
        do i=1,size(definitions)
            collection%groups(i)%definition=definitions(i)
            call update_entry(collection%groups(i),.true.,catalog)
        end do
    end subroutine restore_multipatch_definitions

    integer function count_multipatches(collection) result(n)
        type(multipatch_collection), intent(in) :: collection
        n=0
        if(allocated(collection%groups))n=size(collection%groups)
    end function count_multipatches

    logical function multipatch_contains(collection,id,kind,variant) result(included)
        type(multipatch_collection), intent(in) :: collection
        integer, intent(in) :: id,kind,variant
        included=.false.
        if(id<1.or.id>count_multipatches(collection))return
        associate(definition=>collection%groups(id)%definition)
            included=definition%kind==kind.and.any(definition%members==variant)
        end associate
    end function multipatch_contains

    subroutine update_entry(entry,force,catalog)
        type(multipatch_entry), intent(inout) :: entry
        logical, intent(in) :: force
        type(studio_patch_catalog), intent(in), optional :: catalog
        integer(int64), allocatable :: revisions(:)
        integer :: i

        associate(definition=>entry%definition)
            revisions=[(single_patch_revision(definition%kind,definition%members(i),catalog=catalog),i=1,size(definition%members))]
            if(allocated(entry%revisions).and..not.force)then
                if(size(entry%revisions)==size(revisions))then
                    if(all(entry%revisions==revisions))return
                end if
            end if
            call inspect_multipatch(definition,entry%validation,entry%native,catalog=catalog)
            call move_alloc(revisions,entry%revisions)
        end associate
    end subroutine update_entry

    subroutine refresh_multipatch(collection,id,force,catalog)
        type(multipatch_collection), intent(inout) :: collection
        integer, intent(in) :: id
        logical, intent(in) :: force
        type(studio_patch_catalog), intent(in), optional :: catalog
        if(id<1.or.id>count_multipatches(collection))return
        call update_entry(collection%groups(id),force,catalog)
    end subroutine refresh_multipatch

    subroutine remember(collection,edit)
        type(multipatch_collection), intent(inout) :: collection
        type(multipatch_edit), intent(in) :: edit

        if(.not.allocated(collection%history))allocate(collection%history(0))
        collection%history=[collection%history,edit]
    end subroutine remember

    subroutine save_multipatch(collection,id,definition,error,catalog)
        type(multipatch_collection), intent(inout) :: collection
        integer, intent(inout) :: id
        type(multipatch_draft), intent(in) :: definition
        character(:), allocatable, intent(out) :: error
        type(studio_patch_catalog), intent(inout), optional :: catalog
        type(multipatch_entry) :: entry
        type(multipatch_edit) :: edit
        integer :: n,i

        n=count_multipatches(collection)
        error='Select a valid multipatch.'
        if(id<0.or.id>n)return
        if(.not.allocated(definition%name).or..not.allocated(definition%members))return
        error='Enter a name for the multipatch.'
        if(len_trim(definition%name)==0)return
        do i=1,n
            if(i==id)cycle
            if(collection%groups(i)%definition%name/=trim(definition%name))cycle
            error='A multipatch already uses this name.'
            return
        end do
        entry%definition=definition
        entry%definition%name=trim(definition%name)
        call update_entry(entry,.true.,catalog)
        if(.not.entry%validation%valid)then
            error=entry%validation%message
            return
        end if
        call store_entry(collection,id,entry,edit,error,catalog)
    end subroutine save_multipatch

    subroutine publish_multipatch_volume(collection,name,patches,names,joins,id,error,mesh_samples,mesh,catalog)
        type(multipatch_collection), intent(inout) :: collection
        character(*), intent(in) :: name,names(:)
        type(studio_native_patch), intent(in) :: patches(:)
        type(multipatch_join), intent(in) :: joins(:)
        integer, intent(out) :: id
        character(:), allocatable, intent(out) :: error
        type(studio_patch_catalog), intent(inout), optional :: catalog
        type(multipatch_entry) :: entry
        type(multipatch_edit) :: edit
        integer :: first,last,kind,i
        integer :: resolution(3,size(patches))
        integer, optional, intent(in) :: mesh_samples(3)
        type(studio_mesh_t), optional, intent(out) :: mesh

        id=0
        error='Enter a name for the multipatch.'
        if(len_trim(name)==0)return
        do i=1,count_multipatches(collection)
            if(collection%groups(i)%definition%name/=trim(name))cycle
            error='A multipatch already uses the name '//trim(name)//'. Ask for a different name.'
            return
        end do

        ! Validate the complete native group before publishing any catalog rows.
        entry%definition%kind=3
        entry%definition%name=trim(name)
        entry%definition%members=[(i,i=1,size(patches))]
        entry%definition%joins=joins
        call inspect_multipatch(entry%definition,entry%validation,entry%native,patches)
        if(.not.entry%validation%valid)then
            error=entry%validation%message
            return
        end if

        ! Sample the approved mesh before either object enters the project.
        resolution=4
        if(present(mesh_samples).and.present(mesh))then
            call build_planned_volume_mesh(mesh_samples,joins,entry%native,mesh,error)
            if(len(error)>0)return
            resolution=7
        end if
        call publish_native_patches(patches,names,kind,last,error,resolutions=resolution, &
            operation='multipatch construction',catalog=catalog)
        if(len(error)>0)return
        first=last-size(patches)+1
        entry%definition%members=entry%definition%members+first-1
        entry%definition%joins%patch_a=entry%definition%joins%patch_a+first-1
        entry%definition%joins%patch_b=entry%definition%joins%patch_b+first-1
        entry%revisions=[(single_patch_revision(kind,i,catalog=catalog),i=first,last)]
        edit%batch_stamp=single_patch_undo_stamp(catalog=catalog)
        call store_entry(collection,id,entry,edit,error,catalog)
        if(len(error)>0)return
        if(present(mesh_samples).and.present(mesh))then
            mesh%group=id
            mesh%name=trim(name)//' · Hex mesh'
            mesh%source_group=trim(name)
            mesh%parts%variant=collection%groups(id)%definition%members
            mesh%parts%revision=[(single_patch_mesh_revision(kind,i,catalog=catalog),i=first,last)]
            mesh%construction_stamp=multipatch_undo_stamp(collection,catalog)
        end if
    end subroutine publish_multipatch_volume

    subroutine store_entry(collection,id,entry,edit,error,catalog)
        type(multipatch_collection), intent(inout) :: collection
        integer, intent(inout) :: id
        type(multipatch_entry), intent(inout) :: entry
        type(multipatch_edit), intent(inout) :: edit
        character(:), allocatable, intent(out) :: error
        type(studio_patch_catalog), intent(in), optional :: catalog

        error=''
        if(id==0)then
            id=count_multipatches(collection)+1
            edit%action=created
            call insert_group(collection,id,entry)
        else
            edit%action=revised
            if(.not.save_group_history([collection%groups(id)%definition],edit%previous,patch_history_project(catalog)))then
                error='Could not record the multipatch changes on disk.'
                return
            end if
            call move_group(entry,collection%groups(id))
        end if
        edit%index=id
        edit%stamp=next_patch_revision()
        call remember(collection,edit)
    end subroutine store_entry

    subroutine remove_multipatch(collection,id,catalog,error)
        type(multipatch_collection), intent(inout) :: collection
        integer, intent(in) :: id
        type(studio_patch_catalog), intent(inout), optional :: catalog
        character(:), allocatable, intent(out), optional :: error
        type(multipatch_edit) :: edit

        if(present(error))error='Select a valid multipatch.'
        if(id<1.or.id>count_multipatches(collection))return
        if(.not.save_group_history([collection%groups(id)%definition],edit%previous,patch_history_project(catalog)))then
            if(present(error))error='Could not record the removed multipatch on disk.'
            return
        end if
        edit%action=removed
        edit%index=id
        edit%stamp=next_patch_revision()
        call remember(collection,edit)
        call remove_group(collection,id)
        if(present(error))error=''
    end subroutine remove_multipatch

    integer(int64) function multipatch_undo_stamp(collection,catalog) result(stamp)
        type(multipatch_collection), intent(in) :: collection
        type(studio_patch_catalog), intent(in), optional :: catalog

        stamp=0
        if(.not.allocated(collection%history))return
        if(size(collection%history)==0)return
        associate(edit=>collection%history(size(collection%history)))
            if(edit%batch_stamp>0)then
                if(usable_undo_stamp(edit%batch_stamp,catalog=catalog)==0)return
            end if
            stamp=usable_undo_stamp(edit%stamp,catalog=catalog)
        end associate
    end function multipatch_undo_stamp

    integer(int64) function multipatch_construction_stamp(collection,catalog) result(stamp)
        type(multipatch_collection), intent(in) :: collection
        type(studio_patch_catalog), intent(in), optional :: catalog

        stamp=0
        if(multipatch_undo_stamp(collection,catalog)==0)return
        associate(edit=>collection%history(size(collection%history)))
            if(edit%batch_stamp>0)stamp=edit%stamp
        end associate
    end function multipatch_construction_stamp

    subroutine undo_multipatch(collection,selection,message,creation_stamp,catalog,replay)
        type(multipatch_collection), intent(inout) :: collection
        integer, intent(out) :: selection
        character(:), allocatable, intent(out) :: message
        integer(int64), optional, intent(out) :: creation_stamp
        type(studio_patch_catalog), intent(inout), optional :: catalog
        type(multipatch_replay), intent(out), optional :: replay
        type(multipatch_entry) :: restored
        type(multipatch_draft), allocatable :: definitions(:)
        type(multipatch_edit) :: edit
        type(history_blob) :: displaced
        integer :: at,id,kind
        character(:), allocatable :: error

        selection=0
        message=''
        if(present(creation_stamp))creation_stamp=0
        if(multipatch_undo_stamp(collection,catalog)==0)return
        at=size(collection%history)
        edit=collection%history(at)
        id=edit%index
        if(edit%action/=created)then
            if(.not.load_group_history(edit%previous,definitions))then
                message='Could not read the multipatch history.'
                return
            end if
            if(size(definitions)/=1)return
            restored%definition=definitions(1)
            call update_entry(restored,.true.,catalog)
        end if
        if(present(replay).and.edit%action/=removed)then
            if(.not.save_group_history([collection%groups(id)%definition],displaced,patch_history_project(catalog)))then
                message='Could not write the multipatch redo data.'
                return
            end if
        end if

        select case(edit%action)
        case(created)
            message='Multipatch creation undone'
            if(edit%batch_stamp>0)then
                if(single_patch_undo_stamp(catalog=catalog)/=edit%batch_stamp)then
                    message='Undo the later geometry edit first.'
                    call release_history_blob(displaced)
                    return
                end if
                if(present(replay))then
                    call undo_single_patch(kind,error,catalog=catalog,replay=replay%patches)
                else
                    call undo_single_patch(kind,error,catalog=catalog)
                end if
                if(len(error)>0.or.kind==0)then
                    message=error
                    call release_history_blob(displaced)
                    return
                end if
                message='Multipatch construction undone'
                if(present(creation_stamp))creation_stamp=edit%stamp
            end if
            call remove_group(collection,id)
        case(revised)
            call move_group(restored,collection%groups(id))
            message='Multipatch changes undone'
            selection=id
        case(removed)
            call insert_group(collection,id,restored)
            message='Multipatch restored'
            selection=id
        case default
            message='The multipatch Undo entry has an unknown action.'
            call release_history_blob(displaced)
            return
        end select
        collection%history=collection%history(:at-1)
        if(present(replay))then
            replay%edit=edit
            replay%displaced=displaced
        else
            call release_history_blob(edit%previous)
        end if
    end subroutine undo_multipatch

    subroutine redo_multipatch(collection,replay,changed,message,catalog)
        type(multipatch_collection), intent(inout) :: collection
        type(multipatch_replay), intent(inout) :: replay
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: message
        type(studio_patch_catalog), intent(inout), optional :: catalog
        type(multipatch_draft), allocatable :: definitions(:)
        type(multipatch_entry) :: restored
        integer :: id

        changed=.false.
        message='Could not read the multipatch redo data.'
        id=replay%edit%index
        if(replay%edit%action/=removed)then
            if(.not.load_group_history(replay%displaced,definitions))return
            if(size(definitions)/=1)return
            restored%definition=definitions(1)
        end if
        if(replay%edit%batch_stamp>0)then
            call redo_single_patch(replay%patches,message,catalog)
            if(len(message)>0)return
        end if
        select case(replay%edit%action)
        case(created)
            call update_entry(restored,.true.,catalog)
            call insert_group(collection,id,restored)
        case(revised)
            call update_entry(restored,.true.,catalog)
            call move_group(restored,collection%groups(id))
        case(removed)
            call remove_group(collection,id)
        case default
            message='The multipatch Redo entry has an unknown action.'
            return
        end select
        call remember(collection,replay%edit)
        replay%edit=multipatch_edit()
        call release_history_blob(replay%displaced)
        changed=.true.
        message=''
    end subroutine redo_multipatch

    subroutine release_multipatch_replay(replay)
        type(multipatch_replay), intent(inout) :: replay

        call release_history_blob(replay%edit%previous)
        call release_history_blob(replay%displaced)
        call release_patch_replay(replay%patches)
    end subroutine release_multipatch_replay

    subroutine clear_multipatch_history(collection)
        type(multipatch_collection), intent(inout) :: collection
        integer :: i

        if(.not.allocated(collection%history))return
        do i=1,size(collection%history)
            call release_history_blob(collection%history(i)%previous)
        end do
        deallocate(collection%history)
    end subroutine clear_multipatch_history

    subroutine move_group(source,destination)
        type(multipatch_entry), intent(inout) :: source
        type(multipatch_entry), intent(out) :: destination

        destination%definition%kind=source%definition%kind
        call move_alloc(source%definition%name,destination%definition%name)
        call move_alloc(source%definition%members,destination%definition%members)
        call move_alloc(source%definition%joins,destination%definition%joins)
        destination%validation=source%validation
        call move_alloc(source%native,destination%native)
        call move_alloc(source%revisions,destination%revisions)
    end subroutine move_group

    subroutine insert_group(collection,id,entry)
        type(multipatch_collection), intent(inout) :: collection
        integer, intent(in) :: id
        type(multipatch_entry), intent(inout) :: entry
        type(multipatch_entry), allocatable :: resized(:)
        integer :: n,i

        n=count_multipatches(collection)
        allocate(resized(n+1))
        do i=1,id-1
            call move_group(collection%groups(i),resized(i))
        end do
        call move_group(entry,resized(id))
        do i=id,n
            call move_group(collection%groups(i),resized(i+1))
        end do
        call move_alloc(resized,collection%groups)
    end subroutine insert_group

    subroutine remove_group(collection,id)
        type(multipatch_collection), intent(inout) :: collection
        integer, intent(in) :: id
        type(multipatch_entry), allocatable :: resized(:)
        integer :: n,i

        n=count_multipatches(collection)
        allocate(resized(n-1))
        do i=1,id-1
            call move_group(collection%groups(i),resized(i))
        end do
        do i=id+1,n
            call move_group(collection%groups(i),resized(i-1))
        end do
        call move_alloc(resized,collection%groups)
    end subroutine remove_group

    subroutine truncate_multipatches(collection,n)
        type(multipatch_collection), intent(inout) :: collection
        integer, intent(in) :: n
        type(multipatch_entry), allocatable :: retained(:)
        integer :: i

        if(n==count_multipatches(collection))return
        allocate(retained(n))
        do i=1,n
            call move_group(collection%groups(i),retained(i))
        end do
        call move_alloc(retained,collection%groups)
    end subroutine truncate_multipatches

    subroutine append_multipatch_definitions(collection,definitions,catalog)
        type(multipatch_collection), intent(inout) :: collection
        type(multipatch_draft), intent(in) :: definitions(:)
        type(studio_patch_catalog), intent(in), optional :: catalog
        type(multipatch_entry), allocatable :: groups(:)
        integer :: i,n

        if(size(definitions)==0)return
        n=count_multipatches(collection)
        allocate(groups(n+size(definitions)))
        do i=1,n
            call move_group(collection%groups(i),groups(i))
        end do
        do i=1,size(definitions)
            groups(n+i)%definition=definitions(i)
            call update_entry(groups(n+i),.true.,catalog)
        end do
        call move_alloc(groups,collection%groups)
    end subroutine append_multipatch_definitions

    subroutine erase_multipatch_definitions(collection,positions)
        type(multipatch_collection), intent(inout) :: collection
        integer, intent(in) :: positions(:)
        type(multipatch_entry), allocatable :: retained(:)
        logical, allocatable :: deleted(:)
        integer :: size_before,source,target

        if(size(positions)==0)return
        size_before=count_multipatches(collection)
        allocate(deleted(size_before),source=.false.)
        deleted(positions)=.true.
        allocate(retained(size_before-size(positions)))
        target=0
        do source=1,size_before
            if(deleted(source))cycle
            target=target+1
            call move_group(collection%groups(source),retained(target))
        end do
        call move_alloc(retained,collection%groups)
    end subroutine erase_multipatch_definitions

    subroutine insert_multipatch_definitions(collection,positions,definitions,catalog)
        type(multipatch_collection), intent(inout) :: collection
        integer, intent(in) :: positions(:)
        type(multipatch_draft), intent(in) :: definitions(:)
        type(studio_patch_catalog), intent(in) :: catalog
        type(multipatch_entry), allocatable :: restored(:)
        logical, allocatable :: inserted(:)
        integer :: size_after,source,target,i

        if(size(positions)==0)return
        size_after=count_multipatches(collection)+size(positions)
        allocate(restored(size_after),inserted(size_after))
        inserted=.false.
        inserted(positions)=.true.
        do i=1,size(positions)
            target=positions(i)
            restored(target)%definition=definitions(i)
            call update_entry(restored(target),.true.,catalog)
        end do
        source=0
        do target=1,size_after
            if(inserted(target))cycle
            source=source+1
            call move_group(collection%groups(source),restored(target))
        end do
        call move_alloc(restored,collection%groups)
    end subroutine insert_multipatch_definitions
end module studio_multipatch_collection
